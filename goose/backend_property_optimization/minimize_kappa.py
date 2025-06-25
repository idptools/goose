'''
Code for kappa minimzation. Going to attempt a different strategy
because the current one is insufficient. 
'''
import random
import numpy as np
from goose.backend_property_calculation.calculate_kappa import global_asymmetry, patterning_asymmetry, charge_matrix_to_ternary, kappa, analytical_deltamax
from goose.backend_property_calculation.calculate_properties_vectorized import sequences_to_matrices

def calc_global_sigma(ternarized_seqs):
    """
    Calculate the global sigma value for a set of ternarized sequences.
    
    Parameters
    ----------
    ternarized_seqs : np.ndarray
        2D array of ternarized sequences (shape: n_sequences x sequence_length)

    Returns
    -------
    np.ndarray
        Global asymmetry value for each sequence.
    """
    return global_asymmetry(ternarized_seqs)
    

def calc_sigma_windows(ternarized_seqs, window_size, overall_asymmetry):
    '''
    get sigma value windows for a sequence 
    Parameters
    ----------
    ternarized_seqs : np.ndarray
        2D array of ternarized sequences (shape: n_sequences x sequence_length)
    window_size : int
        Size of the sliding window to calculate local asymmetry
    overall_asymmetry : np.ndarray
        Global asymmetry value for each sequence, used for comparison
    
    Returns
    -------
    np.ndarray
        2D array of squared differences between local window asymmetry and overall asymmetry
        (shape: n_sequences x n_windows)

    '''
    n_sequences, winged_length = ternarized_seqs.shape
    n_windows = winged_length - window_size + 1
    
    # Create all windows at once using vectorized sliding window
    windows_3d = np.lib.stride_tricks.sliding_window_view(
        ternarized_seqs, window_shape=window_size, axis=1
    )    
    # Count charges across all windows
    window_count_pos = np.sum(windows_3d == 1, axis=2)
    window_count_neg = np.sum(windows_3d == -1, axis=2)
    
    # Calculate local fractions
    fa_local = window_count_pos / window_size
    fb_local = window_count_neg / window_size

    denominator = fa_local + fb_local
    window_asymmetry = np.zeros((n_sequences, n_windows))
    valid_mask = denominator > 0
    window_asymmetry[valid_mask] = ((fa_local[valid_mask] - fb_local[valid_mask])**2) / denominator[valid_mask]
    # Calculate squared differences from overall asymmetry
    overall_asymmetry_broadcast = overall_asymmetry[:, np.newaxis]
    squared_diffs = (window_asymmetry - overall_asymmetry_broadcast)**2
    return squared_diffs


def efficient_kappa_minimization(seq_matrices, max_iterations=1000):
    """
    Efficiently minimize kappa by targeting the most problematic windows.
    
    Parameters
    ----------
    seq_matrices : np.ndarray
        2D array of sequence matrices (shape: n_sequences x sequence_length)
    max_iterations : int, default=1000
        Maximum number of optimization iterations
        
    Returns
    -------
    optimized_sequences : np.ndarray
        Sequences with minimized kappa values
    kappa_scores : np.ndarray
        Final kappa scores for each sequence
    """
    # Convert sequence matrices to ternary sequences
    ternary_seqs = seq_matrices

    # get num seqs, seq_length
    n_sequences, seq_length = ternary_seqs.shape
    
    # Calculate initial kappa scores
    current_kappa = kappa(seq_matrices)

    # set best_sequences to None. 
    best_sequences=seq_matrices
    # set best_kappas to current values since they are the best so far.
    best_kappas = current_kappa.copy()
    
    # calculate overall asymmetry, which is necessary for the optimization
    overall_asymmetry = calc_global_sigma(ternary_seqs)
    # track without improvement
    without_improvement = 0
    # set window_size to 5, we leave it at 5 until there is no improvement and then change
    # to 6 or vice versa. This is because kappa is the mean value of windows 5 and 6.
    window_size=6
    
    # Pre-compute charged positions for all sequences (done once)
    all_charged_positions = precompute_charged_positions(ternary_seqs)
    
    # main loop
    for iteration in range(max_iterations):
        # Vectorized window analysis
        window_diffs, worst_windows, top_k_windows = vectorized_window_search(
            ternary_seqs, window_size, overall_asymmetry
        )
        
        # For each sequence, try to improve problematic windows
        improved = False
        
        # Process sequences in parallel-friendly batches
        for seq_idx in range(n_sequences):
            seq_improved = False
            current_window_sum = window_diffs[seq_idx].sum()
            
            # Try multiple problematic windows, not just the worst
            problematic_windows = top_k_windows[seq_idx]
            
            for window_idx in reversed(problematic_windows):  # Start with worst
                window_start = window_idx
                window_end = min(window_start + window_size, seq_length)
                
                # Find charged positions in this window using pre-computed data
                charged_in_seq = all_charged_positions[seq_idx]
                window_charged_pos = charged_in_seq[
                    (charged_in_seq >= window_start) & (charged_in_seq < window_end)
                ]
                
                if len(window_charged_pos) == 0:
                    continue
                
                best_improvement = 0
                best_move = None
                
                # Vectorized approach: process all moves for each from_pos
                for from_pos in window_charged_pos:
                    from_charge = ternary_seqs[seq_idx, from_pos]
                    
                    # Find all positions with different charges (vectorized)
                    different_charges = np.where(ternary_seqs[seq_idx] != from_charge)[0]
                    different_charges = different_charges[different_charges != from_pos]
                    
                    if len(different_charges) == 0:
                        continue
                    
                    # Smart sampling: prioritize positions outside current window
                    outside_window = different_charges[
                        (different_charges < window_start) | (different_charges >= window_end)
                    ]
                    inside_window = different_charges[
                        (different_charges >= window_start) & (different_charges < window_end)
                    ]
                    
                    # Prefer outside window moves, but include some inside window moves
                    if len(outside_window) > 0:
                        if len(outside_window) > 30:
                            candidates = np.random.choice(outside_window, 30, replace=False)
                        else:
                            candidates = outside_window
                        
                        # Add some inside window moves if space allows
                        remaining_slots = 50 - len(candidates)
                        if remaining_slots > 0 and len(inside_window) > 0:
                            n_inside = min(remaining_slots, len(inside_window))
                            inside_sample = np.random.choice(inside_window, n_inside, replace=False)
                            candidates = np.concatenate([candidates, inside_sample])
                    else:
                        # Only inside window moves available
                        if len(inside_window) > 50:
                            candidates = np.random.choice(inside_window, 50, replace=False)
                        else:
                            candidates = inside_window
                    
                    # Evaluate moves using fast window calculation
                    for to_pos in candidates:
                        # Try the swap
                        ternary_seqs[seq_idx, from_pos], ternary_seqs[seq_idx, to_pos] = \
                            ternary_seqs[seq_idx, to_pos], ternary_seqs[seq_idx, from_pos]
                        
                        # Fast calculation of improvement (only affected windows)
                        new_affected_sum = calc_affected_windows_fast(
                            ternary_seqs[seq_idx], from_pos, to_pos, window_size, 
                            overall_asymmetry[seq_idx]
                        )
                        
                        # Calculate original sum for affected windows
                        ternary_seqs[seq_idx, from_pos], ternary_seqs[seq_idx, to_pos] = \
                            ternary_seqs[seq_idx, to_pos], ternary_seqs[seq_idx, from_pos]  # Revert
                        
                        old_affected_sum = calc_affected_windows_fast(
                            ternary_seqs[seq_idx], from_pos, to_pos, window_size, 
                            overall_asymmetry[seq_idx]
                        )
                        
                        improvement = old_affected_sum - new_affected_sum
                        
                        if improvement > best_improvement:
                            best_improvement = improvement
                            best_move = (from_pos, to_pos)
                
                # Apply the best move if found and break (only one move per sequence per iteration)
                if best_move is not None:
                    from_pos, to_pos = best_move
                    ternary_seqs[seq_idx, from_pos], ternary_seqs[seq_idx, to_pos] = \
                        ternary_seqs[seq_idx, to_pos], ternary_seqs[seq_idx, from_pos]
                    seq_improved = True
                    break  # Move to next sequence
            
            if seq_improved:
                improved = True
                without_improvement = 0
        
        # if not improved, increment counter and change window size
        if not improved:
            without_improvement += 1
            if window_size == 5:
                window_size = 6
            else:
                window_size = 5
        if without_improvement >= 2:
            without_improvement = 0
            # calculate current kappa
            current_kappa = kappa(ternary_seqs, is_ternarized=True)
            # if best_sequences is None, save current sequences
            if best_sequences is None:
                best_sequences = ternary_seqs.copy()
                best_kappas = current_kappa.copy()
            else:
                # find which sequences have lower kappa than best_sequences
                for seq_idx in range(n_sequences):
                    if current_kappa[seq_idx] < best_kappas[seq_idx]:
                        best_sequences[seq_idx] = ternary_seqs[seq_idx].copy()
                        best_kappas[seq_idx] = current_kappa[seq_idx]
                
            # shuffle the sequence and keep trying. 
            for seq_idx in range(n_sequences):
                # Shuffle the sequence randomly while preserving charge counts
                # Get current sequence
                current_seq = ternary_seqs[seq_idx].copy()
                # Generate random permutation indices
                perm_indices = np.random.permutation(len(current_seq))
                # Apply permutation to preserve charge counts
                ternary_seqs[seq_idx] = current_seq[perm_indices]
            
            # Update charged positions after shuffling
            all_charged_positions = precompute_charged_positions(ternary_seqs)
        
        # Recalculate global asymmetry after changes
        overall_asymmetry = calc_global_sigma(ternary_seqs)

    return best_sequences

def calc_affected_windows_fast(ternary_seq, from_pos, to_pos, window_size, overall_asymmetry):
    """
    Fast calculation of window difference changes after a swap.
    Only recalculates affected windows instead of all windows.
    
    Parameters
    ----------
    ternary_seq : np.ndarray
        1D array representing a single ternary sequence
    from_pos : int
        Position to swap from
    to_pos : int
        Position to swap to
    window_size : int
        Size of the sliding window
    overall_asymmetry : float
        Global asymmetry value for this sequence
        
    Returns
    -------
    float
        Sum of squared differences for affected windows only
    """
    seq_length = len(ternary_seq)
    
    # Find windows affected by the swap
    affected_windows = set()
    for pos in [from_pos, to_pos]:
        start_window = max(0, pos - window_size + 1)
        end_window = min(pos + 1, seq_length - window_size + 1)
        affected_windows.update(range(start_window, end_window))
    
    if not affected_windows:
        return 0
    
    # Convert to sorted list for vectorized processing
    affected_windows = sorted(list(affected_windows))
    
    # Vectorized calculation for affected windows only
    window_starts = np.array(affected_windows)
    
    # Create windows using advanced indexing
    indices = window_starts[:, np.newaxis] + np.arange(window_size)
    windows = ternary_seq[indices]
    
    # Vectorized asymmetry calculation
    pos_counts = np.sum(windows == 1, axis=1)
    neg_counts = np.sum(windows == -1, axis=1)
    
    fa_local = pos_counts / window_size
    fb_local = neg_counts / window_size
    denominator = fa_local + fb_local
    
    window_asymmetry = np.zeros(len(affected_windows))
    valid_mask = denominator > 0
    window_asymmetry[valid_mask] = ((fa_local[valid_mask] - fb_local[valid_mask])**2) / denominator[valid_mask]
    
    # Calculate squared differences from overall asymmetry
    squared_diffs = (window_asymmetry - overall_asymmetry)**2
    return squared_diffs.sum()


def batch_evaluate_moves(ternary_seq, from_pos, to_positions, window_size, 
                         overall_asymmetry, current_window_sum):
    """
    Evaluate multiple moves from a single position in parallel.
    
    Parameters
    ----------
    ternary_seq : np.ndarray
        1D array representing a single ternary sequence
    from_pos : int
        Position to move from
    to_positions : np.ndarray
        Array of positions to evaluate moving to
    window_size : int
        Size of the sliding window
    overall_asymmetry : float
        Global asymmetry value for this sequence
    current_window_sum : float
        Current sum of window differences (baseline)
        
    Returns
    -------
    np.ndarray
        Array of improvements for each move
    """
    improvements = np.zeros(len(to_positions))
    from_charge = ternary_seq[from_pos]
    
    for i, to_pos in enumerate(to_positions):
        if ternary_seq[to_pos] != from_charge:
            # Temporarily perform swap
            ternary_seq[from_pos], ternary_seq[to_pos] = ternary_seq[to_pos], ternary_seq[from_pos]
            
            # Calculate new sum for affected windows only
            new_sum = calc_affected_windows_fast(ternary_seq, from_pos, to_pos, window_size, overall_asymmetry)
            
            # Calculate original sum for affected windows
            ternary_seq[from_pos], ternary_seq[to_pos] = ternary_seq[to_pos], ternary_seq[from_pos]  # Revert
            old_sum = calc_affected_windows_fast(ternary_seq, from_pos, to_pos, window_size, overall_asymmetry)
            
            improvements[i] = old_sum - new_sum
    
    return improvements

def precompute_charged_positions(ternary_seqs):
    """
    Pre-compute charged positions for all sequences to avoid repeated calculations.
    
    Parameters
    ----------
    ternary_seqs : np.ndarray
        2D array of ternary sequences
        
    Returns
    -------
    list
        List of arrays containing charged positions for each sequence
    """
    n_sequences = ternary_seqs.shape[0]
    charged_positions = []
    
    for seq_idx in range(n_sequences):
        # Find all charged positions (not neutral)
        charged_pos = np.where(ternary_seqs[seq_idx] != 0)[0]
        charged_positions.append(charged_pos)
    
    return charged_positions


def vectorized_window_search(ternary_seqs, window_size, overall_asymmetry,
                             num_top_windows=3):
    """
    Vectorized approach to find the top problematic windows for all sequences.
    
    Parameters
    ----------
    ternary_seqs : np.ndarray
        2D array of ternary sequences
    window_size : int
        Size of the sliding window
    overall_asymmetry : np.ndarray
        Global asymmetry values for all sequences
    num_top_windows : int, default=3
        Number of top problematic windows to return for each sequence
        
    Returns
    -------
    tuple
        (window_diffs, worst_windows, top_k_windows)
    """
    # Calculate window differences for all sequences
    window_diffs = calc_sigma_windows(ternary_seqs, window_size, overall_asymmetry)
    
    # Find worst windows
    worst_windows = np.argmax(window_diffs, axis=1)
    
    # Also find top-k problematic windows for more diverse optimization
    k = min(num_top_windows, window_diffs.shape[1])  # Top 3 or fewer if sequence is short
    top_k_windows = np.argsort(window_diffs, axis=1)[:, -k:]
    
    return window_diffs, worst_windows, top_k_windows

