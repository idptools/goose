'''
Streamlined kappa calculation module for charged residues only.

This module is optimized for:
- Charged residue kappa calculation only (K, R as positive; D, E as negative)
- np.ndarray input (from sequences_to_matrices)
- Minimal redundancy and maximum performance
- Dual-window (5 & 6) calculations as the standard
'''

import numpy as np


def charge_matrix_to_ternary(seq_matrix):
    """
    Convert sequence matrix to ternary array for charge-based kappa calculation.
    
    This function is optimized specifically for charged residues:
    - Positive charges: K(8), R(14) → +1
    - Negative charges: D(2), E(3) → -1
    - All others → 0
    
    Parameters
    ----------
    seq_matrix : np.ndarray
        2D array of shape (n_sequences, sequence_length) with integer amino acid codes
        
    Returns
    -------
    np.ndarray
        2D array of shape (n_sequences, sequence_length) with ternary charge values
    """
    # Use direct boolean indexing - faster than np.isin for small sets
    ternary_arrays = np.zeros_like(seq_matrix, dtype=np.int8)
    
    # Positive charges: K(8), R(14) - direct comparisons are faster than np.isin
    pos_mask = (seq_matrix == 8) | (seq_matrix == 14)
    ternary_arrays[pos_mask] = 1
    
    # Negative charges: D(2), E(3)
    neg_mask = (seq_matrix == 2) | (seq_matrix == 3)
    ternary_arrays[neg_mask] = -1
    
    return ternary_arrays


def analytical_deltamax_vectorized(ternary_arrays, window_size, overall_asymmetry):
    """
    Fully vectorized analytical deltamax computation - eliminates all loops.
    
    This ultra-optimized version uses advanced NumPy broadcasting to compute
    all window positions simultaneously, providing maximum performance.
    """
    n_sequences, seq_length = ternary_arrays.shape
    
    # Count positive and negative charges - vectorized
    count_pos = np.sum(ternary_arrays == 1, axis=1)
    count_neg = np.sum(ternary_arrays == -1, axis=1)
    
    # Pre-compute constants
    winged_length = seq_length + 2 * window_size
    n_windows = winged_length - window_size + 1
    inv_window_size = 1.0 / window_size
    
    # Maximally segregated arrangement positions - these are arrays now
    pos_start = np.full(n_sequences, window_size)  # All sequences have same pos_start
    pos_end = window_size + count_pos  # Shape: (n_sequences,)
    neg_start = window_size + seq_length - count_neg
    neg_end = np.full(n_sequences, window_size + seq_length)  # All sequences have same neg_end
    
    # Create all window positions at once - shape: (n_windows,)
    window_starts = np.arange(n_windows)
    window_ends = window_starts + window_size
    
    # Broadcast for all sequences and windows - shapes: (n_sequences, n_windows)
    pos_start_bc = pos_start[:, np.newaxis]  # (n_sequences, 1)
    pos_end_bc = pos_end[:, np.newaxis]      # (n_sequences, 1)
    neg_start_bc = neg_start[:, np.newaxis]  # (n_sequences, 1)
    neg_end_bc = neg_end[:, np.newaxis]      # (n_sequences, 1)
    
    window_starts_bc = window_starts[np.newaxis, :]  # (1, n_windows)
    window_ends_bc = window_ends[np.newaxis, :]      # (1, n_windows)
    
    # Vectorized window charge calculations for ALL windows and sequences
    window_pos = np.maximum(0, 
                           np.minimum(pos_end_bc, window_ends_bc) - 
                           np.maximum(pos_start_bc, window_starts_bc))
    
    window_neg = np.maximum(0,
                           np.minimum(neg_end_bc, window_ends_bc) - 
                           np.maximum(neg_start_bc, window_starts_bc))
    
    # Calculate window asymmetries for ALL windows - shape: (n_sequences, n_windows)
    total_charges = window_pos + window_neg
    valid_windows = total_charges > 0
    
    # Initialize window asymmetries
    window_asymmetry = np.zeros((n_sequences, n_windows))
    
    if np.any(valid_windows):
        fa_local = window_pos * inv_window_size
        fb_local = window_neg * inv_window_size
        fa_plus_fb = fa_local + fb_local
        fa_minus_fb = fa_local - fb_local
        
        # Only compute for valid windows
        window_asymmetry[valid_windows] = (fa_minus_fb[valid_windows] ** 2) / fa_plus_fb[valid_windows]
    
    # Broadcast overall asymmetry and compute squared differences
    overall_asymmetry_bc = overall_asymmetry[:, np.newaxis]  # (n_sequences, 1)
    squared_diffs = (window_asymmetry - overall_asymmetry_bc) ** 2
    
    # Sum and average over all windows
    return np.sum(squared_diffs, axis=1) / n_windows


def analytical_deltamax(ternary_arrays, window_size, overall_asymmetry):
    """
    Compute the analytical deltamax for charged residue kappa calculation.
    
    This is the optimized analytical solution that avoids building sequences
    and provides perfect accuracy with significant speed and memory improvements.
    
    Parameters
    ----------
    ternary_arrays : np.ndarray
        2D array of shape (n_sequences, sequence_length) with ternary charge values
    window_size : int
        Size of the sliding window
    overall_asymmetry : np.ndarray
        Array of global asymmetry values for each sequence
        
    Returns
    -------
    np.ndarray
        Array of deltamax values for each sequence
    """
    n_sequences, seq_length = ternary_arrays.shape
    
    # Use vectorized version for larger datasets, loop version for smaller ones
    # The crossover point depends on memory vs compute trade-offs
    if n_sequences > 50 or seq_length > 200:
        # Large datasets: use vectorized version (more memory but faster)
        return analytical_deltamax_vectorized(ternary_arrays, window_size, overall_asymmetry)
    
    # Small datasets: use loop version to save memory
    # Count positive and negative charges - vectorized
    count_pos = np.sum(ternary_arrays == 1, axis=1)
    count_neg = np.sum(ternary_arrays == -1, axis=1)
    
    # Pre-compute constants
    winged_length = seq_length + 2 * window_size
    n_windows = winged_length - window_size + 1
    inv_n_windows = 1.0 / n_windows  # Pre-compute division
    inv_window_size = 1.0 / window_size  # Pre-compute division
    
    # Maximally segregated arrangement: wing_zeros | pos | neutral | neg | wing_zeros
    pos_start = window_size
    pos_end = window_size + count_pos
    neg_start = window_size + seq_length - count_neg
    neg_end = window_size + seq_length
    
    # Pre-allocate arrays
    delta_sum = np.zeros(n_sequences)
    window_asymmetry = np.zeros(n_sequences)
    
    # Vectorized window calculation - unroll small loops for better performance
    window_positions = np.arange(n_windows)
    
    for i in window_positions:
        window_start = i
        window_end = i + window_size
        
        # Calculate charges in each window for all sequences - vectorized
        window_pos = np.maximum(0, 
                               np.minimum(pos_end, window_end) - 
                               np.maximum(pos_start, window_start))
        window_neg = np.maximum(0,
                               np.minimum(neg_end, window_end) - 
                               np.maximum(neg_start, window_start))
        
        # Calculate window asymmetries - optimized
        total_charges = window_pos + window_neg
        valid_windows = total_charges > 0
        
        # Reset window asymmetry array (faster than creating new)
        window_asymmetry.fill(0.0)
        
        if np.any(valid_windows):
            fa_local = window_pos[valid_windows] * inv_window_size
            fb_local = window_neg[valid_windows] * inv_window_size
            fa_plus_fb = fa_local + fb_local
            fa_minus_fb = fa_local - fb_local
            window_asymmetry[valid_windows] = (fa_minus_fb * fa_minus_fb) / fa_plus_fb
        
        # Accumulate squared differences - vectorized
        diff = window_asymmetry - overall_asymmetry
        delta_sum += diff * diff
    
    return delta_sum * inv_n_windows


def patterning_asymmetry(ternary_winged, window_size, overall_asymmetry):
    """
    Compute patterning asymmetry using fully vectorized sliding windows.
    
    Parameters
    ----------
    ternary_winged : np.ndarray
        2D array with wings added for proper window calculation
    window_size : int
        Size of the sliding window
    overall_asymmetry : np.ndarray
        Array of global asymmetry values for each sequence
        
    Returns
    -------
    np.ndarray
        Array of patterning asymmetry values for each sequence
    """
    n_sequences, winged_length = ternary_winged.shape
    n_windows = winged_length - window_size + 1
    
    # Create all windows at once using vectorized sliding window
    windows_3d = np.lib.stride_tricks.sliding_window_view(
        ternary_winged, window_shape=window_size, axis=1
    )
    
    # Count charges across all windows
    window_count_pos = np.sum(windows_3d == 1, axis=2)
    window_count_neg = np.sum(windows_3d == -1, axis=2)
    
    # Calculate local fractions
    fa_local = window_count_pos / window_size
    fb_local = window_count_neg / window_size
    
    # Compute window asymmetries
    denominator = fa_local + fb_local
    window_asymmetry = np.zeros((n_sequences, n_windows))
    valid_mask = denominator > 0
    window_asymmetry[valid_mask] = ((fa_local[valid_mask] - fb_local[valid_mask])**2) / denominator[valid_mask]
    
    # Calculate squared differences from overall asymmetry
    overall_asymmetry_broadcast = overall_asymmetry[:, np.newaxis]
    squared_diffs = (window_asymmetry - overall_asymmetry_broadcast)**2
    
    # Return average patterning asymmetry
    return np.sum(squared_diffs, axis=1) / n_windows


def global_asymmetry(ternary_winged):
    """
    Compute global compositional asymmetry for multiple sequences.
    
    Parameters
    ----------
    ternary_winged : np.ndarray
        2D array with ternary charge values including wings
        
    Returns
    -------
    np.ndarray
        Array of global asymmetry values for each sequence
    """
    n_sequences, winged_length = ternary_winged.shape
    
    # Count positive and negative charges
    count_pos = np.sum(ternary_winged == 1, axis=1)
    count_neg = np.sum(ternary_winged == -1, axis=1)
    
    # Calculate fractions
    fa = count_pos / winged_length
    fb = count_neg / winged_length
    
    # Compute asymmetry
    denominator = fa + fb
    asymmetry = np.zeros(n_sequences)
    valid_mask = denominator > 0
    asymmetry[valid_mask] = ((fa[valid_mask] - fb[valid_mask])**2) / denominator[valid_mask]
    
    return asymmetry


def kappa(seq_matrix, flatten=False, is_ternarized=False):
    """
    Calculate kappa for charged residues using dual windows (5 and 6) for maximum performance.
    
    This is the primary function for kappa calculation, optimized for the standard
    approach of averaging window sizes 5 and 6.
    
    Parameters
    ----------
    seq_matrix : np.ndarray
        2D array of shape (n_sequences, sequence_length) with integer amino acid codes
    flatten : bool, optional
        If True, kappa values above 1 will be set to 1
    is_ternarized : bool, optional
        If True, seq_matrix is already in ternary format (1 for positive, -1 for negative, 0 for neutral)
        
    Returns
    -------
    np.ndarray
        Array of kappa values (average of window sizes 5 and 6) for each sequence
    """
    n_sequences, seq_length = seq_matrix.shape
    
    # Input validation
    if seq_length < 6:  # Need at least 6 for window_size=6
        return np.full(n_sequences, -1.0)
    
    # Convert to ternary charge array ONCE
    if not is_ternarized:
        ternary_arrays = charge_matrix_to_ternary(seq_matrix)
    else:
        ternary_arrays = seq_matrix.astype(np.int8)
    
    # Check for valid sequences ONCE with faster operations
    pos_counts = np.sum(ternary_arrays == 1, axis=1)
    neg_counts = np.sum(ternary_arrays == -1, axis=1)
    valid_sequences = (pos_counts > 0) & (neg_counts > 0)
    
    if not np.any(valid_sequences):
        return np.full(n_sequences, -1.0)
    
    # Pre-allocate result array
    kappa_avg = np.full(n_sequences, -1.0)
    
    # Only process valid sequences to save computation
    valid_indices = np.where(valid_sequences)[0]
    if len(valid_indices) == 0:
        return kappa_avg
    
    # Extract only valid sequences for processing
    valid_ternary = ternary_arrays[valid_sequences]
    n_valid = valid_ternary.shape[0]
    
    # OPTIMIZATION: Create both winged arrays in a single allocation
    # Use the larger wing size (6) and create views for the smaller one
    max_wing = 6
    winged_length = seq_length + 2 * max_wing  # 2 * 6 = 12
    ternary_winged_6 = np.zeros((n_valid, winged_length), dtype=np.int8)
    ternary_winged_6[:, max_wing:max_wing+seq_length] = valid_ternary
    
    # Create view for window size 5 (trim 1 from each side)
    ternary_winged_5 = ternary_winged_6[:, 1:-1]
    
    # Calculate asymmetries for both window sizes - only for valid sequences
    overall_asymmetry_5 = global_asymmetry(ternary_winged_5)
    overall_asymmetry_6 = global_asymmetry(ternary_winged_6)
    
    # OPTIMIZATION: Calculate deltas and deltamaxes in sequence to maintain cache locality
    delta_5 = patterning_asymmetry(ternary_winged_5, 5, overall_asymmetry_5)
    deltamax_5 = analytical_deltamax(valid_ternary, 5, overall_asymmetry_5)
    
    delta_6 = patterning_asymmetry(ternary_winged_6, 6, overall_asymmetry_6)
    deltamax_6 = analytical_deltamax(valid_ternary, 6, overall_asymmetry_6)
    
    # Calculate kappa values for both window sizes - vectorized operations
    valid_5 = deltamax_5 > 0
    valid_6 = deltamax_6 > 0
    both_valid = valid_5 & valid_6
    
    # Only compute where both window sizes are valid
    if np.any(both_valid):
        # OPTIMIZATION: Pre-compute reciprocals to avoid division in tight loop
        inv_deltamax_5 = 1.0 / deltamax_5[both_valid]
        inv_deltamax_6 = 1.0 / deltamax_6[both_valid]
        
        kappa_5_valid = delta_5[both_valid] * inv_deltamax_5
        kappa_6_valid = delta_6[both_valid] * inv_deltamax_6
        kappa_avg_valid = (kappa_5_valid + kappa_6_valid) * 0.5
        
        if flatten:
            kappa_avg_valid = np.minimum(kappa_avg_valid, 1.0)
        
        # Map back to original indices
        kappa_avg[valid_indices[both_valid]] = kappa_avg_valid
    
    return kappa_avg


