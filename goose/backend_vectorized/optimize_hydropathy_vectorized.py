'''
code specifically for optimizing the hydropathy value after sequences are generated
uses numpy vectorized operations to maximize efficiency.
'''

import numpy as np
from typing import List, Union

# Constants
HYDROPATHY_SCALE = np.array([6.3, 7.0, 1.0, 1.0, 7.3, 4.1, 1.3, 9.0, 
                            0.6, 8.3, 6.4, 1.0, 2.9, 1.0, 0.0, 3.7, 3.8, 8.7, 3.6, 3.2])


def calculate_hydropathy_batch(seq_matrices):
    """
    Calculate hydropathy scores for multiple sequences simultaneously.
    
    Parameters:
    -----------
    seq_matrices : list of numpy.ndarray
        List of sequences represented as integer matrices
        
    Returns:
    --------
    numpy.ndarray
        Array of hydropathy scores for each sequence
    """
    # Calculate mean hydropathy for each sequence using vectorized operations
    hydropathy_scores = np.array([np.mean(HYDROPATHY_SCALE[seq]) for seq in seq_matrices])
    return hydropathy_scores


def optimize_hydropathy_vectorized(
    sequences: Union[List[str], np.ndarray],
    target_hydropathy: float, 
    preserve_charged: bool = True,
    max_iterations: int = 100,
    tolerance: float = 0.05,
    batch_size: int = 10, 
    only_return_within_tolernace: bool = False,
    return_when_num_hit=None,
    exclude_residues=None) -> List[str]:
    """
    Fully vectorized approach to optimize protein sequences to match a target hydropathy.
    Processes all sequences in parallel and modifies them in batches.
    
    Parameters:
    -----------
    sequences : list of str or numpy.ndarray
        Input sequences to optimize
    target_hydropathy : float
        Target mean hydropathy value to achieve
    preserve_charged : bool
        If True, charged residues (D, E, K, R) will not be modified
    max_iterations : int
        Maximum number of optimization iterations per sequence
    tolerance : float
        Acceptable difference between achieved and target hydropathy
    batch_size : int
        Number of positions to modify in each sequence per iteration
    only_return_within_tolernace : bool
        If True, only return sequences that are within the tolerance
        default is False
    return_when_num_hit : int
        If specified, return when this number of sequences are within tolerance
    exclude_residues : list
        List of residues to exclude from optimization
        
    Returns:
    --------
    list
        Optimized sequences with hydropathy values close to target
    """
    # Convert input to proper format
    if isinstance(sequences, str):
        sequences = [sequences]
    if not isinstance(sequences, np.ndarray):
        sequences = np.array(sequences)

    if return_when_num_hit is not None:
        only_return_within_tolernace = True
        if return_when_num_hit > len(sequences):
            return_when_num_hit = len(sequences)        
    else:
        return_when_num_hit = len(sequences)
    
    # Create amino acid mappings
    aa_to_int = {'A': 0, 'C': 1, 'D': 2, 'E': 3, 'F': 4, 'G': 5, 'H': 6, 'I': 7, 
                'K': 8, 'L': 9, 'M': 10, 'N': 11, 'P': 12, 'Q': 13, 'R': 14, 
                'S': 15, 'T': 16, 'V': 17, 'W': 18, 'Y': 19}
    int_to_aa = {v: k for k, v in aa_to_int.items()}
    
    # Process exclude_residues parameter
    excluded_indices = np.array([], dtype=int)
    if exclude_residues is not None:
        excluded_indices = np.array([aa_to_int[aa] for aa in exclude_residues if aa in aa_to_int])
    
    # Check if all sequences have the same length for optimal vectorization
    seq_lens = np.array([len(seq) for seq in sequences])
    same_length = np.all(seq_lens == seq_lens[0])
    
    # Convert sequences to integer matrices
    seq_matrices = []
    for seq in sequences:
        seq_matrices.append(np.array([aa_to_int[aa] for aa in seq]))
    
    # Convert to 2D matrix if all sequences have same length
    if same_length:
        seq_matrices = np.stack(seq_matrices)
    
    # Define charged and hydropathy parameters
    charged_indices = np.array([2, 3, 8, 14])  # D, E, K, R
    
    # Create lookup arrays for amino acid substitutions
    # Precompute all possible substitutions to increase/decrease hydropathy
    hydro_lookup = np.zeros((20, 20), dtype=bool)  # [current_aa, replacement_aa] -> increases hydropathy?
    for i in range(20):
        for j in range(20):
            hydro_lookup[i, j] = HYDROPATHY_SCALE[j] > HYDROPATHY_SCALE[i]
    
    # Main optimization loop
    for _ in range(max_iterations):
        # Calculate current hydropathy for all sequences
        if same_length:
            # For same-length sequences, we can use advanced vectorization
            current_hydro = np.mean(HYDROPATHY_SCALE[seq_matrices], axis=1)
        else:
            # For variable-length sequences
            current_hydro = np.array([np.mean(HYDROPATHY_SCALE[seq]) for seq in seq_matrices])
        
        # Check which sequences are within tolerance
        delta = current_hydro - target_hydropathy
        within_tolerance = np.abs(delta) <= tolerance

        # if return_when_num_hit is specified, check if we have enough sequences within tolerance
        if np.sum(within_tolerance) >= return_when_num_hit:
            break
        
        # Get sequences that need modification
        need_increase = delta < -tolerance
        need_decrease = delta > tolerance
        
        # Process each sequence that needs adjustment
        for i in range(len(seq_matrices)):
            if within_tolerance[i]:
                continue
                
            current_seq = seq_matrices[i] if not same_length else seq_matrices[i].copy()
            
            # Find mutable positions
            if preserve_charged:
                mutable_mask = ~np.isin(current_seq, charged_indices)
            else:
                mutable_mask = np.ones_like(current_seq, dtype=bool)
                
            mutable_positions = np.where(mutable_mask)[0]
            if len(mutable_positions) == 0:
                continue
                
            # Select batch_size positions to modify (or fewer if not enough)
            num_to_modify = min(batch_size, len(mutable_positions))
            positions_to_modify = np.random.choice(mutable_positions, size=num_to_modify, replace=False)
            
            for pos in positions_to_modify:
                current_aa = current_seq[pos]
                
                # Find appropriate replacements
                if need_increase[i]:
                    # Need higher hydropathy
                    candidates = np.where(hydro_lookup[current_aa])[0]
                    if preserve_charged:
                        candidates = candidates[~np.isin(candidates, charged_indices)]
                    # Exclude specified residues
                    if len(excluded_indices) > 0:
                        candidates = candidates[~np.isin(candidates, excluded_indices)]
                else:
                    # Need lower hydropathy
                    candidates = np.where(~hydro_lookup[current_aa])[0]
                    if preserve_charged:
                        candidates = candidates[~np.isin(candidates, charged_indices)]
                    # Exclude specified residues
                    if len(excluded_indices) > 0:
                        candidates = candidates[~np.isin(candidates, excluded_indices)]
                    
                # Remove self from candidates
                candidates = candidates[candidates != current_aa]
                
                # Make replacement if candidates exist
                if len(candidates) > 0:
                    replacement = np.random.choice(candidates)
                    if same_length:
                        seq_matrices[i, pos] = replacement
                    else:
                        seq_matrices[i][pos] = replacement

    # if we only want to return sequences within tolerance, filter them
    if only_return_within_tolernace:
        # calulate hydropathy for all sequences using matrices
        hydropathy_scores = calculate_hydropathy_batch(seq_matrices)
        within_tolerance = np.abs(hydropathy_scores - target_hydropathy) <= tolerance
        seq_matrices = seq_matrices[within_tolerance]

    
    # Convert optimized sequences back to strings
    if same_length and isinstance(seq_matrices, np.ndarray):
        optimized_sequences = [''.join([int_to_aa[idx] for idx in seq]) for seq in seq_matrices]
    else:
        optimized_sequences = [''.join([int_to_aa[idx] for idx in seq]) for seq in seq_matrices]
    
    return optimized_sequences

