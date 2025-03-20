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
    tolerance: float = 0.05) -> List[str]:
    """
    Efficiently optimize protein sequences to match a target mean hydropathy value
    using vectorized operations.
    
    Parameters:
    -----------
    sequences : list of str or numpy.ndarray
        Input sequences to optimize (can be a single sequence or multiple)
    target_hydropathy : float
        Target mean hydropathy value to achieve
    preserve_charged : bool, optional (default=True)
        If True, charged residues (D, E, K, R) will not be modified
    max_iterations : int, optional (default=100)
        Maximum number of optimization iterations per sequence
    tolerance : float, optional (default=0.05)
        Acceptable difference between achieved and target hydropathy
        
    Returns:
    --------
    list
        Optimized sequences with hydropathy values close to target
    """
    # Convert single string to list if needed
    if isinstance(sequences, str):
        sequences = [sequences]
    
    # Convert to numpy array if not already
    if not isinstance(sequences, np.ndarray):
        sequences = np.array(sequences)
    
    # Create amino acid mappings
    aa_to_int = {'A': 0, 'C': 1, 'D': 2, 'E': 3, 'F': 4, 'G': 5, 'H': 6, 'I': 7, 
                'K': 8, 'L': 9, 'M': 10, 'N': 11, 'P': 12, 'Q': 13, 'R': 14, 
                'S': 15, 'T': 16, 'V': 17, 'W': 18, 'Y': 19}
    int_to_aa = {v: k for k, v in aa_to_int.items()}
    
    # Convert all sequences to integer matrices at once
    seq_matrices = []
    for seq in sequences:
        seq_matrices.append(np.array([aa_to_int[aa] for aa in seq]))
    
    # Identify charged and uncharged residues
    charged_indices = np.array([2, 3, 8, 14])  # D, E, K, R
    uncharged_indices = np.array([i for i in range(20) if i not in charged_indices])
    
    # Sort amino acids by hydropathy
    sorted_aa_by_hydro = np.argsort(HYDROPATHY_SCALE)
    
    # Create maps for amino acid substitution based on hydropathy
    # For each amino acid, find replacements that increase or decrease hydropathy
    increase_hydro_map = {}
    decrease_hydro_map = {}
    
    for aa_idx in range(20):
        current_hydro = HYDROPATHY_SCALE[aa_idx]
        
        # Get replacements that would increase hydropathy
        higher_hydro = sorted_aa_by_hydro[HYDROPATHY_SCALE[sorted_aa_by_hydro] > current_hydro]
        if len(higher_hydro) == 0:
            higher_hydro = sorted_aa_by_hydro[HYDROPATHY_SCALE[sorted_aa_by_hydro] >= current_hydro]
        increase_hydro_map[aa_idx] = higher_hydro
        
        # Get replacements that would decrease hydropathy
        lower_hydro = sorted_aa_by_hydro[HYDROPATHY_SCALE[sorted_aa_by_hydro] < current_hydro]
        if len(lower_hydro) == 0:
            lower_hydro = sorted_aa_by_hydro[HYDROPATHY_SCALE[sorted_aa_by_hydro] <= current_hydro]
        decrease_hydro_map[aa_idx] = lower_hydro
    
    # Process sequences in parallel batches
    for iteration in range(max_iterations):
        # Calculate hydropathy for all sequences simultaneously
        hydropathy_scores = calculate_hydropathy_batch(seq_matrices)
        
        # Check which sequences are already within tolerance
        delta = hydropathy_scores - target_hydropathy
        within_tolerance = np.abs(delta) <= tolerance
        
        # If all sequences are within tolerance, we're done
        if np.all(within_tolerance):
            break
        
        # Process each sequence that needs adjustment
        for i, (seq_array, hydro_score) in enumerate(zip(seq_matrices, hydropathy_scores)):
            if within_tolerance[i]:
                continue  # Skip sequences that are already good
            
            # Determine if we need to increase or decrease hydropathy
            increase_needed = hydro_score < target_hydropathy
            
            # Find positions that can be modified
            if preserve_charged:
                mutable_mask = ~np.isin(seq_array, charged_indices)
            else:
                mutable_mask = np.ones_like(seq_array, dtype=bool)
                
            # Get modifiable positions
            mutable_positions = np.where(mutable_mask)[0]
            if len(mutable_positions) == 0:
                continue  # Can't modify this sequence
            
            # Choose a random position to modify
            mod_position = np.random.choice(mutable_positions)
            current_aa = seq_array[mod_position]
            
            # Choose appropriate replacement based on whether we need to increase or decrease
            if increase_needed:
                candidate_aas = increase_hydro_map[current_aa]
            else:
                candidate_aas = decrease_hydro_map[current_aa]
                
            # Filter out charged residues if preserving charge
            if preserve_charged:
                candidate_aas = candidate_aas[~np.isin(candidate_aas, charged_indices)]
                
            if len(candidate_aas) > 0:
                # Replace the amino acid
                replacement = np.random.choice(candidate_aas)
                seq_matrices[i][mod_position] = replacement
    
    # Convert integer arrays back to sequences
    optimized_sequences = [''.join([int_to_aa[idx] for idx in seq]) for seq in seq_matrices]
    
    return optimized_sequences


def optimize_hydropathy_fully_vectorized(
    sequences: Union[List[str], np.ndarray],
    target_hydropathy: float, 
    preserve_charged: bool = True,
    max_iterations: int = 100,
    tolerance: float = 0.05,
    batch_size: int = 10) -> List[str]:
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
    
    # Create amino acid mappings
    aa_to_int = {'A': 0, 'C': 1, 'D': 2, 'E': 3, 'F': 4, 'G': 5, 'H': 6, 'I': 7, 
                'K': 8, 'L': 9, 'M': 10, 'N': 11, 'P': 12, 'Q': 13, 'R': 14, 
                'S': 15, 'T': 16, 'V': 17, 'W': 18, 'Y': 19}
    int_to_aa = {v: k for k, v in aa_to_int.items()}
    
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
    for iteration in range(max_iterations):
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
        
        # If all sequences are within tolerance, we're done
        if np.all(within_tolerance):
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
                else:
                    # Need lower hydropathy
                    candidates = np.where(~hydro_lookup[current_aa])[0]
                    if preserve_charged:
                        candidates = candidates[~np.isin(candidates, charged_indices)]
                    
                # Remove self from candidates
                candidates = candidates[candidates != current_aa]
                
                # Make replacement if candidates exist
                if len(candidates) > 0:
                    replacement = np.random.choice(candidates)
                    if same_length:
                        seq_matrices[i, pos] = replacement
                    else:
                        seq_matrices[i][pos] = replacement
    
    # Convert optimized sequences back to strings
    if same_length and isinstance(seq_matrices, np.ndarray):
        optimized_sequences = [''.join([int_to_aa[idx] for idx in seq]) for seq in seq_matrices]
    else:
        optimized_sequences = [''.join([int_to_aa[idx] for idx in seq]) for seq in seq_matrices]
    
    return optimized_sequences


