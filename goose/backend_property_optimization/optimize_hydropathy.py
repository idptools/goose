'''
code specifically for optimizing the hydropathy value after sequences are generated
uses numpy vectorized operations to maximize efficiency.
'''

import numpy as np
from typing import List, Union
from goose import parameters
from goose.backend_property_calculation.calculate_properties_vectorized import calculate_hydropathy_batch, sequences_to_matrices, matrices_to_sequences


def optimize_hydropathy_vectorized(
    sequences: Union[List[str], np.ndarray],
    target_hydropathy: float, 
    preserve_charged: bool = True,
    max_iterations: int = 1000,
    tolerance: float = parameters.HYDRO_ERROR,
    batch_size: int = 10, 
    only_return_within_tolernace: bool = True,
    return_when_num_hit=None,
    exclude_residues=None,
    need_to_convert_seqs=False,
    convert_back_to_sequences=False) -> List[str]:
    """
    Hydropathy optimization function that modifies sequences to achieve a target hydropathy value.

    Parameters:
    -----------
    sequences : List[str] or np.ndarray
        Input sequences as a list of strings or a numpy array of strings.
    target_hydropathy : float
        Target mean hydropathy value to achieve.
    preserve_charged : bool
        If True, charged residues (D, E, K, R) will not be modified.
    max_iterations : int
        Maximum number of optimization iterations per sequence.
    tolerance : float
        Acceptable difference between achieved and target hydropathy.
    batch_size : int
        Number of mutations to apply in each iteration.
    only_return_within_tolernace : bool
        If True, only return sequences that are within the specified tolerance.
    return_when_num_hit : int or None
        If specified, return when this many sequences are within tolerance.
    exclude_residues : List[str] or None
        Residues to exclude from modification.
    need_to_convert_seqs : bool
        If True, convert sequences to matrices before processing.
    convert_back_to_sequences : bool
        If True, convert the final matrices back to sequences.

    Returns:
    --------
    List[str]
        Optimized sequences with hydropathy values close to the target.
    """
    # Hydropathy scale constants
    HYDROPATHY_SCALE = np.array([6.3, 7.0, 1.0, 1.0, 7.3, 4.1, 1.3, 9.0, 
                            0.6, 8.3, 6.4, 1.0, 2.9, 1.0, 0.0, 3.7, 3.8, 8.7, 3.6, 3.2])

    if need_to_convert_seqs:    # Convert input to proper format
        if isinstance(sequences, str):
            sequences = [sequences]
        # Convert sequences to matrices
        seq_matrices = sequences_to_matrices(sequences)
    else:
        seq_matrices = sequences

    # see if we need to return a specific number. 
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
    
    # Process exclude_residues parameter
    excluded_indices = set()
    if exclude_residues is not None:
        excluded_indices = set(aa_to_int[aa] for aa in exclude_residues if aa in aa_to_int)
    
    # get num and length
    n_sequences, seq_length = seq_matrices.shape
    
    # Define charged residue indices
    charged_indices = {2, 3, 8, 14} if preserve_charged else set()  # D, E, K, R
    
    # Precompute replacement options for each amino acid (using sets for faster lookups)
    forbidden_indices = charged_indices.union(excluded_indices)
    
    better_options = {}  # aa -> (increase_options, decrease_options)
    for aa in range(20):
        current_hydro = HYDROPATHY_SCALE[aa]
        
        # Options that increase hydropathy
        increase_opts = [i for i in range(20) 
                        if HYDROPATHY_SCALE[i] > current_hydro and i not in forbidden_indices]
        
        # Options that decrease hydropathy  
        decrease_opts = [i for i in range(20)
                        if HYDROPATHY_SCALE[i] < current_hydro and i not in forbidden_indices]
        
        better_options[aa] = (np.array(increase_opts), np.array(decrease_opts))
    
    # Calculate initial hydropathy and track which sequences need work
    current_hydro = calculate_hydropathy_batch(seq_matrices)
    completed_sequences = set()
    
    # Precompute mutable positions for each sequence  
    mutable_positions_cache = {}
    for seq_idx in range(n_sequences):
        if preserve_charged:
            mutable_pos = np.where(~np.isin(seq_matrices[seq_idx], list(charged_indices)))[0]
        else:
            mutable_pos = np.arange(seq_length)
        mutable_positions_cache[seq_idx] = mutable_pos
    
    # Main optimization loop with smart early termination
    for iteration in range(max_iterations):
        # Check convergence
        delta = current_hydro - target_hydropathy
        within_tolerance_mask = np.abs(delta) <= tolerance
        
        # Update completed sequences
        newly_completed = set(np.where(within_tolerance_mask)[0]) - completed_sequences
        completed_sequences.update(newly_completed)
        
        # Early termination
        if len(completed_sequences) >= return_when_num_hit:
            break
            
        # Get sequences that still need work
        active_sequences = [i for i in range(n_sequences) if i not in completed_sequences]
        if not active_sequences:
            break
            
        # Process only active sequences in batches
        sequences_to_recalc = set()
        
        for seq_idx in active_sequences:
            current_delta = delta[seq_idx]
            
            # Skip if already within tolerance (shouldn't happen but safety check)
            if abs(current_delta) <= tolerance:
                continue
                
            mutable_pos = mutable_positions_cache[seq_idx]
            if len(mutable_pos) == 0:
                continue
                
            # Determine direction and select positions
            need_increase = current_delta < -tolerance
            num_to_modify = min(batch_size, len(mutable_pos))
            positions_to_modify = np.random.choice(mutable_pos, size=num_to_modify, replace=False)
            
            # Make modifications
            for pos in positions_to_modify:
                current_aa = seq_matrices[seq_idx, pos]
                increase_opts, decrease_opts = better_options[current_aa]
                
                candidates = increase_opts if need_increase else decrease_opts
                
                if len(candidates) > 0:
                    replacement = np.random.choice(candidates)
                    seq_matrices[seq_idx, pos] = replacement
                    sequences_to_recalc.add(seq_idx)
        
        # Recalculate hydropathy only for modified sequences
        if sequences_to_recalc:
            modified_indices = np.array(list(sequences_to_recalc))
            current_hydro[modified_indices] = calculate_hydropathy_batch(seq_matrices[modified_indices])
    
    # Final filtering if needed
    if only_return_within_tolernace:
        final_hydro = calculate_hydropathy_batch(seq_matrices)
        within_tolerance = np.abs(final_hydro - target_hydropathy) <= tolerance
        seq_matrices = seq_matrices[within_tolerance]
    
    # if we have more than one sequence, sort them so the first sequence is the one that is closest to the target hydropathy
    if len(seq_matrices) > 1:
        final_hydro = calculate_hydropathy_batch(seq_matrices)
        sorted_indices = np.argsort(np.abs(final_hydro - target_hydropathy))
        seq_matrices = seq_matrices[sorted_indices]

    # Convert back to sequences
    if convert_back_to_sequences:
        optimized_sequences = matrices_to_sequences(seq_matrices)
    else:
        optimized_sequences = seq_matrices
    
    return optimized_sequences

