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
    
    # determine better options for each amino acid. 
    if preserve_charged==False and len(excluded_indices) == 0:
        better_options = {0: (np.array([ 1, 4, 7, 9, 10, 17]), np.array([ 2, 3, 5, 6, 8, 11, 12, 13, 14, 15, 16, 18, 19])), 1: (np.array([ 4, 7, 9, 17]), np.array([ 0, 2, 3, 5, 6, 8, 10, 11, 12, 13, 14, 15, 16, 18, 19])), 2: (np.array([ 0, 1, 4, 5, 6, 7, 9, 10, 12, 15, 16, 17, 18, 19]), np.array([ 8, 14])), 3: (np.array([ 0, 1, 4, 5, 6, 7, 9, 10, 12, 15, 16, 17, 18, 19]), np.array([ 8, 14])), 4: (np.array([ 7, 9, 17]), np.array([ 0, 1, 2, 3, 5, 6, 8, 10, 11, 12, 13, 14, 15, 16, 18, 19])), 5: (np.array([ 0, 1, 4, 7, 9, 10, 17]), np.array([ 2, 3, 6, 8, 11, 12, 13, 14, 15, 16, 18, 19])), 6: (np.array([ 0, 1, 4, 5, 7, 9, 10, 12, 15, 16, 17, 18, 19]), np.array([ 2, 3, 8, 11, 13, 14])), 7: (np.array([], dtype=np.int8), np.array([ 0, 1, 2, 3, 4, 5, 6, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19])), 8: (np.array([ 0, 1, 2, 3, 4, 5, 6, 7, 9, 10, 11, 12, 13, 15, 16, 17, 18, 19]), np.array([14])), 9: (np.array([ 7, 17]), np.array([ 0, 1, 2, 3, 4, 5, 6, 8, 10, 11, 12, 13, 14, 15, 16, 18, 19])), 10: (np.array([ 1, 4, 7, 9, 17]), np.array([ 0, 2, 3, 5, 6, 8, 11, 12, 13, 14, 15, 16, 18, 19])), 11: (np.array([ 0, 1, 4, 5, 6, 7, 9, 10, 12, 15, 16, 17, 18, 19]), np.array([ 8, 14])), 12: (np.array([ 0, 1, 4, 5, 7, 9, 10, 15, 16, 17, 18, 19]), np.array([ 2, 3, 6, 8, 11, 13, 14])), 13: (np.array([ 0, 1, 4, 5, 6, 7, 9, 10, 12, 15, 16, 17, 18, 19]), np.array([ 8, 14])), 14: (np.array([ 0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 15, 16, 17, 18, 19]), np.array([], dtype=np.int8)), 15: (np.array([ 0, 1, 4, 5, 7, 9, 10, 16, 17]), np.array([ 2, 3, 6, 8, 11, 12, 13, 14, 18, 19])), 16: (np.array([ 0, 1, 4, 5, 7, 9, 10, 17]), np.array([ 2, 3, 6, 8, 11, 12, 13, 14, 15, 18, 19])), 17: (np.array([7]), np.array([ 0, 1, 2, 3, 4, 5, 6, 8, 9, 10, 11, 12, 13, 14, 15, 16, 18, 19])), 18: (np.array([ 0, 1, 4, 5, 7, 9, 10, 15, 16, 17]), np.array([ 2, 3, 6, 8, 11, 12, 13, 14, 19])), 19: (np.array([ 0, 1, 4, 5, 7, 9, 10, 15, 16, 17, 18]), np.array([ 2, 3, 6, 8, 11, 12, 13, 14]))}
    elif preserve_charged==True and len(excluded_indices) == 0:
        better_options = {0: (np.array([ 1,  4,  7,  9, 10, 17]), np.array([ 5,  6, 11, 12, 13, 15, 16, 18, 19])), 1: (np.array([ 4,  7,  9, 17]), np.array([ 0,  5,  6, 10, 11, 12, 13, 15, 16, 18, 19])), 2: (np.array([ 0,  1,  4,  5,  6,  7,  9, 10, 12, 15, 16, 17, 18, 19]), np.array([], dtype=np.int8)), 3: (np.array([ 0,  1,  4,  5,  6,  7,  9, 10, 12, 15, 16, 17, 18, 19]), np.array([], dtype=np.int8)), 4: (np.array([ 7,  9, 17]), np.array([ 0,  1,  5,  6, 10, 11, 12, 13, 15, 16, 18, 19])), 5: (np.array([ 0,  1,  4,  7,  9, 10, 17]), np.array([ 6, 11, 12, 13, 15, 16, 18, 19])), 6: (np.array([ 0,  1,  4,  5,  7,  9, 10, 12, 15, 16, 17, 18, 19]), np.array([11, 13])), 7: (np.array([], dtype=np.int8), np.array([ 0,  1,  4,  5,  6,  9, 10, 11, 12, 13, 15, 16, 17, 18, 19])), 8: (np.array([ 0,  1,  4,  5,  6,  7,  9, 10, 11, 12, 13, 15, 16, 17, 18, 19]), np.array([], dtype=np.int8)), 9: (np.array([ 7, 17]), np.array([ 0,  1,  4,  5,  6, 10, 11, 12, 13, 15, 16, 18, 19])), 10: (np.array([ 1,  4,  7,  9, 17]), np.array([ 0,  5,  6, 11, 12, 13, 15, 16, 18, 19])), 11: (np.array([ 0,  1,  4,  5,  6,  7,  9, 10, 12, 15, 16, 17, 18, 19]), np.array([], dtype=np.int8)), 12: (np.array([ 0,  1,  4,  5,  7,  9, 10, 15, 16, 17, 18, 19]), np.array([ 6, 11, 13])), 13: (np.array([ 0,  1,  4,  5,  6,  7,  9, 10, 12, 15, 16, 17, 18, 19]), np.array([], dtype=np.int8)), 14: (np.array([ 0,  1,  4,  5,  6,  7,  9, 10, 11, 12, 13, 15, 16, 17, 18, 19]), np.array([], dtype=np.int8)), 15: (np.array([ 0,  1,  4,  5,  7,  9, 10, 16, 17]), np.array([ 6, 11, 12, 13, 18, 19])), 16: (np.array([ 0,  1,  4,  5,  7,  9, 10, 17]), np.array([ 6, 11, 12, 13, 15, 18, 19])), 17: (np.array([7]), np.array([ 0,  1,  4,  5,  6,  9, 10, 11, 12, 13, 15, 16, 18, 19])), 18: (np.array([ 0,  1,  4,  5,  7,  9, 10, 15, 16, 17]), np.array([ 6, 11, 12, 13, 19])), 19: (np.array([ 0,  1,  4,  5,  7,  9, 10, 15, 16, 17, 18]), np.array([ 6, 11, 12, 13]))}
    else:
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
    for _ in range(max_iterations):
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
            
        # Process only active sequences in batches - fully vectorized version
        sequences_to_recalc = set()
        
        if active_sequences:
            # Convert to numpy array for vectorized operations
            active_indices = np.array(active_sequences)
            active_deltas = delta[active_indices]
            
            # Filter out sequences already within tolerance
            valid_mask = np.abs(active_deltas) > tolerance
            valid_indices = active_indices[valid_mask]
            
            # Further filter sequences that have mutable positions
            sequences_with_mutable = []
            for seq_idx in valid_indices:
                if len(mutable_positions_cache[seq_idx]) > 0:
                    sequences_with_mutable.append(seq_idx)
            
            if sequences_with_mutable:
                sequences_with_mutable = np.array(sequences_with_mutable)
                need_increase_mask = delta[sequences_with_mutable] < -tolerance
                
                # Determine the number of modifications for each sequence
                num_to_modify_per_seq = [min(batch_size, len(mutable_positions_cache[seq_idx])) for seq_idx in sequences_with_mutable]
                total_modifications = sum(num_to_modify_per_seq)

                if total_modifications > 0:
                    # Pre-allocate arrays
                    all_seq_indices = np.empty(total_modifications, dtype=np.int32)
                    all_positions = np.empty(total_modifications, dtype=np.int32)
                    all_need_increase = np.empty(total_modifications, dtype=bool)

                    # Fill arrays more efficiently
                    all_seq_indices = np.repeat(sequences_with_mutable, num_to_modify_per_seq)
                    all_need_increase = np.repeat(need_increase_mask, num_to_modify_per_seq)
                    
                    # Efficiently choose positions to modify
                    positions_to_modify_list = [np.random.choice(mutable_positions_cache[seq_idx], size=num, replace=False) 
                                                for seq_idx, num in zip(sequences_with_mutable, num_to_modify_per_seq)]
                    all_positions = np.concatenate(positions_to_modify_list)

                    # Get current amino acids at selected positions
                    all_current_aas = seq_matrices[all_seq_indices, all_positions]
                
                    # Vectorized replacement selection
                    all_replacements = np.full(total_modifications, -1, dtype=np.int32)
                    
                    # Vectorized replacement selection without looping through unique amino acids
                    increase_mask = all_need_increase
                    decrease_mask = ~all_need_increase

                    # Process increases
                    if np.any(increase_mask):
                        increase_indices = np.where(increase_mask)[0]
                        increase_aas = all_current_aas[increase_indices]
                        for aa in np.unique(increase_aas):
                            increase_opts = better_options[aa][0]
                            if len(increase_opts) > 0:
                                aa_specific_mask = (increase_aas == aa)
                                n_increase = np.sum(aa_specific_mask)
                                replacements = np.random.choice(increase_opts, size=n_increase)
                                all_replacements[increase_indices[aa_specific_mask]] = replacements

                    # Process decreases
                    if np.any(decrease_mask):
                        decrease_indices = np.where(decrease_mask)[0]
                        decrease_aas = all_current_aas[decrease_indices]
                        for aa in np.unique(decrease_aas):
                            decrease_opts = better_options[aa][1]
                            if len(decrease_opts) > 0:
                                aa_specific_mask = (decrease_aas == aa)
                                n_decrease = np.sum(aa_specific_mask)
                                replacements = np.random.choice(decrease_opts, size=n_decrease)
                                all_replacements[decrease_indices[aa_specific_mask]] = replacements
                    
                    # Apply only valid replacements
                    valid_mask = all_replacements >= 0
                    if np.any(valid_mask):
                        seq_matrices[all_seq_indices[valid_mask], all_positions[valid_mask]] = all_replacements[valid_mask]
                        sequences_to_recalc.update(np.unique(all_seq_indices[valid_mask]))
        
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
    else:
        return None

    # Convert back to sequences
    if convert_back_to_sequences:
        optimized_sequences = matrices_to_sequences(seq_matrices)
    else:
        optimized_sequences = seq_matrices
    
    return optimized_sequences

