'''
code specifically for optimizing the hydropathy value after sequences are generated
uses numpy vectorized operations to maximize efficiency.

Includes:
- General hydropathy optimization function
- Within class hydropathy optimization function
- Minimal changes hydropathy optimization function
- Hydropathy optimization that avoids original residues

'''
from typing import List, Union
import numpy as np
from sparrow.protein import Protein
from goose import parameters
from goose.backend_property_calculation.calculate_properties_batch import calculate_hydropathy_batch, sequences_to_matrices, matrices_to_sequences
from goose.backend_property_calculation.calculate_properties_single_sequence import sequence_to_array, array_to_sequence, calculate_hydropathy_single_sequence
from goose.data.defined_aa_classes import aa_class_indices

#-=-=-=-=- General hydropathy optimization function -=-=-=-=-=-

def optimize_hydropathy(
    sequences: Union[List[str], np.ndarray],
    target_hydropathy: float, 
    preserve_charged: bool = True,
    max_iterations: int = 1000,
    tolerance: float = parameters.MAXIMUM_HYDRO_ERROR,
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


#-=-=-=-=- Within class hydropathy optimization function -=-=-=-=-=-

def optimize_hydropathy_within_class(
    sequences: Union[List[str], np.ndarray],
    target_hydropathy: float, 
    min_batch_size: int = 5,
    max_iterations: int = 5000,
    tolerance: float = parameters.MAXIMUM_HYDRO_ERROR,
    batch_size: int = 10, 
    only_return_within_tolerance: bool = True,
    return_when_num_hit=None) -> List[str]:
    """
    Optimize protein sequences to match a target hydropathy while keeping amino acids within their classes.
    
    Amino acid classes:
    - Hydrophobic: A, I, L, M, V
    - Aromatic: F, W, Y  
    - Polar: S, T, N, Q
    - Positive: D, E
    - Negative: K, R
    - Glycine: G (cannot change)
    - Cysteine: C (cannot change)
    - Proline: P (cannot change)
    - Histidine: H (cannot change)
    
    Parameters:
    -----------
    sequences : List[str] or np.ndarray
        Input sequences to optimize
    target_hydropathy : float
        Target mean hydropathy value to achieve
    min_batch_size : int
        Minimum batch size for processing sequences
        This is used if a single sequence is provided        
    max_iterations : int
        Maximum number of optimization iterations
    tolerance : float
        Acceptable difference between achieved and target hydropathy
    batch_size : int
        Number of positions to modify in each sequence per iteration
    only_return_within_tolerance : bool
        If True, only return sequences that are within the tolerance
    return_when_num_hit : int or None
        If specified, return when this number of sequences are within tolerance
        
    Returns:
    --------
    List[str]
        Optimized sequences with hydropathy values close to target
    """
    # Hydropathy scale constants (same order as AA_TO_INT)
    HYDROPATHY_SCALE = np.array([6.3, 7.0, 1.0, 1.0, 7.3, 4.1, 1.3, 9.0, 
                                0.6, 8.3, 6.4, 1.0, 2.9, 1.0, 0.0, 3.7, 3.8, 8.7, 3.6, 3.2])

    # Convert input to proper format
    if isinstance(sequences, str):
        sequences = [sequences for _ in range(min_batch_size)]
    if not isinstance(sequences, np.ndarray):  
        sequences = np.array(sequences)

    if return_when_num_hit is not None:
        only_return_within_tolerance = True
        if return_when_num_hit > len(sequences):
            return_when_num_hit = len(sequences)        
    else:
        return_when_num_hit = len(sequences)
    
    # Create amino acid mappings
    aa_to_int = {'A': 0, 'C': 1, 'D': 2, 'E': 3, 'F': 4, 'G': 5, 'H': 6, 'I': 7, 
                'K': 8, 'L': 9, 'M': 10, 'N': 11, 'P': 12, 'Q': 13, 'R': 14, 
                'S': 15, 'T': 16, 'V': 17, 'W': 18, 'Y': 19}
    
    # Define amino acid classes (using indices)
    aa_classes = {
        # Hydrophobic: A, I, L, M, V
        0: [0, 7, 9, 10, 17],   # A -> A, I, L, M, V
        7: [0, 7, 9, 10, 17],   # I -> A, I, L, M, V  
        9: [0, 7, 9, 10, 17],   # L -> A, I, L, M, V
        10: [0, 7, 9, 10, 17],  # M -> A, I, L, M, V
        17: [0, 7, 9, 10, 17],  # V -> A, I, L, M, V
        
        # Aromatic: F, W, Y
        4: [4, 18, 19],         # F -> F, W, Y
        18: [4, 18, 19],        # W -> F, W, Y
        19: [4, 18, 19],        # Y -> F, W, Y
        
        # Polar: S, T, N, Q
        15: [15, 16, 11, 13],   # S -> S, T, N, Q
        16: [15, 16, 11, 13],   # T -> S, T, N, Q
        11: [15, 16, 11, 13],   # N -> S, T, N, Q
        13: [15, 16, 11, 13],   # Q -> S, T, N, Q
        
        # Positive: D, E
        2: [2, 3],              # D -> D, E
        3: [2, 3],              # E -> D, E
        
        # Negative: K, R
        8: [8, 14],             # K -> K, R
        14: [8, 14],            # R -> K, R
        
        # Immutable residues
        5: [5],                 # G -> G (glycine)
        1: [1],                 # C -> C (cysteine)
        12: [12],               # P -> P (proline)
        6: [6],                 # H -> H (histidine)
    }
    
    # Convert sequences to matrices
    seq_matrices = sequences_to_matrices(sequences)
    n_sequences, seq_length = seq_matrices.shape
    
    # Precompute replacement options for each amino acid within its class
    within_class_options = {}  # aa -> (increase_options, decrease_options)
    
    for aa in range(20):
        current_hydro = HYDROPATHY_SCALE[aa]
        class_members = aa_classes.get(aa, [aa])  # Default to self if not in any class
        
        # Find class members that increase hydropathy
        increase_opts = [i for i in class_members 
                        if HYDROPATHY_SCALE[i] > current_hydro]
        
        # Find class members that decrease hydropathy  
        decrease_opts = [i for i in class_members
                        if HYDROPATHY_SCALE[i] < current_hydro]
        
        within_class_options[aa] = (np.array(increase_opts), np.array(decrease_opts))
    
    # Calculate initial hydropathy and track completed sequences
    current_hydro = calculate_hydropathy_batch(seq_matrices)
    completed_sequences = set()
    
    # Precompute mutable positions for each sequence (all positions are mutable within class)
    mutable_positions_cache = {}
    for seq_idx in range(n_sequences):
        # All positions are potentially mutable, but effectiveness depends on class options
        mutable_positions_cache[seq_idx] = np.arange(seq_length)
    
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
            
        # Process only active sequences
        sequences_to_recalc = set()
        
        for seq_idx in active_sequences:
            current_delta = delta[seq_idx]
            
            # Skip if already within tolerance
            if abs(current_delta) <= tolerance:
                continue
                
            mutable_pos = mutable_positions_cache[seq_idx]
            if len(mutable_pos) == 0:
                continue
                
            # Determine direction needed
            need_increase = current_delta < -tolerance
            
            # Find positions that can actually be modified in the desired direction
            # The key fix: we need to check if ANY amino acid in the same class
            # can help move toward the target, not just those above/below current AA
            effective_positions = []
            for pos in mutable_pos:
                current_aa = seq_matrices[seq_idx, pos]
                current_aa_hydro = HYDROPATHY_SCALE[current_aa]
                class_members = aa_classes.get(current_aa, [current_aa])
                
                # Check if any class member can help move toward target
                can_help = False
                for class_member in class_members:
                    if class_member == current_aa:
                        continue  # Skip self
                    
                    member_hydro = HYDROPATHY_SCALE[class_member]
                    
                    if need_increase:
                        # Need to increase hydropathy: any class member with higher hydropathy helps
                        if member_hydro > current_aa_hydro:
                            can_help = True
                            break
                    else:
                        # Need to decrease hydropathy: any class member with lower hydropathy helps
                        if member_hydro < current_aa_hydro:
                            can_help = True
                            break
                
                if can_help:
                    effective_positions.append(pos)
            
            if not effective_positions:
                continue  # No effective modifications possible for this sequence
                
            # Select positions to modify
            num_to_modify = min(batch_size, len(effective_positions))
            positions_to_modify = np.random.choice(effective_positions, size=num_to_modify, replace=False)
            # Make modifications
            for pos in positions_to_modify:
                current_aa = seq_matrices[seq_idx, pos]
                increase_opts, decrease_opts = within_class_options[current_aa]
                
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
    if only_return_within_tolerance:
        final_hydro = calculate_hydropathy_batch(seq_matrices)
        within_tolerance = np.abs(final_hydro - target_hydropathy) <= tolerance
        seq_matrices = seq_matrices[within_tolerance]
    
    # sort by how close they are to the target hydropathy
    if len(seq_matrices) > 1:
        final_hydro = calculate_hydropathy_batch(seq_matrices)
        sorted_indices = np.argsort(np.abs(final_hydro - target_hydropathy))
        seq_matrices = seq_matrices[sorted_indices]

    # Convert back to sequences
    if len(seq_matrices) > 0:
        optimized_sequences = matrices_to_sequences(seq_matrices)
    else:
        optimized_sequences = None
    
    return optimized_sequences


#-=-=-=-=- Minimal changes hydropathy optimization function -=-=-=-=-=-

def optimize_hydropathy_minimal_changes(input_sequence, target_hydropathy, max_iterations=100, 
                                        tolerance=parameters.MAXIMUM_HYDRO_ERROR, preserve_charged=True):
    """
    Optimize hydropathy of a sequence by making minimal changes to achieve a target hydropathy value.
    
    This function strategically selects amino acids that provide maximum improvement toward the target
    hydropathy at each iteration, minimizing the total number of changes needed.

    Parameters:
    -----------
    input_sequence : str
        The input protein sequence.
    target_hydropathy : float
        The target mean hydropathy value to achieve.
    max_iterations : int
        Maximum number of optimization iterations.
    tolerance : float
        Acceptable difference between achieved and target hydropathy.
        set by parameters.MAXIMUM_HYDRO_ERROR
    preserve_charged : bool
        If True, charged residues (D, E, K, R) will not be modified.

    Returns:
    --------
    str
        Optimized sequence with hydropathy close to the target.
    """
    # set aa hydro
    AA_hydro = {"A": 6.3,  "R": 0.0,  "N": 1.0,  "D": 1.0,  "C": 7.0,  "Q": 1.0,  "E": 1.0,  "G": 4.1,  "H": 1.3,  "I": 9.0,  "L": 8.3,  "K": 0.6,  "M": 6.4,  "F": 7.3,  "P": 2.9,  "S": 3.7,  "T": 3.8,  "W": 3.6,  "Y": 3.2,  "V": 8.7  }
    
    # Define charged amino acids
    charged_residues = {'D', 'E', 'K', 'R'} if preserve_charged else set()
    
    # All amino acids sorted by hydropathy for efficient selection
    aa_by_hydropathy = sorted(AA_hydro.items(), key=lambda x: x[1])
    
    def find_optimal_substitution(sequence, target_hydro, current_hydro, excluded_positions=None):
        """
        Find all optimal substitutions that provide maximal improvement and randomly select one.
        This ensures less deterministic behavior and better exploration of the sequence space.
        
        Returns:
        --------
        tuple: (position, old_aa, new_aa, improvement, num_optimal_positions) or None
        """
        if excluded_positions is None:
            excluded_positions = set()
            
        all_substitutions = []
        best_final_distance = float('inf')
        
        sequence_length = len(sequence)
        current_distance = abs(current_hydro - target_hydro)
        
        # First pass: Find all possible substitutions and their resulting distances
        for pos in range(sequence_length):
            if pos in excluded_positions:
                continue
                
            current_aa = sequence[pos]
            
            # Skip if amino acid is protected
            if current_aa in charged_residues:
                continue
                
            current_aa_hydro = AA_hydro[current_aa]
            
            # Try all possible amino acid substitutions
            for candidate_aa, candidate_hydro in AA_hydro.items():
                if candidate_aa == current_aa or candidate_aa in charged_residues:
                    continue
                    
                # Calculate resulting hydropathy
                hydro_change = candidate_hydro - current_aa_hydro
                new_hydro = current_hydro + (hydro_change / sequence_length)
                new_distance = abs(new_hydro - target_hydro)
                
                # Only consider if it improves (gets closer to target)
                if new_distance < current_distance:
                    improvement = current_distance - new_distance
                    all_substitutions.append((pos, current_aa, candidate_aa, improvement, new_distance))
                    
                    # Track the best distance achieved
                    if new_distance < best_final_distance:
                        best_final_distance = new_distance
        
        if not all_substitutions:
            return None
            
        # Second pass: Find all substitutions that achieve the maximal improvement (best distance)
        # Allow for small numerical tolerance
        tolerance = 1e-10
        optimal_substitutions = []
        
        for substitution in all_substitutions:
            pos, current_aa, candidate_aa, improvement, final_distance = substitution
            if abs(final_distance - best_final_distance) <= tolerance:
                optimal_substitutions.append((pos, current_aa, candidate_aa, improvement))
        
        # Randomly select one of the optimal substitutions
        if optimal_substitutions:
            selected_substitution = np.random.choice(len(optimal_substitutions))
            pos, old_aa, new_aa, improvement = optimal_substitutions[selected_substitution]
            return (pos, old_aa, new_aa, improvement, len(optimal_substitutions))
        
        return None
    
    # Start optimization
    current_hydropathy = Protein(input_sequence).hydrophobicity
    sequence = list(input_sequence)
    changes_made = []
    
    for iteration in range(max_iterations):
        # Check if we've reached the target
        current_hydropathy = Protein(''.join(sequence)).hydrophobicity
        distance_to_target = abs(current_hydropathy - target_hydropathy)
        
        if distance_to_target <= tolerance:
            break
            
        # Find the best substitution for this iteration
        excluded_positions = {change[0] for change in changes_made}  # Don't re-modify same positions
        substitution = find_optimal_substitution(sequence, target_hydropathy, current_hydropathy, excluded_positions)
        
        if substitution is None:
            break
            
        pos, old_aa, new_aa, improvement, num_optimal = substitution
        
        # Make the substitution
        sequence[pos] = new_aa
        new_hydropathy = Protein(''.join(sequence)).hydrophobicity
        changes_made.append((pos, old_aa, new_aa))        
        current_hydropathy = new_hydropathy
    
    final_sequence = ''.join(sequence)
    return final_sequence



def optimize_hydropathy_within_class_avoid_original_residues(
        original_sequence,
        variant_sequence,
        target_hydropathy,
        max_iterations=1000,
        tolerance=0.05):
    """
    Optimize hydropathy of a sequence while avoiding original residues.

    Parameters:
    -----------
    original_sequence : str
        The original sequence to avoid.
    variant_sequence : str
        The sequence to optimize.
    target_hydropathy : float
        Target mean hydropathy value to achieve.
    max_iterations : int
        Maximum number of optimization iterations.
    tolerance : float
        Acceptable difference between achieved and target hydropathy.

    Returns:
    --------
    str
        Optimized sequence with hydropathy close to the target.
    """
    
    # Convert sequences to numpy arrays for vectorized operations
    original_seq = sequence_to_array(original_sequence)
    variant_seq = sequence_to_array(variant_sequence)
    
    # set mutable residues. Can't change D,E,K,R,H,G,C,P
    mutable_residues=[0, 4, 7, 8, 9, 10, 11, 13, 14, 15, 16, 17, 18, 19]
    
    current_hydropathy = calculate_hydropathy_single_sequence(variant_seq)
    
    for _ in range(max_iterations):
        if abs(current_hydropathy - target_hydropathy) < tolerance:
            break
        
        # Find positions that can be changed
        changeable_positions = [i for i in range(len(variant_seq)) if variant_seq[i] in mutable_residues]
        
        if not changeable_positions:
            return ''.join(array_to_sequence(variant_seq))
        
        # Randomly select a position to change
        pos_to_change = np.random.choice(changeable_positions)
        
        # Find the best replacement amino acid
        current_aa = variant_seq[pos_to_change]
        best_hydro_error = np.inf
        
        for aa in aa_class_indices[current_aa]:
            if aa not in original_seq:  # Avoid original residues
                new_seq = variant_seq.copy()
                new_seq[pos_to_change] = aa
                new_hydro = calculate_hydropathy_single_sequence(new_seq)

                if abs(new_hydro-target_hydropathy) < best_hydro_error:
                    best_hydro_error = abs(new_hydro - target_hydropathy)
                    variant_seq = new_seq.copy()
                    if best_hydro_error < tolerance:
                        # If we found a suitable amino acid, break early
                        break
    return array_to_sequence(variant_seq)
