'''
Within class hydropathy optimization module. 

This is a restricted optimization approach in that every amino acid must
stay within its class.
Classes:
- Hydrophobic: A, I, L, M, V
- Aromatic: F, W, Y
- Polar: S, T, N, Q
- Positive: D, E
- Negative: K, R
- Glycine: G
- Cysteine: C
- Proline: P
- Histidine: H
'''

import numpy as np
from typing import List, Union
from goose import parameters
from goose.backend_property_calculation.calculate_properties_vectorized import calculate_hydropathy_batch, sequences_to_matrices, matrices_to_sequences


def optimize_hydropathy_within_class_vectorized(
    sequences: Union[List[str], np.ndarray],
    target_hydropathy: float, 
    min_batch_size: int = 5,
    max_iterations: int = 5000,
    tolerance: float = parameters.HYDRO_ERROR,
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
    
    # Convert back to sequences
    if len(seq_matrices) > 0:
        optimized_sequences = matrices_to_sequences(seq_matrices)
    else:
        optimized_sequences = None
    
    return optimized_sequences
