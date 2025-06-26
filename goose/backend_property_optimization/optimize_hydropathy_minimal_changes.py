'''
Function that optimizes hydropathy while minimizing the number
of amino acids changed. This is not necessarily the most efficient,
but the objective is to minimize changes not to optimize efficiency. 
'''

import numpy as np
from sparrow.protein import Protein
from goose import parameters

AA_hydro = {"A": 6.3,  "R": 0.0,  "N": 1.0,  "D": 1.0,  "C": 7.0,  "Q": 1.0,  "E": 1.0,  "G": 4.1,  "H": 1.3,  "I": 9.0,  "L": 8.3,  "K": 0.6,  "M": 6.4,  "F": 7.3,  "P": 2.9,  "S": 3.7,  "T": 3.8,  "W": 3.6,  "Y": 3.2,  "V": 8.7  }


def optimize_hydropathy_minimal_changes(input_sequence, target_hydropathy, max_iterations=100, 
                                        tolerance=parameters.HYDRO_ERROR, preserve_charged=True):
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
        set by parameters.HYDRO_ERROR
    preserve_charged : bool
        If True, charged residues (D, E, K, R) will not be modified.

    Returns:
    --------
    str
        Optimized sequence with hydropathy close to the target.
    """
    
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
