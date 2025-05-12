'''
Code for taking an input sequence and modifying it as little as possible
to match user-specified kappa, FCR, NCPR, and hydropathy values.
'''

import numpy as np
import random
from goose.backend.amino_acids import AminoAcid
from goose.backend_vectorized.calculate_properties_vectorized import (
    sequences_to_matrices, calculate_fcr_batch, calculate_ncpr_batch, calculate_hydropathy_batch
)
from goose.backend_vectorized.calculate_kappa_vectorized import batch_calculate_kappa
from goose.backend.parameters import(
        HYDRO_ERROR, MAXIMUM_KAPPA_ERROR
)


def generate_candidates_vectorized(sequence, mutable_indices, aa_list):
    """
    Generate all possible mutation candidates using vectorized operations.
    
    Parameters:
    -----------
    sequence : list
        The current protein sequence as a list of characters
    mutable_indices : np.ndarray
        Indices of positions that can be mutated
    aa_list : list
        List of amino acids to use for mutations
    
    Returns:
    --------
    tuple
        (candidate_seqs, candidate_positions, candidate_aas)
    """
    seq_str = ''.join(sequence)
    
    # Create all combinations of positions and amino acids
    positions = np.repeat(mutable_indices, len(aa_list))
    aas = np.tile(np.array(aa_list), len(mutable_indices))
    
    # Get original amino acids at each position
    original_aas = np.array([sequence[i] for i in mutable_indices])
    original_aas_repeated = np.repeat(original_aas, len(aa_list))
    
    # Filter out combinations where amino acid doesn't change
    mask = original_aas_repeated != aas
    valid_positions = positions[mask]
    valid_aas = aas[mask]
    
    # Generate candidate sequences efficiently
    candidate_seqs = []
    for pos, aa in zip(valid_positions, valid_aas):
        candidate_seqs.append(seq_str[:pos] + aa + seq_str[pos+1:])
    
    return candidate_seqs, valid_positions.tolist(), valid_aas.tolist()


def minimal_sequence_modification(
    sequence,
    target_kappa=None,
    target_FCR=None,
    target_NCPR=None,
    target_hydropathy=None,
    tolerance_kappa=MAXIMUM_KAPPA_ERROR,
    tolerance_FCR=0.001,
    tolerance_NCPR=0.001,
    tolerance_hydropathy=HYDRO_ERROR,
    max_iterations=1000,
    protected_positions=None,
    weights=None,
    use_simulated_annealing=False,
    temperature=10.0,
    cooling_rate=0.95,
    early_stop_patience=50,
    verbose=True
):
    """
    Iteratively modify a protein sequence to match user-specified kappa, FCR, NCPR, and hydropathy values
    with minimal changes. Uses numpy vectorized operations for efficiency.
    
    Parameters:
    -----------
    sequence : str
        The input protein sequence
    target_kappa, target_FCR, target_NCPR, target_hydropathy : float, optional
        Target values for the respective properties
    tolerance_* : float
        Acceptable tolerance for each property
    max_iterations : int
        Maximum number of iterations to run
    protected_positions : list or None
        List of indices (0-based) that should not be mutated
    weights : dict or None
        Weights for each property (e.g., {'kappa': 2.0, 'FCR': 1.0})
    use_simulated_annealing : bool
        Whether to use simulated annealing to escape local minima
    temperature, cooling_rate : float
        Parameters for simulated annealing
    early_stop_patience : int
        Stop if no improvement for this many iterations
    verbose : bool
        Whether to print progress information
    
    Returns:
    --------
    dict
        Contains modified sequence, number of mutations, and final properties
    """
    aa_list = AminoAcid.standard_amino_acids
    seq = list(sequence.upper())
    seq_len = len(seq)
    
    # Initialize protected positions
    if protected_positions is None:
        protected_positions = []
    mutable_indices = np.array([i for i in range(seq_len) if i not in protected_positions])
    
    # Initialize weights
    if weights is None:
        weights = {'kappa': 1.0, 'FCR': 1.0, 'NCPR': 1.0, 'hydropathy': 1.0}
    
    # Vectorized property calculation
    def get_properties_batch(sequence_list):
        # Convert to matrices for vectorized property calculations
        seq_matrices = sequences_to_matrices(sequence_list)
        
        # Calculate all properties in batch
        fcr_values = calculate_fcr_batch(seq_matrices)
        ncpr_values = calculate_ncpr_batch(seq_matrices)
        hydro_values = calculate_hydropathy_batch(seq_matrices)
        kappa_values = batch_calculate_kappa(sequence_list)
        
        # Create dictionary of properties for each sequence
        properties_list = []
        for i in range(len(sequence_list)):
            props = {
                'kappa': kappa_values[i],
                'FCR': fcr_values[i],
                'NCPR': ncpr_values[i],
                'hydropathy': hydro_values[i]
            }
            properties_list.append(props)
        
        return properties_list
    
    def prop_distance_vectorized(props_list):
        """Calculate distance metrics for a batch of property dictionaries"""
        # Extract property arrays
        kappas = np.array([p['kappa'] for p in props_list])
        fcrs = np.array([p['FCR'] for p in props_list])
        ncprs = np.array([p['NCPR'] for p in props_list])
        hydropathies = np.array([p['hydropathy'] for p in props_list])
        
        # Initialize distances array
        distances = np.zeros(len(props_list))
        
        # Add weighted distances for each property if target is specified
        if target_kappa is not None:
            distances += weights.get('kappa', 1.0) * np.abs(kappas - target_kappa) / tolerance_kappa
        if target_FCR is not None:
            distances += weights.get('FCR', 1.0) * np.abs(fcrs - target_FCR) / tolerance_FCR
        if target_NCPR is not None:
            distances += weights.get('NCPR', 1.0) * np.abs(ncprs - target_NCPR) / tolerance_NCPR
        if target_hydropathy is not None:
            distances += weights.get('hydropathy', 1.0) * np.abs(hydropathies - target_hydropathy) / tolerance_hydropathy
            
        return distances

    # Track mutations
    original_seq = seq.copy()
    best_seq = seq.copy()
    best_dist = float('inf')
    best_props = None
    no_improvement_count = 0
    current_temp = temperature
    
    # Get initial properties
    current_seq_str = ''.join(seq)
    props_list = get_properties_batch([current_seq_str])
    current_props = props_list[0]
    
    # Set target props to original props if None
    if target_kappa is None:
        target_kappa = current_props['kappa']
    if target_FCR is None:
        target_FCR = current_props['FCR']
    if target_NCPR is None:
        target_NCPR = current_props['NCPR']
    if target_hydropathy is None:
        target_hydropathy = current_props['hydropathy']
    
    for iteration in range(max_iterations):
        current_seq_str = ''.join(seq)
        current_dist = prop_distance_vectorized([current_props])[0]
        
        # Update best solution if needed
        if current_dist < best_dist:
            best_dist = current_dist
            best_seq = seq.copy()
            best_props = current_props
            no_improvement_count = 0
        else:
            no_improvement_count += 1
        
        # Check early stopping
        if no_improvement_count >= early_stop_patience:
            if verbose:
                print(f"Early stopping after {iteration} iterations with no improvement")
            break
        
        # Check if all targets are within tolerance
        within = True
        if target_kappa is not None and abs(current_props['kappa'] - target_kappa) > tolerance_kappa:
            within = False
        if target_FCR is not None and abs(current_props['FCR'] - target_FCR) > tolerance_FCR:
            within = False
        if target_NCPR is not None and abs(current_props['NCPR'] - target_NCPR) > tolerance_NCPR:
            within = False
        if target_hydropathy is not None and abs(current_props['hydropathy'] - target_hydropathy) > tolerance_hydropathy:
            within = False
        if within:
            if verbose:
                print(f"All targets met after {iteration} iterations")
            break

        # Generate all possible single substitutions in one batch using vectorized operations
        candidate_seqs, candidate_positions, candidate_aas = generate_candidates_vectorized(seq, mutable_indices, aa_list)
        
        if not candidate_seqs:
            break
            
        # Evaluate all candidates at once using vectorized operations
        props_list = get_properties_batch(candidate_seqs)
        distances = prop_distance_vectorized(props_list)
        
        # Find the best substitutions
        min_dist = np.min(distances)
        
        if use_simulated_annealing and random.random() < 0.3:  # Sometimes try simulated annealing
            # Pick based on probability proportional to improvement
            probs = np.exp(-(distances - current_dist) / current_temp)
            probs = probs / np.sum(probs)
            chosen_idx = np.random.choice(len(candidate_seqs), p=probs)
            best_pos = candidate_positions[chosen_idx]
            best_aa = candidate_aas[chosen_idx]
            current_props = props_list[chosen_idx]
        else:
            # Regular greedy selection
            best_idxs = np.where(distances == min_dist)[0]
            chosen_idx = random.choice(best_idxs)
            best_pos = candidate_positions[chosen_idx]
            best_aa = candidate_aas[chosen_idx]
            current_props = props_list[chosen_idx]
        
        seq[best_pos] = best_aa
        current_temp *= cooling_rate
        
        if verbose and iteration % 10 == 0:
            print(f"Iteration {iteration}, distance: {current_dist:.4f}, temperature: {current_temp:.4f}")

    # Use the best solution found
    mutations = sum(1 for i, aa in enumerate(best_seq) if aa != original_seq[i])
    mutation_positions = [i for i, aa in enumerate(best_seq) if aa != original_seq[i]]
    
    return {
        'sequence': ''.join(best_seq),
        'mutations': mutations,
        'mutation_positions': mutation_positions,
        'original_sequence': ''.join(original_seq),
        'properties': best_props,
        'target_met': all([
            target_kappa is None or abs(best_props['kappa'] - target_kappa) <= tolerance_kappa,
            target_FCR is None or abs(best_props['FCR'] - target_FCR) <= tolerance_FCR,
            target_NCPR is None or abs(best_props['NCPR'] - target_NCPR) <= tolerance_NCPR,
            target_hydropathy is None or abs(best_props['hydropathy'] - target_hydropathy) <= tolerance_hydropathy
        ])
    }
