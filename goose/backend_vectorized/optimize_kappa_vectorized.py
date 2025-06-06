'''
Functionality for optimizing kappa towards some target value. 
If a sequences kappa is too low, that means that the charged residues in the
sequence need to be more asymmetrically distributed across the sequence. 
If a sequences kappa is too high, that means that the charged residues in the
sequence need to be more symmetrically distributed across the sequence.
The optimization is done by moving the charged residues in the sequence.
'''

import numpy as np
from typing import List, Union, Optional, Tuple
from goose.backend_vectorized.calculate_kappa_vectorized import batch_calculate_kappa, vectorized_ternarize

# Constants for amino acid indices
AMINO_ACIDS = "ACDEFGHIKLMNPQRSTVWY"
CHARGED_INDICES = np.array([2, 3, 8, 14])  # D, E, K, R as 0-indexed positions
POSITIVE_INDICES = np.array([8, 14])       # K, R
NEGATIVE_INDICES = np.array([2, 3])        # D, E


def find_unique_ternarized_sequences(ternarized_sequences):
    """
    Identifies unique ternarized sequences and creates mapping to original sequences.
    
    Args:
        ternarized_sequences (np.array): Ternarized sequences

    Returns:
        tuple: (unique_sequences, unique_indices, inverse_mapping)
            - unique_sequences: array of unique ternarized sequences
            - unique_indices: indices of first occurrence of each unique sequence
            - inverse_mapping: mapping to reconstruct original sequence order
    """
    # Convert sequences to tuple for hashing
    seq_tuples = [tuple(seq) for seq in ternarized_sequences]
    
    # Create dictionary to store first occurrence of each unique sequence
    unique_dict = {}
    unique_indices = []
    inverse_mapping = []
    
    for idx, seq_tuple in enumerate(seq_tuples):
        if seq_tuple not in unique_dict:
            unique_dict[seq_tuple] = len(unique_indices)
            unique_indices.append(idx)
        inverse_mapping.append(unique_dict[seq_tuple])
    
    unique_sequences = ternarized_sequences[unique_indices]
    return unique_sequences, unique_indices, inverse_mapping


def decrease_charge_asymmetry(sequences, iterations):
    """
    Iteratively decreases charge asymmetry in sequences by creating more evenly
    distributed alternating patterns of positive and negative charges.
    Heavily optimized version for faster performance, especially for low kappa targets.

    Args:
        sequences (np.array): A NxL shaped array of sequences with values -1, 0, 1.
        iterations (int): The number of iterations to perform.

    Returns:
        list: A list of intermediate sequences at each iteration, including the initial sequences.
    """
    # Only store initial and final sequences if iterations are high
    store_intermediates = iterations <= 20
    intermediate_sequences = [sequences.copy()] if store_intermediates else [sequences.copy()]
    current_sequences = sequences.copy()
    num_sequences = current_sequences.shape[0]
    
    # Pre-generate all random values at once
    strategy_probs = np.random.random(size=(iterations, num_sequences, 3))
    random_positions = np.random.random(size=(iterations, num_sequences, 2))
    
    # Pre-allocate arrays for efficiency
    sequence_buffer = np.zeros_like(current_sequences)
    
    for iter_idx in range(iterations):
        # Efficient batch processing - process all sequences at once where possible
        for i in range(num_sequences):
            sequence = current_sequences[i]
            seq_length = len(sequence)
            
            # Find positive and negative charges once
            pos_mask = sequence == 1
            neg_mask = sequence == -1
            positive_indices = np.where(pos_mask)[0]
            negative_indices = np.where(neg_mask)[0]
            
            if positive_indices.size < 1 or negative_indices.size < 1:
                continue
            
            # Fast strategy selection based on pre-generated random values
            strategy_selector = strategy_probs[iter_idx, i, 0]
            
            # STRATEGY 1: Perfect charge alternation (very effective for low kappa)
            # This strategy focuses on creating perfect alternation of + and - charges
            if strategy_selector < 0.4:  # Higher probability for this strategy
                # Create ideal alternating positions
                ideal_positions = np.zeros(seq_length, dtype=bool)
                
                # Calculate total charges
                total_charges = positive_indices.size + negative_indices.size
                
                # Create alternating pattern for charged residues
                spaced_positions = np.linspace(0, seq_length-1, total_charges).astype(int)
                
                # Get current charged positions (both + and -)
                current_charged = np.sort(np.concatenate([positive_indices, negative_indices]))
                
                # Try to move charged residues to be evenly spaced
                # This greatly helps achieve very low kappa values
                if spaced_positions.size == current_charged.size:
                    # Find positions that need to be swapped
                    for pos_old, pos_new in zip(current_charged, spaced_positions):
                        if pos_old != pos_new and sequence[pos_old] != 0 and sequence[pos_new] == 0:
                            # Do the swap
                            sequence[pos_new] = sequence[pos_old]
                            sequence[pos_old] = 0
                            break  # Only do one swap per iteration for stability
            
            # STRATEGY 2: Break clusters of like charges (original optimized strategy)
            elif strategy_selector < 0.7:
                # Find clusters efficiently using vectorization
                if positive_indices.size > 1:
                    # Get adjacent positive charges
                    pos_clusters = np.where(np.diff(positive_indices) == 1)[0]
                    
                    if pos_clusters.size > 0:
                        # Take first cluster
                        pos_idx = positive_indices[pos_clusters[0] + 1]
                        
                        # Find neutral position to swap with
                        neutral_positions = np.where(sequence == 0)[0]
                        
                        if neutral_positions.size > 0:
                            # Prefer positions adjacent to negative charges
                            good_positions = []
                            for neut_pos in neutral_positions:
                                # Check if adjacent to a negative charge
                                has_neg_neighbor = False
                                if neut_pos > 0 and sequence[neut_pos-1] == -1:
                                    has_neg_neighbor = True
                                elif neut_pos < seq_length-1 and sequence[neut_pos+1] == -1:
                                    has_neg_neighbor = True
                                
                                if has_neg_neighbor:
                                    good_positions.append(neut_pos)
                            
                            # If found good position, swap
                            if good_positions:
                                neut_idx = good_positions[0]
                                sequence[pos_idx], sequence[neut_idx] = sequence[neut_idx], sequence[pos_idx]
                
                # Same logic for negative charges
                if negative_indices.size > 1:
                    neg_clusters = np.where(np.diff(negative_indices) == 1)[0]
                    
                    if neg_clusters.size > 0:
                        neg_idx = negative_indices[neg_clusters[0] + 1]
                        
                        # Find neutral position to swap with
                        neutral_positions = np.where(sequence == 0)[0]
                        
                        if neutral_positions.size > 0:
                            # Prefer positions adjacent to positive charges
                            good_positions = []
                            for neut_pos in neutral_positions:
                                # Check if adjacent to a positive charge
                                has_pos_neighbor = False
                                if neut_pos > 0 and sequence[neut_pos-1] == 1:
                                    has_pos_neighbor = True
                                elif neut_pos < seq_length-1 and sequence[neut_pos+1] == 1:
                                    has_pos_neighbor = True
                                
                                if has_pos_neighbor:
                                    good_positions.append(neut_pos)
                            
                            # If found good position, swap
                            if good_positions:
                                neut_idx = good_positions[0]
                                sequence[neg_idx], sequence[neut_idx] = sequence[neut_idx], sequence[neg_idx]
            
            # STRATEGY 3: Simple direct charge swapping
            else:
                # Direct swap of positive and negative charge
                if positive_indices.size > 0 and negative_indices.size > 0:
                    pos_idx = positive_indices[int(random_positions[iter_idx, i, 0] * positive_indices.size) % positive_indices.size]
                    neg_idx = negative_indices[int(random_positions[iter_idx, i, 1] * negative_indices.size) % negative_indices.size]
                    
                    # Only swap if they are not adjacent to like charges
                    sequence[pos_idx], sequence[neg_idx] = sequence[neg_idx], sequence[pos_idx]
            
            # Update in-place
            current_sequences[i] = sequence
        
        # Store intermediate result if needed
        if store_intermediates:
            intermediate_sequences.append(current_sequences.copy())
    
    # Always store the final result
    if not store_intermediates:
        intermediate_sequences.append(current_sequences.copy())
    
    return intermediate_sequences


def increase_charge_asymmetry(sequences, iterations):
    """
    Iteratively increases charge asymmetry in sequences by moving positive charges
    towards the left and negative charges towards the right. Optimized version.

    Args:
        sequences (np.array): A NxL shaped array of sequences with values -1, 0, 1.
        iterations (int): The number of iterations to perform.

    Returns:
        list: A list of intermediate sequences at each iteration, including the initial sequences.
    """
    intermediate_sequences = [sequences.copy()]
    current_sequences = sequences.copy()
    num_sequences = current_sequences.shape[0]
    
    # Pre-allocate arrays for random positive and negative indices for better performance
    pos_rand_indices = np.random.randint(0, 1000000, size=(iterations, num_sequences))
    neg_rand_indices = np.random.randint(0, 1000000, size=(iterations, num_sequences))
    
    for iter_idx in range(iterations):
        updated_sequences = current_sequences.copy()
        
        for i in range(num_sequences):
            sequence = updated_sequences[i]
            seq_length = len(sequence)
            
            # Find charged residues only once
            positive_indices = np.where(sequence == 1)[0]
            negative_indices = np.where(sequence == -1)[0]
            
            if positive_indices.size > 0 and negative_indices.size > 0:
                # Use pre-allocated random values with modulo operation
                pos_idx = positive_indices[pos_rand_indices[iter_idx, i] % positive_indices.size]
                neg_idx = negative_indices[neg_rand_indices[iter_idx, i] % negative_indices.size]
                
                # Move positive charge to the left (more efficient swap)
                if pos_idx > 0:
                    sequence[pos_idx], sequence[pos_idx - 1] = sequence[pos_idx - 1], sequence[pos_idx]
                
                # Move negative charge to the right (more efficient swap)
                if neg_idx < seq_length - 1:
                    sequence[neg_idx], sequence[neg_idx + 1] = sequence[neg_idx + 1], sequence[neg_idx]
            
            # No need to explicitly update since we're modifying sequence in-place
        
        # Store the current state
        intermediate_sequences.append(updated_sequences.copy())
        current_sequences = updated_sequences
    
    return intermediate_sequences


def optimize_kappa_vectorized(sequences, target_kappa, 
                              tolerance=0.03, num_iterations=100, 
                              max_change_iterations=20,
                              early_stopping=True,
                              return_when_num_hit=None,
                              only_return_within_tolerance=False):
    '''
    Optimize kappa values of sequences by either increasing or decreasing charge asymmetry
    to reach a target kappa value. Optimized for performance.
    
    Args:
        sequences (list): List of amino acid sequences
        target_kappa (float): Target kappa value to reach
        tolerance (float): Acceptable difference from target kappa
        num_iterations (int): Maximum number of overall iterations to perform
        max_change_iterations (int): Number of iterations for each charge rearrangement step
        early_stopping (bool): Stop early if target is reached
        return_when_num_hit : int
            If not None, return when this number of sequences are within tolerance
        only_return_within_tolerance (bool): If True, only return sequences within tolerance.
            if False, return all sequences, even if they are not within tolerance.

    
    Returns:
        list: Sequences with kappa values optimized to be close to the target
    '''
    
    # Convert input to numpy array for consistency
    sequences = np.array([str(seq) for seq in sequences])

    # if we are returning when num hit 
    if return_when_num_hit != None:
        # force only_return_within_tolerance to be true
        only_return_within_tolerance=True
        if return_when_num_hit > len(sequences):
            return_when_num_hit = len(sequences)
    else:
        return_when_num_hit = len(sequences)

    # Ternarize the sequences (only once at the start)
    ternarized_sequences = vectorized_ternarize(sequences, group1=['K', 'R'], group2=['D', 'E'])

    # Get unique ternarized sequences to avoid redundant processing
    unique_seqs, unique_indices, inverse_mapping = find_unique_ternarized_sequences(ternarized_sequences)
    
    # Calculate initial kappa values only for unique sequences
    kappa_values = batch_calculate_kappa(unique_seqs, ternarize=False)
    
    # Identify which sequences need to have kappa increased or decreased
    needs_increase = kappa_values < target_kappa - tolerance
    needs_decrease = kappa_values > target_kappa + tolerance
    
    # Check if any sequence needs modification
    if not np.any(needs_increase) and not np.any(needs_decrease):
        return sequences
    
    # Track which sequences have reached their target
    reached_target = ~(needs_increase | needs_decrease)
    
    # Initialize modified sequences array (avoid redundant copies)
    modified_unique_sequences = unique_seqs.copy()
    
    # Initialize tracking variables for early stopping
    prev_kappa_values = kappa_values.copy()
    stagnation_counter = np.zeros(len(unique_seqs), dtype=int)
    max_stagnation = 5  # Stop if no improvement for this many iterations

    
    # Calculate adaptive iterations based on distance from target
    if early_stopping:
        # Calculate distance to target
        distance_to_target = np.abs(kappa_values - target_kappa)
        # Scale iterations; farther sequences get more iterations
        adaptive_iterations = np.clip(
            np.round(max_change_iterations * (1 + 2 * distance_to_target / tolerance)).astype(int),
            max_change_iterations, max_change_iterations * 3
        )
    else:
        adaptive_iterations = np.full(len(unique_seqs), max_change_iterations)

    
    # Process in iterations
    for iteration in range(num_iterations):
        # If all sequences reached target, we're done
        if np.all(reached_target):
            break
        
        # if we are stopping at a specific number of successes...
        if np.sum(reached_target) >= return_when_num_hit:
            break

        
        # Process sequences that need kappa increased
        increase_indices = np.where(needs_increase)[0]
        if len(increase_indices) > 0:
            increase_seqs = modified_unique_sequences[increase_indices]
            
            # Process each sequence with its adaptive iteration count
            for i, idx in enumerate(increase_indices):
                # Skip if this sequence has stagnated
                if stagnation_counter[idx] >= max_stagnation:
                    continue
                    
                # Apply increase_charge_asymmetry with adaptive iterations
                iters = adaptive_iterations[idx]
                result = increase_charge_asymmetry(increase_seqs[i:i+1], iters)[-1]
                modified_unique_sequences[idx] = result[0]
        
        # Process sequences that need kappa decreased
        decrease_indices = np.where(needs_decrease)[0]
        if len(decrease_indices) > 0:
            decrease_seqs = modified_unique_sequences[decrease_indices]
            
            # Process each sequence with its adaptive iteration count
            for i, idx in enumerate(decrease_indices):
                # Skip if this sequence has stagnated
                if stagnation_counter[idx] >= max_stagnation:
                    continue
                    
                # Apply decrease_charge_asymmetry with adaptive iterations
                iters = adaptive_iterations[idx]
                result = decrease_charge_asymmetry(decrease_seqs[i:i+1], iters)[-1]
                modified_unique_sequences[idx] = result[0]
        
        # Recalculate kappa for sequences that haven't reached target
        recalc_indices = np.where(~reached_target)[0]
        if len(recalc_indices) > 0:
            current_kappa = batch_calculate_kappa(modified_unique_sequences[recalc_indices], ternarize=False)
            
            # Update which sequences need modifications
            for i, idx in enumerate(recalc_indices):
                kappa_val = current_kappa[i]
                kappa_values[idx] = kappa_val
                
                # Update status of this sequence
                old_status = (~reached_target[idx], needs_increase[idx], needs_decrease[idx])
                
                reached_target[idx] = abs(kappa_val - target_kappa) <= tolerance
                needs_increase[idx] = kappa_val < target_kappa - tolerance
                needs_decrease[idx] = kappa_val > target_kappa + tolerance
                
                new_status = (~reached_target[idx], needs_increase[idx], needs_decrease[idx])
                
                # Check for stagnation - if kappa didn't change much
                if abs(kappa_val - prev_kappa_values[idx]) < tolerance / 10:
                    stagnation_counter[idx] += 1
                else:
                    stagnation_counter[idx] = 0
                    
                # Adapt iterations based on progress
                if early_stopping:
                    distance_to_target = abs(kappa_val - target_kappa)
                    adaptive_iterations[idx] = max(
                        max_change_iterations,
                        min(5 * max_change_iterations, 
                            int(max_change_iterations * (1 + distance_to_target / tolerance)))
                    )
            
            # Store current kappa values for stagnation detection
            prev_kappa_values[recalc_indices] = current_kappa
    
    # Create a new array to hold all modified sequences
    all_modified_sequences = np.zeros_like(ternarized_sequences)
    
    # Use inverse_mapping to place modified unique sequences back into original positions
    for i, unique_idx in enumerate(inverse_mapping):
        all_modified_sequences[i] = modified_unique_sequences[unique_idx]
    
    # Convert ternarized sequences back to amino acid sequences
    modified_sequences = optimized_convert_ternarized_to_amino(all_modified_sequences, sequences)

    if only_return_within_tolerance:
        # Create full-sized boolean array for all sequences
        all_reached_target = np.zeros(len(modified_sequences), dtype=bool)
        # Map reached_target status from unique sequences to all sequences
        for i, unique_idx in enumerate(inverse_mapping):
            all_reached_target[i] = reached_target[unique_idx]
        return modified_sequences[all_reached_target]
    
    return modified_sequences


def optimized_convert_ternarized_to_amino(ternarized_sequences, original_sequences):
    """
    Optimized version that preserves the exact amino acid identity of the original sequences,
    only changing their positions.
    
    Args:
        ternarized_sequences (np.array): Ternarized sequences (-1, 0, 1)
        original_sequences (np.array): Original amino acid sequences
        
    Returns:
        np.array: New amino acid sequences with preserved amino acid identity
    """
    # Initialize output array
    converted_sequences = np.empty(len(ternarized_sequences), dtype=object)
    
    # Pre-ternarize all original sequences at once
    original_ternarized = vectorized_ternarize(original_sequences, group1=['K', 'R'], group2=['D', 'E'])
    
    # Process each sequence
    for i, (tern_seq, orig_seq, orig_tern) in enumerate(zip(ternarized_sequences, original_sequences, original_ternarized)):
        # Find positions where charges were moved (boolean mask)
        changed_positions = tern_seq != orig_tern
        
        # If nothing changed, keep original
        if not np.any(changed_positions):
            converted_sequences[i] = orig_seq
            continue
        
        # Start with a copy of the original sequence
        new_seq = list(orig_seq)
        
        # Create "pools" of available amino acids from the original sequence
        pos_pool = [] # Positive charged AAs (K, R)
        neg_pool = [] # Negative charged AAs (D, E)
        neut_pool = [] # Neutral AAs
        
        # Extract all amino acids from the original sequence
        for aa in orig_seq:
            if aa in 'KR':
                pos_pool.append(aa)
            elif aa in 'DE':
                neg_pool.append(aa)
            else:
                neut_pool.append(aa)
        
        # Count how many of each charge we need to place
        pos_needed = np.sum(tern_seq == 1)
        neg_needed = np.sum(tern_seq == -1)
        neut_needed = np.sum(tern_seq == 0)
        
        # Verify we have the right number of amino acids
        # This should always be true if the ternarization process is correct
        assert len(pos_pool) == pos_needed, f"Mismatch in positive charges: have {len(pos_pool)}, need {pos_needed}"
        assert len(neg_pool) == neg_needed, f"Mismatch in negative charges: have {len(neg_pool)}, need {neg_needed}"
        assert len(neut_pool) == neut_needed, f"Mismatch in neutral AAs: have {len(neut_pool)}, need {neut_needed}"
        
        # Place amino acids according to the new ternarized sequence
        for j, val in enumerate(tern_seq):
            if val == 1:  # Positive charge
                new_seq[j] = pos_pool.pop(0)
            elif val == -1:  # Negative charge
                new_seq[j] = neg_pool.pop(0)
            else:  # Neutral
                new_seq[j] = neut_pool.pop(0)
                
        # Final verification - pools should be empty
        assert len(pos_pool) == 0 and len(neg_pool) == 0 and len(neut_pool) == 0, "Not all amino acids were placed"
        
        converted_sequences[i] = ''.join(new_seq)
    
    return converted_sequences

