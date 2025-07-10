'''
Functionality for optimizing kappa towards some target value. 
If a sequences kappa is too low, that means that the charged residues in the
sequence need to be more asymmetrically distributed across the sequence. 
If a sequences kappa is too high, that means that the charged residues in the
sequence need to be more symmetrically distributed across the sequence.
The optimization is done by moving the charged residues in the sequence.
'''

import random
import numpy as np
from typing import List
from goose import parameters
from goose.backend_property_calculation.calculate_kappa import charge_matrix_to_ternary, kappa
from goose.backend_property_calculation.calculate_properties_batch import sequences_to_matrices, matrices_to_sequences
from goose.backend_property_optimization.helper_functions import shuffle_sequence_return_matrix
from goose.backend_property_optimization.minimize_kappa import efficient_kappa_minimization
from goose.backend_property_optimization.modify_kappa import increase_kappa, decrease_kappa, increase_kappa_aggressive

def optimize_kappa(sequence, target_kappa, 
                              num_copies: int = 10,
                              max_change_iterations: int = 10,
                              kappa_tolerance=parameters.MAXIMUM_KAPPA_ERROR, 
                              num_iterations=20000, 
                              return_when_num_hit=None,
                              only_return_within_tolerance=True,
                              convert_input_seq_to_matrix=False,
                              inputting_matrix=True,
                              avoid_shuffle=False) -> List[str]:
    '''
    Optimize kappa values of sequences by either increasing or decreasing charge asymmetry
    to reach a target kappa value. Optimized for performance.
    
    Args:
        sequence (str): Single amino acid sequence to optimize
        target_kappa (float): Target kappa value to reach
        num_copies (int): Number of copies of the sequence to create for optimization
        max_change_iterations (int): Maximum iterations for each sequence modification
        tolerance (float): Acceptable difference from target kappa
        num_iterations (int): Maximum number of overall iterations to perform
        early_stopping (bool): Stop early if target is reached
        return_when_num_hit : int
            If not None, return when this number of sequences are within tolerance
        only_return_within_tolerance (bool): If True, only return sequences within tolerance.
            if False, return all sequences, even if they are not within tolerance.
        convert_input_seq_to_matrix (bool): If True, convert input sequence to matrix format before processing.
        inputting_matrix (bool): If True, inputting something ready to go. 
        avoid_shuffle (bool): If True, avoid shuffling sequences.
    
    Returns:
        list: Sequences with kappa values optimized to be close to the target
    '''

    if inputting_matrix==False:
        if convert_input_seq_to_matrix:
            # Ensure we have a single sequence as input
            if isinstance(sequence, list):
                if len(sequence) == 1:
                    sequence = sequence[0]
                else:
                    raise ValueError("Function expects a single sequence, but got a list of multiple sequences")
            # Convert input to numpy array for consistency
            sequence_matrix = sequences_to_matrices([sequence])
        else:
            sequence_matrix=sequence.copy()
            sequence = matrices_to_sequences(sequence)[0]  # Ensure we have a single sequence
       
        # make sure seq matrix has correct shape
        if sequence_matrix.shape[0] != 1:
            raise ValueError("Input sequence matrix must have shape (1, L) where L is the length of the sequence")

        # make copies of the sequence and shuffle all but the original
        if avoid_shuffle==False:
            sequence_matrix = shuffle_sequence_return_matrix(sequence_matrix, num_shuffles=num_copies)
        else:
            # make a matrix with num_copies of the sequence
            sequence_matrix = np.tile(sequence_matrix, (num_copies, 1))
    else:
        sequence_matrix = sequence.copy()
        sequence = matrices_to_sequences(sequence)

    
    if return_when_num_hit != None:
        # force only_return_within_tolerance to be true
        only_return_within_tolerance=True
        if return_when_num_hit > sequence_matrix.shape[0]:
            return_when_num_hit = sequence_matrix.shape[0]
    else:
        return_when_num_hit = sequence_matrix.shape[0]

    # Ternarize the sequences (only once at the start)
    ternarized_sequences = charge_matrix_to_ternary(sequence_matrix)

    # Calculate initial kappa values for all sequences
    kappa_values = kappa(ternarized_sequences, is_ternarized=True)

    prev_kappa_values = kappa_values.copy()  # Store initial kappa values for stagnation detection

    # Identify which sequences need to have kappa increased or decreased
    needs_increase = kappa_values < target_kappa - kappa_tolerance
    needs_decrease = kappa_values > target_kappa + kappa_tolerance

    # Track which sequences have reached their target
    reached_target = ~(needs_increase | needs_decrease)
    
    # Initialize modified sequences array (avoid redundant copies)
    modified_sequences = ternarized_sequences.copy()

    # set window size
    window_size=5
    # set num not imporoved to 0
    num_not_improved=0

    # Process in iterations
    for num_it in range(num_iterations):
        
        # If all sequences reached target, we're done
        if np.all(reached_target):
            break
        
        # if we are stopping at a specific number of successes...
        if np.sum(reached_target) >= return_when_num_hit:
            break
        
        # Process sequences that need kappa increased
        increase_indices = np.where(needs_increase)[0]

        # Process sequences that need kappa increased
        if len(increase_indices) > 0:
            increase_seqs = modified_sequences[increase_indices]
            for i, idx in enumerate(increase_indices):
                # Apply increase_kappa multiple times based on max_change_iterations
                current_seq = increase_seqs[i:i+1].copy()
                iters = random.randint(1, max_change_iterations)
                for _ in range(iters):
                    # Use aggressive strategy when stagnation is detected
                    if num_not_improved > 1:
                        current_seq = increase_kappa_aggressive(current_seq, window_size=window_size)
                    else:
                        current_seq = increase_kappa(current_seq, window_size=window_size)
                modified_sequences[idx] = current_seq[-1]
        # Process sequences that need kappa decreased
        decrease_indices = np.where(needs_decrease)[0]
        if len(decrease_indices) > 0:
            decrease_seqs = modified_sequences[decrease_indices]
            
            for i, idx in enumerate(decrease_indices):
                iters = random.randint(1, max_change_iterations)
                current_seq = decrease_seqs[i:i+1].copy()
                if target_kappa > 0.15:
                    for _ in range(iters):
                        current_seq = decrease_kappa(current_seq)
                    modified_sequences[idx] = current_seq[-1]
                else:
                    result = efficient_kappa_minimization(current_seq, iters)
                    modified_sequences[idx] = result[0]
                
        # Recalculate kappa for sequences that haven't reached target
        recalc_indices = np.where(~reached_target)[0]
        if len(recalc_indices) > 0:
            current_kappa = kappa(modified_sequences[recalc_indices], is_ternarized=True)
            # Check for stagnation with small tolerance to avoid floating point precision issues
            if np.all(np.abs(current_kappa - prev_kappa_values[recalc_indices]) < 1e-6):
                num_not_improved+=1
            else:
                num_not_improved=0
            # Update which sequences need modifications
            for i, idx in enumerate(recalc_indices):
                kappa_val = current_kappa[i]
                kappa_values[idx] = kappa_val
                
                # Update status of this sequence
                old_status = (~reached_target[idx], needs_increase[idx], needs_decrease[idx])

                reached_target[idx] = abs(kappa_val - target_kappa) <= kappa_tolerance
                needs_increase[idx] = kappa_val < target_kappa - kappa_tolerance
                needs_decrease[idx] = kappa_val > target_kappa + kappa_tolerance

                new_status = (~reached_target[idx], needs_increase[idx], needs_decrease[idx])
                #print(f"  Sequence {idx}: kappa {kappa_val:.6f}, reached_target: {reached_target[idx]}")
                    
            # Store current kappa values for stagnation detection
            prev_kappa_values[recalc_indices] = current_kappa
            #print(f"After iteration {num_it}: {np.sum(reached_target)} sequences reached target")

            if num_not_improved > 2:
                if window_size==5:
                    window_size=6
                else:
                    window_size=5
            if num_not_improved > 3:
                window_size=5
                if inputting_matrix==False:
                    # should have single sequence, just make shuffles
                    if avoid_shuffle==False:
                        modified_sequences[recalc_indices] = charge_matrix_to_ternary(
                            shuffle_sequence_return_matrix(
                            sequences_to_matrices([sequence]), num_shuffles=len(recalc_indices)))
                    else:
                        # set the sequences unable to improve to the original sequence
                        modified_sequences[recalc_indices] = charge_matrix_to_ternary(
                            sequences_to_matrices([sequence]))[0]
                        
                else:
                    # we need to get each sequence that is in recalc_indices and shuffle it. 
                    for n, idx in enumerate(recalc_indices):
                        cur_seq = sequence_matrix[idx:idx+1].copy()
                        modified_sequences[idx] = charge_matrix_to_ternary(
                            shuffle_sequence_return_matrix(cur_seq, num_shuffles=1))[0]
                # reset num_not_improved
                num_not_improved=0
                
    # Convert ternarized sequences back to amino acid sequences
    if inputting_matrix==False:
        modified_amino_sequences = optimized_convert_ternarized_to_amino(modified_sequences, sequence)
    else:
        modified_amino_sequences=[]
        for n, s in enumerate(sequence):
            cur_ternarized_seq = np.array([modified_sequences[n]])
            modified_amino_sequences.append(optimized_convert_ternarized_to_amino(cur_ternarized_seq, s)[0])


    if only_return_within_tolerance:
        # sort sequences by closest to target kappa
        sorted_indices = np.argsort(np.abs(kappa_values - target_kappa))
        # return those within tolerance in sorted order
        modified_amino_sequences = [modified_amino_sequences[i] for i in sorted_indices if abs(kappa_values[i] - target_kappa) <= kappa_tolerance]
        #print(modified_amino_sequences)
        if modified_amino_sequences == []:
            return None
    return modified_amino_sequences


def optimized_convert_ternarized_to_amino(ternarized_sequences, original_sequence):
    """
    Convert ternarized sequences back to amino acid sequences by mapping:
    - Negative charges (-1) back to their original negative amino acids
    - Positive charges (+1) back to their original positive amino acids
    - Neutral (0) back to their original neutral amino acids
    
    Args:
        ternarized_sequences (np.array): NxL array of ternarized sequences with values -1, 0, 1
        original_sequence (str): Single original amino acid sequence
        
    Returns:
        list: List of amino acid sequences with charges rearranged according to ternarized pattern
        
    Raises:
        ValueError: If the number of charges in ternarized sequence doesn't match original sequence
    """
    result_sequences = []
    
    # Get original sequence as list for easier manipulation
    original_chars = list(original_sequence)

    # Separate amino acids by charge type from original sequence
    positive_aas = [aa for aa in original_chars if aa in ['K','R']]
    negative_aas = [aa for aa in original_chars if aa in ['D','E']]
    neutral_aas = [aa for aa in original_chars if aa not in ['K','R','D','E']]
    
    # Count expected numbers for validation
    expected_pos = len(positive_aas)
    expected_neg = len(negative_aas)
    expected_neut = len(neutral_aas)
    
    for seq_idx, ternary_seq in enumerate(ternarized_sequences):
        # Validate that the ternarized sequence has the correct composition
        actual_pos = np.sum(ternary_seq == 1)
        actual_neg = np.sum(ternary_seq == -1)
        actual_neut = np.sum(ternary_seq == 0)
        
        if actual_pos != expected_pos:
            raise ValueError(f"Sequence {seq_idx}: Expected {expected_pos} positive charges but found {actual_pos}")
        if actual_neg != expected_neg:
            raise ValueError(f"Sequence {seq_idx}: Expected {expected_neg} negative charges but found {actual_neg}")
        if actual_neut != expected_neut:
            raise ValueError(f"Sequence {seq_idx}: Expected {expected_neut} neutral residues but found {actual_neut}")
        
        # Create copies for each sequence to maintain composition
        pos_pool = positive_aas.copy()
        neg_pool = negative_aas.copy()
        neut_pool = neutral_aas.copy()
        
        # Initialize result sequence
        result_seq = [''] * len(ternary_seq)
        
        # Counters for each amino acid type
        pos_counter = 0
        neg_counter = 0
        neut_counter = 0
        
        # Map ternarized values back to amino acids
        for j, ternary_val in enumerate(ternary_seq):
            if ternary_val == 1:  # Positive charge
                if pos_counter >= len(pos_pool):
                    raise ValueError(f"Sequence {seq_idx}: Ran out of positive amino acids at position {j}")
                result_seq[j] = pos_pool[pos_counter]
                pos_counter += 1
            elif ternary_val == -1:  # Negative charge
                if neg_counter >= len(neg_pool):
                    raise ValueError(f"Sequence {seq_idx}: Ran out of negative amino acids at position {j}")
                result_seq[j] = neg_pool[neg_counter]
                neg_counter += 1
            else:  # Neutral (0)
                if neut_counter >= len(neut_pool):
                    raise ValueError(f"Sequence {seq_idx}: Ran out of neutral amino acids at position {j}")
                result_seq[j] = neut_pool[neut_counter]
                neut_counter += 1
        
        # Final validation: ensure we used all amino acids
        if pos_counter != len(pos_pool):
            raise ValueError(f"Sequence {seq_idx}: Used {pos_counter} positive amino acids but had {len(pos_pool)}")
        if neg_counter != len(neg_pool):
            raise ValueError(f"Sequence {seq_idx}: Used {neg_counter} negative amino acids but had {len(neg_pool)}")
        if neut_counter != len(neut_pool):
            raise ValueError(f"Sequence {seq_idx}: Used {neut_counter} neutral amino acids but had {len(neut_pool)}")
        
        result_sequences.append(''.join(result_seq))
    
    return result_sequences
