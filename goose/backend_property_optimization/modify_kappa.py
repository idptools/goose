import numpy as np
from numpy.lib.stride_tricks import sliding_window_view


def calculate_linear_ncpr(ternarized_seqs, window_size):
    '''
    Pure numpy implementation using stride_tricks for maximum efficiency.
    
    Parameters
    ----------
    ternarized_seqs : list of np.ndarray
        List of ternarized sequences, where each sequence is a 2D numpy array.
    window_size : int
        The size of the sliding window to use for calculating NCPR.
    
    Returns
    -------
    np.ndarray
        A 2D numpy array where each row corresponds to a sequence and each column corresponds to a position in the sequence.
        The values are the linear NCPR calculated for each position in the sequence.
    '''
    # Convert to numpy array
    seqs_array = np.array(ternarized_seqs, dtype=float)
    
    # Pad sequences to handle edge effects
    half_window = window_size // 2
    padded_seqs = np.pad(seqs_array, ((0, 0), (half_window, half_window)), mode='edge')
    
    # Create sliding windows for all sequences at once
    windows = sliding_window_view(padded_seqs, window_shape=window_size, axis=1)
    
    # Calculate mean across the window dimension
    linear_ncpr = np.mean(windows, axis=-1)
    
    return linear_ncpr

def decrease_kappa(ternarized_sequences):
    '''
    Function to decrease the kappa value of a sequence by applying a sliding window approach.
    Purposefully isn't maximimally efficient. 
    Parameters
    ----------
    ternarized_sequences : list of np.ndarray
        List of ternarized sequences, where each sequence is a 2D numpy array.
    Returns
    -------
    list of np.ndarray
        Modified ternarized sequences with charges added to decrease kappa.
    '''
    # Convert to numpy array for easier manipulation
    seqs_array = np.array(ternarized_sequences)
    window_size = 5
    half_window = window_size // 2
    
    # Calculate the linear NCPR for each sequence
    linear_ncpr = calculate_linear_ncpr(ternarized_sequences, window_size=window_size)
    
    # Create a copy of the sequences to modify
    modified_sequences = seqs_array.copy()
    
    # Process each sequence individually
    for seq_idx in range(len(ternarized_sequences)):
        seq = seqs_array[seq_idx]
        seq_linear_ncpr = linear_ncpr[seq_idx]

        # sort the seq_linear_ncpr by absolute value
        sorted_indices = np.argsort(seq_linear_ncpr)
        # use np.random to get a random index from the first 10% of the sorted indices
        num_to_choose = max(1, len(sorted_indices) // 10)  # at least one position
        chosen_indices = sorted_indices[:num_to_choose]
        random_first_index = np.random.choice(chosen_indices)
        # now get from the last 10%
        last_indices = sorted_indices[-num_to_choose:]
        random_last_index = np.random.choice(last_indices)
        
        # get all positions that are -1 or 1 in the sequence
        charged_positions = np.where(np.abs(seq) == 1)[0]
        # get a charged position from the window in random_first_index
        first_window_start = max(0, random_first_index - half_window)
        first_window_end = min(len(seq), random_first_index + half_window + 1)
        first_window_positions = np.arange(first_window_start, first_window_end)
        first_window_charged_positions = np.intersect1d(charged_positions, first_window_positions)
        # choose a random positon from first_window_charged_positions
        if len(first_window_charged_positions) > 0:
            random_first_pos = np.random.choice(first_window_charged_positions)
        else:
            # If no charged positions in window, skip this sequence
            continue

        # get a random position in the random_last_window
        last_window_start = max(0, random_last_index - half_window)
        last_window_end = min(len(seq), random_last_index + half_window + 1)
        last_window_positions = np.arange(last_window_start, last_window_end)
        random_last_position = np.random.choice(last_window_positions)
        
        # Move the charge from random_first_pos to random_last_position
        if random_first_pos != random_last_position:
            # Store what's currently at the target position
            displaced_value = modified_sequences[seq_idx, random_last_position]
            # Move the charge to the target position
            charge_value = modified_sequences[seq_idx, random_first_pos]
            modified_sequences[seq_idx, random_last_position] = charge_value
            # Put the displaced value where the charge was
            modified_sequences[seq_idx, random_first_pos] = displaced_value
    
    # Convert back to list of arrays
    return [modified_sequences[i] for i in range(len(modified_sequences))]



def increase_kappa(ternarized_sequences, window_size=5):
    '''
    Moves a negative residue to the more negative half of the sequence and a positive residue to the more positive half of the sequence.
    Fully vectorized implementation for improved performance - removes the for loop entirely.
    
    Parameters
    ----------
    ternarized_sequences : list of np.ndarray
        List of ternarized sequences, where each sequence is a 2D numpy array.
    window_size : int, optional
        The size of the sliding window to use for calculating NCPR. Default is 5
    Returns
    -------
    list of np.ndarray
        List of modified sequences with increased kappa values.
    '''
    # Calculate the linear NCPR for each sequence
    linear_ncpr = calculate_linear_ncpr(ternarized_sequences, window_size=window_size)
    # get the mean for the first half of the sequence 
    mean_first_half = np.mean(linear_ncpr[:, :linear_ncpr.shape[1] // 2], axis=1)
    # get the mean for the second half of the sequence
    mean_second_half = np.mean(linear_ncpr[:, linear_ncpr.shape[1] // 2:], axis=1)
    
    quarter_length = linear_ncpr.shape[1] // 4
    if np.mean(mean_first_half) > np.mean(mean_second_half):
        # set positive positions to the indices at the first 1/4 of the sequence
        positive_positions = np.arange(quarter_length)
        negative_positions = np.arange(quarter_length*3, linear_ncpr.shape[1])
    else:
        # set positive positions to the indices at the last 1/4 of the sequence
        negative_positions = np.arange(quarter_length)
        positive_positions = np.arange(quarter_length*3, linear_ncpr.shape[1])
    
    # Convert to numpy array for easier manipulation
    seqs_array = np.array(ternarized_sequences)
    n_sequences, seq_length = seqs_array.shape
    
    # Create modified sequences as a copy
    modified_sequences = seqs_array.copy()
    
    # Vectorized analysis: find all charge positions across all sequences
    positive_masks = (seqs_array == 1)  # Shape: (n_sequences, seq_length)
    negative_masks = (seqs_array == -1)  # Shape: (n_sequences, seq_length)
    neutral_masks = (seqs_array == 0)    # Shape: (n_sequences, seq_length)
    
    # Create target region masks for all sequences
    positive_region_mask = np.zeros((n_sequences, seq_length), dtype=bool)
    negative_region_mask = np.zeros((n_sequences, seq_length), dtype=bool)
    
    # Apply the same target regions to all sequences
    valid_positive_positions = positive_positions[positive_positions < seq_length]
    valid_negative_positions = negative_positions[negative_positions < seq_length]
    
    positive_region_mask[:, valid_positive_positions] = True
    negative_region_mask[:, valid_negative_positions] = True
    
    # Find misplaced charges vectorized across all sequences
    misplaced_positive_masks = positive_masks & ~positive_region_mask  # Positive charges not in positive region
    misplaced_negative_masks = negative_masks & ~negative_region_mask  # Negative charges not in negative region
    
    # For each sequence, try to move one misplaced charge of each type
    # We'll use advanced indexing and random selection
    for seq_idx in range(n_sequences):
        # Get misplaced positions for this sequence
        misplaced_pos_positions = np.where(misplaced_positive_masks[seq_idx])[0]
        misplaced_neg_positions = np.where(misplaced_negative_masks[seq_idx])[0]
        
        # Skip if no misplaced charges
        if len(misplaced_pos_positions) == 0 or len(misplaced_neg_positions) == 0:
            # Fallback: try using neutral positions if available
            neutral_positions = np.where(neutral_masks[seq_idx])[0]
            all_positive_pos = np.where(positive_masks[seq_idx])[0]
            all_negative_pos = np.where(negative_masks[seq_idx])[0]
            
            if len(neutral_positions) >= 2 and len(all_positive_pos) > 0 and len(all_negative_pos) > 0:
                # Determine movement direction based on sequence asymmetry
                first_half_more_positive = mean_first_half[seq_idx] > mean_second_half[seq_idx]
                
                if first_half_more_positive:
                    # Move rightmost positive to left, leftmost negative to right
                    furthest_pos = np.max(all_positive_pos)
                    furthest_neg = np.min(all_negative_pos)
                    pos_target_region = neutral_positions[:max(1, len(neutral_positions)//2)]
                    neg_target_region = neutral_positions[len(neutral_positions)//2:]
                else:
                    # Move leftmost positive to right, rightmost negative to left  
                    furthest_pos = np.min(all_positive_pos)
                    furthest_neg = np.max(all_negative_pos)
                    pos_target_region = neutral_positions[len(neutral_positions)//2:]
                    neg_target_region = neutral_positions[:max(1, len(neutral_positions)//2)]
                
                # Execute moves if targets available
                if len(pos_target_region) > 0:
                    pos_dest = np.random.choice(pos_target_region)
                    if furthest_pos != pos_dest:
                        modified_sequences[seq_idx, pos_dest], modified_sequences[seq_idx, furthest_pos] = \
                            modified_sequences[seq_idx, furthest_pos], modified_sequences[seq_idx, pos_dest]
                        # Remove used position from negative targets
                        neg_target_region = neg_target_region[neg_target_region != pos_dest]
                
                if len(neg_target_region) > 0:
                    neg_dest = np.random.choice(neg_target_region)
                    if furthest_neg != neg_dest:
                        modified_sequences[seq_idx, neg_dest], modified_sequences[seq_idx, furthest_neg] = \
                            modified_sequences[seq_idx, furthest_neg], modified_sequences[seq_idx, neg_dest]
            continue
        
        # Find available target positions (in target regions but not already occupied by desired charge)
        available_pos_targets = np.where(positive_region_mask[seq_idx] & ~positive_masks[seq_idx])[0]
        available_neg_targets = np.where(negative_region_mask[seq_idx] & ~negative_masks[seq_idx])[0]
        
        if len(available_pos_targets) > 0 and len(available_neg_targets) > 0:
            # Select random misplaced charges to move
            source_pos = np.random.choice(misplaced_pos_positions)
            source_neg = np.random.choice(misplaced_neg_positions)
            
            # Select random target positions
            target_pos = np.random.choice(available_pos_targets)
            
            # Ensure no conflict between targets
            available_neg_targets = available_neg_targets[available_neg_targets != target_pos]
            
            if len(available_neg_targets) > 0:
                target_neg = np.random.choice(available_neg_targets)
                
                # Perform the swaps
                if source_pos != target_pos:
                    modified_sequences[seq_idx, target_pos], modified_sequences[seq_idx, source_pos] = \
                        modified_sequences[seq_idx, source_pos], modified_sequences[seq_idx, target_pos]
                
                if source_neg != target_neg:
                    modified_sequences[seq_idx, target_neg], modified_sequences[seq_idx, source_neg] = \
                        modified_sequences[seq_idx, source_neg], modified_sequences[seq_idx, target_neg]
    
    # Convert back to list of arrays
    return [modified_sequences[i] for i in range(len(modified_sequences))] 


def increase_kappa_aggressive(ternarized_sequences, window_size=5):
    '''
    Alternative aggressive increase_kappa function that systematically moves charges 
    toward sequence ends to maximize asymmetry. This function is designed to avoid 
    getting stuck in local optima that plague the standard increase_kappa function.
    
    Parameters
    ----------
    ternarized_sequences : list of np.ndarray
        List of ternarized sequences, where each sequence is a 2D numpy array.
    window_size : int, optional
        The size of the sliding window to use for calculating NCPR. Default is 5
    Returns
    -------
    list of np.ndarray
        List of modified sequences with increased kappa values.
    '''
    # Convert to numpy array for easier manipulation
    seqs_array = np.array(ternarized_sequences)
    n_sequences, seq_length = seqs_array.shape
    
    # Create modified sequences as a copy
    modified_sequences = seqs_array.copy()
    
    # Strategy: Direct clustering approach - move charges to extreme ends
    for seq_idx in range(n_sequences):
        seq = seqs_array[seq_idx]
        
        # Find all charge positions
        positive_positions = np.where(seq == 1)[0]
        negative_positions = np.where(seq == -1)[0]
        
        if len(positive_positions) == 0 or len(negative_positions) == 0:
            continue
            
        # Calculate current center of mass for each charge type
        pos_center = np.mean(positive_positions) if len(positive_positions) > 0 else seq_length / 2
        neg_center = np.mean(negative_positions) if len(negative_positions) > 0 else seq_length / 2
        
        # Determine which end each charge type should move towards
        if pos_center > seq_length / 2:
            # Positive charges are already towards the right, push them further right
            pos_target_region = "right"
            neg_target_region = "left"
        else:
            # Positive charges are towards the left, push them further left
            pos_target_region = "left" 
            neg_target_region = "right"
        
        # Move charges toward the designated ends to maximize asymmetry
        if pos_target_region == "right":
            # Move leftmost positive charge towards right
            leftmost_pos = np.min(positive_positions)
            # Find available positions in right half
            right_half_start = seq_length // 2
            available_right = np.where((seq[right_half_start:] != 1) & 
                                     (np.arange(right_half_start, seq_length) > leftmost_pos))[0]
            if len(available_right) > 0:
                target_pos = right_half_start + np.random.choice(available_right)
                # Perform atomic swap to maintain charge conservation
                modified_sequences[seq_idx, target_pos], modified_sequences[seq_idx, leftmost_pos] = \
                    modified_sequences[seq_idx, leftmost_pos], modified_sequences[seq_idx, target_pos]
        else:
            # Move rightmost positive charge towards left
            rightmost_pos = np.max(positive_positions)
            # Find available positions in left half
            left_half_end = seq_length // 2
            available_left = np.where((seq[:left_half_end] != 1) & 
                                    (np.arange(left_half_end) < rightmost_pos))[0]
            if len(available_left) > 0:
                target_pos = np.random.choice(available_left)
                # Perform atomic swap to maintain charge conservation
                modified_sequences[seq_idx, target_pos], modified_sequences[seq_idx, rightmost_pos] = \
                    modified_sequences[seq_idx, rightmost_pos], modified_sequences[seq_idx, target_pos]
        
        # Do the same for negative charges
        if neg_target_region == "right":
            # Move leftmost negative charge towards right
            leftmost_neg = np.min(negative_positions)
            right_half_start = seq_length // 2
            available_right = np.where((seq[right_half_start:] != -1) & 
                                     (np.arange(right_half_start, seq_length) > leftmost_neg))[0]
            if len(available_right) > 0:
                target_pos = right_half_start + np.random.choice(available_right)
                # Perform atomic swap to maintain charge conservation
                modified_sequences[seq_idx, target_pos], modified_sequences[seq_idx, leftmost_neg] = \
                    modified_sequences[seq_idx, leftmost_neg], modified_sequences[seq_idx, target_pos]
        else:
            # Move rightmost negative charge towards left  
            rightmost_neg = np.max(negative_positions)
            left_half_end = seq_length // 2
            available_left = np.where((seq[:left_half_end] != -1) & 
                                    (np.arange(left_half_end) < rightmost_neg))[0]
            if len(available_left) > 0:
                target_pos = np.random.choice(available_left)
                # Perform atomic swap to maintain charge conservation
                modified_sequences[seq_idx, target_pos], modified_sequences[seq_idx, rightmost_neg] = \
                    modified_sequences[seq_idx, rightmost_neg], modified_sequences[seq_idx, target_pos]
    
    # Convert back to list of arrays
    return [modified_sequences[i] for i in range(len(modified_sequences))] 


def test_charge_conservation(ternarized_sequences, function_to_test):
    '''
    Test function to verify that charge counts are preserved.
    
    Parameters
    ----------
    ternarized_sequences : list of np.ndarray
        List of ternarized sequences
    function_to_test : function
        The function to test (e.g., increase_kappa)
    
    Returns
    -------
    bool
        True if charge conservation is maintained, False otherwise
    '''
    # Count charges before
    original_counts = []
    for seq in ternarized_sequences:
        pos_count = np.sum(seq == 1)
        neg_count = np.sum(seq == -1)
        zero_count = np.sum(seq == 0)
        original_counts.append((pos_count, neg_count, zero_count))
    
    # Apply function
    modified_sequences = function_to_test(ternarized_sequences)
    
    # Count charges after
    modified_counts = []
    for seq in modified_sequences:
        pos_count = np.sum(seq == 1)
        neg_count = np.sum(seq == -1)
        zero_count = np.sum(seq == 0)
        modified_counts.append((pos_count, neg_count, zero_count))
    
    # Check if counts match
    for i, (orig, mod) in enumerate(zip(original_counts, modified_counts)):
        if orig != mod:
            print(f"Sequence {i}: Original {orig} -> Modified {mod}")
            return False
    
    print("All sequences maintain charge conservation!")
    return True
