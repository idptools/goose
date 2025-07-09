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
    This function is vectorized for efficiency.
    
    Parameters
    ----------
    ternarized_sequences : list of np.ndarray
        List of ternarized sequences, where each sequence is a 2D numpy array.
    
    Returns
    -------
    list of np.ndarray
        Modified ternarized sequences with charges moved to decrease kappa.
    '''
    # Convert to numpy array for easier manipulation
    seqs_array = np.array(ternarized_sequences)
    n_sequences, seq_length = seqs_array.shape
    window_size = 5
    half_window = window_size // 2
    
    # Calculate the linear NCPR for each sequence
    linear_ncpr = calculate_linear_ncpr(ternarized_sequences, window_size=window_size)
    
    # Create a copy of the sequences to modify
    modified_sequences = seqs_array.copy()
    
    # Vectorized selection of random indices
    sorted_indices = np.argsort(linear_ncpr, axis=1)
    num_to_choose = max(1, seq_length // 10)
    
    # For each sequence, choose a random index from the first 10% of sorted indices
    first_indices_cols = np.random.randint(0, num_to_choose, size=n_sequences)
    random_first_indices = sorted_indices[np.arange(n_sequences), first_indices_cols]

    # For each sequence, choose a random index from the last 10% of sorted indices
    last_indices_cols = np.random.randint(0, num_to_choose, size=n_sequences)
    random_last_indices = sorted_indices[np.arange(n_sequences), -last_indices_cols - 1]

    # Process each sequence individually for the charge moving logic
    for seq_idx in range(n_sequences):
        seq = seqs_array[seq_idx]
        random_first_index = random_first_indices[seq_idx]
        random_last_index = random_last_indices[seq_idx]

        # Get charged positions
        charged_positions = np.where(np.abs(seq) == 1)[0]
        
        # Find a charged position in the first window
        first_window_start = max(0, random_first_index - half_window)
        first_window_end = min(seq_length, random_first_index + half_window + 1)
        first_window_positions = np.arange(first_window_start, first_window_end)
        first_window_charged_positions = np.intersect1d(charged_positions, first_window_positions)
        
        if len(first_window_charged_positions) == 0:
            continue
        
        random_first_pos = np.random.choice(first_window_charged_positions)
        
        # Find a random position in the last window
        last_window_start = max(0, random_last_index - half_window)
        last_window_end = min(seq_length, random_last_index + half_window + 1)
        last_window_positions = np.arange(last_window_start, last_window_end)
        
        if len(last_window_positions) == 0:
            continue

        random_last_position = np.random.choice(last_window_positions)
        
        # Move the charge
        if random_first_pos != random_last_position:
            displaced_value = modified_sequences[seq_idx, random_last_position]
            charge_value = modified_sequences[seq_idx, random_first_pos]
            modified_sequences[seq_idx, random_last_position] = charge_value
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
    
    # Vectorized swap for misplaced charges
    # Process sequences that have misplaced charges
    has_misplaced_pos = np.any(misplaced_positive_masks, axis=1)
    has_misplaced_neg = np.any(misplaced_negative_masks, axis=1)
    main_swap_mask = has_misplaced_pos & has_misplaced_neg

    if np.any(main_swap_mask):
        # Get indices of sequences to process
        swap_indices = np.where(main_swap_mask)[0]

        # Vectorized random choice of source positions
        source_pos_indices = np.argmax(misplaced_positive_masks[swap_indices], axis=1)
        source_neg_indices = np.argmax(misplaced_negative_masks[swap_indices], axis=1)

        # Vectorized random choice of target positions
        available_pos_targets = positive_region_mask & ~positive_masks
        available_neg_targets = negative_region_mask & ~negative_masks
        
        target_pos_indices = np.argmax(available_pos_targets[swap_indices], axis=1)
        target_neg_indices = np.argmax(available_neg_targets[swap_indices], axis=1)

        # Perform swaps
        # Swap positive charges
        temp_val = modified_sequences[swap_indices, target_pos_indices]
        modified_sequences[swap_indices, target_pos_indices] = modified_sequences[swap_indices, source_pos_indices]
        modified_sequences[swap_indices, source_pos_indices] = temp_val

        # Swap negative charges
        temp_val = modified_sequences[swap_indices, target_neg_indices]
        modified_sequences[swap_indices, target_neg_indices] = modified_sequences[swap_indices, source_neg_indices]
        modified_sequences[swap_indices, source_neg_indices] = temp_val

    # Fallback for sequences that didn't have misplaced charges
    fallback_mask = ~main_swap_mask
    if np.any(fallback_mask):
        fallback_indices = np.where(fallback_mask)[0]
        
        for seq_idx in fallback_indices:
            neutral_positions = np.where(neutral_masks[seq_idx])[0]
            all_positive_pos = np.where(positive_masks[seq_idx])[0]
            all_negative_pos = np.where(negative_masks[seq_idx])[0]

            if len(neutral_positions) >= 2 and len(all_positive_pos) > 0 and len(all_negative_pos) > 0:
                first_half_more_positive = mean_first_half[seq_idx] > mean_second_half[seq_idx]
                
                if first_half_more_positive:
                    furthest_pos = np.max(all_positive_pos)
                    furthest_neg = np.min(all_negative_pos)
                    pos_target_region = neutral_positions[:max(1, len(neutral_positions)//2)]
                    neg_target_region = neutral_positions[len(neutral_positions)//2:]
                else:
                    furthest_pos = np.min(all_positive_pos)
                    furthest_neg = np.max(all_negative_pos)
                    pos_target_region = neutral_positions[len(neutral_positions)//2:]
                    neg_target_region = neutral_positions[:max(1, len(neutral_positions)//2)]

                if len(pos_target_region) > 0:
                    pos_dest = np.random.choice(pos_target_region)
                    if furthest_pos != pos_dest:
                        modified_sequences[seq_idx, pos_dest], modified_sequences[seq_idx, furthest_pos] = \
                            modified_sequences[seq_idx, furthest_pos], modified_sequences[seq_idx, pos_dest]
                        neg_target_region = neg_target_region[neg_target_region != pos_dest]
                
                if len(neg_target_region) > 0:
                    neg_dest = np.random.choice(neg_target_region)
                    if furthest_neg != neg_dest:
                        modified_sequences[seq_idx, neg_dest], modified_sequences[seq_idx, furthest_neg] = \
                            modified_sequences[seq_idx, furthest_neg], modified_sequences[seq_idx, neg_dest]
    
    # Convert back to list of arrays
    return [modified_sequences[i] for i in range(len(modified_sequences))]

def increase_kappa_aggressive(ternarized_sequences, window_size=5):
    """
    Alternative aggressive increase_kappa function that systematically moves charges
    toward sequence ends to maximize asymmetry. This function is designed to avoid
    getting stuck in local optima that plague the standard increase_kappa function.
    This is a fully vectorized implementation for maximum speed.

    Parameters
    ----------
    ternarized_sequences : list of np.ndarray
        List of ternarized sequences, where each sequence is a 2D numpy array.
    window_size : int, optional
        The size of the sliding window to use for calculating NCPR. Default is 5.
        (Note: window_size is not used in this aggressive implementation but is kept for API consistency).

    Returns
    -------
    list of np.ndarray
        List of modified sequences with increased kappa values.
    """
    # Convert to numpy array for easier manipulation
    seqs_array = np.array(ternarized_sequences)
    n_sequences, seq_length = seqs_array.shape

    # Create modified sequences as a copy
    modified_sequences = seqs_array.copy()

    # Masks for charge positions
    positive_masks = (seqs_array == 1)
    negative_masks = (seqs_array == -1)

    # Find sequences that have at least one of each charge
    valid_indices = np.where(np.any(positive_masks, axis=1) & np.any(negative_masks, axis=1))[0]

    if len(valid_indices) == 0:
        return [modified_sequences[i] for i in range(len(modified_sequences))]

    # Process only the valid sequences
    valid_seqs = seqs_array[valid_indices]
    valid_pos_masks = positive_masks[valid_indices]
    valid_neg_masks = negative_masks[valid_indices]

    # Calculate centers of mass for positive charges
    positions = np.arange(seq_length)
    pos_centers = np.sum(valid_pos_masks * positions, axis=1) / (np.sum(valid_pos_masks, axis=1) + 1e-9)

    # Determine target direction
    pos_target_is_right = pos_centers > (seq_length - 1) / 2

    # Find source charges
    leftmost_pos = np.argmax(valid_pos_masks, axis=1)
    rightmost_pos = seq_length - 1 - np.argmax(np.fliplr(valid_pos_masks), axis=1)
    leftmost_neg = np.argmax(valid_neg_masks, axis=1)
    rightmost_neg = seq_length - 1 - np.argmax(np.fliplr(valid_neg_masks), axis=1)

    pos_source = np.where(pos_target_is_right, leftmost_pos, rightmost_pos)
    neg_source = np.where(pos_target_is_right, rightmost_neg, leftmost_neg)

    # --- Fully Vectorized Swaps ---
    half_len = seq_length // 2
    
    # Create a temporary array to hold the results of the first (positive) swap
    temp_sequences = modified_sequences.copy()

    # Positive Swaps
    for i, seq_idx in enumerate(valid_indices):
        source_pos = pos_source[i]
        if pos_target_is_right[i]:
            target_range = np.arange(half_len, seq_length)
            available = target_range[(seqs_array[seq_idx, target_range] != 1) & (target_range > source_pos)]
            if len(available) > 0:
                target_pos = np.random.choice(available)
                temp_sequences[seq_idx, target_pos], temp_sequences[seq_idx, source_pos] = \
                    temp_sequences[seq_idx, source_pos], temp_sequences[seq_idx, target_pos]
        else:
            target_range = np.arange(half_len)
            available = target_range[(seqs_array[seq_idx, target_range] != 1) & (target_range < source_pos)]
            if len(available) > 0:
                target_pos = np.random.choice(available)
                temp_sequences[seq_idx, target_pos], temp_sequences[seq_idx, source_pos] = \
                    temp_sequences[seq_idx, source_pos], temp_sequences[seq_idx, target_pos]

    # Negative Swaps (operates on the result of the positive swaps)
    for i, seq_idx in enumerate(valid_indices):
        source_pos = neg_source[i]
        if not pos_target_is_right[i]: # Negative target is right
            target_range = np.arange(half_len, seq_length)
            available = target_range[(temp_sequences[seq_idx, target_range] != -1) & (target_range > source_pos)]
            if len(available) > 0:
                target_pos = np.random.choice(available)
                temp_sequences[seq_idx, target_pos], temp_sequences[seq_idx, source_pos] = \
                    temp_sequences[seq_idx, source_pos], temp_sequences[seq_idx, target_pos]
        else: # Negative target is left
            target_range = np.arange(half_len)
            available = target_range[(temp_sequences[seq_idx, target_range] != -1) & (target_range < source_pos)]
            if len(available) > 0:
                target_pos = np.random.choice(available)
                temp_sequences[seq_idx, target_pos], temp_sequences[seq_idx, source_pos] = \
                    temp_sequences[seq_idx, source_pos], temp_sequences[seq_idx, target_pos]
    
    # Assign the final modified sequences
    modified_sequences = temp_sequences

    # Convert back to list of arrays
    return [modified_sequences[i] for i in range(len(modified_sequences))]
