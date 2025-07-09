'''
Code for optimizing the disorder of a sequence. 

Includes general disorder optimization and a version that only swaps residues within the same class.
This is useful for generating sequences with high disorder while maintaining specific residue classes.
'''
import random
import numpy as np
import metapredict as meta
from goose.data.defined_aa_classes import aa_classes_by_aa

def optimize_disorder(input_sequence, disorder_cutoff=0.5,
                      max_iterations=100, preserve_charge_placement=True,
                      metapredict_version=3,
                      allowed_fraction_below_cutoff=0.01):
    '''
    Function to optimize the disorder of a sequence. 
    Identifies regions that are below the disorder cutoff and moves residues around
    to try to increase the disorder of the sequence. 

    Parameters
    ----------
    input_sequence : str
        The input protein sequence.
    disorder_cutoff : float
        The disorder cutoff value. Default is 0.5.
    max_iterations : int
        Maximum number of optimization iterations. Default is 100.
    preserve_charge_placement : bool
        If True, charged residues (D, E, K, R) will not be modified. Default is True.
    metapredict_version : int
        The version of metapredict to use for disorder prediction. Default is 3.
    allowed_fraction_below_cutoff : float
        Fraction of residues that can be below the disorder cutoff to allow optimization. Default is 0.01.
    
    Returns
    -------
    str
        Optimized sequence with improved disorder score.
    '''
    allowed_res_below_cutoff = int(len(input_sequence) * allowed_fraction_below_cutoff)

    # make sequence a copy of input_sequence
    sequence = input_sequence
    
    # get initial disorder score
    best_disorder_score = np.min(meta.predict_disorder(sequence, version=metapredict_version))

    # convert sequence to numpy array for vectorized operations
    sequence_array = np.array(list(sequence))
    
    # determine which positions can be modified using vectorized operations
    if preserve_charge_placement:
        # create boolean mask for non-charged residues
        charged_mask = np.isin(sequence_array, ['D', 'E', 'K', 'R'])
        potential_targets = np.where(~charged_mask)[0]
    else:
        # all positions are potential targets
        potential_targets = np.arange(len(sequence))


    # iterate through optimization attempts
    for _ in range(max_iterations):
        # find the region of the sequence below the disorder cutoff
        disorder_scores = meta.predict_disorder(sequence, version=metapredict_version)

        # use vectorized operations to find top disordered indices
        # filter potential targets and get their disorder scores
        target_mask = np.isin(np.arange(len(disorder_scores)), potential_targets)
        valid_indices = np.where(target_mask)[0]
        valid_scores = disorder_scores[valid_indices]
        
        # get top 10 most disordered indices using vectorized operations
        if len(valid_indices) >= 10:
            top_k = 10
        else:
            top_k = len(valid_indices)
        
        # use argpartition for efficient top-k selection (faster than full sort)
        top_indices_local = np.argpartition(valid_scores, -top_k)[-top_k:]
        top_disordered_indices = valid_indices[top_indices_local]
        
        # shuffle using numpy for consistency
        np.random.shuffle(top_disordered_indices)

        # find low disorder indices using vectorized operations
        low_disorder_mask = disorder_scores < disorder_cutoff
        low_disorder_all = np.where(low_disorder_mask)[0]
        
        # filter to only include modifiable positions using vectorized operations
        low_disorder_indices = np.intersect1d(low_disorder_all, potential_targets)

        # if all residues are above cutoff, return the sequence
        if len(low_disorder_indices) == 0:
            return sequence
            
        # randomly select how many changes to make in this iteration
        per_round_iteration = random.randint(1, min(len(low_disorder_indices), len(top_disordered_indices)))
        
        # if per_round_iteration is 0, skip this iteration
        if per_round_iteration == 0:
            continue

        # use numpy random choice for efficient sampling without replacement
        indices_to_change = np.random.choice(low_disorder_indices, size=per_round_iteration, replace=False)
        target_swaps = np.random.choice(top_disordered_indices, size=per_round_iteration, replace=False)

        # create sequence as list for efficient modification
        sequence_list = list(sequence)
        
        # perform swaps using vectorized indexing
        for i, (index_to_change, swap_index) in enumerate(zip(indices_to_change, target_swaps)):
            sequence_list[index_to_change], sequence_list[swap_index] = sequence_list[swap_index], sequence_list[index_to_change]
        
        # check the disorder of the new sequence
        new_sequence = ''.join(sequence_list)
        cur_disorder = meta.predict_disorder(new_sequence, version=metapredict_version)
        new_disorder_min = np.min(cur_disorder)
        
        # if minimum disorder is above cutoff, return immediately
        if new_disorder_min >= disorder_cutoff:
            return new_sequence
        
        # use vectorized operation to count residues below cutoff
        residues_below_cutoff = np.sum(cur_disorder < disorder_cutoff)
        
        # update sequence if improvement is found
        if new_disorder_min > best_disorder_score:
            sequence = new_sequence
            best_disorder_score = new_disorder_min
            # if the number of residues below cutoff is within allowed fraction, return
            if residues_below_cutoff <= allowed_res_below_cutoff:
                return new_sequence

    return sequence



def optimize_disorder_within_class(
                      input_sequence, disorder_cutoff=0.5,
                      max_iterations=100,
                      metapredict_version=3,
                      num_iter_without_improvement_stop=20):
    '''
    Function to optimize the disorder of a sequence. 
    Identifies regions that are below the disorder cutoff and moves residues around
    to try to increase the disorder of the sequence. 
    This function will only swap residues in a sequence with residues 
    that are within the same class.

    Parameters
    ----------
    input_sequence : str
        The input protein sequence.
    disorder_cutoff : float
        The disorder cutoff value. Default is 0.5.
    max_iterations : int
        Maximum number of optimization iterations. Default is 100.
    metapredict_version : int
        The version of metapredict to use for disorder prediction. Default is 3.
    num_iter_without_improvement_stop : int
        Number of iterations without improvement before stopping optimization. Default is 20.

    Notes
    -----
    The function will return the optimized sequence with improved disorder score.
    If no improvement is found after the specified number of iterations, it will return the best sequence found so far.
    
    Returns
    -------
    str
        Optimized sequence with improved disorder score.
    '''

    # make sequence a copy of input_sequence
    sequence = input_sequence
    best_sequence = input_sequence
    
    # get initial disorder score and full disorder profile
    initial_disorder = meta.predict_disorder(sequence, version=metapredict_version)
    best_disorder_score = np.min(initial_disorder)

    # convert sequence to numpy array for vectorized operations
    sequence_array = np.array(list(sequence))
    
    # determine which positions can be modified using vectorized operations
    # create boolean mask for non-charged residues
    immutable_mask = np.isin(sequence_array, ['G', 'P', 'C', 'H'])
    # find indices of reisudes that are not in immutable_mask
    potential_targets = np.where(~immutable_mask)[0]

    # if no potential targets, return the original sequence
    if len(potential_targets) == 0:
        return sequence

    no_improvement = 0

    # iterate through optimization attempts
    for num_i in range(max_iterations):
        # find the region of the sequence below the disorder cutoff
        disorder_scores = meta.predict_disorder(sequence, version=metapredict_version)

        # use vectorized operations to find top disordered indices
        # filter potential targets and get their disorder scores
        target_mask = np.isin(np.arange(len(disorder_scores)), potential_targets)
        valid_indices = np.where(target_mask)[0]
        valid_scores = disorder_scores[valid_indices]
        
        # get top 10 most disordered indices using vectorized operations
        if len(valid_indices) >= 10:
            top_k = 10
        else:
            top_k = len(valid_indices)
        
        # use argpartition for efficient top-k selection (faster than full sort)
        top_indices_local = np.argpartition(valid_scores, -top_k)[-top_k:]
        top_disordered_indices = valid_indices[top_indices_local]
        
        # shuffle using numpy for consistency
        np.random.shuffle(top_disordered_indices)

        # find low disorder indices using vectorized operations
        low_disorder_mask = disorder_scores < disorder_cutoff
        low_disorder_all = np.where(low_disorder_mask)[0]
        
        # filter to only include modifiable positions using vectorized operations
        low_disorder_indices = np.intersect1d(low_disorder_all, potential_targets)

        # if all residues are above cutoff, return the sequence
        if len(low_disorder_indices) == 0:
            return sequence
            
        # create sequence as list for efficient modification
        sequence_list = list(sequence)
        
        # collect all potential swaps with their expected benefit
        potential_swaps = []
        
        # for each top disordered residue, find same-class residues with lower disorder to swap
        for high_disorder_idx in top_disordered_indices:
            high_disorder_residue = sequence_list[high_disorder_idx]
            high_disorder_score = disorder_scores[high_disorder_idx]
            
            # find all positions with the same amino acid class as the high disorder residue
            same_class_residues = aa_classes_by_aa[high_disorder_residue]
            
            # find positions in sequence with same class residues and lower disorder scores
            for i in potential_targets:
                if (sequence_list[i] in same_class_residues and 
                    disorder_scores[i] < high_disorder_score and
                    i != high_disorder_idx):
                    # calculate potential benefit of this swap
                    benefit = high_disorder_score - disorder_scores[i]
                    potential_swaps.append((high_disorder_idx, i, benefit))
        
        # sort potential swaps by benefit (descending) and take the best ones
        potential_swaps.sort(key=lambda x: x[2], reverse=True)
        
        # perform the top swaps, but avoid conflicting swaps
        used_indices = set()
        swaps_made = 0
        max_swaps_per_iteration = min(len(top_disordered_indices), 10)
        
        for high_idx, low_idx, benefit in potential_swaps:
            if swaps_made >= max_swaps_per_iteration:
                break
            if high_idx not in used_indices and low_idx not in used_indices:
                # perform the swap
                sequence_list[high_idx], sequence_list[low_idx] = sequence_list[low_idx], sequence_list[high_idx]
                used_indices.add(high_idx)
                used_indices.add(low_idx)
                swaps_made += 1
        
        # check the disorder of the new sequence
        new_sequence = ''.join(sequence_list)
        cur_disorder = meta.predict_disorder(new_sequence, version=metapredict_version)
        new_disorder_min = np.min(cur_disorder)
        
        # if minimum disorder is above cutoff, return immediately
        if new_disorder_min >= disorder_cutoff:
            return new_sequence
        
        # improved acceptance criteria: consider min
        improvement_found = False
        if new_disorder_min > best_disorder_score:
            improvement_found = True
        elif new_disorder_min == best_disorder_score:
            improvement_found = True
        
        # update sequence if improvement is found
        if improvement_found:
            sequence = new_sequence
            best_sequence = new_sequence
            best_disorder_score = new_disorder_min
            no_improvement = 0
            
        else:
            no_improvement += 1
            
        # if no improvement after specified iterations, return the best sequence found
        if no_improvement >= num_iter_without_improvement_stop:
            break

    return best_sequence

