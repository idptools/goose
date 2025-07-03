'''
Code for optimizing the disorder of a sequence. 
'''
import random
import numpy as np
import metapredict as meta

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


