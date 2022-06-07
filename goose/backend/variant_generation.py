'''
code for variant generation, actually predicts disorder.
'''
from goose.backend.sequence_generation_backend import identify_residue_positions, get_optimal_residue, optimal_residue_key, random_amino_acid, create_seq_by_props, fast_predict_disorder
import metapredict as meta
from amino_acids import AminoAcid

import random

def optimize_disorder_within_class_once(sequence):
    '''
    function to move around residues within a sequence
    within individual classes to maximize disorder in variants 
    where residues for each class must be constant in location
    '''
    # define which residues in each class
    aromatics = ['F', 'W', 'Y']
    polar = ['Q', 'N', 'S', 'T']
    hydrophobics = ['I', 'V', 'L', 'A', 'M']
    positive = ['K', 'R']
    negative = ['D', 'E']

    # get disorder
    disorder = meta.predict_disorder(sequence)
    
    # lowest val for disorder
    lowest_val = min(disorder)

    # take in the lowest val + 0.1 as possible targets
    lowest_val_upper_range = lowest_val + 0.1

    # now get max disorder
    highest_val = max(disorder)

    # take in the max val - 0.1 as possible targets
    highest_val_lower_range = highest_val - 0.1

    # get positions to target
    low_disorder_pos = []
    high_disorder_pos = []

    # also use those ranges to identify target residues classes
    low_disorder_pos_classes = []
    high_disorder_pos_classes = []

    # add positions to list
    for aa_ind in range(0, len(sequence)):
        cur_val = disorder[aa_ind]
        if cur_val >= lowest_val and cur_val <= lowest_val_upper_range:
            low_disorder_pos.append(aa_ind)
            curaa = sequence[aa_ind]
            curclass = AminoAcid.return_AA_class(curaa)
            low_disorder_pos_classes.append(curclass)
        if cur_val <= highest_val and cur_val >= highest_val_lower_range:
            high_disorder_pos.append(aa_ind)
            curaa = sequence[aa_ind]
            curclass = AminoAcid.return_AA_class(curaa)
            high_disorder_pos_classes.append(curclass)

    # now get possible classes to swap
    usable_swap_classes = ['polar', 'hydrophobic', 'positive', 'negative', 'aromatic']
    possible_swap_classes = []
    for cur_class in low_disorder_pos_classes:
        if cur_class in high_disorder_pos_classes:
            if cur_class not in possible_swap_classes:
                if cur_class in usable_swap_classes:
                    possible_swap_classes.append(cur_class)

    # if nothing possible to swap, return the seq
    if possible_swap_classes == []:
        return sequence

    # now narrow down possible residues
    final_residue_low_disorder_loc = []
    final_residue_high_disorder_loc = []

    for aa_loc in range(0, len(low_disorder_pos)):
        if low_disorder_pos_classes[aa_loc] in possible_swap_classes:
            final_residue_low_disorder_loc.append(low_disorder_pos[aa_loc])

    for aa_loc in range(0, len(high_disorder_pos)):
        if high_disorder_pos_classes[aa_loc] in possible_swap_classes:
            final_residue_high_disorder_loc.append(high_disorder_pos[aa_loc])    


    # choose a random loc to swap
    chosen_lower_ind = random.choice(final_residue_low_disorder_loc)

    # get class of that residue
    chosen_lower_aa_class = AminoAcid.return_AA_class(sequence[chosen_lower_ind])

    # get candidate residues from higher residues
    final_chosen_higher_ind_choices = []
    for higher_res_ind in final_residue_high_disorder_loc:
        if AminoAcid.return_AA_class(sequence[higher_res_ind]) == chosen_lower_aa_class:
            final_chosen_higher_ind_choices.append(higher_res_ind)

    if len(final_chosen_higher_ind_choices) > 1:
        chosen_higher_ind = random.choice(final_chosen_higher_ind_choices)
    else:
        chosen_higher_ind = final_chosen_higher_ind_choices[0]

    # get the identity of the higher residue and lower residue
    final_higher_res_identity = sequence[chosen_higher_ind]
    final_lower_res_identity = sequence[chosen_lower_ind]

    # rebuild the sequence
    rebuilt_sequence = ''
    for i in range(0, len(sequence)):
        if i == chosen_lower_ind:
            rebuilt_sequence += final_higher_res_identity
        elif i == chosen_higher_ind:
            rebuilt_sequence += final_lower_res_identity
        else:
            rebuilt_sequence += sequence[i]

    return rebuilt_sequence


def optimize_disorder_within_class(sequence, num_iterations=500):
    '''
    function that uses optimize_disorder_within_class
    but follows multiple iterations to try to get a better
    sequence as far as disorder values. Will keep the residue
    order and type as far as class of amino acid the same. DOES 
    NOT CHANGE SEQUENCE COMPOSITION
    '''
    opt_seq = optimize_disorder_within_class_once(sequence)
    for i in range(0, num_iterations):
        new_seq = optimize_disorder_within_class_once(opt_seq)
        opt_seq = new_seq
    return opt_seq



def optimize_disorder_once_constant_residues(sequence, constant_residues = []):
    '''
    function to move around residues within a sequence
    while keeping desired constant residues.. well, constant.
    '''
    # forbidden positions for changes. forbidden_positions ... probably
    # could have come up with a better name for that list.
    forbidden_positions = []
    for i in constant_residues:
        forbidden_positions.extend(identify_residue_positions(sequence, i))

    # get disorder
    disorder = meta.predict_disorder(sequence)
    
    # lowest val for disorder
    lowest_val = min(disorder)

    # take in the lowest val + 0.1 as possible targets
    lowest_val_upper_range = lowest_val + 0.1

    # now get max disorder
    highest_val = max(disorder)

    # take in the max val - 0.1 as possible targets
    highest_val_lower_range = highest_val - 0.1

    # get positions to target
    low_disorder_pos = []
    high_disorder_pos = []

    # add positions to list
    for aa_ind in range(0, len(sequence)):
        if aa_ind not in forbidden_positions:
            cur_val = disorder[aa_ind]
            if cur_val >= lowest_val and cur_val <= lowest_val_upper_range:
                low_disorder_pos.append(aa_ind)
                curaa = sequence[aa_ind]
            if cur_val <= highest_val and cur_val >= highest_val_lower_range:
                high_disorder_pos.append(aa_ind)
                curaa = sequence[aa_ind]

    if low_disorder_pos == [] or high_disorder_pos == []:
        # add positions to list
        best_low_val = 1
        best_high_val = -1
        for aa_ind in range(0, len(sequence)):
            if aa_ind not in forbidden_positions:
                cur_val = disorder[aa_ind]
                if cur_val < best_low_val:
                    best_low_val = cur_val
                    low_disorder_pos.append(aa_ind)
                    curaa = sequence[aa_ind]
                if cur_val > best_high_val:
                    best_high_val = cur_val
                    high_disorder_pos.append(aa_ind)
                    curaa = sequence[aa_ind]
    
    # choose a random loc to swap
    chosen_lower_ind = random.choice(low_disorder_pos)

    # get class of that residue
    chosen_higher_ind = random.choice(high_disorder_pos)

    # get the identity of the higher residue and lower residue
    final_higher_res_identity = sequence[chosen_higher_ind]
    final_lower_res_identity = sequence[chosen_lower_ind]

    # rebuild the sequence
    rebuilt_sequence = ''
    for i in range(0, len(sequence)):
        if i == chosen_lower_ind:
            rebuilt_sequence += final_higher_res_identity
        elif i == chosen_higher_ind:
            rebuilt_sequence += final_lower_res_identity
        else:
            rebuilt_sequence += sequence[i]

    return rebuilt_sequence


def optimize_disorder_constant_residues(sequence, constant_residues=[], num_iterations=500):
    '''
    function that uses optimize_disorder_within_class
    but follows multiple iterations to try to get a better
    sequence as far as disorder values. Will keep the residue
    order and type as far as class of amino acid the same. DOES 
    NOT CHANGE SEQUENCE COMPOSITION
    '''
    opt_seq = optimize_disorder_once_constant_residues(sequence, constant_residues=constant_residues)
    for i in range(0, num_iterations):
        new_seq = optimize_disorder_once_constant_residues(opt_seq, constant_residues=constant_residues)
        opt_seq = new_seq
    return opt_seq





