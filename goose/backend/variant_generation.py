'''
code for variant generation, actually predicts disorder.
'''

from goose.backend.sequence_generation_backend import identify_residue_positions, fast_predict_disorder, gen_sequence, create_seq_by_props
from goose.backend.amino_acids import AminoAcid
from goose.backend import parameters
from goose.backend.protein import Protein
from goose.goose_exceptions import GooseError, GooseInputError, GooseFail, GooseException
from goose.backend.variant_generation_backend import create_kappa_variant, create_shuffle_variant, create_constant_residue_variant, create_hydropathy_class_variant, create_new_variant, create_constant_class_variant, create_new_var_constant_class_nums

import metapredict as meta

import random


def sequence_variant_disorder(current_sequence, original_disorder_list, cutoff_val=parameters.DISORDER_THRESHOLD, strict=False):
    """
    Function for determining if a sequence variant is disordered. The
    purpose of this over the typical check_disorder function for generated
    sequences is that the input sequence for a disordered variant may have
    regions that are below the cutoff value (more ordered). This function
    gets around that by allowing those regions that were originally not as disordered
    in the variant sequence to remain less disordered. However, the region will
    NEVER BE LESS DISORDERED than the input sequence, but it may become more disordered.

    This function (if strict is set to False) will allow an occassional residue
    to dip below the cutoff value.

    Parameters
    -------------

    current_sequence : String
        The variant sequence that is being checked for disorder

    original_disorder_list : List
        A list of the disorder values for the sequence input into the
        sequence variant generation function
    
    cutoff_value : float
        The cutoff value for disorder. A value of at least 0.6 is
        typically used by GOOSE.

    strict : bool
        if set to true, will not count a sequence as disordered even if a single amino
        acid falls below the cutoff value.

    Returns
    ---------

    Bool
        Returns True if the sequence is disordered
        
        Returns False if the sequence has too many residues below 
        the cutoff_val.

    """

    # first get the list of disordered residues for the current sequence
    variant_disorder_list = meta.predict_disorder(current_sequence)

    # first make a 'modified' disorder value list so that
    # residues can go below cutoff.
    adjusted_disorder_values = []

    # for the disorder values in the disorder list
    for disorder_vals in range(0, len(variant_disorder_list)):

        # get the original disorder value at the current position
        input_disorder = original_disorder_list[disorder_vals]
        
        # if the input disorder value is greater than the specified cutoff value
        if input_disorder > cutoff_val:
            # reduce the necessary input disorder value to the cutoff value
            input_disorder = cutoff_val

        # otherwise, the input disorder value is simply equal to the disorder
        # value for the original sequence residue at the current position
        else:
            input_disorder = input_disorder

        adjusted_disorder_values.append(input_disorder)

    # keep track of consecutively ordered residues and total ordered
    cur_order_streak = 0
    total_ordered_residues = 0

    # allow up to 5 % of residues to be 'ordered' provided they aren't consecutive
    allowed_order_residues = round(0.05*len(current_sequence))
    if allowed_order_residues < 1:
        allowed_order_residues = 1

    # calculate number of conseuctively allowed residues
    # going to allow up to the length of the sequence over 25 but 
    # not greater than 10
    consec_allowed = round(len(current_sequence)/25)
    if consec_allowed < 1:
        consec_allowed = 1
    if consec_allowed > 10:
        consec_allowed = 10

    for disorder_val in range(0, len(adjusted_disorder_values)):
        # get current values for the original and the variant
        variant_disorder = variant_disorder_list[disorder_val]
        original_disorder = adjusted_disorder_values[disorder_val]

        # see if strict set to true
        if strict == True:
            if variant_disorder < original_disorder:
                return False
        else:
            if variant_disorder < original_disorder:
                cur_order_streak += 1
                total_ordered_residues += 1
            else:
                cur_order_streak = 0
            # check to see if order streak is too high
            # or total num ordered residues too high
            if cur_order_streak > consec_allowed:
                return False
            if total_ordered_residues > allowed_order_residues:
                return False

    # if all residues in the sequence variant are greater than either
    # the original value or at least the cutoff value, return True
    return True





def optimize_disorder_within_class_once(sequence):
    '''
    function to move around residues within a sequence
    within individual classes to maximize disorder in variants 
    where residues for each class must be constant in location

    parameters 
    -----------
    sequence : string
        the amino acid sequence as a string

    returns 
    --------
    rebuilt_sequence : string
        returns the rebuilt (optimized) sequence

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

    paramters
    ----------
    sequence : string
        the amino acid sequence as a string

    num_iterations : int
        the number of times to run the optimizer.

    returns
    -------
    rebuilt_sequence : string
        returns the rebuilt (optimized) sequence as a string

    '''
    # make an optimized seq before puting it in for iterations
    opt_seq = optimize_disorder_within_class_once(sequence)
    for i in range(0, num_iterations):
        new_seq = optimize_disorder_within_class_once(opt_seq)
        opt_seq = new_seq
    return opt_seq



def optimize_disorder_once_constant_residues(sequence, constant_residues = []):
    '''
    function to move around residues within a sequence
    while keeping desired constant residues.. well, constant.

    paramters
    ---------
    sequence : string
        the amino acid sequence as a string

    constant_residues : list
        the list of residues that you want to be held constant
        in the input sequence
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

    paramters
    ----------
    sequence : string
        the amino acid sequence as a string

    num_iterations : int
        the number of times to run the optimizer.

    returns
    -------
    rebuilt_sequence : string
        returns the rebuilt (optimized) sequence as a string

    '''

    # optimize the sequence once
    opt_seq = optimize_disorder_once_constant_residues(sequence, constant_residues=constant_residues)
    # now go through and continue optimizations
    for i in range(0, num_iterations):
        new_seq = optimize_disorder_once_constant_residues(opt_seq, constant_residues=constant_residues)
        opt_seq = new_seq
    return opt_seq


#/-/-/-/-/-/-/-/-/-/-/-/-/-/-/-#
#       *       *       *      #
# Sequence Variant Generators  #
#       *       *       *      #
#/-/-/-/-/-/-/-/-/-/-/-/-/-/-/-#

def gen_new_var_constant_class(sequence, attempts=5, 
    disorder_threshold = parameters.DISORDER_THRESHOLD, strict_disorder=False):
    '''
    A function to generate a variant where the sequence composition is new but
    the numbers of each residue from each class is the same. The overall properties
    of the generated sequence will also be constant.

    parameters
    ----------
    sequence : string
        the amino acid sequence as a string
    attempts : int
        the number of times to attempt to make the variant
    disorder_threshold : float
        threshold value for disorder value between 0 and 1. 
    strict_disorder : bool
        whether to use a strict disorder criteria where
        if set to true, any value below the cutoff or that
        disorder position for the variant will not be 
        considered disordered

    returns
    -------
    disordered_seq : string
        returns a disordered sequence as a string
    '''
    # get original sequence disorder
    starting_disorder = meta.predict_disorder(sequence)

    for attempt_num in range(0, attempts):
        disordered_seq = create_new_var_constant_class_nums(sequence)
        if sequence_variant_disorder(disordered_seq, starting_disorder, 
            cutoff_val=disorder_threshold, strict=strict_disorder) == True:
            # if passes the 'disorder test', return the seq
            return disordered_seq

        else:
            newsequence = optimize_disorder_constant_residues(disordered_seq, num_iterations=50)
            if sequence_variant_disorder(newsequence, starting_disorder, 
                cutoff_val=disorder_threshold, strict=strict_disorder) == True:
                return newsequence   

    raise GooseFail('Unable to generate sequence.')





def gen_constant_class_variant(sequence, attempts=5, 
    disorder_threshold = parameters.DISORDER_THRESHOLD, strict_disorder=False):
    '''
    function to generate a variant with the same properties as the 
    input variant as well as the same order of amino acids as
    far as class and the same number in each class

    parameters
    ----------
    sequence : string
        the amino acid sequence as a string
    attempts : int
        the number of times to attempt to make the variant
    disorder_threshold : float
        threshold value for disorder value between 0 and 1. 
    strict_disorder : bool
        whether to use a strict disorder criteria where
        if set to true, any value below the cutoff or that
        disorder position for the variant will not be 
        considered disordered

    returns
    -------
    disordered_seq : string
        returns a disordered sequence as a string
    '''

    # get original sequence disorder
    starting_disorder = meta.predict_disorder(sequence)

    for attempt_num in range(0, attempts):
        disordered_seq = create_constant_class_variant(sequence)
        if sequence_variant_disorder(disordered_seq, starting_disorder, 
            cutoff_val=disorder_threshold, strict=strict_disorder) == True:
            # if passes the 'disorder test', return the seq
            return disordered_seq
        
        else:
            newsequence = optimize_disorder_within_class(disordered_seq, num_iterations=50)
            if sequence_variant_disorder(newsequence, starting_disorder, 
                cutoff_val=disorder_threshold, strict=strict_disorder) == True:
                return newsequence   

    raise GooseFail('Unable to generate sequence.')




def gen_new_variant(sequence, attempts=5, 
    disorder_threshold = parameters.DISORDER_THRESHOLD, strict_disorder=False):
    '''
    function to generate a variant that is completely different
    in sequence to the input but has all the same overall parameters.
    Does not account for specific classes of residues.

    parameters
    ----------
    sequence : string
        the amino acid sequence as a string
    attempts : int
        the number of times to attempt to make the variant
    disorder_threshold : float
        threshold value for disorder value between 0 and 1. 
    strict_disorder : bool
        whether to use a strict disorder criteria where
        if set to true, any value below the cutoff or that
        disorder position for the variant will not be 
        considered disordered

    returns
    -------
    disordered_seq : string
        returns a disordered sequence as a string
    '''

    # get original sequence disorder
    starting_disorder = meta.predict_disorder(sequence)

    for attempt_num in range(0, attempts):
        disordered_seq = create_new_variant(sequence)
        if sequence_variant_disorder(disordered_seq, starting_disorder, 
            cutoff_val=disorder_threshold, strict=strict_disorder) == True:
            # if passes the 'disorder test', return the seq
            return disordered_seq
        else:
            newsequence = optimize_disorder_constant_residues(disordered_seq, num_iterations=50)
            if sequence_variant_disorder(newsequence, starting_disorder, 
                cutoff_val=disorder_threshold, strict=strict_disorder) == True:
                return newsequence   
    raise GooseFail('Unable to generate sequence.')




def gen_hydropathy_class_variant(sequence, hydropathy, allowed_hydro_error = parameters.HYDRO_ERROR,
    attempts=5, disorder_threshold = parameters.DISORDER_THRESHOLD, strict_disorder=False):
    '''
    function to take in a sequence and make a variant that adjusts the
    hydropathy while keeping the position and nuimber of amino acids the
    same by class of amino acid

    parameters
    ----------
    sequence : string
        the amino acid sequence as a string
    hydropathy : float
        the hydropathy value desired for the variant between 0 and 9
    allowed_hydro_error : float
        how far off the final sequence hydropathy can be from the 
        desired hydropathy value
    attempts : int
        the number of times to attempt to make the variant
    disorder_threshold : float
        threshold value for disorder value between 0 and 1. 
    strict_disorder : bool
        whether to use a strict disorder criteria where
        if set to true, any value below the cutoff or that
        disorder position for the variant will not be 
        considered disordered

    returns
    -------
    disordered_seq : string
        returns a disordered sequence as a string
    '''

    # get original sequence disorder
    starting_disorder = meta.predict_disorder(sequence)

    for attempt_num in range(0, attempts):
        disordered_seq = create_hydropathy_class_variant(sequence, hydro=hydropathy, 
            allowed_hydro_error = allowed_hydro_error)
        if sequence_variant_disorder(disordered_seq, starting_disorder, 
            cutoff_val=disorder_threshold, strict=strict_disorder) == True:
            # if passes the 'disorder test', return the seq
            return disordered_seq
        
        else:
            newsequence = optimize_disorder_within_class(disordered_seq, num_iterations=50)
            if sequence_variant_disorder(newsequence, starting_disorder, 
                cutoff_val=disorder_threshold, strict=strict_disorder) == True:
                return newsequence

    raise GooseFail('Unable to generate sequence.')




def gen_constant_residue_variant(sequence, constant_residues = [],
    attempts=5, disorder_threshold = parameters.DISORDER_THRESHOLD, strict_disorder=False):
    '''
    function that will generate a new sequence variant
    where specific residues are held constant. The 
    variant will have the same aggregate properties
    as the original sequence.

    parameters
    ----------
    sequence : string
        The amino acid sequence to make a variant of as a string

    constant_residues : list
        A list of residues to hold constant in the sequence variant
    attempts : int
        the number of times to attempt to make the variant
    disorder_threshold : float
        threshold value for disorder value between 0 and 1. 
    strict_disorder : bool
        whether to use a strict disorder criteria where
        if set to true, any value below the cutoff or that
        disorder position for the variant will not be 
        considered disordered

    returns
    -------
    disordered_seq : string
        returns a disordered sequence as a string
    '''

    # get original sequence disorder
    starting_disorder = meta.predict_disorder(sequence)

    for attempt_num in range(0, attempts):
        disordered_seq = create_constant_residue_variant(sequence, 
            constant_residues=constant_residues)
        if sequence_variant_disorder(disordered_seq, starting_disorder, 
            cutoff_val=disorder_threshold, strict=strict_disorder) == True:
            # if passes the 'disorder test', return the seq
            return disordered_seq
        else:
            #optimize_disorder_constant_residues(sequence, constant_residues=[], num_iterations=500):
            newsequence = optimize_disorder_constant_residues(disordered_seq, 
                constant_residues = constant_residues, num_iterations=50)
            if sequence_variant_disorder(newsequence, starting_disorder, 
                cutoff_val=disorder_threshold, strict=strict_disorder) == True:
                return newsequence

    raise GooseFail('Unable to generate sequence.')



def gen_shuffle_variant(sequence, shuffle_regions = [], use_index=False,
    attempts=5, disorder_threshold = parameters.DISORDER_THRESHOLD, strict_disorder=False):
    '''
    Function that will shuffle specific regions of an IDR.
    Multiple regions can be specified simultaneously.

    parameters
    ----------
    sequences : string
        the string of the amino acid sequence to shuffle

    shuffle_regions : list of lists
        A nested list containing the regions of the
        protein to shuffle. If a single list is input, 
        the function just does that region (no need
        for unnecessary listed nests.)

    use_index : Bool
        Whether the first amino acid should be considered '1'
        or '0'. By default will start at 1 (for the first amino acid.), 
        which doesn't follow standard coding formatting but
        will be easier for most end user. If it drives you bananas,
        feel free to set use_index to True.

    attempts : int
        the number of times to attempt to make the variant
    disorder_threshold : float
        threshold value for disorder value between 0 and 1. 
    strict_disorder : bool
        whether to use a strict disorder criteria where
        if set to true, any value below the cutoff or that
        disorder position for the variant will not be 
        considered disordered

    returns
    -------
    disordered_seq : string
        returns a disordered sequence as a string
    '''

    # get original sequence disorder
    starting_disorder = meta.predict_disorder(sequence)

    for attempt_num in range(0, attempts):
        disordered_seq = create_shuffle_variant(sequence, 
            shuffle_regions=shuffle_regions, use_index=use_index)
        if sequence_variant_disorder(disordered_seq, starting_disorder, 
            cutoff_val=disorder_threshold, strict=strict_disorder) == True:
            # if passes the 'disorder test', return the seq
            return disordered_seq
    raise GooseFail('Unable to generate sequence.')



def gen_kappa_variant(sequence, kappa, allowed_kappa_error = parameters.MAXIMUM_KAPPA_ERROR,
    attempts=5, disorder_threshold = parameters.DISORDER_THRESHOLD, strict_disorder=False):
    '''
    Function to generate a sequence with a user-defined
    kappa value. Requires kappa calculation using 
    SPARROW. Kappa is a function of charge asymmetry, larger
    kappa values have more asymmetrically distributed
    charged residues.

    parameters
    ----------
    sequence : str
        the amino acid sequence as a sting

    kappa : float
        the desired kappa value between 0 and 1 as a float
        1 = max asymmetry
        0 = minimum charge asymmetry
    allowed_kappa_error : float
        How much difference there can be between the 
        specified kappa value and the kappa value
        calculated for the returned sequence
    attempts : int
        the number of times to attempt to make the variant
    disorder_threshold : float
        threshold value for disorder value between 0 and 1. 
    strict_disorder : bool
        whether to use a strict disorder criteria where
        if set to true, any value below the cutoff or that
        disorder position for the variant will not be 
        considered disordered

    returns
    -------
    disordered_seq : string
        returns a disordered sequence as a string
    '''

    # get original sequence disorder
    starting_disorder = meta.predict_disorder(sequence)

    for attempt_num in range(0, attempts):
        disordered_seq = create_kappa_variant(sequence, kappa=kappa, 
            allowed_kappa_error = allowed_kappa_error)

        if sequence_variant_disorder(disordered_seq, starting_disorder, 
            cutoff_val=disorder_threshold, strict=strict_disorder) == True:
            # if passes the 'disorder test', return the seq
            return disordered_seq
        else:
            newsequence = optimize_disorder_within_class(disordered_seq, num_iterations=50)
            if sequence_variant_disorder(newsequence, starting_disorder, 
                cutoff_val=disorder_threshold, strict=strict_disorder) == True:
                return newsequence
    raise GooseFail('Unable to generate sequence.')
