'''
code for variant generation, actually predicts disorder.
'''

from goose.backend.sequence_generation_backend import identify_residue_positions, fast_predict_disorder, gen_sequence, create_seq_by_props
from goose.backend.amino_acids import AminoAcid
from goose.backend import parameters
from goose.backend.protein import Protein
from goose.goose_exceptions import GooseError, GooseInputError, GooseFail, GooseException
from goose.backend.variant_generation_backend import create_kappa_variant, create_region_shuffle_variant, create_constant_residue_variant, create_hydropathy_class_variant, create_new_variant, create_constant_class_variant, create_new_var_constant_class_nums, create_asymmetry_variant, create_fcr_class_variant, create_ncpr_class_variant, create_all_props_class_variant
from goose.backend.seq_by_dimension_backend import make_rg_re_variant, predict_rg, predict_re

import metapredict as meta

import random

def sequence_variant_disorder(current_seq_or_disorder, original_disorder_list, 
    cutoff_val=parameters.DISORDER_THRESHOLD, strict=False, input_disorder_val=False):
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

    current_seq_or_disorder : String
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

    input_disorder_val : bool  
        whether the disorder value is being input

    Returns
    ---------

    Bool
        Returns True if the sequence is disordered
        
        Returns False if the sequence has too many residues below 
        the cutoff_val.

    """
    # first get the list of disordered residues if not input
    if input_disorder_val == False:
        variant_disorder_list = meta.predict_disorder(current_seq_or_disorder)
    else:
        variant_disorder_list = current_seq_or_disorder

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
    allowed_order_residues = round(0.05*len(current_seq_or_disorder))
    if allowed_order_residues < 1:
        allowed_order_residues = 1

    # calculate number of conseuctively allowed residues
    # going to allow up to the length of the sequence over 25 but 
    # not greater than 10
    consec_allowed = round(len(current_seq_or_disorder)/25)
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


def optimize_disorder_within_class(sequence, num_iterations=500, cutoff=parameters.DISORDER_THRESHOLD):
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

    cutoff : float
        The min value required for an amino acid to be considered disordered

    returns
    -------
    rebuilt_sequence : string
        returns the rebuilt (optimized) sequence as a string

    '''
    # make an optimized seq before puting it in for iterations
    opt_seq = optimize_disorder_within_class_once(sequence)
    min_disorder = min(meta.predict_disorder(opt_seq))
    if min_disorder >= cutoff:
        return opt_seq
    else:
        for i in range(0, num_iterations):
            new_seq = optimize_disorder_within_class_once(opt_seq)
            new_dis = meta.predict_disorder(new_seq)
            # if the min disorder of the new seq is greater than the min
            # of the previous seq, update min_disorder and opt_seq values
            if min(new_dis) > min_disorder:
                min_disorder = min(new_dis)
                opt_seq=new_seq
            # if disordered enough, kill the loop and return the seq
            if min(new_dis) >= cutoff:
                return opt_seq
    # if haven't goptten to objective, return seq as is.
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



def optimize_disorder_constant_residues(sequence, constant_residues=[], num_iterations=500, cutoff=parameters.DISORDER_THRESHOLD):
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

    cutoff : float
        The min value required for an amino acid to be considered disordered

    returns
    -------
    rebuilt_sequence : string
        returns the rebuilt (optimized) sequence as a string

    '''

    # optimize the sequence once
    opt_seq = optimize_disorder_once_constant_residues(sequence, constant_residues=constant_residues)
    min_disorder = min(meta.predict_disorder(opt_seq))
    if min_disorder >= cutoff:
        return opt_seq
    else:
        for i in range(0, num_iterations):
            new_seq = optimize_disorder_once_constant_residues(opt_seq, constant_residues=constant_residues)
            new_dis = meta.predict_disorder(new_seq)
            # if the min disorder of the new seq is greater than the min
            # of the previous seq, update min_disorder and opt_seq values
            if min(new_dis) > min_disorder:
                min_disorder = min(new_dis)
                opt_seq=new_seq
            # if disordered enough, kill the loop and return the seq
            if min(new_dis) >= cutoff:
                return opt_seq
    # return the best sequence
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



def gen_region_shuffle_variant(sequence, shuffle_regions = [], use_index=False,
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
        disordered_seq = create_region_shuffle_variant(sequence, 
            shuffle_regions=shuffle_regions, use_index=use_index)
        if sequence_variant_disorder(disordered_seq, starting_disorder, 
            cutoff_val=disorder_threshold, strict=strict_disorder) == True:
            # if passes the 'disorder test', return the seq
            return disordered_seq
    raise GooseFail('Unable to generate sequence.')



def gen_kappa_variant(sequence, kappa, allowed_kappa_error = parameters.MAXIMUM_KAPPA_ERROR,
    attempts=10, disorder_threshold = parameters.DISORDER_THRESHOLD, strict_disorder=False):
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
    try:
        disordered_seq = create_kappa_variant(sequence, kappa=kappa, 
            allowed_kappa_error = allowed_kappa_error)
    except:
        raise GooseFail('Unable to generate kappa variant.')
    if sequence_variant_disorder(disordered_seq, starting_disorder, 
        cutoff_val=disorder_threshold, strict=strict_disorder) == True:
        # if passes the 'disorder test', return the seq
        return disordered_seq
    else:
        newsequence = optimize_disorder_within_class(disordered_seq, num_iterations=50)
        if sequence_variant_disorder(newsequence, starting_disorder, 
            cutoff_val=disorder_threshold, strict=strict_disorder) == True:
            return newsequence
    raise GooseFail('Unable to generate kappa variant.')


def gen_asymmetry_variant(sequence, increase_decrease, aa_class, num_change=None, window_size = 6,
    attempts=10, disorder_threshold = parameters.DISORDER_THRESHOLD, strict_disorder=False):
    '''
    function to change the asymmetry of a residue / residue class in a sequence

    paramters
    ---------
    sequence : str
        the amino acid sequence as a string
    increase_decrease : str
        whether to increase or decrease the charge asymmetry
            set to 'decrease' to decrease the asymmetry
            set to 'increase' to increase the charge asymmetry
    aa_class : string or list
        the residues or parameter to increase or decrease asymmetry of
            string options include 
                'negative' - negative residues
                'positive' - positive residues    
                'proline' - prolines
                'aromatic' - W Y F
                'aliphatic' - I V L A M
                'polar' - Q N S T
            you can also specify a list of amino acids to change.
            To do this, simply input a list. Ex. 
    window_size : Int
        The size of the window to break the sequence down into as 
        far as regions calcualted for the specific res class.
    num_change : Int
        The number of times to change the asymmetry of the sequence for the variant
    attempts : int
        the number of times ot try to make the sequence

    disorder_threshold : float
        the threshold value required for an amino acid
        to be considered disordered

    strict_disorder : Bool
        whether or not to require all disorder values to be 
        over threshold or if it si okay to use the values
        from the input sequence
    returns
    -------
    disorder_seq : string
        returns the amino acid sequence of teh generated sequence as a string
    '''
    # get original sequence disorder
    starting_disorder = meta.predict_disorder(sequence)

    # define classes
    class_dict = {'aromatic':['W', 'Y', 'F'], 'aliphatic':['I', 'V', 'L', 'A', 'M'], 'polar':['Q', 'N', 'S', 'T'], 'proline': ['P'], 'negative':['D', 'E'], 'positive':['K', 'R']}
    # if user inputs a list, use that list
    if type(aa_class) == list:
        target_amino_acids = aa_class
    else:
        try:
            target_amino_acids = class_dict[aa_class]
        except:
            raise GooseInputError('The specified aa_class is not a list of amino acids or a class. Classes allowed are: 1. aromatic, 2. aliphatic, 3. polar, 4. proline')


    for attempt_num in range(0, attempts):
        disordered_seq = create_asymmetry_variant(sequence, increase_decrease=increase_decrease, 
                            aa_class = aa_class, window_size=window_size, num_change=num_change)

        if sequence_variant_disorder(disordered_seq, starting_disorder, 
            cutoff_val=disorder_threshold, strict=strict_disorder) == True:
            # if passes the 'disorder test', return the seq
            return disordered_seq
        else:
            newsequence = optimize_disorder_constant_residues(disordered_seq, constant_residues=target_amino_acids, num_iterations=500)
            if sequence_variant_disorder(newsequence, starting_disorder, 
                cutoff_val=disorder_threshold, strict=strict_disorder) == True:
                return newsequence
    raise GooseFail('Unable to generate sequence.')


def gen_fcr_class_variant(sequence, fcr, constant_ncpr=True, use_closest=True, attempts=10, 
    disorder_threshold=parameters.DISORDER_THRESHOLD, strict_disorder=False):
    '''
    function that will alter the FCR of a sequence
    and keep everything else the same while also 
    minimizing change to aromatics and 'special' amino acids
    ie. W, F, Y, P, C.

    **Will change between W, F, and Y as needed, but will
    try to keep it aromatic.

    parameters
    -----------
    sequence : str
        The amino acid sequence as a string

    fcr : float
        the FCR value between 0 and 1

    constant_ncpr : Bool
        whether to allow changes to the NCPR when changing the FCR

    use_closest : Bool
        whether to just use the closest FCR val to the input value.
    
    attempts : int
        the number of times ot try to make the sequence

    disorder_threshold : float
        the threshold value required for an amino acid
        to be considered disordered

    strict_disorder : Bool
        whether or not to require all disorder values to be 
        over threshold or if it si okay to use the values
        from the input sequence
    '''    
    # get original sequence disorder
    starting_disorder = meta.predict_disorder(sequence)    
    
    # attempt to build sequence
    for attempt_num in range(0, attempts):
        disordered_seq = create_fcr_class_variant(sequence, fraction=fcr,
         constant_ncpr=constant_ncpr, use_closest=use_closest)

        if sequence_variant_disorder(disordered_seq, starting_disorder, 
            cutoff_val=disorder_threshold, strict=strict_disorder) == True:
            # if passes the 'disorder test', return the seq
            return disordered_seq
        else:
            newsequence = optimize_disorder_within_class(disordered_seq, num_iterations=500)
            if sequence_variant_disorder(newsequence, starting_disorder, 
                cutoff_val=disorder_threshold, strict=strict_disorder) == True:
                return newsequence
    raise GooseFail('Unable to generate sequence.')



def gen_ncpr_class_variant(sequence, ncpr, constant_fcr=True, use_closest=True, attempts=10, 
    disorder_threshold=parameters.DISORDER_THRESHOLD, strict_disorder=False):
    '''
    function that will alter the FCR of a sequence
    and keep everything else the same while also 
    minimizing change to aromatics and 'special' amino acids
    ie. W, F, Y, P, C.

    **Will change between W, F, and Y as needed, but will
    try to keep it aromatic.

    parameters
    -----------
    sequence : str
        The amino acid sequence as a string

    ncpr : float
        the ncpr value between -1 and 1

    constant_fcr : Bool
        whether to allow changes to the fcr when changing the ncpr

    use_closest : Bool
        whether to just use the closest ncpr val to the input value.
    
    attempts : int
        the number of times ot try to make the sequence

    disorder_threshold : float
        the threshold value required for an amino acid
        to be considered disordered

    strict_disorder : Bool
        whether or not to require all disorder values to be 
        over threshold or if it si okay to use the values
        from the input sequence
    '''    
    # get original sequence disorder
    starting_disorder = meta.predict_disorder(sequence)    
    
    # attempt to build sequence
    for attempt_num in range(0, attempts):
        disordered_seq = create_ncpr_class_variant(sequence, net_charge=ncpr,
         constant_fcr=constant_fcr, use_closest=use_closest)

        if sequence_variant_disorder(disordered_seq, starting_disorder, 
            cutoff_val=disorder_threshold, strict=strict_disorder) == True:
            # if passes the 'disorder test', return the seq
            return disordered_seq
        else:
            newsequence = optimize_disorder_within_class(disordered_seq, num_iterations=500)
            if sequence_variant_disorder(newsequence, starting_disorder, 
                cutoff_val=disorder_threshold, strict=strict_disorder) == True:
                return newsequence
    raise GooseFail('Unable to generate sequence.')


def gen_all_props_class_variant(sequence, hydropathy=None, fcr=None, ncpr=None, kappa=None, 
    attempts=10, disorder_threshold=parameters.DISORDER_THRESHOLD, strict_disorder=False):
    '''
    function to make a variant where you can change hydropathy, fraction charged
    residues, net charge, and kappa all at once while minimizing the changes to
    residues by class in the sequence. As you change the sequence more, you will
    further alter the sequence, even outside of the classes of amino acids. This
    function simply attempts to minimize those changes.

    parameters
    ----------
    sequence : string
        the amino acid sequence to use to make the variant

    hydropathy : float
        a value for mean hydropathy between 0 6.1

    fcr : float
        the FCR value between 0 and 1
    
    ncpr : float
        the ncpr value between -1 and 1

    kappa : float
        the kappa vlaue for the variant between 0 and 1

    attempts : int
        the number of times ot try to make the sequence

    disorder_threshold : float
        the threshold value required for an amino acid
        to be considered disordered

    strict_disorder : Bool
        whether or not to require all disorder values to be 
        over threshold or if it si okay to use the values
        from the input sequence
    '''
    # get original sequence disorder
    starting_disorder = meta.predict_disorder(sequence)    
    
    # attempt to build sequence
    for attempt_num in range(0, attempts):
        disordered_seq = create_all_props_class_variant(sequence, hydropathy=hydropathy,
        fraction=fcr, net_charge=ncpr,kappa=kappa)

        if sequence_variant_disorder(disordered_seq, starting_disorder, 
            cutoff_val=disorder_threshold, strict=strict_disorder) == True:
            # if passes the 'disorder test', return the seq
            return disordered_seq
        else:
            newsequence = optimize_disorder_within_class(disordered_seq, num_iterations=500)
            if sequence_variant_disorder(newsequence, starting_disorder, 
                cutoff_val=disorder_threshold, strict=strict_disorder) == True:
                return newsequence
    raise GooseFail('Unable to generate sequence.')


def gen_targeted_shuffle_variant(sequence, target_aas, attempts=10, 
    disorder_threshold=parameters.DISORDER_THRESHOLD, strict_disorder=False):
    '''
    function that will let you shuffle a sequence by specifying residues
    or classes of residues to shuffle. This is the opposite behavior of
    the create excluded shuffle variant where you specify which
    residues or classes of residues to not target.

    parameters
    ----------
    sequence : str
        the amino acid sequence as a string

    target_aas : str or list
        a list of amino acids to target for shuffling
        or a class of amino acids to target for shuffling
        Possible target classes:
            charged : DEKR
            polar : QNST
            aromatic : FYW
            aliphatic : IVLAM
            negative: DE
            positive : KR

    attempts : int
        the number of times to try to make the sequence

    disorder_threshold : float
        the threshold value required for an amino acid
        to be considered disordered

    strict_disorder : Bool
        whether or not to require all disorder values to be 
        over threshold or if it is okay to use the values
        from the input sequence
    '''
    # dict of classes that are possible to choose
    classdict={'charged':['D', 'E', 'K', 'R'], 'polar':['Q', 'N', 'S', 'T'], 'aromatic':
    ['F', 'W', 'Y'], 'aliphatic': ['I', 'V', 'L', 'A', 'M'], 'negative':['D', 'E'], 'positive':['K', 'R']}
    
    # possible amino acids
    amino_acids = ['A', 'C', 'D', 'E', 'F', 'G', 'H', 'I', 'K', 'L', 'M', 'N', 'P', 'Q', 'R', 'S', 'T', 'V', 'W', 'Y']

    # verify target aas
    if type(target_aas)==str:
        if target_aas in amino_acids:
            raise GooseInputError('You only specified a single amino acid. This will not change your sequence because the amino acids will just change places with itself.')
        elif target_aas in classdict.keys():
            target_aas = classdict[target_aas]
        else:
            raise GooseInputError('The specified target_aas is not a valid amino acid or class of amino acids.')
    elif type(target_aas)==list:
        for i in target_aas:
            if i not in amino_acids:
                raise GooseInputError('The specified target_aas is not a valid amino acid.')
    else:
        raise GooseInputError('The specified target_aas must be type list or string.')

    # get original sequence disorder
    starting_disorder = meta.predict_disorder(sequence)    
    
    # identify target aas.
    target_aa_list=[]
    building_seq=''

    for curind, aa in enumerate(sequence):
        if aa in target_aas:
            target_aa_list.append(aa)
            building_seq+='0'
        else:
            building_seq+=aa    


    # attempt to build sequence
    for attempt_num in range(0, attempts):
        # copy list
        cur_attempt_list = target_aa_list.copy()
        # shuffle list 
        random.shuffle(cur_attempt_list)

        # build final seq
        final_seq=''
        for i in building_seq:
            if i == '0':
                final_seq+=cur_attempt_list.pop()
            else:
                final_seq+=i
        
        # check disorder.
        if sequence_variant_disorder(final_seq, starting_disorder, 
            cutoff_val=disorder_threshold, strict=strict_disorder) == True:
            # if passes the 'disorder test', return the seq
            return final_seq

    # if it doesn't work, raise an error
    raise GooseFail('Unable to generate sequence.')
    
    
def gen_targeted_reposition_variant(sequence, target_aas, attempts=10, 
    disorder_threshold=parameters.DISORDER_THRESHOLD, strict_disorder=False):
    '''
    function that will let you alter a sequence by specifying residues
    or classes of residues to reposition. This will change the positions
    of the specified residues without changing the order of non-specified
    residues within the sequence.

    parameters
    ----------
    sequence : str
        the amino acid sequence as a string

    target_aas : str or list
        a list of amino acids to target for repositioning
        or a class of amino acids to target for repositioning
        Possible target classes:
            charged : DEKR
            polar : QNST
            aromatic : FYW
            aliphatic : IVLAM
            negative: DE
            positive : KR

    attempts : int
        the number of times to try to make the sequence

    disorder_threshold : float
        the threshold value required for an amino acid
        to be considered disordered

    strict_disorder : Bool
        whether or not to require all disorder values to be 
        over threshold or if it is okay to use the values
        from the input sequence
    '''
    # dict of classes that are possible to choose
    classdict={'charged':['D', 'E', 'K', 'R'], 'polar':['Q', 'N', 'S', 'T'], 'aromatic':
    ['F', 'W', 'Y'], 'aliphatic': ['I', 'V', 'L', 'A', 'M'], 'negative':['D', 'E'], 'positive':['K', 'R']}
    
    # possible amino acids
    amino_acids = ['A', 'C', 'D', 'E', 'F', 'G', 'H', 'I', 'K', 'L', 'M', 'N', 'P', 'Q', 'R', 'S', 'T', 'V', 'W', 'Y']

    # verify target aas
    if type(target_aas)==str:
        if target_aas in classdict.keys():
            target_aas = classdict[target_aas]
        else:
            raise GooseInputError('The specified target_aas is not a valid amino acid or class of amino acids.')
    elif type(target_aas)==list:
        for i in target_aas:
            if i not in amino_acids:
                raise GooseInputError('The specified target_aas is not a valid amino acid.')
    else:
        raise GooseInputError('The specified target_aas must be type list or string.')

    # get original sequence disorder
    starting_disorder = meta.predict_disorder(sequence)    
    
    # identify target aas, nontarget aas, and positions.
    target_aa_list=[]
    nontarget_aa_list=[]
    target_aa_positions=[]
    all_positions=[]
    for i, aa in enumerate(sequence):
        if aa in target_aas:
            target_aa_list.append(aa)
            target_aa_positions.append(i)
        else:
            nontarget_aa_list.append(aa)
        all_positions.append(i)

    # attempt to build sequence
    for attempt_num in range(0, attempts):
        # # copy lists
        cur_attempt_target_aas = target_aa_list.copy()
        cur_attempt_nontarget_aas = nontarget_aa_list.copy()
        
        # randomly generate new positions
        new_positions=random.sample(all_positions, k=len(target_aa_list))

        # build final seq
        final_seq=''
        for i in range(len(sequence)):
            if i in new_positions:
                final_seq+=cur_attempt_target_aas.pop(0)
            else:
                final_seq+=cur_attempt_nontarget_aas.pop(0)

        # check disorder.
        if sequence_variant_disorder(final_seq, starting_disorder, 
            cutoff_val=disorder_threshold, strict=strict_disorder) == True:
            # if passes the 'disorder test', return the seq
            return final_seq

    # if it doesn't work, raise an error
    raise GooseFail('Unable to generate sequence.')

def gen_weighted_shuffle_variant(sequence, target_aas, shuffle_weight, attempts=10, 
    disorder_threshold=parameters.DISORDER_THRESHOLD, strict_disorder=False):
    '''
    function that will let you perform a weighted shuffle a sequence by 
    specifying residues or classes of residues to shuffle and a weight that
    corresponds to the degree of shuffling that you want to perform. 
    The weight is a number between 0.0-1.0 and corresponds to the probability
    of moving a residue during shuffling. If you specify target amino acids, only
    those amino acids are included in the shuffling and weighting can still be 
    applied to only those target amino acids.
    parameters
    ----------
    sequence : str
        the amino acid sequence as a string
    target_aas : str or list
        a list of amino acids to target for shuffling
        or a class of amino acids to target for shuffling
        Possible target classes:
            charged : DEKR
            polar : QNST
            aromatic : FYW
            aliphatic : IVLAM
            negative: DE
            positive : KR
            
    shuffle_weight : float
        a weight between 0.0-1.0 representing the probability of 
        moving a residue during shuffling
    attempts : int
        the number of times to try to make the sequence
    disorder_threshold : float
        the threshold value required for an amino acid
        to be considered disordered
    strict_disorder : Bool
        whether or not to require all disorder values to be 
        over threshold or if it is okay to use the values
        from the input sequence
    '''
    # dict of classes that are possible to choose
    classdict={'charged':['D', 'E', 'K', 'R'], 'polar':['Q', 'N', 'S', 'T'], 'aromatic':
    ['F', 'W', 'Y'], 'aliphatic': ['I', 'V', 'L', 'A', 'M'], 'negative':['D', 'E'], 'positive':['K', 'R']}

    # possible amino acids
    amino_acids = ['A', 'C', 'D', 'E', 'F', 'G', 'H', 'I', 'K', 'L', 'M', 'N', 'P', 'Q', 'R', 'S', 'T', 'V', 'W', 'Y']

    # verify target aas
    if type(target_aas)==str:
        if target_aas in amino_acids:
            raise GooseInputError('You only specified a single amino acid. This will not change your sequence because the amino acids will just change places with itself.')
        elif target_aas in classdict.keys():
            target_aas = classdict[target_aas]
        else:
            raise GooseInputError('The specified target_aas is not a valid amino acid or class of amino acids.')
    elif type(target_aas)==list:
        for i in target_aas:
            if i not in amino_acids:
                raise GooseInputError('The specified target_aas is not a valid amino acid.')
    else:
        raise GooseInputError('The specified target_aas must be type list or string.')

    # get original sequence disorder
    starting_disorder = meta.predict_disorder(sequence)    

    # define probabilities of not relocating a residue and relocating a residue, respectively
    shuffle_weights=[1-shuffle_weight, shuffle_weight]

    # get list of target amino acids from the sequence
    target_aa_list = [aa for aa in sequence if aa in target_aas]

    # attempt to build sequence
    for attempt_num in range(0, attempts):
        # perform weighted sample to determine residues to shuffle
        target_mask = random.choices([False, True], weights=shuffle_weights, k=len(target_aa_list))
        mask = [target_mask.pop(0) if aa in target_aas else False for aa in sequence]

        # gather positions and identities of residues to shuffle
        orig_scramble_positions = [i for i, val in enumerate(mask) if val == True]
        orig_aas = [sequence[i] for i in orig_scramble_positions]

        # perform Fisher-Yates shuffle only with target_aas marked for shuffling
        for i in range(len(orig_aas) - 1, 0, -1):
            remaining_reposition_sites = orig_scramble_positions[:i]

            # randomly select new position for relocation of the target aa
            new_position = random.choice(remaining_reposition_sites)
            new_position_index = remaining_reposition_sites.index(new_position)

            # swap positions with another target aa marked for shuffling
            orig_aas[i], orig_aas[new_position_index] = orig_aas[new_position_index], orig_aas[i]

        # build final seq. Uses shuffled residues at sites marked for repositioning. Otherwise, uses the original residue at that site
        final_seq = ''.join( [orig_aas.pop(0) if i in orig_scramble_positions else aa for i, aa in enumerate(sequence)] )

        # check disorder
        if sequence_variant_disorder(final_seq, starting_disorder, 
            cutoff_val=disorder_threshold, strict=strict_disorder) == True:
            # if passes the 'disorder test', return the seq
            return final_seq

    # if it doesn't work, raise an error
    raise GooseFail('Unable to generate sequence.')


def gen_excluded_shuffle_variant(sequence, exclude_aas, attempts=10, 
    disorder_threshold=parameters.DISORDER_THRESHOLD, strict_disorder=False):
    '''
    function that will let you shuffle a sequence by specifying residues
    or classes of residues to NOT shuffle. This is the opposite behavior of
    the create targeted shuffle variant where you specify which
    residues or classes of residues to target.

    parameters
    ----------
    sequence : str
        the amino acid sequence as a string

    target_aas : str or list
        a list of amino acids to target for shuffling
        or a class of amino acids to target for shuffling
        Possible target classes:
            charged : DEKR
            polar : QNST
            aromatic : FYW
            aliphatic : IVLAM
            negative: DE
            positive : KR

    attempts : int
        the number of times ot try to make the sequence

    disorder_threshold : float
        the threshold value required for an amino acid
        to be considered disordered

    strict_disorder : Bool
        whether or not to require all disorder values to be 
        over threshold or if it si okay to use the values
        from the input sequence
    '''
    # dict of classes that are possible to choose
    classdict={'charged':['D', 'E', 'K', 'R'], 'polar':['Q', 'N', 'S', 'T'], 'aromatic':
    ['F', 'W', 'Y'], 'aliphatic': ['I', 'V', 'L', 'A', 'M'], 'negative':['D', 'E'], 'positive':['K', 'R']}
    
    # possible amino acids
    amino_acids = ['A', 'C', 'D', 'E', 'F', 'G', 'H', 'I', 'K', 'L', 'M', 'N', 'P', 'Q', 'R', 'S', 'T', 'V', 'W', 'Y']

    # verify target aas
    if type(exclude_aas)==str:
        if exclude_aas in amino_acids:
            exclude_aas=[exclude_aas]
        elif exclude_aas in classdict.keys():
            exclude_aas = classdict[exclude_aas]
        else:
            raise GooseInputError('The specified exclude_aas is not a valid amino acid or class of amino acids.')
    elif type(exclude_aas)==list:
        for i in exclude_aas:
            if i not in amino_acids:
                raise GooseInputError('The specified exclude_aas is not a valid amino acid.')
    else:
        raise GooseInputError('The specified exclude_aas must be type list or string.')

    # get original sequence disorder
    starting_disorder = meta.predict_disorder(sequence)  

    # identify exclude aas.
    exclude_aa_list=[]
    building_seq=''

    for curind, aa in enumerate(sequence):
        if aa in exclude_aas:
            building_seq+=aa
        else:
            exclude_aa_list.append(aa)
            building_seq+='0'

    # attempt to build sequence
    for attempt_num in range(0, attempts):
        # copy list
        cur_attempt_list = exclude_aa_list.copy()
        # shuffle list 
        random.shuffle(cur_attempt_list)

        # build final seq
        final_seq=''
        for i in building_seq:
            if i == '0':
                final_seq+=cur_attempt_list.pop()
            else:
                final_seq+=i
        
        # check disorder.
        if sequence_variant_disorder(final_seq, starting_disorder, 
            cutoff_val=disorder_threshold, strict=strict_disorder) == True:
            # if passes the 'disorder test', return the seq
            return final_seq

    # if it doesn't work, raise an error
    raise GooseFail('Unable to generate sequence.')


def gen_dimensions_variant(sequence, increase_or_decrease, rg_or_re, return_all=False, 
    return_all_interval=0.2, disorder_threshold=parameters.DISORDER_THRESHOLD, 
    strict_disorder=False, include_original=False, num_attempts=None):
    '''
    Function to make a sequence more compact or expanded
    based on predicted Rg or Re values. 
    Keeps sequence composition constant. 

    parameters
    -----------
    sequence : str
        The amino acid sequence as a string

    increase_or_decrease : str
        Whether to make the sequence more compact or expanded
        Options are 'increase' or 'decrease'

    rg_or_re :  str
        Whether to alter rg or re
        Options are 'rg' or 're'
    
    return_all : bool
        whether to return all sequences. 

    return_all_interval : float
        the minimal difference between each sequence required 
        to be included in the return_all list. 

    disorder_threshold : float
        The minimum disorder value required for a position

    strict_disorder : Bool
        whether or not to require all disorder values to be 
        over threshold or if it is okay to use the values
        from the input sequence

    include_original : Bool
        whether to inlcude the orignal sequence.
        If set to True, you will get a dictionary with 2 keys:
            'original' and 'variants' where each one will have a dictionary
            of sequences as keys that correspond to values that are the Rg 
            or Re for that sequence

    num_attempts : int
        how many times to attempt to increase or decrease the Rg. 

    Returns 
    --------
    dictionary or nested dictionary:
        Depending on if you want the original sequence or not, will
        return a dictionary or a nested dictionary.

    '''    
    # get starting dimensions
    if rg_or_re=='rg':
        # get the starting dimensions
        starting_dimensions=predict_rg(sequence)
    else:
        starting_dimensions=predict_re(sequence)

    # get original sequence disorder
    starting_disorder = meta.predict_disorder(sequence)
    
    if num_attempts==None:
        num_attempts=len(sequence)
    else:
        if num_attempts<1:
            raise GooseException('cannot have number of attempts be below 1.')

    # iterate over attempts
    for attempt in range(0, num_attempts):
        
        success=False
        try:
            # get all the sequences. 
            seqs_to_dims=make_rg_re_variant(sequence, 
                increase_or_decrease,
                rg_or_re, numseqs=256)
            success=True
        except:
            continue

        # if se made sequences with varying Rg / Re, continue to check for disorder and increase or decrease
        if success==True:

            # predict disorder
            nest_seq_disorder_vals=meta.predict_disorder_batch(list(seqs_to_dims.keys()), 
                                                                show_progress_bar=False)
            # get confirmed disordered seqs
            confirmed_dims_to_seq = {}
            for seq, disorder in nest_seq_disorder_vals:
                if strict_disorder==False:
                    if sequence_variant_disorder(disorder, starting_disorder, 
                                                    cutoff_val=disorder_threshold, 
                                                    strict=False, input_disorder_val=True):
                        confirmed_dims_to_seq[seqs_to_dims[seq]]=seq
                else:
                    if min(disorder) > disorder_threshold:
                        confirmed_dims_to_seq[seqs_to_dims[seq]]=seq

            # make sure we have disordered sequences
            if confirmed_dims_to_seq!={}:

                # get seqs to return. 
                final_dim_vals_sorted = sorted(list(confirmed_dims_to_seq.keys()))

                # set continue_to_return_seqs ==False
                continue_to_return_seqs=False

                if increase_or_decrease=='increase':
                    if final_dim_vals_sorted[-1]>starting_dimensions:
                        continue_to_return_seqs=True
                else:
                    if final_dim_vals_sorted[0]<starting_dimensions:
                        continue_to_return_seqs=True

                if continue_to_return_seqs==True:

                    # dict for final variants
                    final_seqs={}
                    if return_all==True:
                        final_seqs[confirmed_dims_to_seq[final_dim_vals_sorted[0]]]=final_dim_vals_sorted[0]
                        cur_val=final_dim_vals_sorted[0]
                        for dim_val in range(1, len(final_dim_vals_sorted)):
                            cur_dim = final_dim_vals_sorted[dim_val]
                            if abs(cur_dim-cur_val)> return_all_interval:
                                final_seqs[confirmed_dims_to_seq[final_dim_vals_sorted[dim_val]]]=final_dim_vals_sorted[dim_val]
                                cur_val=final_dim_vals_sorted[dim_val]

                        # make sure we get the whole range. 
                        if final_dim_vals_sorted[-1] not in final_seqs:
                            final_seqs[confirmed_dims_to_seq[final_dim_vals_sorted[-1]]]=final_dim_vals_sorted[-1]

                    else:
                        if increase_or_decrease=='increase':
                            final_seqs[confirmed_dims_to_seq[final_dim_vals_sorted[-1]]]=final_dim_vals_sorted[-1]
                        else:
                            final_seqs[confirmed_dims_to_seq[final_dim_vals_sorted[0]]]=final_dim_vals_sorted[0]

                    # add original if needed. 
                    if include_original==True:
                        return{'original':{sequence:starting_dimensions}, 'variants':final_seqs}
                    else:
                        return final_seqs

    raise GooseFail('Unable to generate sequence.')
