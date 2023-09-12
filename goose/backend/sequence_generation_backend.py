import random
import math
from random import randint

from goose.backend import lists
from goose.backend.protein import Protein
from goose.backend.amino_acids import AminoAcid
from goose.goose_exceptions import GooseError, GooseInputError
from goose.backend import parameters



def optimal_residue_key(four_amino_acids):

    """

    Function to translate amino acids in a sequence to the proper key for
    the dis_value dict. The dis_value dict holds the predicted disorder 
    values for thousands of distinct 5 amino acid combinations. To keep
    this dict from being so massive that it slows down GOOSE, the dict
    did not use every possible amino acid in the 4 amino acids that are
    used for prediting potential disorder values in the 5th amino acid.
    For example, all aromatics were treated as equivalent. This reduced
    the total lenth of the dict by several orders of magnitude.


    Parameters
    -------------

    four_amino_acids : String
        The 4 amino acids used from the input sequence that 
        are to be used for prdicting disorder


    Returns
    ---------

    key : String
        Returns 4 amino acids that coorespond to the amino acids used
        for generating the dis_value dict.

    """

    # make an empty string to hold the 4 amino acids
    key = ""

    # for the amino acids in the input amino acids, change each one
    # to the corresponding amino acid used for the dis_val_dict
    for i in four_amino_acids:
        if i == "F" or i == "W" or i == "Y":
            key += "W"
        elif i == "C":
            key += "C"
        elif i == "L" or i == "V" or i == "I":
            key += "L"
        elif i == "M":
            key += "M"
        elif i == "A":
            key += "A"
        elif i == "H":
            key += "H"
        elif i == "K" or i == "R":
            key += "K"
        elif i == "Q" or i == "N":
            key += "Q"
        elif i == "D" or i == "E":
            key += "D"
        elif i == "G" or i == "S":
            key += "G"
        elif i == "T":
            key += "T"
        elif i == "P":
            key += "P"
        else:
            continue
    
    # return the final 4 amino acids that are the key to the 
    # dis_val_dict
    return key


def get_optimal_residue(four_amino_acids, exclude_residues = [], cutoff_disorder = None, return_all=False):

    """

    Clever shortcut that allows bypassing predicting disorder
    of each amino acid for ever amino acid to be added to a sequence
    when generating a sequence. 

    Basically, a giant dict where 4 amino acids at the end of an IDR 
    were varied and then followed by a 5th amino acid. The 5th 
    amino acid tested was every possible standard amino acid. 
    The disorder value of the 5th amino acid was added to a
    list with known orders of amino acids (alphabetical).
    Now, that dict can be accessed to get approximate
    disorder values for potential amino acids to follow any 
    4 amino acids from a sequence being generated. This is 
    much faster than iteratively predicting sequence disorder 
    and is surprisingly accurate despite the massive assumptions that 
    had to be taken when making the dict.

    Importantly, this function will return a random amino acid from the
    list as opposed to the one with the highest value which gives
    GOOSE some stochasticity in sequence generation (generatiing the
    same sequence every time isn't all that useful...).

    Parameters
    -------------

    four_amino_acids : String
        The 4 amino acids used from the input sequence that 
        are to be used for prdicting disorder.

    exclude_residues : List
        List of residues to be excluded from possible residues to be
        returned.
    
    cutoff_dis_val : Float
        The cutoff value to be used for considering something as disordered

    return_all : Bool
        Whether to return all candidate residues

    Returns
    ---------

    String
        A single amino acid as a string.

    """

    # Need this if statement to avoid accidental infinite
    # loops due to my while statement downstream.
    if len(exclude_residues) >= 20:
        raise GooseInputError("You cannot exclude all amino acids.")

    # set order of amino acids for potential_AA_vals 
    # This order matches what was used for the input 5th
    # amino acid during generation of the aa_dis_val_4_V2 dict
    ordered_amino_acids = ['A', 'C', 'D', 'E', 'F', 'G', 'H', 'I', 'K', 'L', 'M', 'N', 'P', 'Q', 'R', 'S', 'T', 'V', 'W', 'Y']
    
    # make empty list to hold potential resiudes if more than
    # one would in theory work.
    potential_residue_numbers = []
    
    # translate the sequence to keys that are in the dict
    four_amino_acids_key = optimal_residue_key(four_amino_acids)
    
    # get vals from dict
    potential_AA_vals = lists.aa_dis_val_4_v3[four_amino_acids_key]
    
    # sort the values from highest to lowest
    potential_residue_numbers = sorted(potential_AA_vals, reverse=True)
    
    # make empty list to hold potential amino_acids
    candidate_amino_acids = [] 

    # adjust the cutoff value for the dis_val list
    # Note - this was empiracally determined based on disorder vals
    # in the dis_val list.

    # using the cutoff value equal to the average predicted score in the dict + 1 standard deviation.
    # average disorder in dict = 0.6397844593942738, stdev of scores in dict = 0.07225962272855471
    # THIS IS DIFFERENT THAN parameters.DISORDER_THRESHOLD due to how the precomputed dict was made.
    if cutoff_disorder == None:
        cutoff_disorder = 0.7120440821228284


    # make sure cutoff dis val doesn't get too high
    if cutoff_disorder > 0.95:
        cutoff_disorder = 0.95

    # iterate through potential amino acids
    for i in range(0, len(potential_AA_vals)):
        # get the predicted disorder value for that amino acid
        predicted_value = potential_AA_vals[i]
        # if that value is greater than the cutoff
        if predicted_value > cutoff_disorder:
            # Get the corresponding amino acid
            corresponding_amino_acid = ordered_amino_acids[i]
            # if the amino acid is not supposed to be excluded...
            if corresponding_amino_acid not in exclude_residues:
                # add that amino acid to the candidate amino acid list
                candidate_amino_acids.append(corresponding_amino_acid)

    # if the candidate amino acids list is still empty, need to do something else...
    if candidate_amino_acids == []:
        # setting arbitrary index value to iterate through the potential residue numbers
        amino_acid_index = 0
        # while we don't yet have the best possible amino acid under 0.7...
        while candidate_amino_acids == []:
            # figure out what index value corresponds to the highest current value
            current_index_value = potential_residue_numbers[amino_acid_index]
            # figure out where the residue is in the original unsorted list
            current_residue_position = potential_AA_vals.index(current_index_value)
            
            # get the corresponding amino acid. Amino acids always in order
            # for list ordered_amino_acids, so just call index as is
            corresponding_amino_acid = ordered_amino_acids[current_residue_position]        
            # if that amino acid is not to be excluded...
            if corresponding_amino_acid not in exclude_residues:
                # add it to the candidate list
                candidate_amino_acids.append(corresponding_amino_acid)
            # go to next amino acid index
            amino_acid_index += 1
            if amino_acid_index==20:
                for amino_acid in ordered_amino_acids:
                    if amino_acid not in exclude_residues:
                        candidate_amino_acids.append(amino_acid)


    # choose a random amino acid from the list to return
    if return_all == False:
        return candidate_amino_acids[randint(0, len(candidate_amino_acids)-1)]
    else:
        return candidate_amino_acids



#function that returns a random amino acid from a specified list.
def random_amino_acid(seq_list):

    """

    Function that returns a random amino acid from
    a list of amino acids.
    

    Parameters
    -------------

    seq_list : List
        List of amino acids from which to select a random residue.


    Returns
    ---------
    
    String
        Returns a single amino acid from the selected list as a string.

    """

    return seq_list[randint(0, len(seq_list)-1)]


def fast_predict_disorder(sequence):
    '''
    a fast disorder predictor that uses
    pre-computed values to guestimate the disorder
    across the sequence

    Parameters
    ----------
    sequence: String
        The amino acid sequence as a string

    Returns
    -------
    disorder_values : List
        Returns the disorder values on a residue-by-residue basis as a list
        of float values.
    '''
    # make list to hold disorder values
    disorder_values = []
    # amino acids in the order used to generate aa_dis_val_4_v3_dict
    amino_acids = ['A', 'C', 'D', 'E', 'F', 'G', 'H', 'I', 'K', 'L', 'M', 'N', 'P', 'Q', 'R', 'S', 'T', 'V', 'W', 'Y']
    # iterate through sequence
    for i in range(0, len(sequence)):
        # first 3 values need to be tweeked such that we can get disorder values back
        # this is because the 'pseudo predictor' requries 4 residues to get the predicted
        # disorder for the following residue.
        beginning = [0, 1, 2, 3]

        # if the value is in the first 4 amino acids, just use the first 4 as the 'base'
        if i in beginning:
            cur_sequence = sequence[0:4]
        else:
            # if not first 4 amino acids, take the 4 amino acids
            # preceeding the next amino acid to get the optimal
            # residue key to use in the aa_dis_val_4_v3_dict dictionary
            cur_sequence = sequence[i-4:i]
        
        # now input the base (4 amino acicds) as cur_sequence to get the 
        # corresponding key for the dict holding pseudo disorder values
        cur_sequence = optimal_residue_key(cur_sequence)
        # get the disorder values using the key from the aa_dis_val_4_v3 dict
        all_disorder_values=(lists.aa_dis_val_4_v3[cur_sequence])
        # get the next residue in the sequence
        next_residue = sequence[i]
        # figure out the index value for that residue
        index_value = amino_acids.index(next_residue)
        # get the corresponding disorder value
        disorder_value = all_disorder_values[index_value]
        # append the disorder value to the growing list of disorder values
        disorder_values.append(disorder_value)
    # return a list of disorder values corresponding to each amino acid
    return disorder_values



def identify_residue_positions(sequence, residue):

    """

    Function to identify the positions of any given residue
    in an input sequence. Used for generating sequence variants
    where a specific residue is to be held constant.


    Parameters
    -------------

    sequence : String
        The sequence that is being examined for the specified residue.

    Returns
    ---------

    Residue_cooordinates : list 
        Returns a list holding the coordinates of the specified residue
        within the input sequence.

    """

    # make sure residue that was input is uppercase since this
    # equivalence is case-sensitive in Python
    residue = residue.upper()

    # create empty list to store coordinates of residue in seuqence
    residue_coordinates = []

    # iterate over the amino acids in the sequence
    for AA in range(0, len(sequence)):
        # if the current residue is the specified residue
        if sequence[AA] == residue:
            # append the corresponding location to the
            # residue_coordinates list
            residue_coordinates.append(AA)

    # return the reisude coordinates list            
    return residue_coordinates



def optimize_once(sequence):
    '''

    takes in a sequence and re-orders the amino acids to try to 
    maximize the likelihood that the resulting sequence is disordered

    Parameters
    -------------
    sequence : String
        The sequence that is being examined for the specified residue.

    Returns
    ---------
    best_sequence : String 
        The best sequence as far as optimization

    '''
    # get pseudo disorder
    sequence_disorder = fast_predict_disorder(sequence)
    # find worst residue based on pseudo disorder predicotor
    worst_residue_index = sequence_disorder.index(min(sequence_disorder))
    # get worst residue
    worst_residue = sequence[worst_residue_index]
    # figure out all amino acids present in the sequence
    all_possible_amino_acids = set([char for char in sequence])
    # now make list of sequences where the amino acids are switched between the worst and other amino acids
    test_sequence = []

    # now for all the possible amino acids
    for i in all_possible_amino_acids:
        # find their locations
        possible_indices = identify_residue_positions(sequence, i)
        if len(possible_indices)>1:
            cur_index = possible_indices[randint(0, len(possible_indices)-1)]
        else:
            cur_index = sequence.index(i)
        # figure out if the cur res is not the residue trying to be swapped
        cur_res = sequence[cur_index]
        if cur_res != worst_residue:
            if worst_residue_index < cur_index:
                cur_sequence = sequence[0:worst_residue_index] + cur_res + sequence[worst_residue_index+1:cur_index] + worst_residue + sequence[cur_index+1:] 
            else:
                cur_sequence = sequence[0:cur_index] + worst_residue + sequence[cur_index+1:worst_residue_index] + cur_res + sequence[worst_residue_index+1:]
            # add the sequence to the test sequence list
            test_sequence.append(cur_sequence)

    # set the best disorder value to the minimume of seq disorder
    # because we only want stuff better than that.
    best_disorder_value = min(sequence_disorder)
    
    # set the best sequence to the input because if we cnan't make anyhting
    # better we will return the input sequence as is
    best_sequence = sequence

    # for each sequence to be tested
    for i in test_sequence:
        # get it's fast disorder
        cur_disorder = fast_predict_disorder(i)
        # figure out the current disorder at the worst residue position
        cur_dis_at_worst = cur_disorder[worst_residue_index]
        # if it's better than the best disorder value replace it
        if cur_dis_at_worst > best_disorder_value:
            best_sequence = i
            best_disorder_value = cur_dis_at_worst

    # return the best sequenc
    return best_sequence



def optimize_sequence(sequence, iterations=100):
    '''
    Uses optimize_once function iteratively to fully optimize a sequence.
    Iterations is the number of optimization attempts before just calling it
    so that this doesn't take too much time.

    Parameters
    ----------
    sequence : String
        The sequence to be optimized as a string

    iterations :int
        The number of iterations of optimize_once to try before just calling it
    
    Returns
    -------
    new_sequence : String
        A new sequence that should have a better chance of being disordered

    '''
    # keep track of list of sequences. That way if we start getting the same sequence
    # again we can kill the optimization
    already_used = [sequence]

    # keep track if repeat sequence is made. If the sequence is repeated than we have
    # more or less optimized as much as we can.
    repeat_sequence = False

    # keep track of iterations
    cur_iter = 0

    # iteratively optimized sequecnce
    while repeat_sequence == False:
        new_sequence = optimize_once(sequence)
        if already_used.count(new_sequence) > 1 or cur_iter == iterations:
            # if we are at the max number of iterations or we have 
            # the same sequence in the list twice, return it and kill
            # the optimization
            return new_sequence
            repeat_sequence = True
        # add the generated sequence to the list of sequences already generated
        # to keep track of possible duplicates
        already_used.append(new_sequence)
        sequence = new_sequence
        # update cur_iter
        cur_iter += 1



def shuffle_seq(seq):

    """

    Function to shuffle an input sequnce


    Parameters
    -------------

    seq : String
        The sequence that is being shuffled

    Returns
    ---------

    String 
        Returns a shuffled version of seq

    """

    return "".join(random.sample(seq, len(seq)))


def random_optimization(sequence, min_random_iterations = 100, min_is_max = False):
    '''
    function that I'm testing to see if randomly choosing 2 amino acids
    to swap and then checking the 'fast disorder' score to see if it is
    greater than the original is faster than doing anything rationally.
    
    Parameters
    ----------
    sequence : String
        The sequence to be optimized as a string

    min_random_iterations :int
        minimum number of iterations of optimize_once to try before just calling it

    min_is_max: bool
        if set to True, then min_random_iterations is the max number of iterations.


    Returns
    -------
    new_sequence : String
        A new sequence that should have a better chance of being disordered
    
    '''
    # keep track of best sequence
    best_min_disorder = min(fast_predict_disorder(sequence))
    # set best sequences = sequence in case nothing better is made
    best_sequence = sequence
    # adjust random iterations based on seq length
    input_random_iterations = int(len(sequence)/3)
    # adjust number itartions based on seq size if min_is_max is not true
    if min_is_max == False:
        if input_random_iterations > min_random_iterations:
            random_iterations = input_random_iterations
        else:
            random_iterations = min_random_iterations
    else:
        random_iterations = min_random_iterations

    # generate some random seqs
    for i in range(0, random_iterations):
        generated_sequence = shuffle_seq(sequence)
        cur_min_disorder = min(fast_predict_disorder(generated_sequence))
        if cur_min_disorder > best_min_disorder:
            best_sequence = generated_sequence
            best_min_disorder = cur_min_disorder

    # return the best sequence
    return best_sequence

def check_hydropathy(sequence, objective_hydropathy, hydro_error = parameters.HYDRO_ERROR):
    '''
    function to check the hydropathy of a sequence and see if it is 
    within the appropriate error (hydro_error) or not

    parameters
    -----------
    sequence : String
        The sequence to be optimized as a string
    objective_hydropathy : float
        objective hydropathy value as a float
    hydro_error : float
        the max allowed error in the hydropathy between the sequence and the 
        objective hydrpoatyh

    '''
    cur_hydro = round(Protein(sequence).hydropathy, 5)
    error = abs(cur_hydro-objective_hydropathy)
    if error <= hydro_error:
        return True
    return False


def partially_random_optimization(sequence, nonrandom_iterations = None, random_iterations = None):
    '''
    Function combining the random optimizer (which worked remarkably well)
    and the original optimizer (that uses the pseudo disorder predictor)
    
    Parameters
    ----------
    sequence : String
        The sequence to be optimized as a string

    nonrandom_iterations :int
        number of iterations of optimize_once to try before just calling it

    random_iterations: int
        Number of iterations of random_optimizer to try before just calling it


    Returns
    -------
    new_sequence : String
        A new sequence that should have a better chance of being disordered
    
    '''

    # carry out random optimzations
    if random_iterations == None:
        input_random_seq = random_optimization(sequence)
    else:
        input_random_seq = random_optimization(sequence, min_random_iterations=random_iterations, min_is_max=True)

    # now carry out non random optimzations
    if nonrandom_iterations == None:
        number_nonrandom_iterations = int(len(sequence)/10)           
        if number_nonrandom_iterations < 50:
            number_nonrandom_iterations = 50
        optimized_sequence = optimize_sequence(input_random_seq, iterations=number_nonrandom_iterations)
    else:
        optimized_sequence = optimize_sequence(input_random_seq, iterations=nonrandom_iterations)

    # return the optimized sequence
    return optimized_sequence



def gen_charged_positions(length, num_charged):
    '''
    function to generate random positoning for charged residues.
    Used in sequence generation

    Parameters
    ----------
    length : Int
        length of the sequence

    Returns
    -------
    num_charged : Int
        number of charged residues to be added to the sequence being
        generated
    '''

    pseudo_seq = ''
    for i in range(0, length-num_charged):
        pseudo_seq += '0'
    for i in range(0, num_charged):
        pseudo_seq += '+'
    pseudo_seq = shuffle_seq(pseudo_seq)
    charged_positions = identify_residue_positions(pseudo_seq, '+')
    return charged_positions


def gen_sequence(length, usedlist=[]):
    '''
    Function to generate a sequence from a list of possible
    amino acids (weighted or otherwise). 

    Parameters
    ----------
    length : int
        the length of the sequence being generated

    used_list : list
        A list of amino acids that is to
        be used when generating the sequence

    Returns
    -------
    final_sequence : String
        Returns the final sequence as a string

    '''

    final_sequence = ''
    if usedlist == []:
        usedlist = ['A', 'C', 'D', 'E', 'F', 'G', 'H', 'I', 'K', 'L', 'M', 'N', 'P', 'Q', 'R', 'S', 'T', 'V', 'W', 'Y']
    for i in range(0, length):
        final_sequence += random_amino_acid(usedlist)
    return final_sequence




def all_excluded_residues_hydro(sequence, objective_hydropathy, no_charge=False, input_exclusion=[]):

    """
    
    Function to identify which residues to exclude when trying to
    optimize hydropathy. Basically looks at the residues in a sequence
    as well as the objective hydropathy and then based on the current
    hydropathy of the input sequence will tell which residues should not
    be considered for use downstream.

    Parameters
    -------------

    sequence : String
        The input sequence examine for residues to exclude

    objective_hydropathy : Float
        The objective hydropatyh value of the final sequence

    no_charge : Bool
        Whether or not charged residues should be excluded

    input_exclusion : List
        A list of additional residues to exclude from potential residues


    Returns
    ---------
        
    exclude_amino_acids : List
        A list of amino acids to exclude from downstream use

    """

    # make an empty list to hold residues to exclude as potential residues to use
    exclude_amino_acids = []

    # add any input excluded residues to the excluded list
    for i in input_exclusion:
        if i not in exclude_amino_acids:
            exclude_amino_acids.append(i)

    # if the objective hydro is less than the current hydro add residues greater than objective to list
    if Protein(sequence).hydropathy > objective_hydropathy:
        for i in lists.amino_acids:
            if AminoAcid.hydro(i) > objective_hydropathy:
                if i not in exclude_amino_acids:
                    exclude_amino_acids.append(i)

    # if the objective hydro is greater than the current hydro add residues less than objective to list
    if Protein(sequence).hydropathy < objective_hydropathy:
        for i in lists.amino_acids:
            if AminoAcid.hydro(i) < objective_hydropathy:
                if i not in exclude_amino_acids:
                    exclude_amino_acids.append(i)

    # if the list must include charged residues, make sure they are in the list 
    if no_charge == True:
        charged_list = ['D', 'E', 'K', 'R']
        for i in charged_list:
            if i not in exclude_amino_acids:
                exclude_amino_acids.append(i) 

    # if the list includes every single residue, find the least bad residue 
    if len(exclude_amino_acids) == 20:
        best_possible_value = 10
        best_amino_acid = ""
        for i in lists.amino_acids:
            if no_charge==False:
                current_amino_acid = i
            else:
                if i not in charged_list:
                    current_amino_acid = i
            aa_value_difference = abs(AminoAcid.hydro(current_amino_acid) - objective_hydropathy)
            if aa_value_difference < best_possible_value:
                best_possible_value = aa_value_difference
                best_amino_acid = current_amino_acid
        # remove the least bad residue from the list
        exclude_amino_acids.pop(exclude_amino_acids.index(best_amino_acid))

    # return the list of amino acids to be excluded    
    return exclude_amino_acids



def optimize_hydro(sequence, final_hydropathy, use_charged_residues=False, cutoff_dis = 0.7120440821228284, excluded_residues=[]):
    """
    
    Function to optimize the hydropathy of an input sequence to be
    closer to the objective hydropathy. This function is very fast
    and prioritizes getting hydropathy correct. However, it isn't the
    best at making sure the sequence is disordered. It is intended as a
    first step to get a sequence that is being generated close to the
    final sequence but the sequence output by this function typically
    needs some disorder optimizations.

    Parameters
    -------------

    sequence : String
        The input sequence to optimize

    final_hydropathy : Float
        The objective hydropathy value of the final sequence

    use_charged_residues : Bool
        Whether or not charged residues should be used. For
        the functions that generate sequences with combined 
        FCR / NCPR and hydro, this allows charged residues to 
        be avoided to keep the sequence from having incorrect
        FCR / NCPR due to changes during optimization.

    cutoff_dis : Float
        The disorder value to be used as a cutoff. 0.82 was empiraclly
        determined to give best results.

    excluded_residues : List
        List of residues to exclude from possible residues to use
        during the optimization.


    Returns
    ---------
        
    final_sequence : String
        A sequence with hydropathy value closer to the objective
        hydropathy.

    """

    #set current_hydro equal to the current hydropathy of the sequence
    current_hydro = Protein(sequence).hydropathy
    
    #set worst value to 0
    value_coordinate = 0
    
    # determine whether or not to stop the optimzation
    if abs(current_hydro - final_hydropathy) < 0.01:
        stop=True
    else:
        stop=False
    
    # make initial list of residues to exclude based on those input
    exclude_these_residues = excluded_residues
    
    #set values for charged residues
    charged_list = ['D', 'E', 'K', 'R']
    
    # add charged residues to exclude_these_residues if necessary
    if use_charged_residues == False:
        for residue in charged_list:
            if residue not in exclude_these_residues:
                exclude_these_residues.append(residue)

    # if current hydropathy is too high, find amino acids
    # to change that will lower it
    if current_hydro > final_hydropathy:
        worst_value = -100
        for amino_acid_index in range(3, len(sequence)):
            if sequence[amino_acid_index] not in exclude_these_residues:
                cur_hydro_value = Protein(sequence[amino_acid_index]).hydropathy
                if cur_hydro_value > worst_value:
                    worst_value = cur_hydro_value
                    value_coordinate = amino_acid_index

    # if current hydropathy is too high, find amino acids
    # to change that will raise it
    elif current_hydro < final_hydropathy:
        worst_value = 100
        for amino_acid_index in range(3, len(sequence)):
            if sequence[amino_acid_index] not in exclude_these_residues:
                cur_hydro_value = float(Protein(sequence[amino_acid_index]).hydropathy)
                if cur_hydro_value < float(worst_value):
                    worst_value = cur_hydro_value
                    value_coordinate = amino_acid_index

    # just another check to make sure optimization does not occur if hydropatyh
    # of the sequence is not greater than or less than the final_hydropathy
    else:
        stop == True

    # if optimization is not to be stopped..
    if stop == False:
        # if the value coordinate is 0, we won't be able to get an optimal residue.
        # so we will have to simply return the sequence.
        if value_coordinate == 0:
            final_sequence = sequence

        else:
            #figure out what residues need to be excluded
            exclude_vals = all_excluded_residues_hydro(sequence=sequence,
                        objective_hydropathy=final_hydropathy, no_charge=True,
                         input_exclusion=excluded_residues)
            
            # make sure all excluded_residues are in exclude_vals
            for AA in excluded_residues:
                if AA not in exclude_vals:
                    exclude_vals.append(AA)

            #figure out what amino acids precede the worst value coordinate
            optimal_key = sequence[int(value_coordinate-3):int(value_coordinate + 1)]
            
            # get best residue based on the amino acid chosen to change and the 
            # residues that are to be excluded
            best_residue = get_optimal_residue(optimal_key, exclude_vals)
            
            # put the sequence back together
            seq_part1 = sequence[0:value_coordinate]

            # if value coordinate is the very end of the sequence, there is
            # nothing to add after the best_residue
            if value_coordinate == len(sequence)-1:
                seq_part_2 = ""
            # figure out seq_part_2 such that it avoids the residue to be removed
            else:
                seq_part_2 = sequence[value_coordinate+1:]

            # make the final sequence
            final_sequence = seq_part1 + best_residue + seq_part_2
    
    # if stop == true, set the final sequence to the input sequence
    else:
        final_sequence = sequence

    # return the final sequence
    return final_sequence



def hydro_seq(length, mean_hydro, just_neutral=False, allowed_error=None, return_best_seq = False, exclude_residues = []):
    """
    This will return a protein sequence with a specified length and
    mean hydropathy. Separate from the disordered generation because it
    can be tricky to generate sequences with specific hydropathy values,
    so it's nifty to keep this separate.

    Parameters
    ----------
    length : Int
        The length of the sequence to be generated as an integer value

    mean_hydro : Float
        The meany hydropathy value wanted for the output sequence

    just_neutral : Bool
        Whether or not the generated sequence is allowed to contain charged values.
        Useful for generating sequences with specified hydropathies that have specified
        NCPR or FCR values.

    allowed_error : Float
        The allowed error between the input value for mean_hydro and the resulting 
        mean_hydro of the final sequence.

    return_best_seq : Bool
        Whether or not to just return the best sequence (closest hydropathy)
        that was managed regardless if it is within the allowed error. This
        allows for downstream optimization.

    exclude_residues : List
        A list of residues not allowed to be used for sequence generation

    Returns
    -------
    final_sequence : String
        Returns a string that is the final amino acid sequence
        having the specified length and mean_hydro value that has
        a mean_hydro value within the value of allowed_error.

    """ 

    #  Choose weighted list for seq generation
    #--------------------------------------------#
    # if no customized allowed error, set to 0.05
    if allowed_error == None:
        allowed_error = 0.05

    # set best_error to stupidly high error
    best_error = 10000

    #convert mean_hydro into a string.
    rounded_value = str(round(mean_hydro, 1))
    len_hydro = len(rounded_value)

    # use value to get proper dict key for weighted list
    if len_hydro > 1:
        final_value = "{}_{}".format(rounded_value[0], rounded_value[2])
    else:
        final_value = "{}_0".format(rounded_value[0])
    
    #if mean hydro is any value greater than the max that I can make disordered (6.1).
    if just_neutral == False:
        if mean_hydro > 6.1:
            used_list = lists.HydroDict['Hydro_dis_6_1']
        else:
            dict_key = "Hydro_dis_"+str(final_value)
            used_list = lists.HydroDict[dict_key]
    else:
        '''
        if just_neutral == True then could be input for 
        generating specific hydro, fcr, ncpr sequence
        '''

        if mean_hydro > 9.0:
            used_list = lists.NeutralHydroDict['Neutral_hydro_dis_9_0']
        elif mean_hydro <= 1.0:
            used_list = lists.NeutralHydroDict['Neutral_hydro_dis_1_0']
        else:
            dict_key = "Neutral_hydro_dis_"+str(final_value)
            used_list = lists.NeutralHydroDict[dict_key]        

    #  Start attempts to build the sequence
    #--------------------------------------------#
    # set correct hydro to false
    correct_hydro = False
    
    # keep track of iterations
    iters = 0

    if exclude_residues != []:
        excluded_list = []
        for i in used_list:
            if i not in exclude_residues:
                excluded_list.append(i)
        used_list = excluded_list
        if used_list == []:
            raise GooseInputError('The function hydro_seq in /backend/sequence_generation_backend.py is attempting to build a sequence with an empty list due to specified excluded residues.')

    # while the correct hydro is false, try to generate sequence with correct hydro
    while correct_hydro==False:
        # make a sequence
        current_sequence = gen_sequence(length, usedlist=used_list)
        # figure out hydropathy
        current_hydropathy = round(Protein(current_sequence).hydropathy, 4)
        # figure out current error
        cur_error = abs(mean_hydro - current_hydropathy)
        # see if it matches mean_hydro within allowed_error
        if cur_error <= allowed_error:
            # set final sequence
            final_sequence = current_sequence
            # return the sequence
            return final_sequence
            # set correct_hydro to True
            correct_hydro == True
        else:
            if cur_error < best_error:
                best_error = cur_error
                best_sequence = current_sequence

        iters += 1

        if iters == 30000:
            if return_best_seq == True:
                return best_sequence
            else:
                raise GooseError('Unable to generate sequence with correct hydropathy value.')
                # set correct_hydro to True
                correct_hydro == True


def generate_charged_residues(length, FCR, objective_hydropathy):
    '''
    function to try to help the FCR function 
    make a charged list of residues that can actually
    result in a hydropathy value that is possible to
    generate to balance the final objective hydropathy

    parameters
    ----------
    length : int
        the length of the sequence

    FCR : float
        FCR fraction value as a float between 0 and 1

    objective_hydropathy : float
        the objective hydropathy value for the sequence.
        helps determine number of each charged residue if the
        hydropathy of the charged residues differ.

    returns
    -------
    final_charged_res : string
        returns the final_charged_res as a string
    '''
    # empty string to hold charged residues
    final_charged_res = ''

    # hydro values for the charged residues
    charged_dict = {'R':0.0, 'K':0.6, 'D':1.0, 'E':1.0}

    # first need to figure out the ideal charge value for the charged residues
    total_hydro = length * objective_hydropathy
    number_charged = round(FCR*length)
    number_non_charged = length-number_charged
    # because the lowest non-charged hydro is 1, we can figure out what value we need
    # the charged residues to get to by subtracting 1*number_non_charged
    total_hydro_charged = total_hydro - number_non_charged 
    average_hydro_val = total_hydro_charged / number_charged
    D_E = ['D', 'E','D', 'E']
    # if total_hydro_charged < 1, just add 'R'
    if average_hydro_val < 0:
        for i in range(0, number_charged):
            final_charged_res += 'R'
        return final_charged_res
    # if total_hydro_charged > 1, just add 'D' or 'E'
    elif average_hydro_val > 1:
        for i in range(0, number_charged):
            final_charged_res += random_amino_acid(D_E)
        return final_charged_res
    # otherwise we need to figure out how to get as close as possilbe.
    else:
        # if the value is over 0.6, need a combination of D, E, and K
        if average_hydro_val > 0.6:
            final_charged_res = 'K'
            remaining_hydro = total_hydro_charged - 0.6
            for i in range(0, number_charged-1):
                cur_average_hydro = Protein(final_charged_res).hydropathy
                if remaining_hydro == 0:
                    final_charged_res += 'R'
                elif remaining_hydro == 0.6:
                    final_charged_res += 'K'
                elif remaining_hydro == 1:
                    final_charged_res += random_amino_acid(D_E)
                else:
                    if cur_average_hydro < average_hydro_val:
                        final_charged_res += random_amino_acid(D_E)
                        remaining_hydro -= 1
                    else:
                        final_charged_res += 'K'
                        remaining_hydro -= 0.6
        else:
            final_charged_res = 'R'
            remaining_hydro = total_hydro_charged
            for i in range(0, number_charged-1):
                cur_average_hydro = Protein(final_charged_res).hydropathy
                if remaining_hydro == 0:
                    final_charged_res += 'R'
                elif remaining_hydro == 0.6:
                    final_charged_res += 'K'
                elif remaining_hydro == 1:
                    final_charged_res += random_amino_acid(D_E)
                else:
                    if cur_average_hydro < average_hydro_val:
                        final_charged_res += random_amino_acid('K')
                        remaining_hydro -= 0.6
                    else:
                        final_charged_res += 'R'
        return final_charged_res



def sigma_FCR_NCPR(length, sigma_value):
    """
    Returns NCPR and FCR values for a specific sigma value

    Parameters
    ------------
    length : Int
        length of desired disordered sequence

    sigma_value : Float
        The sigma_value as a decimal

    Returns
    -----------
    dict
        returns a dict where:
        {'FCR': FCR_value, 'NCPR': NCPR_value}
    """      
    # calculate interval for potential FCRs
    FCR_fraction = 1/length
    # figure out all possible FCR values based off of length of seq, FCR can't be zero.
    all_FCR_values = []
    start_val = FCR_fraction
    while start_val < 1:
        all_FCR_values.append(start_val)
        start_val = start_val + FCR_fraction

    # figure out all possible FCR NCPR combinations based off of sigma_value
    # first use more stringent criteria first go. If nothing comes back, add 
    # all possible values.
    # make empty dict to hold FCR_NCPR combinations
    possible_FCR_combos = {}
    # make empty list to hold values that can be called randomly later
    possible_FCR_values = []
    #make empty dict to hold best value pair
    best_FCR_NCPR_pair = {}

    # make arbitrary best value to overwrite later
    best_difference = 10
    # for all the FCR values possible
    for i in all_FCR_values:
        # determine possible NCPR values
        possible_NCPR_Value = math.sqrt(i * sigma_value)
        # make sure the possible NCPR value is compatible with the FCr value
        if i >= possible_NCPR_Value:
            potential_error = abs((int(possible_NCPR_Value * length)) - (possible_NCPR_Value * length))
            # if the potential error is the lowest so far...
            if potential_error < best_difference:
                # First clear the dict
                best_FCR_NCPR_pair = {}
                # Then add the pair
                best_FCR_NCPR_pair[i] = possible_NCPR_Value
                # make the best value the new value
                best_difference = potential_error
                # save the best_possible_FCR value for later
                best_possible_FCR = i
            # Make sure the value pair will return a close sigma value
            if abs((int(possible_NCPR_Value * length)) - (possible_NCPR_Value * length))< 0.1:
                # add entry to dict and append to list
                possible_FCR_combos[i] = possible_NCPR_Value
                possible_FCR_values.append(i)

    # if there were no possible combos within the error threshhold
    if possible_FCR_combos == {}:
        # use the best combo
        possible_FCR_combos = best_FCR_NCPR_pair

    # if there are values in the possible values list
    if possible_FCR_values != []:
        #choose a random possilbe FCR value as the FCR value to use
        FCR_value = possible_FCR_values[random.randint(0, len(possible_FCR_values)-1)]

        # calculate the FCR_value based off of the chosen random NCPR Value and sigma_value
        NCPR_value = possible_FCR_combos[FCR_value]

    #if there were no vlaues within the error threshhold..
    else:
        NCPR_value = best_FCR_NCPR_pair[best_possible_FCR]
        FCR_value = best_possible_FCR

    #randomly decide if NCPR value will be negative or positive
    negative_or_not = random.randint(0, 1)

    if negative_or_not == 1:
        NCPR_value = NCPR_value * -1

    # return ncpr and FCR as dict
    return {'FCR': FCR_value, 'NCPR': NCPR_value}


def K_R_optimization(sequence, objective_hydropathy):
    '''
    changes K to R or R to K in sequence to get closer to objective hydropathy

    Parameters
    ----------
    sequence : String
        The input amino acid sequence as a string
    
    objective_hydropathy : Float
        The hydropathy value intended for the final sequence


    Returns
    -------
    final_seq : String
        The final amino acid sequence as a string

    '''
    seq_hydropathy = Protein(sequence).hydropathy
    
    if objective_hydropathy < seq_hydropathy:
        KtoR = True
    else:
        KtoR = False

    hydro_error = abs(seq_hydropathy - objective_hydropathy)
    total_hydro_error = hydro_error*len(sequence)
    num_to_change = round(total_hydro_error / 0.6)

    # figure out if changing K to R or R to K
    if KtoR == True:
        possible_changes = sequence.count('K')
    else:
        possible_changes = sequence.count('R')

    # figure out number residues to change        
    if possible_changes >= num_to_change:
        num_changed_residues = num_to_change
    else:
        num_changed_residues = possible_changes
    
    # figure out locations of residues to change
    if KtoR == True:    
        locations = identify_residue_positions(sequence, 'K')
    else:
        locations = identify_residue_positions(sequence, 'R')

    # shuffle the locations
    random.shuffle(locations)

    # make empty list to hold chosen locations to change
    chosen_locations=[]

    # choose locations to use
    if len(locations) > num_changed_residues:
        for i in range(0, num_changed_residues):
            chosen_locations.append(locations[i])
    else:
        chosen_locations = locations

    # build the final sequence
    final_seq = ''
    for i in range(0, len(sequence)):
        if i in chosen_locations:
            if KtoR == True:
                final_seq += 'R'
            else:
                final_seq += 'K'
        else:
            final_seq += sequence[i]
    # return the final sequence
    return final_seq


def fraction_net_charge(length, fraction, net_charge):

    """

    Function for figuring out NCPR / FCR values for a sequence based 
    on the sequence length, the wanted fraction, and the wanted
    net_charge. The reason this is a bit convoluted is that the
    function was written to get the fraction and net_charge as
    close as possible, and if the FCR / NCPR / length 
    combination is not possible, to prioritize
    getting the NCPR correct over the FCR. 


    Parameters
    -------------

    length : Int
        How long the sequence is.

    fraction : Float
        The objective FCR as a decimal value.

    net_charge : Float
        The objective NCPR as a decimal value.


    Returns
    ---------

    Dict
        Returns a dictionary holding the NCPR and FCR
        values based on the length, fraction, and net_charge.

    """

    # figure out how many charged residues will be needed
    added_NCPR_residues = round(length * abs(net_charge))
    
    #figure out the remaining FCR fraction
    remaining_FCR_fraction = fraction - abs(net_charge)

    # based on the remamining_FCR_fraction, figure out how many residues to add
    added_FCR_residues = round((length * remaining_FCR_fraction)/2)
    remaining_residues_over_two = round(length - added_NCPR_residues/2)

    # figure ot the number of possible residues to add to the sequence
    remaining_residues = (length - added_NCPR_residues)

    # if the number of added_FCR_residues are less than or equal to 
    # remaining residues over 2, added_residues = added_FCR_residues        
    if added_FCR_residues <= remaining_residues_over_two:
        added_residues = added_FCR_residues 

    # otherwise, added_residues are qualt to remaining_residues_over_two       
    else:
        added_residues = remaining_residues_over_two

    # if the number of residues to add is greater than the remaining number of residues...
    if (added_residues * 2) > remaining_residues:
        added_residues = added_residues - 1

    # otherwise, don't do anything to added_residues       
    else:
        added_residues = added_residues 

    # return the final dict
    return {'NCPR_residues': added_NCPR_residues, 'FCR_residues': added_residues}   


def calculate_max_charge(hydropathy):
    '''    
    Function to determine the maximum charge value depending
    on the objective hydropathy of a sequence. Empirically 
    determined by having GOOSE try to generate tons of different
    FCR / NCPR / hydropathy sequences and determining the slope of
    the line that is the cutoff between values that can generate
    disordered sequences vs. those that can't. 

    Parameters
    ----------
    hydropathy : Float
        The objective hydropathy for the final sequence

    Returns
    -------
        The maximum possible charge that a sequence can have
        for a specific hydropathy value. Two different equations
        were determined that differ slightly depending on the version
        of metapredict used. Takes the minimum of both in order to maximize
        chances that the user stays in a regime where a disordered sequence
        can be generated.

    '''
    # calculate the maximum charge values for a given hydropathy
    MAXIMUM_CHARGE_WITH_HYDRO_1 = 1.1907 + (-0.2050 * hydropathy)
    MAXIMUM_CHARGE_WITH_HYDRO_2 = 1.2756 + (-0.2289 * hydropathy)
    # return the lower value between the 2 possibilities.
    return min([MAXIMUM_CHARGE_WITH_HYDRO_1, MAXIMUM_CHARGE_WITH_HYDRO_2])


def hydropathy_optimization(sequence, objective_hydropathy, allowed_error = parameters.HYDRO_ERROR):
    '''
    function to optimize hydropathy of a sequence to bring it closer to
    the objective hydropathy

    parameters
    ----------
    sequence : string
        the amino acid sequence as as string

    objective_hydropathy : float
        the objective hydropathy of the sequence as a float

    allowed_error : float
        the allowed error between the hydropathy of the
        returned sequence and the objective_hydropathy

    returns
    -------
    current_sequence : string
        returns the final sequence that I for some really
        dumb reason named 'current_sequence'.
    '''

    optimizer=0
    current_sequence = sequence
    # set number of possible optimizations to legnth * 3 to limit how long it does this.
    # objective is to just get the sequence reasonably close here.
    while optimizer < len(sequence)*2:
        # try optimizing the sequence
        non_optimized_seq = optimize_hydro(current_sequence, objective_hydropathy, use_charged_residues=False)
        # if optimization did nothing, kill the loop
        if non_optimized_seq == current_sequence:
            current_sequence = non_optimized_seq
            return current_sequence
            optimizer = optimizer + length*40
        # set sequence equal to non_optimized sequence for next optimization
        current_sequence = non_optimized_seq
        # add one to optimizer value
        optimizer = optimizer + 1

        current_hydropathy = round(Protein(current_sequence).hydropathy, 4)
        # see if it matches mean_hydro within allowed_error
        if abs(objective_hydropathy - current_hydropathy) <= allowed_error:
            # set final sequence
            final_sequence = current_sequence
            # return the sequence
            return final_sequence
            optimizer = optimizer + length*40
    return current_sequence


def replace_residues(sequence, residue, replacement):
    '''
    function to replace the first instance of a residue with a different residue

    parameters
    ----------
    sequence : string
        the amino acid sequence as a string

    residue : string
        the residue to identify as a string

    replacement : string
        the amino acid to replace the first instance of the
        residue 'residue' with

    returns
    --------
    sequence : string
        returns the sequnce with the altered residue as a string
    '''
    if residue in sequence:
        cur_res_loc = sequence.index(residue)
        if cur_res_loc == 0:
            final_seq = replacement + sequence[1:]
        elif cur_res_loc == len(sequence)-1:
            final_seq = sequence[:len(sequence)-1] + replacement
        else:
            final_seq = sequence[:cur_res_loc] + replacement + sequence[cur_res_loc+1:]
        return final_seq
    return sequence



def FCR_optimization(sequence, objective_hydropathy, allowed_error=parameters.HYDRO_ERROR):
    '''
    function to modify residues in the hydropathy and FCR 
    function for when the charged residues interfere with
    a correct hydropathy value

    parameters
    ----------
    sequence : string
        the amino acid sequence as a string

    objective_hydropathy : float
        the objective hydropathy of the sequence as a float

    allowed error : float
        the allowed error in hydropathy value between the final 
        sequence hydropathy and the objective hydropathy

    returns
    -------
    final_seq : string
        returns the final sequence as a string
    '''
    
    cur_hydro = Protein(sequence).hydropathy
    
    if abs(objective_hydropathy-cur_hydro) > allowed_error:
        new_sequence = sequence
        if cur_hydro > objective_hydropathy:
            # replace D's with K
            count_D = sequence.count('D')
            for i in range(0, count_D):
                cur_hydro = Protein(new_sequence).hydropathy
                if cur_hydro > objective_hydropathy:
                    new_sequence = replace_residues(new_sequence, 'D', 'K')
            # replace E's with K
            if cur_hydro > objective_hydropathy:
                count_E = sequence.count('E')
                for i in range(0, count_E):
                    cur_hydro = Protein(new_sequence).hydropathy
                    if cur_hydro > objective_hydropathy:
                        new_sequence = replace_residues(new_sequence, 'E', 'K')
                # return the sequence after K_R optimizations
                final_seq = K_R_optimization(new_sequence, objective_hydropathy)
        else:
            # replace K's with D
            count_K = sequence.count('K')
            for i in range(0, count_K):
                cur_hydro = Protein(new_sequence).hydropathy
                if cur_hydro < objective_hydropathy:
                    new_sequence = replace_residues(new_sequence, 'K', 'D')
            # replace R's with E
            if cur_hydro < objective_hydropathy:
                count_R = sequence.count('R')
                for i in range(0, count_R):
                    cur_hydro = Protein(new_sequence).hydropathy
                    if cur_hydro > objective_hydropathy:
                        new_sequence = replace_residues(new_sequence, 'R', 'E')
            # return the sequence after K_R optimizations
            final_seq = K_R_optimization(new_sequence, objective_hydropathy)
    else:
        final_seq = sequence
    return final_seq   




def create_seq_by_props(length, FCR=None, NCPR=None, hydropathy=None, attempts=1, allowed_hydro_error = parameters.HYDRO_ERROR, exclude = []):

    '''
    Function that allows the user to generate a sequence with specific
    properties. Sequence generated is very likely to be disordered but
    is not actually checked for disorder itself.

    Parameters
    ----------
    length : Int
        The length of the wanted protein sequence as an integer value

    FCR : Float
        The fraction of charged residues wanted for the sequence as a 
        decimal value.

    NCPR : Float
        The wanted net charge of the sequence given as a decimal value

    hydropathy : Float 
        The wanted mean hydropathy value of the sequence. Uses ajusted
        Kyte-doolittle hydropathy scale that goes from 0 to 9 where
        0 is least hydrophobic and 9 is most hydrophobic.

    attempts : Int
        The number of times to attempt to build a sequence before throwing
        in the towel

    allowed_hydro_error : Float
        The allowed error for hydropathy between the value of hydropathy and
        the final hydropathy value of the generated sequence

    exclude : list
        A list of residues that are not allowed to be used in the generation
        the specific sequecne


    Returns
    -------
    final_seq : String
        Returns the final sequence that was specified as a string.
        Once again, this sequence is typically disordered but IS
        NOT GAURANTEED TO BE DISORDERED!
    '''

    for num_attempts in range(0, attempts):

        # make an empty string to hold 4 amino acids that are 'starter'
        # amino acids for the get_optimal_residue function.
        seq = ''
        for i in range(0, 4):
            seq += random_amino_acid(lists.disordered_list)

        # empty string to hold amino acids for the final sequence
        final_seq = ''

        # just length
        #===============================================#
        if FCR==None and NCPR==None and hydropathy==None:
            for i in range(0, length):
                if len(final_seq) < 4:
                    cur_input = seq[len(seq)-4:]
                    residue = get_optimal_residue(cur_input, exclude_residues = exclude)
                    final_seq += residue
                    seq += residue
                else:
                    cur_input = final_seq[len(final_seq)-4:]
                    residue = get_optimal_residue(cur_input, exclude_residues = exclude)
                    final_seq += residue
            return final_seq


        # Just hydropathy
        #===============================================#
        elif FCR == None and NCPR == None and hydropathy != None:
            try:
                final_seq = hydro_seq(length, hydropathy, exclude_residues=exclude)
                return final_seq
            except:
                raise GooseError('Failed to make sequence.')


        # Just FCR
        #===============================================#
        elif FCR != None and NCPR == None and hydropathy == None:
            # Generate random list of charged residue positions
            number_charged = round(length * FCR)
            charged_positions = gen_charged_positions(length, number_charged)
            # account for necessary residue exclusions
            if exclude != []:
                final_exclusion = ['D', 'E', 'K', 'R']
                for aa in exclude:
                    if aa not in final_exclusion:
                        final_exclusion.append(aa)
            else:
                final_exclusion = ['D', 'E', 'K', 'R']

            for i in range(0, length):
                if i not in charged_positions:                
                    if len(final_seq) < 4:
                        cur_input = seq[len(seq)-4:]
                        residue = get_optimal_residue(cur_input, exclude_residues = final_exclusion)
                        final_seq += residue
                        seq += residue
                    else:
                        cur_input = final_seq[len(final_seq)-4:]
                        residue = get_optimal_residue(cur_input, exclude_residues = final_exclusion)
                        final_seq += residue
                else:
                    residue = random_amino_acid(lists.charged_list)
                    seq += residue
                    final_seq += residue
            return final_seq



        # Just NCPR
        #===============================================#
        elif FCR == None and NCPR != None and hydropathy == None:

            #  Set up necessary residues to match net_charge
            #--------------------------------------------#
            # figure out how many charged residues to add
            num_charged = round(abs(NCPR)*length)

            # Make a list to hold charged residues
            charged_residues = ''
            # if the net_charge is negative...
            if NCPR < 0:
                for i in range(0, num_charged):
                    # add negative residues
                    charged_residues += random_amino_acid(lists.D_E)

            # if the net charge is positive
            elif NCPR > 0:
                for i in range(0, num_charged):
                    # add positive residues
                    charged_residues += random_amino_acid(lists.K_R)

            #  Randomly adjust FCR while matching net_Charge
            #-----------------------------------------------#
            # figure out the possible additional charges that could be added to the seq
            possible_additional_charges = round((0.9 - abs(NCPR)) * length)

            # if we can still add charged residues
            if possible_additional_charges > 0:
                # figure out a random number of charged residues to add
                added_additional_charges = round(random.randint(0, possible_additional_charges)/2)

                # for the number of added additional charges 
                for i in range(0, added_additional_charges):
                    # add one positive and one negative to keep the sequence NCPR from changing
                    charged_residues += random_amino_acid(lists.D_E)
                    charged_residues += random_amino_acid(lists.K_R)


            # Generate random list of charged residue positions
            charged_positions = gen_charged_positions(length, len(charged_residues))

            # turn charged_residues into shuffled list
            charged_residues = shuffle_seq(charged_residues)
            charged_residue_list = list([char for char in charged_residues])

            # account for necessary residue exclusions
            if exclude != []:
                final_exclusion = ['D', 'E', 'K', 'R']
                for aa in exclude:
                    if aa not in final_exclusion:
                        final_exclusion.append(aa)
            else:
                final_exclusion = ['D', 'E', 'K', 'R']

            for i in range(0, length):
                if i not in charged_positions:                
                    if len(final_seq) < 4:
                        cur_input = seq[len(seq)-4:]
                        residue = get_optimal_residue(cur_input, exclude_residues = final_exclusion)
                        final_seq += residue
                        seq += residue
                    else:
                        cur_input = final_seq[len(final_seq)-4:]
                        residue = get_optimal_residue(cur_input, exclude_residues = final_exclusion)
                        final_seq += residue
                else:
                    residue = charged_residue_list.pop()
                    seq += residue
                    final_seq += residue
            return final_seq


        # FCR and NPCR
        #===============================================#
        elif FCR != None and NCPR != None and hydropathy == None:
            
            # figure out how many residues for FCR and how many for NCPR
            charged_residue_dict = fraction_net_charge(length, fraction=FCR, net_charge=NCPR)

            # now use that to make charged_residues str
            charged_residues = ''
            # if the net_charge is negative...
            if NCPR < 0:
                for i in range(0, charged_residue_dict['NCPR_residues']):
                    # add negative residues
                    charged_residues += random_amino_acid(lists.D_E)

            # if the net charge is positive
            elif NCPR > 0:
                for i in range(0, charged_residue_dict['NCPR_residues']):
                    # add positive residues
                    charged_residues += random_amino_acid(lists.K_R)


            # now accomodate for the FCR residues if necessary
            if charged_residue_dict['FCR_residues'] != 0:
                for fcr_residues in range(0, charged_residue_dict['FCR_residues']):
                    charged_residues += random_amino_acid(lists.K_R)
                    charged_residues += random_amino_acid(lists.D_E)

            # Generate random list of charged residue positions
            charged_positions = gen_charged_positions(length, len(charged_residues))

            # turn charged_residues into shuffled list
            charged_residues = shuffle_seq(charged_residues)
            charged_residue_list = list([char for char in charged_residues])
            
            # account for necessary residue exclusions
            if exclude != []:
                final_exclusion = ['D', 'E', 'K', 'R']
                for aa in exclude:
                    if aa not in final_exclusion:
                        final_exclusion.append(aa)
            else:
                final_exclusion = ['D', 'E', 'K', 'R']

            for i in range(0, length):
                if i not in charged_positions:                
                    if len(final_seq) < 4:
                        cur_input = seq[len(seq)-4:]
                        residue = get_optimal_residue(cur_input, exclude_residues = final_exclusion)
                        final_seq += residue
                        seq += residue
                    else:
                        cur_input = final_seq[len(final_seq)-4:]
                        residue = get_optimal_residue(cur_input, exclude_residues = final_exclusion)
                        final_seq += residue
                else:
                    residue = charged_residue_list.pop()
                    seq += residue
                    final_seq += residue
            return final_seq


        # FCR and hydropathy
        #===============================================#
        elif FCR != None and NCPR == None and hydropathy !=None:
            # if FCR set to 0, just make a neutral hydropathy sequence
            if FCR == 0:
                final_seq = hydro_seq(length, hydropathy, just_neutral=True)
            
            else:
                # first figure out total hydropathy needed
                total_hydropathy = hydropathy*length
                # figure out how many residues are needed
                number_charged_res = round(FCR*length)
                # now figure out how many res you get for rest of seq
                number_hydro_res = length - number_charged_res
                # make string to hold charged residues
                charged_residues = ''


                # add necessary charged residues to the string
                for fcr_residues in range(0, number_charged_res):    
                    charged_residues += random_amino_acid(lists.charged_list)

                # figure out how much hydropathy is taken by charged residues
                charged_residue_hydropathy = Protein(charged_residues).hydropathy


                # now figure out hydropathy needed to balance the charged residues
                total_charged_hydro = charged_residue_hydropathy * len(charged_residues)
                remaining_needed_hydro = total_hydropathy - total_charged_hydro
                
                if number_hydro_res > 0:
                    # figure out hydropathy for subseq to make that baalnces out the charged residues
                    # to get the total objective hydropathy
                    new_mean_hydro = round(remaining_needed_hydro/number_hydro_res, 4)           

                    if new_mean_hydro < 1:
                        charged_residues = generate_charged_residues(length, FCR, hydropathy)

                    # figure out how much hydropathy is taken by charged residues
                    charged_residue_hydropathy = Protein(charged_residues).hydropathy


                    # now figure out hydropathy needed to balance the charged residues
                    total_charged_hydro = charged_residue_hydropathy * len(charged_residues)
                    remaining_needed_hydro = total_hydropathy - total_charged_hydro

                    new_mean_hydro = round(remaining_needed_hydro/number_hydro_res, 4) 

                    # can't make something with hydropathy over 9.
                    if new_mean_hydro > 9:
                        new_mean_hydro = 9

                    # figure out allowed error. This needs to be adjusted because the hydro_seq function will 
                    # simply account for the non-charged residues, but because we have additional residues
                    # we can actually have more error on hydro_seq thus making this move a lot faster.
                    standard_hydro_error = allowed_hydro_error
                    if FCR != 0:
                        adjusted_error = standard_hydro_error / (1-FCR)
                    else:
                        adjusted_error = standard_hydro_error
                    if adjusted_error < allowed_hydro_error:
                        adjusted_error = allowed_hydro_error

                    # now make the sequence that will be used
                    input_hydro_sequence = hydro_seq(number_hydro_res, new_mean_hydro, just_neutral=True, allowed_error = adjusted_error, return_best_seq = True, exclude_residues=exclude)
                    
                    # Generate random list of charged residue positions
                    charged_positions = gen_charged_positions(length, len(charged_residues))

                    # turn charged_residues into shuffled list
                    charged_residues = shuffle_seq(charged_residues)
                    charged_residue_list = list([char for char in charged_residues])

                    # now use input_hydro_sequence and charged_residues to build the final sequence
                    final_seq = ''
                    # keep track of where you are for the input_hydro_sequence
                    input_hydro_sequence_index = 0
                    for i in range(0, length):
                        if i not in charged_positions:                
                            final_seq += input_hydro_sequence[input_hydro_sequence_index]
                            input_hydro_sequence_index += 1
                        else:
                            residue = charged_residue_list.pop()
                            final_seq += residue
                    # check hydropathy
                    if check_hydropathy(final_seq, hydropathy, allowed_hydro_error):
                        final_seq = final_seq
                    else:
                        final_seq = hydropathy_optimization(final_seq, hydropathy, allowed_error = allowed_hydro_error)    
                    # check again
                    #check_hydropathy(sequence, objective_hydropathy, cutoff_val):
                    if check_hydropathy(final_seq, hydropathy, allowed_hydro_error):
                        final_seq = final_seq
                    # if still off, try K_R optimizations.
                    else:
                        final_seq = K_R_optimization(final_seq, hydropathy)

                else:
                    final_seq = shuffle_seq(charged_residues)
                    # adjust the K and R values to match objective hydropathy.
                    final_seq = K_R_optimization(final_seq, hydropathy)
            # final check on the hydropathy
            if check_hydropathy(final_seq, hydropathy, allowed_hydro_error):       
                return final_seq
            else:
                # one last possible save.
                final_seq = FCR_optimization(final_seq, hydropathy, allowed_hydro_error)
                if check_hydropathy(final_seq, hydropathy, allowed_hydro_error):
                    return final_seq


        # NCPR and hydropathy
        #===============================================#
        elif FCR == None and NCPR != None and hydropathy != None:

            # figure out total hydropathy
            total_hydropathy = length*hydropathy

            # now use that to make charged_residues str
            charged_residues = ''

            # figure out max hydropathy 
            max_FCR = calculate_max_charge(hydropathy)

            # figure out if there is room for additional charge
            # with a little wiggle room 
            if max_FCR - 0.1 - (abs(NCPR)) > 0:
                possible_additional_charge = max_FCR - 0.1 - (abs(NCPR))
                # figure out min possible unit based on length
                min_unit = 1/length
                # use to populate list of potential FCR values
                potential_FCR_values = []
                start_charge = 0
                while start_charge < possible_additional_charge:
                    potential_FCR_values.append(start_charge)
                    start_charge += (min_unit*2)
                # now pick a random charge value
                if potential_FCR_values != []:
                    if len(potential_FCR_values) > 1:
                        FCR_value = potential_FCR_values[randint(0, len(potential_FCR_values)-1)]
                    else:
                        FCR_value = potential_FCR_values[0]
                if FCR_value != 0:
                    num_residues = length*FCR_value
                    if num_residues >= 2:
                        if num_residues % 2 != 0:
                            num_residues = num_residues - 1
                        # figure out if you need a specific charge value
                        charge_iterations = round(num_residues / 2)
                        for i in range(0, charge_iterations):
                            charged_residues += random_amino_acid(lists.D_E)
                            if hydropathy > 1:
                                charged_residues += random_amino_acid(lists.K_R)
                            else:
                                charged_residues += 'R'

            # now take care of FCR
            ncpr_residues_to_add = round(length*abs(NCPR))
            # if the net_charge is negative...
            if NCPR < 0:
                for i in range(0, ncpr_residues_to_add):
                    # add negative residues
                    charged_residues += random_amino_acid(lists.D_E)

            # if the net charge is positive
            elif NCPR > 0:
                for i in range(0, ncpr_residues_to_add):
                    # add positive residues
                    if hydropathy < 1:
                        charged_residues += 'R'
                    else:                
                        charged_residues += random_amino_acid(lists.K_R)


            # if NCPR was set to 0 and the FCR ends at 0, just make the sequence
            if charged_residues == '':
                final_seq = hydro_seq(length, hydropathy, just_neutral=True, exclude_residues = exclude)

            else:
                # calculate FCR
                final_FCR = round(len(charged_residues) / length, 4)


                # now need to take care of the hydropathy part.
                number_hydro_res = length - len(charged_residues)

                # figure out how much hydropathy is taken by charged residues
                charged_residue_hydropathy = Protein(charged_residues).hydropathy

                # now figure out hydropathy needed to balance the charged residues
                total_charged_hydro = charged_residue_hydropathy * len(charged_residues)
                remaining_needed_hydro = total_hydropathy - total_charged_hydro

                if number_hydro_res > 0:

                    new_mean_hydro = round(remaining_needed_hydro/number_hydro_res, 4)

                    # make sure new_mean_hydro not too high.
                    if new_mean_hydro > 9:
                        new_mean_hydro = 9

                    # figure out allowed error
                    standard_hydro_error = allowed_hydro_error
                    if final_FCR != 0:
                        adjusted_error = standard_hydro_error / (1-final_FCR)
                    else:
                        adjusted_error = allowed_hydro_error
                    if adjusted_error < allowed_hydro_error:
                        adjusted_error = allowed_hydro_error


                    # now make the sequence that will be used
                    input_hydro_sequence = hydro_seq(number_hydro_res, new_mean_hydro, just_neutral=True, allowed_error =  adjusted_error, return_best_seq=True, exclude_residues = exclude)
                    
                    # Generate random list of charged residue positions
                    charged_positions = gen_charged_positions(length, len(charged_residues))

                    # turn charged_residues into shuffled list
                    charged_residues = shuffle_seq(charged_residues)
                    charged_residue_list = list([char for char in charged_residues])

                    # now use input_hydro_sequence and charged_residues to build the final sequence
                    final_seq = ''
                    # keep track of where you are for the input_hydro_sequence
                    input_hydro_sequence_index = 0
                    for i in range(0, length):
                        if i not in charged_positions:                
                            final_seq += input_hydro_sequence[input_hydro_sequence_index]
                            input_hydro_sequence_index += 1
                        else:
                            residue = charged_residue_list.pop()
                            final_seq += residue
                    
                    # check to see if hydropathy optimization needed
                    if check_hydropathy(final_seq, hydropathy, allowed_hydro_error):
                        final_seq = final_seq
                    else:
                        final_seq = hydropathy_optimization(final_seq, hydropathy, allowed_error = allowed_hydro_error)
                    
                    # check to see of K_R optimization needed
                    if check_hydropathy(final_seq, hydropathy, allowed_hydro_error):
                        final_seq = final_seq
                    else:
                        final_seq = K_R_optimization(final_seq, hydropathy)        
                else:
                    final_seq = shuffle_seq(charged_residues)
                    final_seq = K_R_optimization(final_seq, hydropathy)

            # check hydropathy, return the final sequence
            if check_hydropathy(final_seq, hydropathy, allowed_hydro_error):    
                return final_seq 


        # NCPR and hydropathy and FCR
        #===============================================#
        elif FCR != None and NCPR != None and hydropathy != None:
            
            if FCR == 0 and NCPR == 0:
                final_seq = hydro_seq(length, hydropathy, just_neutral=True, exclude_residues=exclude)

            else:
                # figure out how many residues for FCR and how many for NCPR
                charged_residue_dict = fraction_net_charge(length, fraction=FCR, net_charge=NCPR)

                # now use that to make charged_residues str
                charged_residues = ''

                # if the net_charge is negative...
                if NCPR < 0:
                    for i in range(0, charged_residue_dict['NCPR_residues']):
                        # add negative residues
                        charged_residues += random_amino_acid(lists.D_E)

                # if the net charge is positive
                elif NCPR > 0:
                    for i in range(0, charged_residue_dict['NCPR_residues']):
                        # add positive residues
                        if hydropathy < 1:
                            charged_resides += 'R'
                        else:
                            charged_residues += random_amino_acid(lists.K_R)


                # now accomodate for the FCR residues if necessary
                if charged_residue_dict['FCR_residues'] != 0:
                    for fcr_residues in range(0, charged_residue_dict['FCR_residues']):

                        if hydropathy > 1:
                            charged_residues += random_amino_acid(lists.K_R)
                        else:
                            charged_residues += 'R'

                        charged_residues += random_amino_acid(lists.D_E)

                # first figure out total hydropathy needed
                total_hydropathy = hydropathy*length
                # now figure out how many res you get for rest of seq
                number_hydro_res = length - len(charged_residues)

                # figure out how much hydropathy is taken by charged residues
                charged_residue_hydropathy = Protein(charged_residues).hydropathy

                # now figure out hydropathy needed to balance the charged residues
                total_charged_hydro = charged_residue_hydropathy * len(charged_residues)
                remaining_needed_hydro = total_hydropathy - total_charged_hydro

                if number_hydro_res > 0:

                    new_mean_hydro = round(remaining_needed_hydro/number_hydro_res, 4)

                    # can't make something with hydropathy over 9.
                    if new_mean_hydro > 9:
                        new_mean_hydro = 9
                    
                    # figure out allowed error
                    standard_hydro_error = allowed_hydro_error
                    if FCR != 0:
                        adjusted_error = standard_hydro_error / (1-FCR)
                    else:
                        adjusted_error = allowed_hydro_error
                    if adjusted_error < allowed_hydro_error:
                        adjusted_error = allowed_hydro_error

                    # now make the sequence that will be used
                    input_hydro_sequence = hydro_seq(number_hydro_res, new_mean_hydro, just_neutral=True, allowed_error = adjusted_error, return_best_seq = True, exclude_residues=exclude)

                    # Generate random list of charged residue positions
                    charged_positions = gen_charged_positions(length, len(charged_residues))

                    # turn charged_residues into shuffled list
                    charged_residues = shuffle_seq(charged_residues)
                    charged_residue_list = list([char for char in charged_residues])

                    # now use input_hydro_sequence and charged_residues to build the final sequence
                    final_seq = ''
                    # keep track of where you are for the input_hydro_sequence
                    input_hydro_sequence_index = 0
                    for i in range(0, length):
                        if i not in charged_positions:                
                            final_seq += input_hydro_sequence[input_hydro_sequence_index]
                            input_hydro_sequence_index += 1
                        else:
                            residue = charged_residue_list.pop()
                            final_seq += residue

                    # check hydropathy
                    if check_hydropathy(final_seq, hydropathy, allowed_hydro_error):
                        final_seq = final_seq
                    else:
                       final_seq = hydropathy_optimization(final_seq, hydropathy, allowed_error = allowed_hydro_error)

                    # check again
                    if check_hydropathy(final_seq, hydropathy, allowed_hydro_error):
                        final_seq = final_seq
                    # if still off, try K_R optimizations.
                    else:
                        final_seq = K_R_optimization(final_seq, hydropathy)
                # if no non-charged residues, shuffle the sequence and optimize for hydro with K/R if needed
                else:
                    final_seq = shuffle_seq(charged_residues)
                    final_seq = K_R_optimization(final_seq, hydropathy)

            if check_hydropathy(final_seq, hydropathy, allowed_hydro_error):       
                return final_seq

    # if nothing works, raise error that GOOSE failed to make the sequence.
    raise GooseError('Failed to make sequence.')


def test_seq_by_props(length, FCR=None, NCPR=None, hydropathy=None):
    """
    TO DO - add docstring

    """
    final_sequence = create_seq_by_props(length, FCR=FCR, NCPR=NCPR, hydropathy=hydropathy)

    tmp = Protein(final_sequence)
    seq_NCPR = tmp.NCPR
    seq_FCR = tmp.FCR
    seq_hydropathy = tmp.hydropathy
    if FCR != None:
        print(f'Objective FCR = {FCR} : Actual FCR = {seq_FCR}')
    if NCPR != None:
        print(f'Objective NCPR = {NCPR} : Actual NCPR = {seq_NCPR}')
    if hydropathy != None:
        print(f'Objective hydropathy = {hydropathy} : Actual hydropathy = {seq_hydropathy}')
    print(final_sequence)

#/-/-/-/-/-/-/-/-/-/-/-/-/-/-/-/-/-/-/-/-/-/-/-/-/-/-/-/-/-/-/-/-/-/-/-
#/-/-/-/-/-/-/-/-/-/-/- Amino acid Fractions /-/-/-/-/-/-/-/-/-/-/-/-/-
#/-/-/-/-/-/-/-/-/-/-/-/-/-/-/-/-/-/-/-/-/-/-/-/-/-/-/-/-/-/-/-/-/-/-/-

def create_seq_by_fracs(length, max_aa_fractions={}, choose_optimized_residue=True, **kwargs):
    """
    This will return a sequence with the specified fraction of
    amino acids. To use simply specify the amino acid followed by
    the fraction you want as a decimal.
    
    Parameters
    ------------
    length : int
        length of desired disordered sequence

    max_aa_fractions : dict 
        Dictionary which, if provided, allows the user to over-ride the 
        fraction of a sequence which can be made up of any given amino
        acid. The passed dictionary should contain key/value pairs, where
        keys are one of the twenty standard amino acids and values is a
        float between 0 and 1. If amino acids are missing then the default
        thresholds set by GOOSE are used.

    choose_optimized_residue : bool
        Whether to use the weighted probabiltiy dictionary to select residues
        that have the highest likelihood of being disordered in the sequence. This
        does generally work and improve sequence generation speed but comes at the
        drawback of reducing possible sequence space. For some sequence compositions
        this can make it difficult to generate a disordered sequence, so this can be
        set to False in an attempt to generate those sequecnces

    <each of the 20 amino acids> : float
        Specify the fraction of the sequence that should be made up of one or more
        of the 20 natural amino acids (e.g. A=0.2, Y=0.05) etc.


    Returns
    -----------
    string
       A string of the amino acid sequence

    """

    # dict holding the max fractions that each amino acid can be individually specified as
    max_fraction = {"A": parameters.MAX_FRACTION_A,
    "A": parameters.MAX_FRACTION_A,
    "R": parameters.MAX_FRACTION_R,
    "N": parameters.MAX_FRACTION_N,
    "D": parameters.MAX_FRACTION_D,
    "C": parameters.MAX_FRACTION_C,
    "Q": parameters.MAX_FRACTION_Q,
    "E": parameters.MAX_FRACTION_E,
    "G": parameters.MAX_FRACTION_G,
    "H": parameters.MAX_FRACTION_H,
    "I": parameters.MAX_FRACTION_I,
    "L": parameters.MAX_FRACTION_L,
    "K": parameters.MAX_FRACTION_K,
    "M": parameters.MAX_FRACTION_M,
    "F": parameters.MAX_FRACTION_F,
    "P": parameters.MAX_FRACTION_P,
    "S": parameters.MAX_FRACTION_S,
    "T": parameters.MAX_FRACTION_T,
    "W": parameters.MAX_FRACTION_W,
    "Y": parameters.MAX_FRACTION_Y,
    "V": parameters.MAX_FRACTION_V
    }

    # if we passed an over-ride dictionary for the max amino acid fraction...
    for aa in max_aa_fractions:
        max_fraction[aa] = max_aa_fractions[aa]
            

    # before we do any actual design we can extract the unchanging sequence composition
    # features. The sequence_list we build here is essentially the compositional
    # constraint we're putting on the design

    #           Checks on input parameters
    #=-=#=-=#=-=#=-=#=-=#=-=#=-=#=-=#=-=#=-=#=-=#=-=#
    # make sure no fractions specified as greater than they can be generated as in testing
    for amino in kwargs.keys():
        if kwargs[amino] > max_fraction[amino]:
            exception_value = f'Specified fraction for {amino} of {kwargs[amino]} is greater than max allowed fraction for {amino} of {max_fraction[amino]}!'
            raise GooseInputError(exception_value)

    total_fraction = sum(kwargs.values())
    if round(total_fraction, 5) > 1:
        raise GooseInputError('Cannot specify a total fraction of residues greater than 1!')


    #  Start attempts to build the sequence
    #=-=#=-=#=-=#=-=#=-=#=-=#=-=#=-=#=-=#=-=#=-=#=-=#
    # all amino acids
    amino_acids = ['A', 'C', 'D', 'E', 'F', 'G', 'H', 'I', 'K', 'L', 'M', 'N', 'P', 'Q', 'R', 'S', 'T', 'V', 'W', 'Y']
    # make empty final sequence
    final_sequence = ""

    sequence_list=[]

    used_AAs = []

    # for each amino acid for which we passed fractional information...
    for cur_AA in kwargs.keys():

        # keep track of residues can't use down the line.
        used_AAs.append(cur_AA)

        # remove the amino acids from the list for 
        # sequence generation of fraction != 1.
        amino_acids.remove(cur_AA)

        # res_count is the actual number of cur_AA expected in a sequence
        # of $length residues with $cur_frac fraction. Rounds until it can't. 
        res_count = round(kwargs[cur_AA]*length)
        if len(sequence_list) + res_count > length:
            res_count=int(kwargs[cur_AA]*length)

        # generate a list of the right number copies of the current
        # amino acid
        AA_homopolymer = [cur_AA] * res_count
        
        # extend the ever-growing sequence_list to include these
        # residues
        sequence_list.extend(AA_homopolymer)


    # shuffle here so no risk of introducing bias based on order in which
    # residues are introduced
    random.shuffle(sequence_list)

    # At this point sequence_list is a list that is equal to or shorter than
    # $length, and contains the residues expected in the design sequence at the
    # correct proportions once the sequence is actually $length residues long
    # if all the residues are accounted for...
    if total_fraction == 1:
        # make sure length is reached. This can be an issue if the user doesn't
        # carefully choose compatible lengths and fractional values.
        # only need to do if len sequence lenght != length && total fraction == 1
        if len(sequence_list)<length:
            # need to decide which residue(s) to add to finish sequence.
            # basically figure out which residues are off by the most as far as 
            # what their inteded value was and their actual residue count.
            off_by_dict={}
            for res in used_AAs:
                off_by_dict[res]=abs((kwargs[res]*length)-sequence_list.count(res))
            # now make a sorted list of off by vals...
            add_in_order=[]
            for off_by_val in set(sorted(off_by_dict.values(), reverse=True)):
                for res in off_by_dict:
                    if off_by_dict[res]==off_by_val:
                        add_in_order.append(res)
            # add residues in order of off_by list until seq length met. 
            res_order_added=0
            off_by=length-len(sequence_list)
            for val in range(0, off_by):
                sequence_list.append(add_in_order[res_order_added])
                res_order_added+=1
                # if past max list index, reset.
                if res_order_added >= len(sequence_list):
                    res_order_added=0
        elif len(sequence_list)>length:
            raise GooseError('Error in sequence_generation_backend.py causing create_seq_by_fracs() function to make sequences too long. Please contact the developers or post an issue on Github!')

        else:
            # create a string from the sequence_list
            all_fraction_sequence = "".join(sequence_list)
            # shuffle the seq before returning    
            sequence = shuffle_seq(all_fraction_sequence)
            return sequence

    # ELSE we still need some extra residues
    else:
        # how many residues are we missing?            
        number_of_additional_res = length - len(sequence_list)

        # build a disordered starter sequence which we'll use to help us select 
        # randomly generated disordered residues 
        starter_seq = ''
        for i in range(0, 4):
            starter_seq += random_amino_acid(lists.disordered_list)

        # create a copy of sequence_list. 
        local_sequence_list = sequence_list.copy()

        # add '0' to the end of the local_sequence_list such that local_sequence_list
        # becomes a list where len(local_sequence_list) == length and there are 
        # $number_of_additional_res copies of '0' in the list. Then, shuffle the list,
        # which acts to randomly distribute the '0' elements across the list
        local_sequence_list.extend(['0']*number_of_additional_res)
        random.shuffle(local_sequence_list)


        # finally build the actual sequence, whereby all the '0' elements are replaced with
        # an optimized residue based on the preceding 4 residues
        sequence = ''
        for i in range(0, length):

            # if found a position where we need a new amino acid...
            if local_sequence_list[i] == '0':
                if choose_optimized_residue==True:
                    # if at N-terminus of sequence use the starter_seq to provide 'disordered' context
                    if len(sequence) < 4:
                        # figure out where in the sequence we are
                        start_seq_number = len(starter_seq)-4
                        # use that info to make a slice of the sequence to be used as the input for get_optimal_residue()
                        next_residue_key = starter_seq[start_seq_number:]
                        # choose next residue
                        chosen_residue = get_optimal_residue(next_residue_key, exclude_residues=used_AAs)
                        starter_seq += chosen_residue


                    # else use the actual sequence and select the optimimum residue
                    else:
                        # figure out where in the sequence we are
                        start_seq_number = len(sequence)-4
                        # use that info to make a slice of the sequence to be used as the input for get_optimal_residue()
                        next_residue_key = sequence[start_seq_number:]
                        # choose next residue
                        chosen_residue = get_optimal_residue(next_residue_key, exclude_residues=used_AAs)
                else:
                    # otherwise choose a random amino acid from the list that contains all 
                    # amino acids except for those that have been specified.
                    chosen_residue = amino_acids[random.randint(0, len(amino_acids)-1)]

            # ... else just use the amino acid from the shuffled local_sequence_list
            else:
                chosen_residue = local_sequence_list[i]

            # add chosen residue to sequence
            sequence = sequence + chosen_residue
        # return the sequence
        return sequence





