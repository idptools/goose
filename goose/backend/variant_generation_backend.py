'''
backend functionality for variant generation
'''
import random
from random import randint
import numpy as np

from sparrow import Protein as pr
from goose.goose_exceptions import GooseInputError, GooseBackendBug
from goose.backend.protein import Protein
from goose.backend.sequence_generation_backend import identify_residue_positions, get_optimal_residue, optimal_residue_key, create_seq_by_props, calculate_max_charge, fraction_net_charge
from goose.backend import lists
from goose.backend.amino_acids import AminoAcid
from goose.backend import parameters
from goose.backend.ginell_clustering_parameter import calculate_average_inverse_distance_from_sequence as clustering_param
from goose.backend_vectorized.optimize_kappa_vectorized import optimize_kappa_vectorized

def return_num_for_class(sequence):
    '''
    function to return the number of residues in 
    each class for an input sequence
    
    classes broken down as follows:
    Aromatics : F W Y
    Polar : Q N S T
    Hydrophobic : I V L A M
        Positive charge : K R
        Negative Charge : D E
        Special cases (not in any class)
            Glycine (G)
            Proline (P)
            Histidine (H)    
            Cystine (C)

    parameters
    ----------
    sequence : string
        the amino acid sequence as a string

    returns : dict
        returns dict with number of amino acids in each
        class as defined above

    '''
    # make sure sequence is uppercase
    sequence=sequence.upper()

    # count number of residues in each class.
    num_aromatics = sequence.count('F')+sequence.count('W')+sequence.count('Y')
    num_polar = sequence.count('Q')+sequence.count('N')+sequence.count('S')+sequence.count('T')
    num_hydro = sequence.count('I')+sequence.count('V')+sequence.count('L')+sequence.count('A')+sequence.count('M')
    num_positive = sequence.count('K')+sequence.count('R')
    num_negative = sequence.count('D')+sequence.count('E')
    num_G = sequence.count('G')
    num_P = sequence.count('P')
    num_H = sequence.count('H')
    num_C = sequence.count('C')

    return {'aromatic': num_aromatics, 'polar': num_polar, 'hydrophobic': num_hydro, 'positive': num_positive,
    'negative': num_negative, 'C': num_C, 'H': num_H, 'P': num_P, 'G': num_G}


def get_optimal_res_within_class(residue_class, sequence, additional_exclusion=[], return_all=False):
    '''
    function to get the optimal residues within a class 
    for a sequence where the residue returned will
    be the residue to add to the end of the sequence
    
    parameters  
    ----------
    residue_class : string
        the residue class as an amino acid

    sequence : string
        the amino acid sequence as a string

    additional_exclusion : list
        a list of residues to not allow to be returned

    return_all : Bool
        whether to return all possible values or not

    '''

    # if the residue class is a single AA, just return it.
    if residue_class == 'G':
        return 'G'
    elif residue_class == 'P':
        return 'P'
    elif residue_class == 'H':
        return 'H'
    elif residue_class == 'C':
        return 'C'
    else:
        # first get residue 'key'
        if len(sequence) >= 4:
            four_residues = sequence[len(sequence)-4:]
        else:
            needed_residues = 4-len(sequence)
            four_residues = ''
            for i in range(0, needed_residues):
                four_residues += random.choice(lists.disordered_list)
            four_residues += sequence

        # classes of AAs
        aromatics = ['F', 'W', 'Y']
        polar = ['Q', 'N', 'S', 'T']
        hydrophobic = ['I', 'V', 'L', 'A', 'M']
        positive = ['K', 'R']
        negative = ['D', 'E']
        cannot_change = ['G', 'P', 'H', 'C']

        # residues to exclude for each class of AAs for get_optimal_residue function
        exclude_for_aromatics = ['A', 'C', 'D', 'E', 'G', 'H', 'I', 'K', 'L', 'M', 'N', 'P', 'Q', 'R', 'S', 'T', 'V']
        
        exclude_for_polar = ['A', 'C', 'D', 'E', 'F', 'G', 'H', 'I', 'K', 'L', 'M', 'P', 'R', 'V', 'W', 'Y']
        
        exclude_for_hydrophobic = [ 'C', 'D', 'E', 'F', 'G', 'H', 'K', 'N', 'P', 'Q', 'R', 'S', 'T', 'W', 'Y']
        
        exclude_for_positive = ['A', 'C', 'D', 'E', 'F', 'G', 'H', 'I', 'L', 'M', 'N', 'P', 'Q', 'S', 'T', 'V', 'W', 'Y']
        
        exclude_for_negative = ['A', 'C', 'F', 'G', 'H', 'I', 'K', 'L', 'M', 'N', 'P', 'Q', 'R', 'S', 'T', 'V', 'W', 'Y']
        
        
        if additional_exclusion != []:
            exclude_for_aromatics.extend(additional_exclusion)
            exclude_for_polar.extend(additional_exclusion)
            exclude_for_hydrophobic.extend(additional_exclusion)
            exclude_for_positive.extend(additional_exclusion)
            exclude_for_negative.extend(additional_exclusion)

        # now use get key for optimal residue function
        key = optimal_residue_key(four_residues)
        # now get optimal residue, exclude necessary residues from choice
        if residue_class == 'aromatic':
            chosen_residues = get_optimal_residue(four_residues, exclude_residues=exclude_for_aromatics, return_all=return_all)
        elif residue_class == 'polar':
            chosen_residues = get_optimal_residue(four_residues, exclude_residues=exclude_for_polar, return_all=return_all)
        elif residue_class == 'hydrophobic':
            chosen_residues = get_optimal_residue(four_residues, exclude_residues=exclude_for_hydrophobic, return_all=return_all)        
        elif residue_class == 'positive':
            chosen_residues = get_optimal_residue(four_residues, exclude_residues=exclude_for_positive, return_all=return_all)
        elif residue_class == 'negative':
            chosen_residues = get_optimal_residue(four_residues, exclude_residues=exclude_for_negative, return_all=return_all)
        else:
            raise GooseInputError('chosen residue class is not valid. Options are aromatic, polar, hydrophobic, positive, or negative.')
        
        # return chosen residue
        return chosen_residues



def hydro_range_constant_classes(sequence):
    """

    Function to calculate the minimum and maximum possible hydropathy values
    that a sequence can have if a variant is to be made that keeps the
    number of amino acids in each class the same. The amino acid classes
    are as follows
        Aromatics : F W Y
        Polar : Q N S T
        Hydrophobic : I V L A M
        Positive charge : K R
        Negative Charge : D E
        Special cases (not in any class)
            Glycine (G)
            Proline (P)
            Histidine (H)    
            Cystine (C)

    Parameters
    -------------
    sequence : String
        Amino acid sequence as a string


    Returns
    ---------
    List
        Returns a list where the first value is the lowest possible hydropathy
        for the variant and the second value is the highest possible hydropathy

    """    

    # set min and max values
    min_possible = 0
    max_possible = 0

    # min possible hydro values from each class
    AA_hydro_min = {"A": 6.3,
    "R": 0.0,
    "N": 1.0,
    "D": 1.0,
    "C": 7.0,
    "Q": 1.0,
    "E": 1.0,
    "G": 4.1,
    "H": 1.3,
    "I": 6.3,
    "L": 6.3,
    "K": 0.0,
    "M": 6.3,
    "F": 3.2,
    "P": 2.9,
    "S": 1.0,
    "T": 1.0,
    "W": 3.2,
    "Y": 3.2,
    "V": 6.3
    }

    # max possible hydro values for each class
    AA_hydro_max = {"A": 9.0,
    "R": 0.6,
    "N": 3.8,
    "D": 1.0,
    "C": 7.0,
    "Q": 3.8,
    "E": 1.0,
    "G": 4.1,
    "H": 1.3,
    "I": 9.0,
    "L": 9.0,
    "K": 0.6,
    "M": 9.0,
    "F": 7.3,
    "P": 2.9,
    "S": 3.8,
    "T": 3.8,
    "W": 7.3,
    "Y": 7.3,
    "V": 9.0
    }

    # for each of the amino acids, figure out the value to add for each class
    for aa in sequence:
        min_possible += AA_hydro_min[aa]
        max_possible += AA_hydro_max[aa]

    # return the max and min values dividued by sequence length in a list 
    return[round(min_possible/len(sequence), 10), round(max_possible/len(sequence), 10)]


def residue_optimize_hydropathy_within_class(sequence, objective_hydropathy, target_residue_index):
    '''
    function to get hydropathy value closer to objective hydropathy
    while keeping the residues within the sequence within the same
    class

    This function only changes one residue at a time, which corresponds to 'target_residue_index'!

    parameters
    ----------
    sequence : string 
        amino acid sequence as a string

    objective_hydropathy : Float
        objective hydropathy as a float value. Between 0 and 9.

    target_residue_index : int
        the index number of the target residue

    returns
    -------
    build_sequence : string
        returns the optimized sequence as a string

    '''

    # define which residues in each class
    aromatics = ['F', 'W', 'Y']
    polar = ['Q', 'N', 'S', 'T']
    hydrophobics = ['I', 'V', 'L', 'A', 'M']
    positive = ['K', 'R']

    # make an aa class dict
    aa_class_dict = {'aromatic' : ['F', 'W', 'Y'], 'polar' : ['Q', 'N', 'S', 'T'], 'positive' : ['K', 'R'], 'negative' : ['D', 'E'], 'hydrophobic' : ['I', 'V', 'L', 'A', 'M']}

    #residues that can't be changed or are useless to change
    dont_change = ['D', 'E', 'G', 'P', 'H', 'C']

    # empty string to hold built sequence
    build_sequence = ''

    # get starting hydropathy
    starting_hydro = Protein(sequence).hydropathy

    # decide whether to increase or decrease hydropathy
    if starting_hydro > objective_hydropathy:
        change_hydro = 'decrease_hydropathy'
    else:
        change_hydro = 'increase_hydropathy'

    # figure out the ideal difference between the residues to change.
    total_hydro = starting_hydro * len(sequence)
    total_objective_hydro = objective_hydropathy * len(sequence)
    ideal_residue_difference = abs(total_hydro - total_objective_hydro)
    
    # iterate through sequence
    for cur_aa_num in range(0, len(sequence)):
        cur_aa = sequence[cur_aa_num]
        potential_res=[]
        if cur_aa_num != target_residue_index:
            best_residue = cur_aa
            potential_res.append(best_residue)
        else:
            if cur_aa in dont_change:
                best_residue = cur_aa
                potential_res.append(best_residue)
            else:
                amino_acid_class = AminoAcid.return_AA_class(cur_aa)
                cur_AA_hydro = AminoAcid.hydro(cur_aa)
                if change_hydro == 'decrease_hydropathy':

                    amino_acid_class = AminoAcid.return_AA_class(cur_aa)
                    best_hydropathy = cur_AA_hydro
                    best_residue = cur_aa
                    best_difference = abs(ideal_residue_difference - cur_AA_hydro)

                    for possible_residues in aa_class_dict[amino_acid_class]:
                        cur_poss_hydro = AminoAcid.hydro(possible_residues)
                        if cur_poss_hydro < cur_AA_hydro:
                            potential_res.append(possible_residues)
                else:
                    amino_acid_class = AminoAcid.return_AA_class(cur_aa)
                    best_hydropathy = cur_AA_hydro
                    best_residue = cur_aa
                    best_difference = abs(ideal_residue_difference - cur_AA_hydro)
                    for possible_residues in aa_class_dict[amino_acid_class]:
                        cur_poss_hydro = AminoAcid.hydro(possible_residues)
                        if cur_poss_hydro > cur_AA_hydro:
                            potential_res.append(possible_residues)
        if potential_res == []:
            potential_res.append(cur_aa)
        # randomly select residue.
        build_sequence += random.choice(potential_res)
    # return the sequence
    return build_sequence  


def optimize_hydropathy_within_class(sequence, objective_hydropathy, allowed_hydro_error = parameters.HYDRO_ERROR):
    '''
    This function will optimie the hydropathy of a sequence such that it is
    within the allowed_hydro_error value.
    Will not change the classes of residues in the sequence or their position.

    parameters
    ----------
    sequence : string 
        amino acid sequence as a string

    objective_hydropathy : Float
        objective hydropathy as a float value. Between 0 and 9.

    allowed_error : float
        the allowed amount of error between the objective hydropathy and the final hydropathy

    returns 
    -------
    new_sequence : string
        returns the new (optimzed) sequence as a string

    '''
    # first figure out if possible
    possible_hydro_range = hydro_range_constant_classes(sequence)
    # make sure objective hydropatyh within range of what is possible
    # as far as changing residues within the same class.
    if objective_hydropathy < possible_hydro_range[0] or objective_hydropathy > possible_hydro_range[1]:
        error_message = (f'\n\nUnable to get to objective hydropathy without changing classes of residues.\nFor this sequence the lowest possible hydrpathy is {possible_hydro_range[0]}.\nFor this sequence the highest possible hydropathy is {possible_hydro_range[1]}.\n')
        raise GooseInputError(error_message)
    else:
        # iterate over every residue in the sequence as necessary
        for amino_acid in range(0, len(sequence)):
            new_sequence = residue_optimize_hydropathy_within_class(sequence, objective_hydropathy, amino_acid)
            cur_hydro = Protein(new_sequence).hydropathy
            if abs(cur_hydro - objective_hydropathy) <= allowed_hydro_error:
                return new_sequence
            sequence = new_sequence
    # if iterations didn't get within the ideal value, return the sequence
    return new_sequence


def get_charge_locs(sequence):
    '''
    function that returns the locations of residues
    in a sequence that are charegd. 
    Returns dict of negative locations and positive locationcs

    parameter
    ----------
    sequence : string
        amino acid sequence as a string

    returns
    -------
    dict
        returns a dictionary with the key 'negative'
        corresponding to the location of negatively
        charged residues a 'positive' for the location
        of positively charged residues from the input
        sequence.

    '''

    negatives = []
    positives = []

    for aa_ind, aa in enumerate(sequence):
        if aa == 'K' or aa =='R':
            positives.append(aa_ind)
        elif aa == 'D' or aa == 'E':
            negatives.append(aa_ind)
    return {'negative':negatives, 'positive':positives}


#=/=/=/=/=/=/=/=/=/=/=/=/=/=/
#VARIANT SEQUENCE GENERATORS
#=/=/=/=/=/=/=/=/=/=/=/=/=/=/
'''
Below are the actual sequence generators.
These just generate the sequence.
They do not check for disorder.
'''

def create_new_var_constant_class_nums(sequence):
    '''
    function that takes an input sequence and returns 
    a sequence with the same properties generated from 
    residues that match the number for the number in that
    specific class in the input sequnece

    parameters
    ----------
    sequence : string
        amino acid sequence as a string

    returns
    -------
    build_sequence : string
        returns the final build sequence variant
    '''

    # define which residues in each class
    aromatics = ['F', 'W', 'Y']
    polar = ['Q', 'N', 'S', 'T']
    hydrophobics = ['I', 'V', 'L', 'A', 'M']
    positive = ['K', 'R']
    negative = ['D', 'E']

    # get starting hydropathy of sequence
    starting_hydropathy = Protein(sequence).hydropathy

    # get starting kappa
    starting_kappa = Protein(sequence).kappa

    # get the amount of each class
    residues_by_class = return_num_for_class(sequence)

    # empty string to build sequence on
    build_sequence = ''

    # make list to go through as generating sequence to choose class at random
    classes_list = []
    # add classes to list
    for classval in residues_by_class.keys():
        num_times = residues_by_class[classval]
        if num_times != 0:
            for i in range(0, num_times):
                classes_list.append(classval)
    # shuffle list
    random.shuffle(classes_list)
    
    # build seq
    for cur_class in classes_list:
        build_sequence += get_optimal_res_within_class(cur_class, build_sequence)
    
    # fix hydropathy
    build_sequence = optimize_hydropathy_within_class(build_sequence, objective_hydropathy=starting_hydropathy)

    # fix kappa
    build_sequence = create_kappa_variant(build_sequence, starting_kappa)

    return build_sequence


def create_constant_class_variant(sequence):
    '''
    function to generate a variant with the same properties as the 
    input variant as well as the same order of amino acids as
    far as class and the same number in each class

    parameters
    ----------
    sequence : string
        amino acid sequence as a string

    returns
    -------
    final_seq : string
        returns the final build sequence variant

    '''
    # build sequence to optimize disorder
    build_sequence = ''
    sequence_placeholder=''
    place_neg_charged=[]
    place_pos_charged=[]
    place_aro=[]
    aromatics=['Y', 'W', 'F']
    neg_charged=['D', 'E']
    pos_charged=['K', 'R']
    for amino_acid in sequence:
        cur_class = AminoAcid.return_AA_class(amino_acid)
        if len(build_sequence)>1:
            if build_sequence.count('I')/len(build_sequence) > 0.1:
                potential_residues = get_optimal_res_within_class(cur_class, build_sequence, additional_exclusion=['I'], return_all = True)
                exclude_I=True
            else:
                potential_residues = get_optimal_res_within_class(cur_class, build_sequence, additional_exclusion=[], return_all = True)
                exclude_I=False
        else:
            potential_residues = get_optimal_res_within_class(cur_class, build_sequence, additional_exclusion=[], return_all = True)
            exclude_I=False
        if len(potential_residues) == 1: 
            chosen_residue=potential_residues[0]
        else:
            if amino_acid in potential_residues:
                if amino_acid == 'I':
                    if exclude_I == False:
                        potential_residues.remove(amino_acid)
                else:
                    potential_residues.remove(amino_acid)
            chosen_residue = random.choice(potential_residues)
        
        # add to the sequence we are using to build the base sequence and choose optimal residues
        build_sequence+=chosen_residue
        
        # if negative positive or aromatic add to list, update placeholder
        # sequence we will use later.
        if cur_class=='negative':
            sequence_placeholder+='-'
            place_neg_charged.append(chosen_residue)
        elif cur_class=='positive':
            sequence_placeholder+='+'
            place_pos_charged.append(chosen_residue)
        elif cur_class=='aromatic':
            sequence_placeholder+='0'
            place_aro.append(chosen_residue)
        else:
            sequence_placeholder+='_'
    
    # shuffle the residues to place if aromatic, negative, or positive
    random.shuffle(place_aro)
    random.shuffle(place_pos_charged)
    random.shuffle(place_neg_charged)
    
    # place the residues
    final_sequence=''
    for aa_ind, aa in enumerate(sequence_placeholder):
        if aa == '0':
            final_sequence+=place_aro.pop()
        elif aa == '-':
            final_sequence+=place_neg_charged.pop()
        elif aa == '+':
            final_sequence+=place_pos_charged.pop()
        else:
            final_sequence+=build_sequence[aa_ind]
    
    # now add back in +/- charged residues and aromatics
    # now correct hydropathy
    starting_hydro = Protein(sequence).hydropathy
    final_seq = optimize_hydropathy_within_class(final_sequence, starting_hydro)
    return final_seq


def create_new_variant(sequence):
    '''
    function to generate a variant that is completely different
    in sequence to the input but has all the same overall parameters.
    Does not account for specific classes of residues.

    Getting the SCD right is a little hacky, but this works fast and does a good
    job on the disorder.


    parameters
    ----------
    sequence : string
        amino acid sequence as a string

    returns
    -------
    final_seq : string
        returns the final build sequence variant    
    '''
    input_length = len(sequence)
    input_FCR = Protein(sequence).FCR
    input_NCPR = Protein(sequence).NCPR
    input_hydro = Protein(sequence).hydropathy
    new_sequence = create_seq_by_props(input_length, FCR=input_FCR, NCPR=input_NCPR, hydropathy=input_hydro)
    
    # now fix charged locations to keep kappa / scd constant
    charge_locs = get_charge_locs(sequence)
    neg_charged = []
    pos_charged = []
    charged_res = ['D', 'E', 'K', 'R']
    no_charges = ''
    for aa in new_sequence:
        if aa not in charged_res:
            no_charges += aa
        else:
            if aa == 'D' or aa == 'E':
                neg_charged.append(aa)
            else:
                pos_charged.append(aa)
    # now rebuild seq
    final_seq = ''
    precursor_seq_ind = 0
    for aa_ind in range(0, len(sequence)):
        if aa_ind in charge_locs['negative']:
            final_seq += neg_charged.pop()
        elif aa_ind in charge_locs['positive']:
            final_seq += pos_charged.pop()
        else:
            final_seq += no_charges[precursor_seq_ind]
            precursor_seq_ind += 1

    return final_seq



def create_hydropathy_class_variant(sequence, hydro, allowed_hydro_error = parameters.HYDRO_ERROR):
    '''
    function to take in a sequence and make a variant that adjusts the
    hydropathy while keeping the position and nuimber of amino acids the
    same by class of amino acid

    parameters
    ----------
    sequence : string
        amino acid sequence as a string
    hydro : float
        the hydropathy value you want to change to as a float
    allowed_hydro_error : float
        the allowed error between the objective hydropathy (hydro)
        and the hydropathy of the returned sequence

    returns
    -------
    final_seq : string
        returns the final build sequence variant
    '''
    # change the hydropathy using the optimization function
    final_seq = optimize_hydropathy_within_class(sequence, hydro, allowed_hydro_error)
    return final_seq


def create_constant_residue_variant(sequence, constant_residues = []):
    '''
    function that will generate a new sequence variant
    where specific residues are held constant. The 
    variant will have the same aggregate properties
    as the original sequence.

    Paramters
    ----------
        sequence : string
            The amino acid sequence to make a variant of as a string

        constant_residues : list
            A list of residues to hold constant in the sequence variant

    returns
    -------
    rebuilt_sequence : string
        returns the final sequence as a string

    '''
    # make sure sequence is uppercase (this would otherwise wreck this function)
    sequence = sequence.upper()

    # if no constant residues specified, just generate a sequence and return it
    if constant_residues == []:
        raise GooseInputError('Please specify residues to hold constant. constant_residues cannot be equal to an empty list.')

    # otherwise, get to work!
    # first make sure all residues are capitalized and no overlapping residues.
    upper_case_residues = []
    for conres in constant_residues:
        if conres.upper() not in upper_case_residues:
            upper_case_residues.append(conres.upper())

    # now overwrite constant residues with uppercase version
    constant_residues = upper_case_residues

    # if K or R are in the constant variant, need to pull out both to make sure
    # they aren't altered during variant generation. Also need to do for D and E.
    if 'K' in constant_residues and 'R' not in constant_residues:
        if 'R' not in constant_residues:
            constant_residues.append('R')

    if 'R' in constant_residues:
        if 'K' not in constant_residues:
            constant_residues.append('K')

    if 'D' in constant_residues:
        if 'E' not in constant_residues:
            constant_residues.append('E')

    if 'E' in constant_residues:
        if 'D' not in constant_residues:
            constant_residues.append('D')


    # next identify the residues and their positions in sequence
    # makes a list that has position then amino acid such that positions can
    # be searched in the list and amino acid added to seq later easily
    positions_then_res = []
    for res in constant_residues:
        if res in sequence:
            curpositions = identify_residue_positions(sequence, res)
            for pos in curpositions:
                positions_then_res.append(pos)
                positions_then_res.append(res)

    # get original positions of positive and negative residues
    positives_and_negatives_pos_then_sign = []
    charged = ['E', 'D', 'K', 'R']
    for res in charged:
        if res in sequence:
            if res not in constant_residues:
                if res == 'E' or res == 'D':
                    sign = '-'
                else:
                    sign = '+'
                positions_of_res = identify_residue_positions(sequence, res)
                for pos in positions_of_res:
                    positives_and_negatives_pos_then_sign.append(pos)
                    positives_and_negatives_pos_then_sign.append(sign)

    # now strip the residues that must be held in the same position
    seq_variant = ''
    for res in sequence:
        if res not in constant_residues:
            seq_variant += res

    # Now calculate sequence backbone paramters
    input_FCR    = Protein(seq_variant).FCR
    input_NCPR   = Protein(seq_variant).NCPR 
    input_hydro  = Protein(seq_variant).hydropathy 
    input_length = len(seq_variant)


    # now make the sequence
    seq_variant = create_seq_by_props(length = input_length, FCR=input_FCR, NCPR = input_NCPR, hydropathy = input_hydro, exclude=constant_residues)

    # now strip the charged residues from seq_variant
    input_backbone = ''
    negative_residues_for_var = []
    positive_residues_for_var = []
    
    for res in seq_variant:
        if res in charged:
            if res == 'D' or res == 'E':
                negative_residues_for_var.append(res)
            else:
                positive_residues_for_var.append(res)
        else:
            input_backbone += res

    # now rebuild the sequence
    rebuilt_sequence = ''

    seq_variant_pos = 0

    for i in range(0, len(sequence)):
        if i in positives_and_negatives_pos_then_sign:
            cur_pos = positives_and_negatives_pos_then_sign.index(i)
            if positives_and_negatives_pos_then_sign[cur_pos+1] == '-':
                cur_res = negative_residues_for_var.pop()
            else:
                cur_res = positive_residues_for_var.pop()
        elif i in positions_then_res:
            cur_pos = positions_then_res.index(i)
            cur_res = positions_then_res[cur_pos+1]
        else:
            cur_res = input_backbone[seq_variant_pos]
            seq_variant_pos += 1
        # add the necessary residue
        rebuilt_sequence += cur_res

    return rebuilt_sequence



def seq_chunks_from_regions_list(sequence, regions=[]):
    '''
    function that takes in specified regions and returns
    a list of sequences

    parameters
    ----------
    sequence : str
        amino acid sequenc as a string

    regions : list of lists
        list of lists where sublists specify regions to break up 
        into chunks

    returns 
    -------
    complete_list : list
        returns a list of lists that cover all reagions of the
        input sequence as opposed to just the regions that 
        you want to change.
    '''

    # list to add additional lists to to close gaps
    complete_list = []

    # make sure list is list of lists
    if type(regions[0]) == int:
        regions = [regions]
    

    # iterate through lists to make sublists
    for sublist in range(0, len(regions)):
        curlist = regions[sublist]
        if sublist == len(regions)-1:
            if len(regions) == 1:
                if curlist[0] != 0:
                    complete_list.append([0, curlist[0]])
            complete_list.append(curlist)
            if curlist[1] != len(sequence):
                complete_list.append([curlist[1], len(sequence)])
        else:
            if sublist == 0:
                if curlist[0] != 0:
                    complete_list.append([0,curlist[0]])
            nextlist = regions[sublist+1]
            complete_list.append(curlist)
            if curlist[1]!=nextlist[0]:
                complete_list.append([curlist[1],nextlist[0]])

    return complete_list



def create_region_shuffle_variant(sequence, shuffle_regions=[], use_index=False):
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

    Returns
    -------
    shuffled_seq : String
        Returns your sequence with the specified regions shuffled.
    '''

    # now make a nested list to keep things from getting annoying
    if type(shuffle_regions[0])==int:
        shuffle_regions = [shuffle_regions]

    # now convert to index values for the computering peeps
    # also modifies values such that the range includes specified residues
    # unlike default behavior where the final residue number is not included
    if use_index == False:
        full_list = []
        for i in shuffle_regions:
            full_list.append([i[0]-1, i[1]])
        shuffle_regions = full_list

    # get all regions
    all_regions = seq_chunks_from_regions_list(sequence, regions = shuffle_regions)

    #build the seq
    shuffled_sequence = ''

    for subregion in all_regions:
        curblock = sequence[subregion[0]:subregion[1]]
        if subregion in shuffle_regions:
            shuffled_sequence += ''.join(random.sample(curblock, len(curblock)))
        else:
            shuffled_sequence += curblock

    # return the sequence
    return shuffled_sequence


# function to identify residue closest tot eh middle (from list or otherwise)
def find_middle_res(sequence, residues):
    '''
    function to find the residue from residue(s)
    of interest that are closest to the middle of a sequence
    
    Parameters
    -----------
    sequence : string 
        amino acid sequence as a string

    residues : list or string
        can be a list of residues or a single residue as a string

    Returns
    -------
    approx_middle : int
        returns thte index of the resiude approximately in the middle

    '''

    # make sure consistnetly working with a list
    if type(residues) == str:
        residues = [residues]

    # find approximate middle of sequence
    approx_middle = int(len(sequence)/2)
    # if the approx middle residue is the res of interest, great
    if sequence[approx_middle] in residues:
        return approx_middle
    else:
        for i in range(0, approx_middle):
            cur_residue_1 = sequence[approx_middle+i]
            cur_residue_2 = sequence[approx_middle-i]
            if cur_residue_1 in residues:
                return approx_middle+i
            if cur_residue_2 in residues:
                return approx_middle-i

        # if we still have nothing, we have to return an end value

        if sequence[0] in residues:
            return 0
        elif sequence[len(sequence)] in residues:
            return len(sequence)
        else:
            raise GooseBackendBug('Error in finding middle residue in find_middle_res in variant_generation_backend.py')
    raise GooseBackendBug('Error in finding middle residue, made to end of function in find_middle_res in variant_generation_backend.py')


def random_decrease_kappa(sequence):
    '''
    Function to quickly decrease kappa. 

    Parameters
    ----------
    sequence : str
        A string of amino acids.
    Returns
    -------
    sequence : str
        A string of amino acids with kappa increased (usually).
    '''
    # calculate linear ncpr
    lin_ncpr=[]
    seq_coords=[]
    if len(sequence)<20:
        blob_len=random.choice([3,4,5,6,7])
    else:
        blob_len=random.choice([3,4,5,6,7,8,9])
    for i in range(0, len(sequence)-(blob_len+1)):
        blob = sequence[i:i+blob_len]
        seq_coords.append([i, i+blob_len])
        lin_ncpr.append(pr(blob).NCPR)
    # get the min and max ncpr regions
    min_ncpr = min(lin_ncpr)
    max_ncpr = max(lin_ncpr)
    possible_min_ncprs=[]
    possible_max_ncprs=[]
    taken_range_min=[]
    taken_range_max=[]
    for ncpr_ind, ncpr in enumerate(lin_ncpr):
        if ncpr==min_ncpr:
            if ncpr_ind not in taken_range_min:
                possible_min_ncprs.append(ncpr_ind)
                taken_range_min.extend([ncpr_ind, ncpr_ind+blob_len])
        if ncpr==max_ncpr:
            if ncpr_ind not in taken_range_max:
                possible_max_ncprs.append(ncpr_ind)
                taken_range_max.extend([ncpr_ind, ncpr_ind+blob_len])

    min_ncpr_coords=seq_coords[random.choice(possible_min_ncprs)]
    max_ncpr_corods=seq_coords[random.choice(possible_max_ncprs)]

    # get residiues of interest for thoses regions. 
    min_ncpr_residues = sequence[min_ncpr_coords[0]:min_ncpr_coords[1]]
    max_ncpr_residues = sequence[max_ncpr_corods[0]:max_ncpr_corods[1]]
    
    # swap a res between them or move one res. 
    neg_res_in_min_neg_ncpr=[]
    for num, res in enumerate(min_ncpr_residues):
        if res in ['D', 'E']:
            neg_res_in_min_neg_ncpr.append(num)
    
    pos_res_in_max_pos_ncpr=[]
    for num, res in enumerate(max_ncpr_residues):
        if res in ['K', 'R']:
            pos_res_in_max_pos_ncpr.append(num)

    if neg_res_in_min_neg_ncpr==[]:
        swap_from_neg=random.randint(0, len(min_ncpr_residues)-1)
    else:
        swap_from_neg=random.choice(neg_res_in_min_neg_ncpr)

    if pos_res_in_max_pos_ncpr==[]:
        swap_from_pos=random.randint(0, len(max_ncpr_residues)-1)
    else:
        swap_from_pos=random.choice(pos_res_in_max_pos_ncpr)

    # get index of resiudes to swap based on regions. 
    swap_from_neg_seq_index = swap_from_neg+min_ncpr_coords[0]
    swap_from_pos_seq_index = swap_from_pos+max_ncpr_corods[0]
    # make sure we don't swap the same res for the same res. 
    # if this happens, we could accidentally duplicate a res. 
    if swap_from_pos_seq_index==swap_from_neg_seq_index:
        swap_from_pos_seq_index=random.choice([aa for aa in range(0, len(sequence)-blob_len) if aa != swap_from_neg])
    #
    swap_from_neg_seq_to_pos = sequence[swap_from_neg_seq_index]
    swap_from_pos_seq_to_neg = sequence[swap_from_pos_seq_index]

    # rebuild sequence
    final_sequence=''
    for resnum in range(0, len(sequence)):
        if resnum==swap_from_neg_seq_index:
            final_sequence+=swap_from_pos_seq_to_neg
        elif resnum==swap_from_pos_seq_index:
            final_sequence+=swap_from_neg_seq_to_pos
        else:
            final_sequence+=sequence[resnum]
    return final_sequence


def random_increase_kappa(sequence, objective_kappa):
    '''
    Function to slowly increase kappa. 

    Parameters
    ----------
    sequence : str
        A string of amino acids.
    Returns
    -------
    sequence : str
        A string of amino acids with kappa increased (usually).
    '''
    # function to slowly increase kappa
    pos_res_of_interest = ('K', 'R')
    neg_res_of_interest = ('D', 'E')
    charged_res = pos_res_of_interest + neg_res_of_interest
    other_res_ind = [aa for aa in range(0, len(sequence)) if sequence[aa] not in charged_res]

    gen_seqs = [sequence]
    gen_kappa = [pr(sequence).kappa]

    # figure out non-charged residues to move. 
    approx_center = int(len(sequence)/2)
    furthest_from_center = 0
    for aa_ind in other_res_ind:
        if abs(approx_center-aa_ind) > furthest_from_center:
            furthest_from_center = abs(approx_center-aa_ind)
            furthest_from_center_ind = aa_ind
    
    if furthest_from_center != 0:
        target_other_res_id = sequence[furthest_from_center_ind]
        target_other_res_ind = furthest_from_center_ind
    else:
        target_other_res_id = 0
        target_other_res_ind = -1

    swap_with_ind = random.choice([aa for aa in range(int(0.3*len(sequence)), int(0.7*len(sequence))) if aa != target_other_res_ind])
    swap_with_res=sequence[swap_with_ind]

    if target_other_res_ind != -1:
        mod_seq=''
        for aa_ind in range(0, len(sequence)):
            if aa_ind == swap_with_ind:
                mod_seq+=target_other_res_id
            elif aa_ind == target_other_res_ind:
                mod_seq+=swap_with_res
            else:
                mod_seq+=sequence[aa_ind]
        gen_seqs.append(mod_seq)
        gen_kappa.append(abs(objective_kappa-pr(mod_seq).kappa))
        sequence = mod_seq
        

    pos_res_ind=[aa for aa in range(0, len(sequence)) if sequence[aa] in pos_res_of_interest]
    neg_res_ind=[aa for aa in range(0, len(sequence)) if sequence[aa] in neg_res_of_interest]
    if pos_res_ind == [] or neg_res_ind == []:
        return sequence

    # get location of where pos res are and neg res are.
    pos_res_weight = sum(pos_res_ind)/len(pos_res_ind)
    neg_res_weight = sum(neg_res_ind)/len(neg_res_ind)

    if pos_res_weight==neg_res_weight:
        # randomly choose.
        if random.randint(0,1)==0:
            pos_res_weight+=0.0001
        else:
            neg_res_weight+=0.0001


    # get start loc for neg and pos
    if len(pos_res_ind) > int(len(sequence)/2):
        start_pos_res = len(pos_res_ind)
        start_neg_res = len(sequence)-start_pos_res
    else:
        start_pos_res = int(len(sequence)/2)

    if len(neg_res_ind) > int(len(sequence)/2):
        start_neg_res = len(neg_res_ind)
        start_pos_res = len(sequence)-start_neg_res
    else:
        start_neg_res = int(len(sequence)/2)

    if pos_res_weight > neg_res_weight:
        target_pos_res_ind = min(pos_res_ind)
        target_neg_res_ind = max(neg_res_ind)
        avoid_res = [target_pos_res_ind, target_neg_res_ind]        
        pos_res_target_locs = [aa for aa in range(len(sequence)-start_pos_res, len(sequence)) if aa not in avoid_res]
        avoid_res.extend(pos_res_target_locs)
        neg_res_target_locs = [aa for aa in range(0, start_neg_res) if aa not in avoid_res]

    else:
        target_pos_res_ind = max(pos_res_ind)
        target_neg_res_ind = min(neg_res_ind)
        avoid_res = [target_pos_res_ind, target_neg_res_ind]     
        pos_res_target_locs = [aa for aa in range(0, start_pos_res) if aa not in avoid_res]
        avoid_res.extend(pos_res_target_locs)
        neg_res_target_locs = [aa for aa in range(len(sequence)-start_neg_res,len(sequence)) if aa not in avoid_res]


    target_pos_res = sequence[target_pos_res_ind]
    target_neg_res = sequence[target_neg_res_ind]

    # if either res_target_locs == [], swap pos and neg charges
    if pos_res_target_locs == [] or neg_res_target_locs == []:
        avoid_res = [target_pos_res_ind, target_neg_res_ind]
        if pos_res_weight > neg_res_weight:
            final_pos_res_loc = min([aa for aa in range(0, int(len(sequence)/2)) if aa not in avoid_res])
            final_neg_res_loc = max([aa for aa in range(int(len(sequence)/2), len(sequence)) if aa not in avoid_res])
        else:
            final_neg_res_loc = min([aa for aa in range(0, int(len(sequence)/2)) if aa not in avoid_res])
            final_pos_res_loc = max([aa for aa in range(int(len(sequence)/2), len(sequence)) if aa not in avoid_res])
    else:

        final_pos_res_loc = random.choice(pos_res_target_locs)
        final_neg_res_loc = random.choice(neg_res_target_locs)

    charge_mod=''
    for aa_ind in range(0, len(sequence)):
        if aa_ind == final_pos_res_loc:
            charge_mod += target_pos_res
        elif aa_ind == final_neg_res_loc:
            charge_mod += target_neg_res
        elif aa_ind == target_pos_res_ind:
            charge_mod += sequence[final_pos_res_loc]
        elif aa_ind == target_neg_res_ind:
            charge_mod += sequence[final_neg_res_loc]
        else:
            charge_mod += sequence[aa_ind]

    return charge_mod


def increase_kappa_to_val(sequence, objective_kappa, sub_attempts=5000, objective_error=0.01):
    '''
    Function to increase kappa to a value. 

    Parameters
    ----------
    sequence : str
        A string of amino acids.
    Returns
    -------
    sequence : str
        A string of amino acids with kappa increased (usually).
    '''
    best_kappa=pr(sequence).kappa
    best_diff = abs(objective_kappa-best_kappa)
    if best_diff < objective_error:
        return sequence
    best_seq=sequence
    max_consec_fails=100
    consec_fails=0
    unique_kappas={}
    for i in range(0, sub_attempts):
        new_seq = random_increase_kappa(best_seq, objective_kappa)
        new_kappa=pr(new_seq).kappa
        cur_diff = abs(objective_kappa-new_kappa)
        if cur_diff<best_diff:
            if cur_diff < objective_error:
                return new_seq
            unique_kappas[cur_diff]=new_seq
            best_kappa = new_kappa
            best_seq = new_seq 
            best_diff=cur_diff
            consec_fails=0
        else:
            consec_fails+=1
        if consec_fails > max_consec_fails:
            temp=list(best_seq)
            random.shuffle(temp)
            best_seq=''.join(temp)
            consec_fails=0
            best_kappa=pr(best_seq).kappa
            best_diff=abs(objective_kappa-best_kappa)

    # return the best thing we made the whole time...
    return unique_kappas[min(unique_kappas.keys())]


def decrease_kappa_to_val(sequence, objective_kappa, sub_attempts=5000, objective_error=0.01):
    '''
    Function to decrease kappa to a value. 

    Parameters
    ----------
    sequence : str
        A string of amino acids.
    Returns
    -------
    sequence : str
        A string of amino acids with kappa increased (usually).
    '''
    best_kappa=pr(sequence).kappa
    best_diff = abs(objective_kappa-best_kappa)
    if best_diff < objective_error:
        return sequence
    best_seq=sequence
    max_consec_fails=50
    consec_fails=0
    unique_kappas={}
    for i in range(0, sub_attempts):
        new_seq = random_decrease_kappa(best_seq)
        new_kappa=pr(new_seq).kappa
        cur_diff = abs(objective_kappa-new_kappa)
        if cur_diff<best_diff:
            if cur_diff < objective_error:
                return new_seq
            unique_kappas[cur_diff]=new_seq
            best_kappa = new_kappa
            best_seq = new_seq 
            best_diff=cur_diff
            consec_fails=0
        else:
            consec_fails+=1
        if consec_fails > max_consec_fails:
            temp=list(best_seq)
            random.shuffle(temp)
            best_seq=''.join(temp)
            consec_fails=0
            best_kappa=pr(best_seq).kappa
            best_diff=abs(objective_kappa-best_kappa)

    # return the best thing we made the whole time...
    return unique_kappas[min(unique_kappas.keys())]

def check_seq_comp(seq1,seq2):
    amino_acids=['A', 'C', 'D', 'E', 'F', 'G', 'H', 'I', 'K', 'L', 'M', 'N', 'P', 'Q', 'R', 'S', 'T', 'V', 'W', 'Y']
    for aa in amino_acids:
        if seq1.count(aa)!=seq2.count(aa):
            raise Exception('Different sequence compositions!')


def create_kappa_variant(sequence, kappa, allowed_kappa_error = parameters.MAXIMUM_KAPPA_ERROR, attempts=5):
    '''
    Function to generate a sequence with a user-defined
    kappa value. Requires kappa calculation using 
    SPARROW.

    parameters
    -----------
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
    '''

    num_iterations = 20*attempts
    new_seq = optimize_kappa_vectorized([sequence], target_kappa=kappa, 
                              tolerance=allowed_kappa_error, num_iterations=num_iterations, 
                              max_change_iterations=20,
                              early_stopping=True,
                              return_when_num_hit=20,
                              only_return_within_tolerance=True)
    
    # if we have any sequences returned, return them. 
    if len(new_seq)>1:
        return new_seq[0]
    elif len(new_seq)==1:
        if isinstance(new_seq, list):
            return new_seq[0]
    else:
        return new_seq

    # set number sub attempts
    sub_attempts=5000

    # if user overrides the minimium allowed error, use that instead of the objective error. 
    # otherwise we try to get within 0.01 but will tolerate up to 0.03.
    if kappa<0.6:
        objective_error = min([allowed_kappa_error, 0.02])
    else:
        objective_error = min([allowed_kappa_error, 0.01])

    # make the best sequence the current sequence
    best_sequence=sequence
    current_kappa = pr(best_sequence).kappa
    best_kappa = current_kappa
    best_error = abs(best_kappa - kappa)
    if best_error < allowed_kappa_error:
        return sequence

    # start trying to make sequence
    for attempt in range(0, attempts):
        if current_kappa > kappa:
            new_seq = decrease_kappa_to_val(best_sequence, objective_kappa=kappa, sub_attempts=sub_attempts, 
                objective_error=objective_error)
            current_kappa = pr(new_seq).kappa
            if abs(current_kappa - kappa) < allowed_kappa_error:
                try:
                    check_seq_comp(new_seq, sequence)
                except:
                    raise GooseBackendBug('Sequence composition changed when increasing kappa using decrease_kappa_to_val')                
                return new_seq
            if abs(current_kappa - kappa) < best_error:
                best_sequence = new_seq
                best_kappa = current_kappa
                best_error = abs(best_kappa-kappa)
        if current_kappa < kappa:
            new_seq = increase_kappa_to_val(best_sequence, objective_kappa=kappa, sub_attempts=sub_attempts, 
                objective_error=objective_error)
            current_kappa = pr(new_seq).kappa
            if abs(current_kappa - kappa) < allowed_kappa_error:
                try:
                    check_seq_comp(new_seq, sequence)
                except:
                    raise GooseBackendBug('Sequence composition changed when increasing kappa using increase_kappa_to_val')
                return new_seq
            if abs(current_kappa - kappa) < best_error:
                best_sequence = new_seq
                best_kappa = current_kappa
                best_error = abs(best_kappa-kappa)            

    fail_message = 'Unable to make kappa sequence in variant_generation_backend.py'        
    raise GooseBackendBug(fail_message)


def radiate_out(starting_value):
    # radiates outward from starting value to get value as close to original as possible
    closest_vals = []
    for i in np.arange(0.02, 1, 0.005):
        i = float(i)
        val_1 = round(starting_value - i, 5)
        val_2 = round(starting_value + i, 5)
        if val_1 >= 0 and val_1 <= 1:
            if val_1 not in closest_vals:
                closest_vals.append(val_1)
        if val_2 >= 0 and val_2 <= 1:
            if val_2 not in closest_vals:
                closest_vals.append(val_2)
    return closest_vals


def create_closest_kappa_variant(sequence, ideal_value):
    if ideal_value == -1:
        return sequence
    else:
        # makes the closest kappa variant possible
        all_values = radiate_out(ideal_value)
        for i in all_values:
            try:
                seq = create_kappa_variant(sequence, kappa=i, attempts=2000)
                return seq
            except:
                continue
    raise GooseBackendBug('Unable to use create_closest_kappa_variant in variant_generation_backend.py')


def possible_fcr_vals(sequence):
    # figure out possible FCR vals without changing NCPR
    all_fcr_vals = []
    fractional_val = 1/len(sequence)
    minimal_fcr = abs(Protein(sequence).NCPR)
    num_vals = round(1/fractional_val)+2
    for i in range(0, num_vals):
        curval = round(minimal_fcr + (i*fractional_val), 6)
        if curval < 1:
            all_fcr_vals.append(curval)
    return all_fcr_vals

def get_closest_fcr(sequence, input_fcr_val):
    all_fcr_vals = possible_fcr_vals(sequence)
    closest = 1000
    for val in all_fcr_vals:
        if abs(val-input_fcr_val) < closest:
            closest = abs(val-input_fcr_val)
            best_val = val
    return best_val



def count_charged(sequence):
    '''
    mini fucntion to count charged residues
    '''
    num_pos = sequence.count('K') + sequence.count('R')
    num_neg = sequence.count('D') + sequence.count('E')
    return{'negative':num_neg, 'positive':num_pos}


def create_fcr_class_variant(sequence, fraction, constant_ncpr=True, use_closest = True, ignore_kappa=False):
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

    fraction : float
        the FCR value between 0 and 1

    constant_ncpr : Bool
        whether to allow changes to the NCPR when changing the FCR

    use_closest : Bool
        whether to just use the closest FCR val to the input value.

    ignore_kappa : bool
        whether to ignore getting the kappa value correct.

    '''
    # get starting hydropathy
    original_hydropathy = Protein(sequence).hydropathy
    original_kappa = pr(sequence).kappa
    original_fcr = Protein(sequence).FCR
    original_ncpr = Protein(sequence).NCPR

    # figure out max possible FCR that can maintain disorder
    max_FCR = calculate_max_charge(original_hydropathy)
    if fraction > max_FCR:
        error_message = f'\n\nThe input FCR of {fraction} is greater than the max possible FCR to keep hydropathy the same and maintain disorder, which is {max_FCR}\n'
        raise GooseInputError(error_message)

    # first figure out best changes to modify charge.
    first_changes = ['N', 'Q']
    second_changes = ['S', 'T', 'G']
    third_changes = ['I', 'A', 'L', 'V', 'M']
    cant_change = ['W', 'Y', 'F', 'P', 'C', 'D', 'E', 'K', 'R', 'H']
    change_in_order = ['N', 'Q', 'S', 'T', 'G', 'A','L','M', 'V', 'I']
    
    # now figure out how many charged residues
    original_charged_dict = count_charged(sequence)
    num_negative = original_charged_dict['negative']
    num_positive = original_charged_dict['positive']
    original_total_charged = num_negative+num_positive

    # figure out closest number of fcr to use without changing ncpr
    if fraction != 0:
        new_closest_fraction = get_closest_fcr(sequence, fraction)
    else:
        new_closest_fraction=0


    # if FCR won't change just return it
    if new_closest_fraction == original_fcr:
        return sequence

    # otherwise get to it
    neg_res_locs = []
    if 'D' in sequence:
        neg_res_locs.extend(identify_residue_positions(sequence, 'D'))
    if 'E' in sequence:
        neg_res_locs.extend(identify_residue_positions(sequence, 'E'))
    if len(neg_res_locs) > 1:
        random.shuffle(neg_res_locs)

    pos_res_locs =[]
    if 'K' in sequence:
        pos_res_locs.extend(identify_residue_positions(sequence, 'K'))
    if 'R' in sequence:
        pos_res_locs.extend(identify_residue_positions(sequence, 'R'))
    if len(pos_res_locs) > 1:
        random.shuffle(pos_res_locs)

    # figure otu total number of residus to change
    num_residues_to_change = abs(original_total_charged - round(new_closest_fraction*len(sequence)))

    # if reducing fcr ...
    if new_closest_fraction < original_fcr:
        # get number of times to take out 1 positive and 1 negative residue
        target_residues=[]
        if new_closest_fraction != 0:
            num_changes = round(num_residues_to_change/2)
            for i in range(0, num_changes):
                target_residues.append(neg_res_locs.pop())
                target_residues.append(pos_res_locs.pop())
        else:
            target_residues = []
            target_residues.extend(pos_res_locs)
            target_residues.extend(neg_res_locs)
        # build the starting sequence
        starter_seq=''
        for aa in range(0, len(sequence)):
            if aa not in target_residues:
                starter_seq += sequence[aa]
            else:
                starter_seq += random.choice(first_changes)
    else:
        residues_to_change = []
        for res in change_in_order:
            if len(residues_to_change) < num_residues_to_change:
                if res in sequence:
                    possible_residues = identify_residue_positions(sequence, res)
                    random.shuffle(possible_residues)
                    for aa in possible_residues:
                        if len(residues_to_change) < num_residues_to_change:
                            residues_to_change.append(aa)

        # just a check to make sure we have found enough residues to target
        if len(residues_to_change) < num_residues_to_change:
            raise GooseInputError('Unable to find enough residues that can be changed without altering H, P, W, F, Y, or C.')

        # build the seq
        starter_seq = ''
        add_a_negative=True
        for aa in range(0, len(sequence)):
            if aa not in residues_to_change:
                starter_seq += sequence[aa]
            else:
                if add_a_negative == True:
                    starter_seq += random.choice(['D', 'E', 'D', 'E'])
                    add_a_negative=False
                else:
                    starter_seq += random.choice(['K', 'R', 'K', 'R'])
                    add_a_negative=True

    # now fix kappa, don't try to fix if val is -1!
    if ignore_kappa==False:
        if pr(starter_seq).kappa != -1:
            fixed_kappa_seq = create_kappa_variant(starter_seq, original_kappa)
        else:
            fixed_kappa_seq = starter_seq
    else:
        fixed_kappa_seq=starter_seq

    # now fix hyropathy
    fixed_hydro_seq = optimize_hydropathy_within_class(fixed_kappa_seq, original_hydropathy)

    # now return the seq
    return fixed_hydro_seq




# function to calculate the linear profile of some sequence parameter
def calculate_linear_sequence_profile(sequence, mode, window_size=6):
    # first break seq up into windows
    seq_windows = []
    num_windows = len(sequence)-window_size+1
    for i in range(0, num_windows):
        seq_windows.append(sequence[i:i+window_size])
    # get vals over windows
    seq_vals = []
    # correct some modes
    if mode == 'fcr':
        mode = 'FCR'
    if mode == 'ncpr':
        mode = 'NCPR'
    if mode == 'hydrophobicity':
        mode == 'hydropathy'
    if mode == 'hydro':
        mode = 'hydropathy'
    # if mode is a list, just calculate the fractions of specified amino acids across
    # the sequence
    if type(mode) == list:
        for seq in seq_windows:
            temp_val = 0
            for i in mode:
                temp_val+=seq.count(i)
            seq_vals.append(round(temp_val/window_size, 5))
    elif mode == 'FCR':
        for seq in seq_windows:
            seq_vals.append(Protein(seq).FCR)
    elif mode == 'NCPR':
        for seq in seq_windows:
            seq_vals.append(Protein(seq).NCPR)    
    elif mode == 'hydropathy':
        for seq in seq_windows:
            seq_vals.append(Protein(seq).hydropathy)
    elif mode == 'aromatic':
        for seq in seq_windows:
            temp_val = 0
            for i in ['W', 'Y', 'F']:
                temp_val+=seq.count(i)
            seq_vals.append(round(temp_val/window_size, 5))
    elif mode == 'aliphatic':
        for seq in seq_windows:
            temp_val = 0
            for i in ['I', 'V', 'L', 'A', 'M']:
                temp_val+=seq.count(i)
            seq_vals.append(round(temp_val/window_size, 5))    
    elif mode == 'polar':
        for seq in seq_windows:
            temp_val = 0
            for i in ['Q', 'N', 'S', 'T']:
                temp_val+=seq.count(i)
            seq_vals.append(round(temp_val/window_size, 5))            
    elif mode == 'proline':
        for seq in seq_windows:
            temp_val = 0
            temp_val+=seq.count('P')
            seq_vals.append(round(temp_val/window_size, 5))     
    elif mode == 'positive':
        for seq in seq_windows:
            temp_val = 0
            for i in ['K', 'R']:
                temp_val+=seq.count(i)
            seq_vals.append(round(temp_val/window_size, 5)) 
    elif mode == 'negative':
        for seq in seq_windows:
            temp_val = 0
            for i in ['D', 'E']:
                temp_val+=seq.count(i)
            seq_vals.append(round(temp_val/window_size, 5)) 
    else:
        raise GooseInputError('Specified mode does not exist. Please input a list of amino acids or a property.')
    return seq_vals


def create_asymmetry_variant_once(sequence, increase_decrease, aa_class, window_size=6):
    '''
    function to change the asymmetry of a rseidue in a sequence

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
    '''

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

    # get the profile of the specified class        
    seq_profile = calculate_linear_sequence_profile(sequence, target_amino_acids, window_size=window_size)
    
    # find regions with values for specified val.
    if increase_decrease == 'increase':
        target_frac_region = max(seq_profile)
    else:
        target_frac_region = min(seq_profile)
    window_index = []
    for i in range(0, len(seq_profile)):
        curval = seq_profile[i]
        if curval == target_frac_region:
            window_index.append(i)
    
    # now convert into seq coordinates.
    seq_coords = []
    for i in window_index:
        seq_coords.append([i, i+window_size])
    
    # now get positions of all other parts of sequence with the specified residue
    residue_positions = []
    for i in target_amino_acids:
        residue_positions.extend(identify_residue_positions(sequence, i))

    # now try to get a residue position not within the various regions
    off_limits = []
    for i in seq_coords:
        for j in range(i[0], i[1]+1):
            if j < len(sequence):
                if sequence[j] in target_amino_acids:
                    if j not in off_limits:
                        off_limits.append(j)

    if off_limits == residue_positions:
        target_residue_pos = residue_positions[random.randint(0, len(off_limits)-1)]
    else:
        potential_targets = []
        for i in residue_positions:
            if i not in off_limits:
                potential_targets.append(i)

        if len(potential_targets) > 1:
            target_residue_pos = potential_targets[random.randint(0, len(potential_targets)-1)]
        elif potential_targets == []:
            return sequence
        else:
            target_residue_pos = potential_targets[0]

    # select random region to add the residue to.
    target_placement_region = seq_coords[random.randint(0, len(seq_coords)-1)]
    target_indices = []
    for i in range(target_placement_region[0], target_placement_region[1]+1):
        if i < len(sequence):
            if i != target_residue_pos:
                target_indices.append(i)
    # select final target placement
    final_target_placement = target_indices[random.randint(0, len(target_indices)-1)]
    # build final sequence
    pulled_residue = sequence[target_residue_pos]
    built_sequence=''
    for i in range(0, len(sequence)):
        if i == final_target_placement:
            built_sequence += pulled_residue
            built_sequence += sequence[i]
        elif i == target_residue_pos:
            built_sequence += ''
        else:
            built_sequence += sequence[i]

    return built_sequence

def total_asym(sequence, residues):
    '''
    quick mini function to take in a list of residues for the 
    clustering_param function and return the total
    asymmetry of those residues for the sequence
    '''
    if type(residues) == str:
        residues = list(residues)
    else:
        cur_val = 0
        for i in residues:
            cur_val += clustering_param(sequence, i)
    return cur_val


def create_asymmetry_variant(sequence, increase_decrease, aa_class, window_size=6, num_change=None):
    '''
    function to change the asymmetry of a rseidue in a sequence

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
    num_change : int
        the number of times to increase / decrease the asymmetry of the residue in the sequence
    '''
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

    # get number of residues in the sequence
    total_res = 0
    for aa in target_amino_acids:
        total_res += sequence.count(aa)
    if total_res == 0:
        return sequence
    # if the number of times not listed, do 1/10 the len of the seq
    if num_change == None:
        num_change = int(total_res/10)
    if num_change == 0:
        num_change = 1
    # keep track of the seqs made.
    list_of_seqs =[]
    newseq = create_asymmetry_variant_once(sequence, increase_decrease, aa_class, window_size=window_size)
    if num_change==1:
        return newseq
    else:
        bestseq=newseq
        cur_asym = total_asym(newseq, target_amino_acids)
        increase_window=False
        for i in range(0, num_change-1):
            newseq = create_asymmetry_variant_once(newseq, increase_decrease, aa_class, window_size=window_size)
            new_asym = total_asym(newseq, target_amino_acids)
            if increase_decrease == 'decrease':
                if new_asym < cur_asym:
                    cur_asym = new_asym
                    bestseq=newseq
            else:
                if new_asym > cur_asym:
                    cur_asym = new_asym
                    bestseq=newseq
            # add new seq to growing list of seqs
            list_of_seqs.append(newseq)
            # if a reallllly high number is listed, eventually the asymmetry wont change and
            # you end up with the same sequence. This will kill the function if that happens too much.
            if list_of_seqs.count(newseq) >= 2:
                # reduce window size by 1 and try again
                if window_size >= 3:
                    if window_size == 3:
                        increase_window = True
                    if increase_window == False:
                        if window_size != 3:
                            window_size = window_size-1
                    else:
                        if window_size < 6:
                            window_size = window_size + 1
                        else:
                            increase_window=False
                    

                    #newseq = create_asymmetry_variant_once(newseq, increase_decrease, aa_class, window_size=window_size)
                    #list_of_seqs.append(newseq)

                if list_of_seqs.count(newseq) >=4:
                    # if you now have gotten the sequence AGAIN
                    return bestseq
    # return the final seq
    return bestseq



def needed_charged_residues(length, fraction, net_charge):
    '''
    function to figure out how many positive and 
    how many negative charged residues are need
    to satisfy the ncpr and fcr of a 
    given sequence.
    '''
    total_charged = fraction_net_charge(length, fraction, net_charge)

    if net_charge > 0:
        num_positive = total_charged['NCPR_residues'] + round(total_charged['FCR_residues'])
        num_negative = round(total_charged['FCR_residues'])
    elif net_charge < 0:
        num_negative = total_charged['NCPR_residues'] + round(total_charged['FCR_residues'])
        num_positive = round(total_charged['FCR_residues'])
    else:
        num_positive = round(total_charged['FCR_residues'])
        num_negative = round(total_charged['FCR_residues'])

    return {'positive':num_positive, 'negative':num_negative}



def create_ncpr_class_variant(sequence, net_charge, constant_fcr=True, use_closest = True, ignore_kappa=False):
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

    net_charge : float
        the NCPR value between -1 and 1

    constant_fcr : Bool
        whether to allow changes to the FCR when changing the NCPR

    use_closest : Bool
        whether to just use the closest FCR val to the input value.

    ignore_kappa : bool
        whether to ignore kappa values

    '''
    # get starting hydropathy, kappa, fcr
    original_protein = Protein(sequence)
    original_hydropathy = original_protein.hydropathy
    original_kappa = original_protein.kappa
    original_FCR = original_protein.FCR

    # figure out max possible FCR that can maintain disorder
    max_FCR = calculate_max_charge(original_hydropathy)
    if abs(net_charge) > max_FCR:
        error_message = f'\n\nThe input NCPR of {net_charge} is not within the range of possible NCPR values to keep hydropathy the same and maintain disorder, which is {max_FCR}\n'
        raise GooseInputError(error_message)


    # get numbers for final objective residues
    if constant_fcr == True:
        if abs(net_charge) > original_FCR:
            raise GooseInputError('Cannot have constant FCR with a variant NCPR greater than the input sequence FCR.')
        final_objective = fraction_net_charge(len(sequence), original_FCR, net_charge)
    else:
        if abs(net_charge) > original_FCR:
            new_FCR = abs(net_charge)
            sequence = create_fcr_class_variant(sequence, new_FCR, ignore_kappa=ignore_kappa)
            final_objective = fraction_net_charge(len(sequence), new_FCR, net_charge)
        else:
            final_objective = fraction_net_charge(len(sequence), original_FCR, net_charge)

    # figure out how many negative and positively charged res needed for
    # new sequence
    needed_charged = needed_charged_residues(len(sequence), Protein(sequence).FCR, net_charge)
    num_negative_needed = needed_charged['negative']
    num_positive_needed = needed_charged['positive']

    # figure out how many residues to change
    current_negative = sequence.count('D') + sequence.count('E')
    current_positive = sequence.count('K') + sequence.count('R')

    # get original ncpr
    original_ncpr = Protein(sequence).NCPR

    # figure out what direction changing
    if net_charge == original_ncpr:
        return sequence
    elif net_charge > original_ncpr:
        change_residues = num_positive_needed - current_positive
        target_residues = ['D', 'E']
        swap_residues = ['K', 'R']
    else:
        change_residues = num_negative_needed - current_negative
        target_residues = ['K', 'R']
        swap_residues = ['D', 'E']

    # get target residues locs
    target_residue_locs = []
    for aa in target_residues:
        if aa in sequence:
            target_residue_locs.extend(identify_residue_positions(sequence, aa))
    
    # shuffle the locs
    random.shuffle(target_residue_locs)

    # get final targets
    final_target_locs = []
    for i in range(0, change_residues):
        final_target_locs.append(target_residue_locs.pop())

    # build initial seq
    built_sequence=''
    for i in range(0, len(sequence)):
        if i in final_target_locs:
            built_sequence += random.choice(swap_residues)
        else:
            built_sequence += sequence[i]

    # now fix kappa
    if ignore_kappa==False:
        if pr(built_sequence).kappa != -1:
            fixed_kappa_seq = create_kappa_variant(built_sequence, original_kappa)
        else:
            fixed_kappa_seq = built_sequence
    else:
        fixed_kappa_seq = built_sequence
    
    # now fix hyropathy
    fixed_hydro_seq = optimize_hydropathy_within_class(fixed_kappa_seq, original_hydropathy)

    # now return the seq
    return fixed_hydro_seq


def create_all_props_class_variant(sequence, hydropathy=None, fraction=None, net_charge=None, kappa=None):
    '''
    function to make a variant where you can change hydropathy, fraction charged
    residues, net charge, and kappa all at once while minimizing the changes to
    residues by class in the sequence. As you change the sequence more, you will
    further alter the sequence, even outside of the classes of amino acids. This
    function simply attempts to minimize those changes.
    '''
    original_kappa = pr(sequence).kappa

    if fraction != None:
        if net_charge!= None:
            sequence = create_fcr_class_variant(sequence, fraction, constant_ncpr=False, use_closest = True, ignore_kappa=True)
        else:
            sequence = create_fcr_class_variant(sequence, fraction, constant_ncpr=False, use_closest = True, ignore_kappa=True)
    if net_charge != None:
        if fraction == None:
            sequence = create_ncpr_class_variant(sequence, net_charge, constant_fcr=False, ignore_kappa=True)
        else:
            sequence = create_ncpr_class_variant(sequence, net_charge, constant_fcr=False, ignore_kappa=True)
    if hydropathy != None:
        sequence = create_hydropathy_class_variant(sequence, hydro=hydropathy)
    if kappa == None:
        sequence = create_kappa_variant(sequence, original_kappa)
    else:
        sequence = create_kappa_variant(sequence, kappa)
    return sequence


