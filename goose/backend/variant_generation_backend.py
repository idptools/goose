'''
backend functionality for variant generation
'''
import random
from random import randint

from goose.goose_exceptions import GooseInputError
from goose.backend.protein import Protein
from goose.backend.sequence_generation_backend import identify_residue_positions, get_optimal_residue, optimal_residue_key, random_amino_acid, create_seq_by_props, fast_predict_disorder
from goose.backend import lists
from goose.backend.amino_acids import AminoAcid
from goose.backend import parameters


# might have to remove
from sparrow import Protein as pr




'''
=/=/=/=/=/=/=/=/=/=/=
Backend functionality
=/=/=/=/=/=/=/=/=/=/=
Below is the backend functionality necessary for the variant
generators.
'''

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


    '''
    # define which residues in each class
    aromatics = ['F', 'W', 'Y']
    polar = ['Q', 'N', 'S', 'T']
    hydrophobic = ['I', 'V', 'L', 'A', 'M']
    positive = ['K', 'R']
    negative = ['D', 'E']


    # first count the number of each class in the sequence
    num_aromatics = 0
    num_polar = 0
    num_hydro = 0
    num_positive = 0
    num_negative = 0
    num_G = 0
    num_P = 0
    num_H = 0
    num_C = 0

    # now iterate through the sequence adding to each class as necessary
    for amino_acid in sequence:
        if amino_acid == 'G':
            num_G += 1
        elif amino_acid == 'P':
            num_P += 1
        elif amino_acid == 'H':
            num_H += 1
        elif amino_acid == 'C':
            num_C += 1
        elif amino_acid in negative:
            num_negative += 1
        elif amino_acid in positive:
            num_positive += 1
        elif amino_acid in hydrophobic:
            num_hydro += 1
        elif amino_acid in polar:
            num_polar += 1
        elif amino_acid in aromatics:
            num_aromatics += 1
        else:
            raise GooseInputError('Invalid amino acid detected!')

    return {'aromatic': num_aromatics, 'polar': num_polar, 'hydrophobic': num_hydro, 'positive': num_positive,
    'negative': num_negative, 'C': num_C, 'H': num_H, 'P': num_P, 'G': num_G}


def get_optimal_res_within_class(residue_class, sequence, additional_exclusion=[], return_all=False):
    '''
    function to get the optimal residues within a class 
    for a sequence where the residue returned will
    be the residue to add to the end of the sequence
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
                four_residues += random_amino_acid(lists.disordered_list)
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
    return[round(min_possible/len(sequence), 6), round(max_possible/len(sequence), 6)]



def residue_optimize_hydropathy_within_class(sequence, objective_hydropathy, target_residue_index):
    '''
    function to get hydropathy value closer to objective hydropathy
    while keeping the residues within the sequence within the same
    class

    This function only changes one residue at a time, which corresponds to 'target_residue_index'!
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
    starting_hydro = Protein.calc_mean_hydro(sequence)

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
        if cur_aa_num != target_residue_index:
            best_residue = cur_aa
        else:
            if cur_aa in dont_change:
                best_residue = cur_aa
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
                            if abs(ideal_residue_difference - cur_poss_hydro) >= best_difference:
                                best_residue = possible_residues
                                best_hydropathy = cur_poss_hydro
                                best_difference = abs(ideal_residue_difference - cur_poss_hydro)
                else:
                    amino_acid_class = AminoAcid.return_AA_class(cur_aa)
                    best_hydropathy = cur_AA_hydro
                    best_residue = cur_aa
                    best_difference = abs(ideal_residue_difference - cur_AA_hydro)
                    for possible_residues in aa_class_dict[amino_acid_class]:
                        cur_poss_hydro = AminoAcid.hydro(possible_residues)
                        if cur_poss_hydro > cur_AA_hydro:
                            if abs(ideal_residue_difference - cur_poss_hydro) <= best_difference:
                                best_residue = possible_residues
                                best_hydropathy = cur_poss_hydro
                                best_difference = abs(ideal_residue_difference - cur_poss_hydro)                       
        build_sequence += best_residue
    # return the sequence
    return build_sequence  


def optimize_hydropathy_within_class(sequence, objective_hydropathy, allowed_hydro_error = parameters.HYDRO_ERROR):
    '''
    This function will optimie the hydropathy of a sequence such that it is
    within the allowed_hydro_error value.
    Will not change the classes of residues in the sequence or their position.
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
            cur_hydro = Protein.calc_mean_hydro(new_sequence)
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
    '''

    negatives = []
    positives = []

    for aa_ind in range(0, len(sequence)):
        aa = sequence[aa_ind]
        if aa == 'K':
            positives.append(aa_ind)
        elif aa == 'R':
            positives.append(aa_ind)
        elif aa == 'D':
            negatives.append(aa_ind)
        elif aa == 'E':
            negatives.append(aa_ind)

    return {'negative':negatives, 'positive':positives}


'''
=/=/=/=/=/=/=/=/=/=/=/=/=/=/
VARIANT SEQUENCE GENERATORS
=/=/=/=/=/=/=/=/=/=/=/=/=/=/
'''

'''
Below are the actual sequence generators.
Once again these just generate the sequence.
They do not check for disorder.
'''

def gen_new_var_constant_class_nums(sequence):
    '''
    function that takes an input sequence and returns 
    a sequence with the same properties generated from 
    residues that match the number for the number in that
    specific class in the input sequnece
    '''

    # define which residues in each class
    aromatics = ['F', 'W', 'Y']
    polar = ['Q', 'N', 'S', 'T']
    hydrophobics = ['I', 'V', 'L', 'A', 'M']
    positive = ['K', 'R']
    negative = ['D', 'E']

    # get starting hydropathy of sequence
    starting_hydropathy = Protein.calc_mean_hydro(sequence)

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

    return build_sequence



def gen_constant_class_variant(sequence):
    '''
    function to generate a variant with the same properties as the 
    input variant as well as the same order of amino acids as
    far as class and the same number in each class
    '''
    # build sequence to optimize disorder
    build_sequence = ''
    for amino_acid in sequence:
        cur_class = AminoAcid.return_AA_class(amino_acid)
        potential_residues = get_optimal_res_within_class(cur_class, build_sequence, additional_exclusion=[], return_all = True)
        if len(potential_residues) == 1:
            build_sequence += potential_residues[0]
        else:
            if amino_acid in potential_residues:
                potential_residues.remove(amino_acid)
            if len(potential_residues) == 1:
                build_sequence += potential_residues[0]
            else:
                build_sequence += potential_residues[randint(0, len(potential_residues)-1)]
    # now correct hydropathy
    starting_disorder = Protein.calc_mean_hydro(sequence)
    final_seq = optimize_hydropathy_within_class(build_sequence, starting_disorder)
    return final_seq




def gen_new_variant(sequence):
    '''
    function to generate a variant that is completely different
    in sequence to the input but has all the same overall parameters.
    Does not account for specific classes of residues.

    Getting the SCD right is a little hacky, but this works fast and does a good
    job on the disorder.
    '''
    input_length = len(sequence)
    input_FCR = Protein.calc_FCR(sequence)
    input_NCPR = Protein.calc_NCPR(sequence)
    input_hydro = Protein.calc_mean_hydro(sequence)
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



def hydropathy_class_variant(sequence, hydro, allowed_hydro_error = parameters.HYDRO_ERROR):
    '''
    function to take in a sequence and make a variant that adjusts the
    hydropathy while keeping the position and nuimber of amino acids the
    same by class of amino acid

    MAKES THE SAME SEQUENCE EVERY TIME! MIGHT HAVE TO MODIFY THIS!

    '''
    # change the hydropathy using the optimization function
    final_seq = optimize_hydropathy_within_class(sequence, hydro, allowed_hydro_error)
    return final_seq
    


def gen_constant_residue_variant(sequence, constant_residues = []):
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


    '''
    # make sure sequence is uppercase (this would otherwise wreck this function)
    sequence = sequence.upper()

    # if no constant residues specified, just generate a sequence and return it
    if constant_residues == []:
        return gen_new_variant(sequence)

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
    if 'K' in constant_residues:
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
    input_FCR = Protein.calc_FCR(seq_variant)
    input_NCPR = Protein.calc_NCPR(seq_variant)
    input_hydro = Protein.calc_mean_hydro(seq_variant)
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
            cur_res = seq_variant[seq_variant_pos]
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



def gen_shuffle_variant(sequence, shuffle_regions=[], use_index=False):
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
            raise exception('Error in finding middle residue')
    raise exception('Error in finding middle residue, made to end of function')




def increase_charge_asymmetry_once(sequence, exclude = []):
    
    """
    
    Function to increase charge asymmetry of an input sequence.
    Uses some randomness during generation in case the generated 
    does not reach the cutoff value (then can retry a different
    sequence)

    Parameters
    -------------

    sequence : String
        The amino acid sequence as a string.

    exclude : List
        A list of residues to exclude moving when altering asymmetry


    Returns
    ---------
        
    sequence : String
        Returns an amino acid sequence with an increased charge asymmetry

    """

    # identify positions of negative and positive residues in a sequence
    if 'E' not in exclude:
        E_coordinates = identify_residue_positions(sequence, 'E')
    else:
        E_coordinates = []

    if 'D' not in exclude:
        D_coordinates = identify_residue_positions(sequence, 'D')
    else:
        D_coordinates = []

    if 'K' not in exclude:
        K_coordinates = identify_residue_positions(sequence, 'K')
    else:
        K_coordinates = []

    if 'R' not in exclude:
        R_coordinates = identify_residue_positions(sequence, 'R')
    else:
        R_coordinates = []

    negative_coordinates = (E_coordinates + D_coordinates)
    positive_coordinates = (K_coordinates + R_coordinates)

    # figure out what residues can be targeted and where to put them based on charge distribution
    if len(positive_coordinates) > 0:
        positive_weight = sum(positive_coordinates)/len(positive_coordinates)
    else:
        positive_weight = 0
    
    if len(negative_coordinates) > 0:
        negative_weight = sum(negative_coordinates)/len(negative_coordinates)
    else:
        negative_weight = 0

    N=0


    # if no charged residues just return the sequence
    if positive_weight == 0 and negative_weight == 0:
        return sequence

    # if no positive charges but have negative charges
    elif positive_weight == 0 and negative_weight > 0:
        
        # make empty string to hold final sequence
        nonmoved_res = ""
        target_negative_coordinates = []
        target_negative_residue = max(negative_coordinates)

        # identify residues to pull from final sequence
        for residue in range(0, len(sequence)):
            if residue == target_negative_residue:
                moved_negative_residue = sequence[residue]
            else:
                nonmoved_res += sequence[residue]

        for i in range(0, int(0.25*len(nonmoved_res))):
            if i < target_negative_residue:
                target_negative_coordinates.append(i)  
        # if no negative_charges but have positive charges
        #elif positive_weight > 0 and negative weight == 0:

        if len(target_negative_coordinates) > 1:
            chosen_negative_position = target_negative_coordinates[random.randint(0, len(target_negative_coordinates)-1)]
        else:
            chosen_negative_position = 0        

        # build the final sequence
        final_sequence = ""

        for i in range(0, len(nonmoved_res)):
            if i == chosen_negative_position:
                final_sequence += moved_negative_residue
                final_sequence += nonmoved_res[i]
            else:
                final_sequence += nonmoved_res[i]


    # if no negative charges but have positive charges
    elif positive_weight > 0 and negative_weight == 0:
        # make empty string to hold final sequence
        nonmoved_res = ""
        target_positive_coordinates = []
        target_positive_residue = max(positive_coordinates)

        # identify residues to pull from final sequence
        for residue in range(0, len(sequence)):
            if residue == target_positive_residue:
                moved_positive_residue = sequence[residue]
            else:
                nonmoved_res += sequence[residue]

        for i in range(0, int(0.25*len(nonmoved_res))):
            if i < target_positive_residue:
                target_positive_coordinates.append(i)  
        # if no negative_charges but have positive charges
        #elif positive_weight > 0 and negative weight == 0:

        if len(target_positive_coordinates) > 1:
            chosen_positive_position = target_positive_coordinates[random.randint(0, len(target_positive_coordinates)-1)]
        else:
            chosen_positive_position = 0        

        # build the final sequence
        final_sequence = ""

        for i in range(0, len(nonmoved_res)):
            if i == chosen_positive_position:
                final_sequence += moved_positive_residue
                final_sequence += nonmoved_res[i]
            else:
                final_sequence += nonmoved_res[i]


    else:
        # randomly decide whether to target a positive or negative residue
        chosen_val = random.randint(2, 10)
        if chosen_val % 2 == 0:
            change_positive = True
            change_negative = False
        else:
            change_positive = False
            change_negative = True

        # make empty string to hold final sequence
        nonmoved_res = ""
        target_positive_coordinates = []
        target_negative_coordinates = []
        # figure out where to put residues
        if positive_weight > negative_weight:
            # decide if targeting postiive or negative
            if change_positive == True:
                target_positive_residue = min(positive_coordinates)
            else:
                target_negative_residue = max(negative_coordinates)

            # identify residues to pull from final sequence
            for residue in range(0, len(sequence)):
                if change_positive == True:
                    if residue == target_positive_residue:
                        moved_positive_residue = sequence[residue]
                    else:
                        nonmoved_res += sequence[residue]

                else:
                    if residue == target_negative_residue:
                        moved_negative_residue = sequence[residue]
                    else:
                        nonmoved_res += sequence[residue]

            if change_positive == True:
                for i in range(int(0.75*len(nonmoved_res)), len(nonmoved_res)):
                    if i > target_positive_residue:
                        target_positive_coordinates.append(i)
            else:
                for i in range(0, int(0.25*len(nonmoved_res))):
                    if i < target_negative_residue:
                        target_negative_coordinates.append(i)  

        else:
            if change_positive == True:
                target_positive_residue = max(positive_coordinates)
            else:
                target_negative_residue = min(negative_coordinates)

            # identify residues to pull from final sequence
            for residue in range(0, len(sequence)):
                if change_positive == True:
                    if residue == target_positive_residue:
                        moved_positive_residue = sequence[residue]
                    else: 
                        nonmoved_res += sequence[residue]
                
                else:
                    if residue == target_negative_residue:
                        moved_negative_residue = sequence[residue]
                    else:
                        nonmoved_res += sequence[residue]

            if change_positive == True:
                for i in range(0, int(0.75*len(nonmoved_res))):
                    if i < target_positive_residue:
                        target_positive_coordinates.append(i)
            else:
                for i in range(int(0.25*len(nonmoved_res)), len(nonmoved_res)):
                    if i > target_negative_residue:
                        target_negative_coordinates.append(i)            

        # choose a random target residue
        if change_positive == True:
            if len(target_positive_coordinates) > 1:
                chosen_positive_position = target_positive_coordinates[random.randint(0, len(target_positive_coordinates)-1)]
            else:
                if positive_weight > negative_weight:
                    chosen_positive_position = len(nonmoved_res)-1
                else:
                    chosen_positive_position = 0
        else:
            if len(target_negative_coordinates) > 1:
                chosen_negative_position = target_negative_coordinates[random.randint(0, len(target_negative_coordinates)-1)]
            else:
                if negative_weight > positive_weight:
                    chosen_negative_position = len(nonmoved_res)-1
                else:
                    chosen_negative_position = 0        


        # build the final sequence
        final_sequence = ""

        for i in range(0, len(nonmoved_res)):
            if change_positive == True:
                if i == chosen_positive_position:
                    final_sequence += moved_positive_residue
                    final_sequence += nonmoved_res[i]
                else:
                    final_sequence += nonmoved_res[i]

            if change_negative == True:
                if i == chosen_negative_position:
                    final_sequence += moved_negative_residue
                    final_sequence += nonmoved_res[i]
                else:
                    final_sequence += nonmoved_res[i]

    return final_sequence




def decrease_charge_asymmetry_once(sequence):
    """
    
    Function to decrease charge asymmetry of an input sequence.
    Uses some randomness during generation in case the generated 
    does not reach the cutoff value (then can retry a different
    sequence)

    Parameters
    -------------

    sequence : String
        The amino acid sequence as a string.

    Returns
    ---------
        
    sequence : String
        Returns an amino acid sequence with an increased charge asymmetry

    """


    # set arbitrary lowest and highest ncpr areas
    lowest_NCPR = 100
    highest_NCPR = -100
    N=0
    bloblen=5

    # make lists to hold possible coords to select from randomly
    possible_lowest_NCPR_coords = []
    possible_highest_NCPR_coords = []

    for i in range(0, len(sequence)-bloblen):
        cur_blob = sequence[i: i+bloblen]
        cur_NCPR = Protein.calc_NCPR(cur_blob)
        if cur_NCPR >= highest_NCPR:
            if cur_NCPR > highest_NCPR:
                possible_highest_NCPR_coords = [[i, i+bloblen]]
            else:
                possible_highest_NCPR_coords.append([i, i+bloblen])
            highest_NCPR = cur_NCPR
        if cur_NCPR <= lowest_NCPR:
            if cur_NCPR < lowest_NCPR:
                possible_lowest_NCPR_coords = [[i, i+bloblen]]
            else:
                possible_lowest_NCPR_coords.append([i, i+bloblen])
            lowest_NCPR = cur_NCPR


    # select random lowest and highest NCPR intervals
    if len(possible_lowest_NCPR_coords) > 1:
        lowest_NCPR_blob = possible_lowest_NCPR_coords[random.randint(0, len(possible_lowest_NCPR_coords)-1)]
    else:
        lowest_NCPR_blob = possible_lowest_NCPR_coords[0]

    if len(possible_highest_NCPR_coords) > 1:
        highest_NCPR_blob = possible_highest_NCPR_coords[random.randint(0, len(possible_highest_NCPR_coords)-1)]
    else:
        highest_NCPR_blob = possible_highest_NCPR_coords[0]


    # figure out possible targets
    possible_positive_targets = []
    possible_negative_targets = []

    for aa in range(highest_NCPR_blob[0], highest_NCPR_blob[1]):
        if sequence[aa] == 'K' or sequence[aa] == 'R':
            possible_positive_targets.append(aa)

    for aa in range(lowest_NCPR_blob[0], lowest_NCPR_blob[1]):
        if sequence[aa] == 'D' or sequence[aa] == 'E':
            possible_negative_targets.append(aa)    


    # select random residue to change
    if len(possible_negative_targets) > 1:
        selected_negative_residue = possible_negative_targets[random.randint(0, len(possible_negative_targets)-1)]
    else:
        if possible_negative_targets != []:
            selected_negative_residue = possible_negative_targets[0]
        else:
            selected_negative_residue = ""

    if len(possible_positive_targets) > 1:
        selected_positive_residue = possible_positive_targets[random.randint(0, len(possible_positive_targets)-1)]
    else:
        if possible_positive_targets != []:
            selected_positive_residue = possible_positive_targets[0]
        else:
            selected_positive_residue = ""

    # make sure that there is a residue chosen to swap no matter what!

    if selected_negative_residue == "":
        selected_negative_residue = random.randint(lowest_NCPR_blob[0], lowest_NCPR_blob[1])

    if selected_positive_residue == "":
        selected_positive_residue = random.randint(highest_NCPR_blob[0], highest_NCPR_blob[1])

    #build the final sequence
    final_sequence = ""

    for i in range(0, len(sequence)):
        if i == selected_negative_residue:
            final_sequence += sequence[selected_positive_residue]
        elif i == selected_positive_residue:
            final_sequence += sequence[selected_negative_residue]
        else:
            final_sequence += sequence[i]

    return final_sequence





def decrease_kappa_below_value(sequence, max_kappa, attempts=None):
    '''
    function to decrease the kappa value below
    a target value.

    Parameters
    ----------
    sequence : str
        the amino acid sequence as a string

    max_kappa : float
        the max allowed kappa value as a float

    Returns
    -------
    final_sequence : str
        returns the sequence with a reduced kappa value as a string
    '''
    original_kappa = pr(sequence).kappa
    if original_kappa < max_kappa:
        return sequence

    else:
        # if user doesn't set number of attempts, set number
        if attempts == None:
            # set attempts to 20 * length of the seqeunece
            attempts = 20*len(sequence)
        new_sequence = decrease_charge_asymmetry_once(sequence)
        curkappa = pr(new_sequence).kappa
        # reduce kappa number ot times there are 
        # residues in the sequence
        for i in range(0, attempts):

            if curkappa < max_kappa:
                return new_sequence
            else:
                new_sequence = decrease_charge_asymmetry_once(new_sequence)
                curkappa = pr(new_sequence).kappa
    raise Exception('Unable to return sequence with kappa value below specified value')



def gen_kappa_variant(sequence, kappa, allowed_kappa_error = parameters.MAXIMUM_KAPPA_ERROR, attempts=20000):
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
    if attempts==None:
        attempts = len(sequence)*50


    # first decrease kappa 
    kappa_below = kappa - 0.01
    if kappa_below < 0.005:
        kappa_below = 0.005

    decreased_kappa_seq = decrease_kappa_below_value(sequence, kappa_below)
    starting_kappa = pr(decreased_kappa_seq).kappa
    if abs(starting_kappa-kappa) < allowed_kappa_error:
        return decreased_kappa_seq
    else:
        new_sequence = increase_charge_asymmetry_once(decreased_kappa_seq)
        curkappa = pr(new_sequence).kappa
        for i in range(0, attempts):
            if abs(curkappa - kappa) < allowed_kappa_error:
                return new_sequence
            else:
                new_sequence = increase_charge_asymmetry_once(new_sequence)
                curkappa = pr(new_sequence).kappa
                # if you go too high in value reset and try again.
                # should be stochastic enough to work....
                if curkappa > kappa + (allowed_kappa_error*3):
                    new_sequence = decreased_kappa_seq
                    curkappa = starting_kappa
    raise Exception('unable to generate sequence with desired kappa')


test = 'IKLANATKKVGTKPAESDKKEEEKSAETKE'
print(Protein.calc_all_properties(test))
print(Protein.calc_all_properties((test)))

'''
need to do the following:

*when have multiple monitors, should definitely
move common code between variant and sequence generation
to a common backend location to import....


1. Minimal variant
3. Increase asymmetry (specific res or charge, polar, etc.)
4. increase charge asymmetry

then add all to docs and to the variant_generation.py stuff...

variants to add to main generation point
----------------------------------------

gen_new_var_constant_class_nums => returned sequence has same properties and same number of residues by class as the input sequence but the order is not the same. POSITIONS OF RESIDUES WITHIN EACH CLASS CAN MOVE RELATIVE TO ONE ANOTHER!

gen_constant_class_variant => returned sequence has same number AND POSITION of amino acis as input sequence BY CLASS. Properties will be same in returend sequence too.

gen_new_variant => returned sequence will have the same properies as the input sequence. That's it. Number of amino acids by class is allowed to change.

hydropathy_class_variant => will keep the position and number of all amino acids by class the same in the returned sequence. Allows adjustments to hydropathy while keeping everyhting else constant.

gen_constant_residue_variant => returns sequence that has specified residues held constant. The returned sequence will have the same properties as the input sequence

gen_minimal_variant -> change the properties of your input sequence while minimizing changes to the sequence

gen_shuffle_variant -> choose indiviudal regions of your sequence to shuffle

increase_asymmetry - > choose residue(s) or class of residues to increase asymmetry of in the sequence.

increase_charge_asymmetry -> increase kappa

***
also need to add in the charge distribution functionality
***



#def_shuffle_variant:

#def gen_minimal_variant:



seq_var = 'IKLANATKKVGTKPAESDKKEEEKSAETKEPTKEPTKVEE'
newseq = gen_new_variant(seq_var)

print(Protein.calc_all_properties(seq_var))
print(Protein.calc_all_properties(newseq))
'''



