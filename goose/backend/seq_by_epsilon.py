'''
Code for desigining sequences based on epsilon. 

Epsilon is a value that is a summation of pairwise interaction values between
amino acids. It is used to predict the interaction between two proteins.
This allowes us to make sequences with epsilon values similar to a starting
sequence as well as epsilon values that are simply predicted to be overall
attractive or repulsive.  

Because epsilon is composed of the summation of an attractive and repulsive
interaction vector between two sequences, we can make sequences that have
similar interaction vectors. This can be useful for making sequences that
interact with a specific protein in a similar way to another protein. 
'''

'''
TO DO: (if time)
    • Also need seq by chunked epsilon.
    • Add in 'make more attractive' or 'make more repulsive' optimizations. 
'''

from collections import Counter
import itertools
import random

import matplotlib.pyplot as plt
import numpy as np

from finches import epsilon_calculation
from finches.forcefields.mpipi import Mpipi_model
from goose.backend import lists
from goose import goose_exceptions


# -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
#        -=-=-=-=-=-= Code for epsilon related stuffs =-=-=-=-=-=-=
# -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=

def get_interaction_vectors(seq1, seq2, model='Mpipi_GGv1', approach='mean'):
    '''
    Function that returns the interaction vectors between two sequences. 
    The vector will return a list where the first list is seq1 vector
    and the second item in the list is seq2 vector. 

    Should be able to use this to map reigions in our protein that need modifications
    to comply with an objective interaction vector. 

    Parameters
    -----------
    seq1 : string
        the amino acid sequence for seq1 as a string

    seq2 : string
        the amino acid sequence for seq2 as a string

    model : string
        The specific model parameters we are using
        default = Mpipi_GGv1

    approach : string
        How the vectors should be calculated. Options are sum and mean.
        Default = sum. 
        Mean (should) give approximate values per amino acid in the sequence which
        in theory will identify interaction propensities that are closer to values
        in the pairwise interaction dict. 

    Returns
    -------
    list of np.arrays
        Returns a list of np.arrays. First item in the list is the
        interaction vector for seq1, the second is for seq2. 
    '''
    # temporary hold over unitl I implement other forcefields
    if model != 'Mpipi_GGv1':
        raise Exception('Other forcefields have not been implemented. Please set to Mpipi_GGv1.')
    
    # initialize forcefield parameters
    Mpipi_GGv1_model = Mpipi_model('Mpipi_GGv1')
    
    # make IMC Object
    IMC_object = epsilon_calculation.InteractionMatrixConstructor(Mpipi_GGv1_model)
    
    # use IMC_object to calculate heterotypic matrix
    interaction_matrix=IMC_object.calculate_pairwise_heterotypic_matrix(seq1,seq2)
    # now get the mean or sum values
    if approach=='mean':
        t1_vector=(np.mean(interaction_matrix, axis=1))
        t2_vector=(np.mean(interaction_matrix, axis=0))
    elif approach == 'sum':
        t1_vector=(np.sum(interaction_matrix, axis=1))
        t2_vector=(np.sum(interaction_matrix, axis=0))
    else:
        raise Exception('approach can only be set to mean or sum at the moment.')

    return [t1_vector, t2_vector]


def get_epsilon_vectors(seq1, seq2, model='Mpipi_GGv1'):
    '''
    Function that returns the interaction vectors between two sequences. 
    The vector will return a list where the first list is seq1 vector
    and the second item in the list is seq2 vector. 

    Should be able to use this to map reigions in our protein that need modifications
    to comply with an objective interaction vector. 

    Parameters
    -----------
    seq1 : string
        the amino acid sequence for seq1 as a string

    seq2 : string
        the amino acid sequence for seq2 as a string

    model : string
        The specific model parameters we are using
        default = Mpipi_GGv1

    Returns
    -------
    list of np.arrays
        Returns a list of np.arrays. First item in the list is the
        interaction vector for seq1, the second is for seq2. 
    '''
    # temporary hold over unitl I implement other forcefields
    if model != 'Mpipi_GGv1':
        raise Exception('Other forcefields have not been implemented. Please set to Mpipi_GGv1.')
    
    # initialize forcefield parameters
    Mpipi_GGv1_model = Mpipi_model('Mpipi_GGv1')
    
    # make IMC Object
    IMC_object = epsilon_calculation.InteractionMatrixConstructor(Mpipi_GGv1_model)
    
    # use IMC_object to calculate heterotypic matrix
    interaction_matrix=IMC_object.calculate_epsilon_vectors(seq1,seq2)

    return [interaction_matrix[0], interaction_matrix[1]]

def get_epsilon_value(seq1, seq2, model='Mpipi_GGv1'):
    '''
    Function to get the epsilon value using the Mpipi_GGv1 model. 

    Parameters
    -----------
    seq1 : string
        the amino acid sequence for seq1 as a string

    seq2 : string
        the amino acid sequence for seq2 as a string

    model : string
        The specific model parameters we are using
        default = Mpipi_GGv1

    Returns
    -------
    list of np.arrays
        Returns a list of np.arrays. First item in the list is the
        interaction vector for seq1, the second is for seq2. 
    '''
    # temporary hold over unitl I implement other forcefields
    if model != 'Mpipi_GGv1':
        raise Exception('Other forcefields have not been implemented. Please set to Mpipi_GGv1.')
    
    # initialize forcefield parameters
    Mpipi_GGv1_model = Mpipi_model('Mpipi_GGv1')
    
    # make IMC Object
    IMC_object = epsilon_calculation.InteractionMatrixConstructor(Mpipi_GGv1_model)
    
    # use IMC_object to calculate heterotypic matrix
    epslon_value=IMC_object.calculate_epsilon_value(seq1,seq2)

    return epslon_value


# -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
#    -=-=-=-=-=-=-=-=- Code for modifiying sequences =-=-=-=-=-=-=-=-=-
# -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=



def optimize_to_epsilon_value(starting_sequence, interacting_sequence, objective_epsilon,
    allowed_error=0.1, optimization_iterations=None, exclude_aas=[]):
    '''
    Function to optimze a sequence to have a specific epsilon value relative to another sequence. 

    Parameters
    -----------
    starting_sequence : str
        The sequence that you want to optimize to have a specific epsilon value relative to
        the interacting sequence.

    interacting_sequence : str
        The sequence that the starting sequence is interacting with.

    objective_epsilon : float
        The epsilon value that you want the starting sequence to have relative to the interacting

    allowed_error : float
        The allowed error in the epsilon value. Default is 0.1

    optimization_iterations : int
        The number of iterations to run the optimization. Default is length of the sequence

    exclude_aas=[]
        A list of amino acids to exclude from being added to the sequence during optimization.


    Returns
    -------
    str
        The optimized sequence that has the closest epsilon value to the objective epsilon value.
    '''
    # pairwise interaction values for each amino acid and U (RNA) against all
    # amino acids and U. Uses mPiPi-GGv1.
    epsilon_aas={'A': {'A': 0.1729, 'C': 0.1479, 'D': 0.1122, 'E': 0.0991, 'F': -0.1814, 'G': 0.1355, 'H': -0.1472, 'I': 0.1184, 'K': 0.1819, 'L': 0.1302, 'M': 0.1617, 'N': 0.0295, 'P': 0.0992, 'Q': 0.0159, 'R': -0.1, 'S': 0.1578, 'T': 0.1822, 'V': 0.1273, 'W': -0.3609, 'Y': -0.2128}, 'C': {'A': 0.1479, 'C': 0.1244, 'D': 0.0871, 'E': 0.0731, 'F': -0.2188, 'G': 0.1127, 'H': -0.1831, 'I': 0.1859, 'K': 0.1586, 'L': 0.1751, 'M': 0.1461, 'N': 0.0009, 'P': 0.0673, 'Q': -0.0136, 'R': -0.1343, 'S': 0.1352, 'T': 0.16, 'V': 0.1827, 'W': -0.4051, 'Y': -0.2514}, 'D': {'A': 0.1122, 'C': 0.0871, 'D': 2.574, 'E': 2.4772, 'F': 0.0495, 'G': 0.0783, 'H': -1.1931, 'I': 0.1469, 'K': -2.5449, 'L': 0.136, 'M': 0.1068, 'N': -0.0382, 'P': 0.0172, 'Q': -0.0538, 'R': -2.4973, 'S': 0.0989, 'T': 0.1225, 'V': 0.1447, 'W': -0.0261, 'Y': 0.0362}, 'E': {'A': 0.0991, 'C': 0.0731, 'D': 2.4772, 'E': 2.3835, 'F': 0.0398, 'G': 0.0643, 'H': -1.1415, 'I': 0.1344, 'K': -2.4518, 'L': 0.1232, 'M': 0.0931, 'N': -0.056, 'P': -0.0018, 'Q': -0.0721, 'R': -2.3995, 'S': 0.0854, 'T': 0.1096, 'V': 0.1323, 'W': -0.0378, 'Y': 0.0262}, 'F': {'A': -0.1814, 'C': -0.2188, 'D': 0.0495, 'E': 0.0398, 'F': -0.6114, 'G': -0.2051, 'H': -0.566, 'I': -0.1716, 'K': 0.0014, 'L': -0.1832, 'M': -0.2144, 'N': -0.3574, 'P': -0.393, 'Q': -0.3819, 'R': -0.2339, 'S': -0.199, 'T': -0.1842, 'V': -0.1662, 'W': -0.8217, 'Y': -0.6489}, 'G': {'A': 0.1355, 'C': 0.1127, 'D': 0.0783, 'E': 0.0643, 'F': -0.2051, 'G': -0.0747, 'H': -0.1708, 'I': 0.1667, 'K': 0.1411, 'L': 0.1568, 'M': 0.1302, 'N': -0.0007, 'P': 0.0525, 'Q': -0.0154, 'R': -0.1286, 'S': -0.0354, 'T': 0.1446, 'V': 0.1647, 'W': -0.3788, 'Y': -0.2356}, 'H': {'A': -0.1472, 'C': -0.1831, 'D': -1.1931, 'E': -1.1415, 'F': -0.566, 'G': -0.1708, 'H': 0.8857, 'I': -0.1362, 'K': 1.3994, 'L': -0.1475, 'M': -0.178, 'N': -0.3183, 'P': -0.3451, 'Q': -0.3419, 'R': 1.289, 'S': -0.1641, 'T': -0.1492, 'V': -0.1313, 'W': -0.7715, 'Y': -0.6027}, 'I': {'A': 0.1184, 'C': 0.1859, 'D': 0.1469, 'E': 0.1344, 'F': -0.1716, 'G': 0.1667, 'H': -0.1362, 'I': 0.0416, 'K': 0.2281, 'L': 0.0644, 'M': 0.1023, 'N': 0.056, 'P': 0.149, 'Q': 0.0432, 'R': -0.0806, 'S': 0.1954, 'T': 0.2247, 'V': 0.0536, 'W': -0.3651, 'Y': -0.2054}, 'K': {'A': 0.1819, 'C': 0.1586, 'D': -2.5449, 'E': -2.4518, 'F': 0.0014, 'G': 0.1411, 'H': 1.3994, 'I': 0.2281, 'K': 2.2789, 'L': 0.2165, 'M': 0.1853, 'N': 0.0263, 'P': 0.1123, 'Q': 0.0126, 'R': 2.0644, 'S': 0.1688, 'T': 0.1977, 'V': 0.2233, 'W': 0.017, 'Y': 0.0275}, 'L': {'A': 0.1302, 'C': 0.1751, 'D': 0.136, 'E': 0.1232, 'F': -0.1832, 'G': 0.1568, 'H': -0.1475, 'I': 0.0644, 'K': 0.2165, 'L': 0.0794, 'M': 0.1153, 'N': 0.045, 'P': 0.1346, 'Q': 0.0319, 'R': -0.0924, 'S': 0.1849, 'T': 0.2138, 'V': 0.0761, 'W': -0.3771, 'Y': -0.2171}, 'M': {'A': 0.1617, 'C': 0.1461, 'D': 0.1068, 'E': 0.0931, 'F': -0.2144, 'G': 0.1302, 'H': -0.178, 'I': 0.1023, 'K': 0.1853, 'L': 0.1153, 'M': 0.1501, 'N': 0.0155, 'P': 0.0957, 'Q': 0.0016, 'R': -0.1241, 'S': 0.1566, 'T': 0.1844, 'V': 0.1118, 'W': -0.4093, 'Y': -0.2485}, 'N': {'A': 0.0295, 'C': 0.0009, 'D': -0.0382, 'E': -0.056, 'F': -0.3574, 'G': -0.0007, 'H': -0.3183, 'I': 0.056, 'K': 0.0263, 'L': 0.045, 'M': 0.0155, 'N': -0.127, 'P': -0.0983, 'Q': -0.1452, 'R': -0.2736, 'S': 0.0151, 'T': 0.0356, 'V': 0.0562, 'W': -0.5511, 'Y': -0.3917}, 'P': {'A': 0.0992, 'C': 0.0673, 'D': 0.0172, 'E': -0.0018, 'F': -0.393, 'G': 0.0525, 'H': -0.3451, 'I': 0.149, 'K': 0.1123, 'L': 0.1346, 'M': 0.0957, 'N': -0.0983, 'P': 0.0539, 'Q': -0.1179, 'R': -0.2801, 'S': 0.0821, 'T': 0.1148, 'V': 0.1451, 'W': -0.6428, 'Y': -0.4368}, 'Q': {'A': 0.0159, 'C': -0.0136, 'D': -0.0538, 'E': -0.0721, 'F': -0.3819, 'G': -0.0154, 'H': -0.3419, 'I': 0.0432, 'K': 0.0126, 'L': 0.0319, 'M': 0.0016, 'N': -0.1452, 'P': -0.1179, 'Q': -0.1638, 'R': -0.2956, 'S': 0.001, 'T': 0.0222, 'V': 0.0435, 'W': -0.5805, 'Y': -0.417}, 'R': {'A': -0.1, 'C': -0.1343, 'D': -2.4973, 'E': -2.3995, 'F': -0.2339, 'G': -0.1286, 'H': 1.289, 'I': -0.0806, 'K': 2.0644, 'L': -0.0924, 'M': -0.1241, 'N': -0.2736, 'P': -0.2801, 'Q': -0.2956, 'R': 2.08, 'S': -0.1167, 'T': -0.098, 'V': -0.0777, 'W': -0.3863, 'Y': -0.2948}, 'S': {'A': 0.1578, 'C': 0.1352, 'D': 0.0989, 'E': 0.0854, 'F': -0.199, 'G': -0.0354, 'H': -0.1641, 'I': 0.1954, 'K': 0.1688, 'L': 0.1849, 'M': 0.1566, 'N': 0.0151, 'P': 0.0821, 'Q': 0.001, 'R': -0.1167, 'S': 0.1455, 'T': 0.1698, 'V': 0.1921, 'W': -0.3809, 'Y': -0.2308}, 'T': {'A': 0.1822, 'C': 0.16, 'D': 0.1225, 'E': 0.1096, 'F': -0.1842, 'G': 0.1446, 'H': -0.1492, 'I': 0.2247, 'K': 0.1977, 'L': 0.2138, 'M': 0.1844, 'N': 0.0356, 'P': 0.1148, 'Q': 0.0222, 'R': -0.098, 'S': 0.1698, 'T': 0.1965, 'V': 0.2204, 'W': -0.371, 'Y': -0.2168}, 'V': {'A': 0.1273, 'C': 0.1827, 'D': 0.1447, 'E': 0.1323, 'F': -0.1662, 'G': 0.1647, 'H': -0.1313, 'I': 0.0536, 'K': 0.2233, 'L': 0.0761, 'M': 0.1118, 'N': 0.0562, 'P': 0.1451, 'Q': 0.0435, 'R': -0.0777, 'S': 0.1921, 'T': 0.2204, 'V': 0.0726, 'W': -0.3555, 'Y': -0.1992}, 'W': {'A': -0.3609, 'C': -0.4051, 'D': -0.0261, 'E': -0.0378, 'F': -0.8217, 'G': -0.3788, 'H': -0.7715, 'I': -0.3651, 'K': 0.017, 'L': -0.3771, 'M': -0.4093, 'N': -0.5511, 'P': -0.6428, 'Q': -0.5805, 'R': -0.3863, 'S': -0.3809, 'T': -0.371, 'V': -0.3555, 'W': -1.0436, 'Y': -0.8617}, 'Y': {'A': -0.2128, 'C': -0.2514, 'D': 0.0362, 'E': 0.0262, 'F': -0.6489, 'G': -0.2356, 'H': -0.6027, 'I': -0.2054, 'K': 0.0275, 'L': -0.2171, 'M': -0.2485, 'N': -0.3917, 'P': -0.4368, 'Q': -0.417, 'R': -0.2948, 'S': -0.2308, 'T': -0.2168, 'V': -0.1992, 'W': -0.8617, 'Y': -0.687}}    
    
    # list of all amino acids
    all_aas=['A', 'C', 'D', 'E', 'F', 'G', 'H', 'I', 'K', 'L', 'M', 'N', 'P', 'Q', 'R', 'S', 'T', 'V', 'W', 'Y']
        
    # make list of exclude aas non redundant
    exclude_aas=list(set(exclude_aas))

    # make sure we aren't removing all of the amino acids
    if len(exclude_aas)==20:
        raise goose_exceptions.GooseInputError('Cannot exclude all amino acids.')

    # get rid of amino acids as needed
    for aa in exclude_aas:
        all_aas.remove(aa)

    # if optimization_iterations=-None, set to lenght of sequence
    if optimization_iterations is None:
        optimization_iterations=len(starting_sequence)

    # calculate the epsilon value of the starting sequence
    eps_vectors=get_epsilon_vectors(starting_sequence, interacting_sequence)
    starting_epsilon=np.sum(eps_vectors)
    # calculate the difference between the starting epsilon and the objective epsilon
    epsilon_diff=objective_epsilon-starting_epsilon
    # set the epsilon difference to be the current difference
    current_epsilon_diff=epsilon_diff
    # set the current sequence to be the starting sequence
    current_sequence=starting_sequence

    # now iterate over sequence to get towards the objective value. 
    for i in range(0, optimization_iterations):
        # if i > 1, recalculate the eps vectors. 
        if i > 1:
            eps_vectors=get_epsilon_vectors(current_sequence, interacting_sequence)
        # choose a random amino acid to change. 
        aa_to_change=random.randint(0, len(starting_sequence)-1)
        
        # dict to hold possible aas
        possible_aas={}

        # iterate through all aas
        for aa in all_aas:
            curv=0
            for a in interacting_sequence:
                curv=curv+epsilon_aas[aa][a]
            possible_aas[aa]=curv/len(interacting_sequence)
        # sort the possible aas by the difference between the current epsilon value and the objective epsilon value
        possible_aas=sorted(possible_aas, key=lambda x:abs(possible_aas[x]-current_epsilon_diff))
        # get the new amino acid
        new_aa=possible_aas[0]
        # change the amino acid
        current_sequence=current_sequence[:aa_to_change]+new_aa+current_sequence[aa_to_change+1:]
        # recalculate the epsilon value
        eps_vectors=get_epsilon_vectors(current_sequence, interacting_sequence)
        current_epsilon=np.sum(eps_vectors)
        # calculate the new difference
        current_epsilon_diff=objective_epsilon-current_epsilon
        # if the difference is less than the allowed error, break
        if abs(current_epsilon_diff) < allowed_error:
            return current_sequence
    raise goose_exceptions.GooseFail('Failed to optimize to epsilon value.')



def return_constrained_aa_list(input_sequence, exclude_aas=[]):
    '''
    returns a constrained aa list depending on if we get too many of
    any specific amino acid. 
    
    Parameters
    -----------
    input_sequence : str
        The sequence that we are checking for amino acid counts

    exclude_aas : list
        A list of amino acids to exclude from the generated sequence.

    Returns
    -------
    list
        A list of amino acids that are allowed to be used in the sequence

    '''
    # list all aas
    all_aas=['A', 'C', 'D', 'E', 'F', 'G', 'H', 'I', 'K', 'L', 'M', 'N', 'P', 'Q', 'R', 'S', 'T', 'V', 'W', 'Y']
    
    # lists we will watch out for 
    aro=['W', 'F', 'Y']
    ali=['I', 'L', 'V', 'A','M']

    # make list of exclude aas non redundant
    exclude_aas=list(set(exclude_aas))

    # make sure not all aas are excluded
    if len(exclude_aas)==20:
        raise goose_exceptions.GooseInputError('Cannot exclude all amino acids.')

    # exclude some aa if needed
    for aa in exclude_aas:
        all_aas.remove(aa)

    # get the count of each amino acid
    aa_counts=Counter(input_sequence)
    
    total_aro = sum([aa_counts[aa] for aa in aro])
    total_ali = sum([aa_counts[aa] for aa in ali])
    
    # make sure we don't get a weird amount of C
    if 'C' in input_sequence:
        if aa_counts['C']/len(input_sequence) > 0.15:
            all_aas.remove('C')

    # make sure we don't get too many aro
    if total_aro/len(input_sequence) > 0.2:
        for aa in aro:
            # check in case we already excluded it
            if aa in all_aas:
                all_aas.remove(aa)

    # modulate maximum allowable aliphatics dependent on aromaitcs. 
    max_ali = 0.3 - total_aro/len(input_sequence)

    # make sure we don't get too many ali
    if total_ali/len(input_sequence) > max_ali:
        for aa in ali:
            if aa in all_aas:
                all_aas.remove(aa)

    # make sure all aas not equal to []
    if all_aas == []:
        # if list is empty, we will add anything the user didn't explicitly exclude. 
        # this is very unlikely, but we should priortize that over the fractions.
        for aa in ['A', 'C', 'D', 'E', 'F', 'G', 'H', 'I', 'K', 'L', 'M', 'N', 'P', 'Q', 'R', 'S', 'T', 'V', 'W', 'Y']:
            if aa not in exclude_aas:
                all_aas.append(aa)

    return all_aas

# -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
#    -=-=-=-=-=-=-=-=-=-=- Code for sequence creation =-=-=-=-=-=-=-=-=-
# -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=


def create_starting_sequence(length, exclude_aas=[]):
    '''
    Function to create a random amino acid sequence of a specific length.
    The randomness is contstrained at least to increase probability
    of making something disordered but otherwise generally doesn't make
    the same sequence twice.

    Parameters
    -----------
    length : int
        The length of the sequence you want to create

    exclude_aas : list
        A list of amino acids to exclude from the generated sequence.
        Default = []

    Returns
    -------
    str
        The random amino acid sequence of the specified length
    '''
    if length <= 0:
        raise goose_exceptions.GooseInputError('Length must be greater than 0.')

    # make list of exclude aas non redundant
    exclude_aas=list(set(exclude_aas))

    if exclude_aas != []:
        if len(exclude_aas)==20:
            raise goose_exceptions.GooseInputError('Cannot exclude all amino acids.')
        fl=[]
        for aa in lists.disordered_list:
            if aa not in exclude_aas:
                fl.append(aa)
    else:
        fl=lists.disordered_list

    return ''.join(np.random.choice(fl, length))


def create_seq_by_epsilon_vectors(sequence_of_interest, interacting_sequence=None, 
    exclude_aas=[]):
    '''
    Function to create a sequence with approximately similar 
    interaction vectors to a sequence of interest. This function
    has some tolerance in that it will randomly choose from the
    best 2 amino acids at each position, so the total epsilon value
    will differ. However, the patterna cross the sequence should be quite close. 

    Parameters
    -----------
    sequence_of_interest : str
        The sequence that you want your generated sequence to
        have an epsilon value relative to.

    interacting_sequence : str
        The seqence that the sequence of interest is interacting with.
        If None, then it will be the same as sequence of interest (
        so will make something with the same homotypic interactions)

    exclude_aas : list
        A list of amino acids to exclude from the generated sequence. 
        Default = []
    '''
    # pairwise interaction values for each amino acid and U (RNA) against all
    # amino acids and U. Uses mPiPi-GGv1.
    epsilon_aas={'A': {'A': 0.1729, 'C': 0.1479, 'D': 0.1122, 'E': 0.0991, 'F': -0.1814, 'G': 0.1355, 'H': -0.1472, 'I': 0.1184, 'K': 0.1819, 'L': 0.1302, 'M': 0.1617, 'N': 0.0295, 'P': 0.0992, 'Q': 0.0159, 'R': -0.1, 'S': 0.1578, 'T': 0.1822, 'V': 0.1273, 'W': -0.3609, 'Y': -0.2128}, 'C': {'A': 0.1479, 'C': 0.1244, 'D': 0.0871, 'E': 0.0731, 'F': -0.2188, 'G': 0.1127, 'H': -0.1831, 'I': 0.1859, 'K': 0.1586, 'L': 0.1751, 'M': 0.1461, 'N': 0.0009, 'P': 0.0673, 'Q': -0.0136, 'R': -0.1343, 'S': 0.1352, 'T': 0.16, 'V': 0.1827, 'W': -0.4051, 'Y': -0.2514}, 'D': {'A': 0.1122, 'C': 0.0871, 'D': 2.574, 'E': 2.4772, 'F': 0.0495, 'G': 0.0783, 'H': -1.1931, 'I': 0.1469, 'K': -2.5449, 'L': 0.136, 'M': 0.1068, 'N': -0.0382, 'P': 0.0172, 'Q': -0.0538, 'R': -2.4973, 'S': 0.0989, 'T': 0.1225, 'V': 0.1447, 'W': -0.0261, 'Y': 0.0362}, 'E': {'A': 0.0991, 'C': 0.0731, 'D': 2.4772, 'E': 2.3835, 'F': 0.0398, 'G': 0.0643, 'H': -1.1415, 'I': 0.1344, 'K': -2.4518, 'L': 0.1232, 'M': 0.0931, 'N': -0.056, 'P': -0.0018, 'Q': -0.0721, 'R': -2.3995, 'S': 0.0854, 'T': 0.1096, 'V': 0.1323, 'W': -0.0378, 'Y': 0.0262}, 'F': {'A': -0.1814, 'C': -0.2188, 'D': 0.0495, 'E': 0.0398, 'F': -0.6114, 'G': -0.2051, 'H': -0.566, 'I': -0.1716, 'K': 0.0014, 'L': -0.1832, 'M': -0.2144, 'N': -0.3574, 'P': -0.393, 'Q': -0.3819, 'R': -0.2339, 'S': -0.199, 'T': -0.1842, 'V': -0.1662, 'W': -0.8217, 'Y': -0.6489}, 'G': {'A': 0.1355, 'C': 0.1127, 'D': 0.0783, 'E': 0.0643, 'F': -0.2051, 'G': -0.0747, 'H': -0.1708, 'I': 0.1667, 'K': 0.1411, 'L': 0.1568, 'M': 0.1302, 'N': -0.0007, 'P': 0.0525, 'Q': -0.0154, 'R': -0.1286, 'S': -0.0354, 'T': 0.1446, 'V': 0.1647, 'W': -0.3788, 'Y': -0.2356}, 'H': {'A': -0.1472, 'C': -0.1831, 'D': -1.1931, 'E': -1.1415, 'F': -0.566, 'G': -0.1708, 'H': 0.8857, 'I': -0.1362, 'K': 1.3994, 'L': -0.1475, 'M': -0.178, 'N': -0.3183, 'P': -0.3451, 'Q': -0.3419, 'R': 1.289, 'S': -0.1641, 'T': -0.1492, 'V': -0.1313, 'W': -0.7715, 'Y': -0.6027}, 'I': {'A': 0.1184, 'C': 0.1859, 'D': 0.1469, 'E': 0.1344, 'F': -0.1716, 'G': 0.1667, 'H': -0.1362, 'I': 0.0416, 'K': 0.2281, 'L': 0.0644, 'M': 0.1023, 'N': 0.056, 'P': 0.149, 'Q': 0.0432, 'R': -0.0806, 'S': 0.1954, 'T': 0.2247, 'V': 0.0536, 'W': -0.3651, 'Y': -0.2054}, 'K': {'A': 0.1819, 'C': 0.1586, 'D': -2.5449, 'E': -2.4518, 'F': 0.0014, 'G': 0.1411, 'H': 1.3994, 'I': 0.2281, 'K': 2.2789, 'L': 0.2165, 'M': 0.1853, 'N': 0.0263, 'P': 0.1123, 'Q': 0.0126, 'R': 2.0644, 'S': 0.1688, 'T': 0.1977, 'V': 0.2233, 'W': 0.017, 'Y': 0.0275}, 'L': {'A': 0.1302, 'C': 0.1751, 'D': 0.136, 'E': 0.1232, 'F': -0.1832, 'G': 0.1568, 'H': -0.1475, 'I': 0.0644, 'K': 0.2165, 'L': 0.0794, 'M': 0.1153, 'N': 0.045, 'P': 0.1346, 'Q': 0.0319, 'R': -0.0924, 'S': 0.1849, 'T': 0.2138, 'V': 0.0761, 'W': -0.3771, 'Y': -0.2171}, 'M': {'A': 0.1617, 'C': 0.1461, 'D': 0.1068, 'E': 0.0931, 'F': -0.2144, 'G': 0.1302, 'H': -0.178, 'I': 0.1023, 'K': 0.1853, 'L': 0.1153, 'M': 0.1501, 'N': 0.0155, 'P': 0.0957, 'Q': 0.0016, 'R': -0.1241, 'S': 0.1566, 'T': 0.1844, 'V': 0.1118, 'W': -0.4093, 'Y': -0.2485}, 'N': {'A': 0.0295, 'C': 0.0009, 'D': -0.0382, 'E': -0.056, 'F': -0.3574, 'G': -0.0007, 'H': -0.3183, 'I': 0.056, 'K': 0.0263, 'L': 0.045, 'M': 0.0155, 'N': -0.127, 'P': -0.0983, 'Q': -0.1452, 'R': -0.2736, 'S': 0.0151, 'T': 0.0356, 'V': 0.0562, 'W': -0.5511, 'Y': -0.3917}, 'P': {'A': 0.0992, 'C': 0.0673, 'D': 0.0172, 'E': -0.0018, 'F': -0.393, 'G': 0.0525, 'H': -0.3451, 'I': 0.149, 'K': 0.1123, 'L': 0.1346, 'M': 0.0957, 'N': -0.0983, 'P': 0.0539, 'Q': -0.1179, 'R': -0.2801, 'S': 0.0821, 'T': 0.1148, 'V': 0.1451, 'W': -0.6428, 'Y': -0.4368}, 'Q': {'A': 0.0159, 'C': -0.0136, 'D': -0.0538, 'E': -0.0721, 'F': -0.3819, 'G': -0.0154, 'H': -0.3419, 'I': 0.0432, 'K': 0.0126, 'L': 0.0319, 'M': 0.0016, 'N': -0.1452, 'P': -0.1179, 'Q': -0.1638, 'R': -0.2956, 'S': 0.001, 'T': 0.0222, 'V': 0.0435, 'W': -0.5805, 'Y': -0.417}, 'R': {'A': -0.1, 'C': -0.1343, 'D': -2.4973, 'E': -2.3995, 'F': -0.2339, 'G': -0.1286, 'H': 1.289, 'I': -0.0806, 'K': 2.0644, 'L': -0.0924, 'M': -0.1241, 'N': -0.2736, 'P': -0.2801, 'Q': -0.2956, 'R': 2.08, 'S': -0.1167, 'T': -0.098, 'V': -0.0777, 'W': -0.3863, 'Y': -0.2948}, 'S': {'A': 0.1578, 'C': 0.1352, 'D': 0.0989, 'E': 0.0854, 'F': -0.199, 'G': -0.0354, 'H': -0.1641, 'I': 0.1954, 'K': 0.1688, 'L': 0.1849, 'M': 0.1566, 'N': 0.0151, 'P': 0.0821, 'Q': 0.001, 'R': -0.1167, 'S': 0.1455, 'T': 0.1698, 'V': 0.1921, 'W': -0.3809, 'Y': -0.2308}, 'T': {'A': 0.1822, 'C': 0.16, 'D': 0.1225, 'E': 0.1096, 'F': -0.1842, 'G': 0.1446, 'H': -0.1492, 'I': 0.2247, 'K': 0.1977, 'L': 0.2138, 'M': 0.1844, 'N': 0.0356, 'P': 0.1148, 'Q': 0.0222, 'R': -0.098, 'S': 0.1698, 'T': 0.1965, 'V': 0.2204, 'W': -0.371, 'Y': -0.2168}, 'V': {'A': 0.1273, 'C': 0.1827, 'D': 0.1447, 'E': 0.1323, 'F': -0.1662, 'G': 0.1647, 'H': -0.1313, 'I': 0.0536, 'K': 0.2233, 'L': 0.0761, 'M': 0.1118, 'N': 0.0562, 'P': 0.1451, 'Q': 0.0435, 'R': -0.0777, 'S': 0.1921, 'T': 0.2204, 'V': 0.0726, 'W': -0.3555, 'Y': -0.1992}, 'W': {'A': -0.3609, 'C': -0.4051, 'D': -0.0261, 'E': -0.0378, 'F': -0.8217, 'G': -0.3788, 'H': -0.7715, 'I': -0.3651, 'K': 0.017, 'L': -0.3771, 'M': -0.4093, 'N': -0.5511, 'P': -0.6428, 'Q': -0.5805, 'R': -0.3863, 'S': -0.3809, 'T': -0.371, 'V': -0.3555, 'W': -1.0436, 'Y': -0.8617}, 'Y': {'A': -0.2128, 'C': -0.2514, 'D': 0.0362, 'E': 0.0262, 'F': -0.6489, 'G': -0.2356, 'H': -0.6027, 'I': -0.2054, 'K': 0.0275, 'L': -0.2171, 'M': -0.2485, 'N': -0.3917, 'P': -0.4368, 'Q': -0.417, 'R': -0.2948, 'S': -0.2308, 'T': -0.2168, 'V': -0.1992, 'W': -0.8617, 'Y': -0.687}}
    
    # list of all amino acids
    all_aas=['A', 'C', 'D', 'E', 'F', 'G', 'H', 'I', 'K', 'L', 'M', 'N', 'P', 'Q', 'R', 'S', 'T', 'V', 'W', 'Y']
    
    if len(exclude_aas)==20:
        raise goose_exceptions.GooseInputError('Cannot exclude all amino acids.')

    # remove AAs we are excluding
    for a in exclude_aas:
        all_aas.remove(a)

    # see if we need to make interacting sequence the same as sequence of interest
    if interacting_sequence is None:
        interacting_sequence = sequence_of_interest

    # calculate eps vectors
    eps_vectors=get_epsilon_vectors(sequence_of_interest, interacting_sequence)
    # sum them to get overall interaction vectors
    interaction_vectors=eps_vectors[0]+eps_vectors[1]

    # make empty string to hold new sequence
    new_sequence = ''

    # now we need to choose amino acids that match the interaction vector of our 
    # sequence of interest. 
    for i in range(0, len(sequence_of_interest)):
        # get the interaction value for the amino acid
        aa_interaction=interaction_vectors[i]
        # get the amino acids that are closest to the interaction value
        possible_vals={}
        for aa in all_aas:
            cur_tot=0
            for aa_interacting in interacting_sequence:
                cur_tot=cur_tot+epsilon_aas[aa][aa_interacting]
            possible_vals[aa]=cur_tot/len(interacting_sequence)
        possible_aas=sorted(possible_vals, key=lambda x:abs(possible_vals[x]-aa_interaction))
        # based on the class, allow variable amounts of 'randomness'
        # this lets us keep the chemical interaction specificity close
        # but explore maximal sequecne space.
        if possible_aas[0] in ['D', 'E', 'K', 'R']:
            randomness=2
        elif possible_aas[0] in ['W', 'Y', 'F']:
            randomness=3
        elif possible_aas[0] in ['A', 'I', 'L', 'V', 'M']:
            randomness=5
        elif possible_aas[0] in ['N', 'Q', 'S', 'T', 'G']:
            randomness=4
        else:
            randomness=3
        # choose amino acid to add
        new_sequence=new_sequence+random.choice(possible_aas[:randomness])
    return new_sequence



def create_seq_by_epsilon(sequence_of_interest, interacting_sequence=None, 
    exclude_aas=[]):
    '''
    Function that creates a sequence with the same total matching epsilon. Doesn't
    Care about linear space as much as previous function. 

    Parameters
    -----------
    sequence_of_interest : str
        The sequence that you want your generated sequence to
        have an epsilon value relative to.

    interacting_sequence : str
        The seqence that the sequence of interest is interacting with.
        If None, then it will be the same as sequence of interest (
        so will make something with the same homotypic interactions)

    exclude_aas : list
        A list of amino acids to exclude from the generated sequence. 
        Default = []        
    '''
    # pairwise interaction values for each amino acid and U (RNA) against all
    # amino acids and U. Uses mPiPi-GGv1.
    epsilon_aas={'A': {'A': 0.1729, 'C': 0.1479, 'D': 0.1122, 'E': 0.0991, 'F': -0.1814, 'G': 0.1355, 'H': -0.1472, 'I': 0.1184, 'K': 0.1819, 'L': 0.1302, 'M': 0.1617, 'N': 0.0295, 'P': 0.0992, 'Q': 0.0159, 'R': -0.1, 'S': 0.1578, 'T': 0.1822, 'V': 0.1273, 'W': -0.3609, 'Y': -0.2128}, 'C': {'A': 0.1479, 'C': 0.1244, 'D': 0.0871, 'E': 0.0731, 'F': -0.2188, 'G': 0.1127, 'H': -0.1831, 'I': 0.1859, 'K': 0.1586, 'L': 0.1751, 'M': 0.1461, 'N': 0.0009, 'P': 0.0673, 'Q': -0.0136, 'R': -0.1343, 'S': 0.1352, 'T': 0.16, 'V': 0.1827, 'W': -0.4051, 'Y': -0.2514}, 'D': {'A': 0.1122, 'C': 0.0871, 'D': 2.574, 'E': 2.4772, 'F': 0.0495, 'G': 0.0783, 'H': -1.1931, 'I': 0.1469, 'K': -2.5449, 'L': 0.136, 'M': 0.1068, 'N': -0.0382, 'P': 0.0172, 'Q': -0.0538, 'R': -2.4973, 'S': 0.0989, 'T': 0.1225, 'V': 0.1447, 'W': -0.0261, 'Y': 0.0362}, 'E': {'A': 0.0991, 'C': 0.0731, 'D': 2.4772, 'E': 2.3835, 'F': 0.0398, 'G': 0.0643, 'H': -1.1415, 'I': 0.1344, 'K': -2.4518, 'L': 0.1232, 'M': 0.0931, 'N': -0.056, 'P': -0.0018, 'Q': -0.0721, 'R': -2.3995, 'S': 0.0854, 'T': 0.1096, 'V': 0.1323, 'W': -0.0378, 'Y': 0.0262}, 'F': {'A': -0.1814, 'C': -0.2188, 'D': 0.0495, 'E': 0.0398, 'F': -0.6114, 'G': -0.2051, 'H': -0.566, 'I': -0.1716, 'K': 0.0014, 'L': -0.1832, 'M': -0.2144, 'N': -0.3574, 'P': -0.393, 'Q': -0.3819, 'R': -0.2339, 'S': -0.199, 'T': -0.1842, 'V': -0.1662, 'W': -0.8217, 'Y': -0.6489}, 'G': {'A': 0.1355, 'C': 0.1127, 'D': 0.0783, 'E': 0.0643, 'F': -0.2051, 'G': -0.0747, 'H': -0.1708, 'I': 0.1667, 'K': 0.1411, 'L': 0.1568, 'M': 0.1302, 'N': -0.0007, 'P': 0.0525, 'Q': -0.0154, 'R': -0.1286, 'S': -0.0354, 'T': 0.1446, 'V': 0.1647, 'W': -0.3788, 'Y': -0.2356}, 'H': {'A': -0.1472, 'C': -0.1831, 'D': -1.1931, 'E': -1.1415, 'F': -0.566, 'G': -0.1708, 'H': 0.8857, 'I': -0.1362, 'K': 1.3994, 'L': -0.1475, 'M': -0.178, 'N': -0.3183, 'P': -0.3451, 'Q': -0.3419, 'R': 1.289, 'S': -0.1641, 'T': -0.1492, 'V': -0.1313, 'W': -0.7715, 'Y': -0.6027}, 'I': {'A': 0.1184, 'C': 0.1859, 'D': 0.1469, 'E': 0.1344, 'F': -0.1716, 'G': 0.1667, 'H': -0.1362, 'I': 0.0416, 'K': 0.2281, 'L': 0.0644, 'M': 0.1023, 'N': 0.056, 'P': 0.149, 'Q': 0.0432, 'R': -0.0806, 'S': 0.1954, 'T': 0.2247, 'V': 0.0536, 'W': -0.3651, 'Y': -0.2054}, 'K': {'A': 0.1819, 'C': 0.1586, 'D': -2.5449, 'E': -2.4518, 'F': 0.0014, 'G': 0.1411, 'H': 1.3994, 'I': 0.2281, 'K': 2.2789, 'L': 0.2165, 'M': 0.1853, 'N': 0.0263, 'P': 0.1123, 'Q': 0.0126, 'R': 2.0644, 'S': 0.1688, 'T': 0.1977, 'V': 0.2233, 'W': 0.017, 'Y': 0.0275}, 'L': {'A': 0.1302, 'C': 0.1751, 'D': 0.136, 'E': 0.1232, 'F': -0.1832, 'G': 0.1568, 'H': -0.1475, 'I': 0.0644, 'K': 0.2165, 'L': 0.0794, 'M': 0.1153, 'N': 0.045, 'P': 0.1346, 'Q': 0.0319, 'R': -0.0924, 'S': 0.1849, 'T': 0.2138, 'V': 0.0761, 'W': -0.3771, 'Y': -0.2171}, 'M': {'A': 0.1617, 'C': 0.1461, 'D': 0.1068, 'E': 0.0931, 'F': -0.2144, 'G': 0.1302, 'H': -0.178, 'I': 0.1023, 'K': 0.1853, 'L': 0.1153, 'M': 0.1501, 'N': 0.0155, 'P': 0.0957, 'Q': 0.0016, 'R': -0.1241, 'S': 0.1566, 'T': 0.1844, 'V': 0.1118, 'W': -0.4093, 'Y': -0.2485}, 'N': {'A': 0.0295, 'C': 0.0009, 'D': -0.0382, 'E': -0.056, 'F': -0.3574, 'G': -0.0007, 'H': -0.3183, 'I': 0.056, 'K': 0.0263, 'L': 0.045, 'M': 0.0155, 'N': -0.127, 'P': -0.0983, 'Q': -0.1452, 'R': -0.2736, 'S': 0.0151, 'T': 0.0356, 'V': 0.0562, 'W': -0.5511, 'Y': -0.3917}, 'P': {'A': 0.0992, 'C': 0.0673, 'D': 0.0172, 'E': -0.0018, 'F': -0.393, 'G': 0.0525, 'H': -0.3451, 'I': 0.149, 'K': 0.1123, 'L': 0.1346, 'M': 0.0957, 'N': -0.0983, 'P': 0.0539, 'Q': -0.1179, 'R': -0.2801, 'S': 0.0821, 'T': 0.1148, 'V': 0.1451, 'W': -0.6428, 'Y': -0.4368}, 'Q': {'A': 0.0159, 'C': -0.0136, 'D': -0.0538, 'E': -0.0721, 'F': -0.3819, 'G': -0.0154, 'H': -0.3419, 'I': 0.0432, 'K': 0.0126, 'L': 0.0319, 'M': 0.0016, 'N': -0.1452, 'P': -0.1179, 'Q': -0.1638, 'R': -0.2956, 'S': 0.001, 'T': 0.0222, 'V': 0.0435, 'W': -0.5805, 'Y': -0.417}, 'R': {'A': -0.1, 'C': -0.1343, 'D': -2.4973, 'E': -2.3995, 'F': -0.2339, 'G': -0.1286, 'H': 1.289, 'I': -0.0806, 'K': 2.0644, 'L': -0.0924, 'M': -0.1241, 'N': -0.2736, 'P': -0.2801, 'Q': -0.2956, 'R': 2.08, 'S': -0.1167, 'T': -0.098, 'V': -0.0777, 'W': -0.3863, 'Y': -0.2948}, 'S': {'A': 0.1578, 'C': 0.1352, 'D': 0.0989, 'E': 0.0854, 'F': -0.199, 'G': -0.0354, 'H': -0.1641, 'I': 0.1954, 'K': 0.1688, 'L': 0.1849, 'M': 0.1566, 'N': 0.0151, 'P': 0.0821, 'Q': 0.001, 'R': -0.1167, 'S': 0.1455, 'T': 0.1698, 'V': 0.1921, 'W': -0.3809, 'Y': -0.2308}, 'T': {'A': 0.1822, 'C': 0.16, 'D': 0.1225, 'E': 0.1096, 'F': -0.1842, 'G': 0.1446, 'H': -0.1492, 'I': 0.2247, 'K': 0.1977, 'L': 0.2138, 'M': 0.1844, 'N': 0.0356, 'P': 0.1148, 'Q': 0.0222, 'R': -0.098, 'S': 0.1698, 'T': 0.1965, 'V': 0.2204, 'W': -0.371, 'Y': -0.2168}, 'V': {'A': 0.1273, 'C': 0.1827, 'D': 0.1447, 'E': 0.1323, 'F': -0.1662, 'G': 0.1647, 'H': -0.1313, 'I': 0.0536, 'K': 0.2233, 'L': 0.0761, 'M': 0.1118, 'N': 0.0562, 'P': 0.1451, 'Q': 0.0435, 'R': -0.0777, 'S': 0.1921, 'T': 0.2204, 'V': 0.0726, 'W': -0.3555, 'Y': -0.1992}, 'W': {'A': -0.3609, 'C': -0.4051, 'D': -0.0261, 'E': -0.0378, 'F': -0.8217, 'G': -0.3788, 'H': -0.7715, 'I': -0.3651, 'K': 0.017, 'L': -0.3771, 'M': -0.4093, 'N': -0.5511, 'P': -0.6428, 'Q': -0.5805, 'R': -0.3863, 'S': -0.3809, 'T': -0.371, 'V': -0.3555, 'W': -1.0436, 'Y': -0.8617}, 'Y': {'A': -0.2128, 'C': -0.2514, 'D': 0.0362, 'E': 0.0262, 'F': -0.6489, 'G': -0.2356, 'H': -0.6027, 'I': -0.2054, 'K': 0.0275, 'L': -0.2171, 'M': -0.2485, 'N': -0.3917, 'P': -0.4368, 'Q': -0.417, 'R': -0.2948, 'S': -0.2308, 'T': -0.2168, 'V': -0.1992, 'W': -0.8617, 'Y': -0.687}}
    
    # list of all amino acids
    all_aas=['A', 'C', 'D', 'E', 'F', 'G', 'H', 'I', 'K', 'L', 'M', 'N', 'P', 'Q', 'R', 'S', 'T', 'V', 'W', 'Y']
    
    if len(exclude_aas)==20:
        raise goose_exceptions.GooseInputError('Cannot exclude all aminon acids.')

    # remove excluded amino acids
    for aa in exclude_aas:
        all_aas.remove(aa)


    # make sequence interactor self if None is specified.
    if interacting_sequence==None:
        interacting_sequence=squence_of_interest
    
    # get epsilon value between sequence and interacting_sequence
    epsilon=get_epsilon_value(sequence_of_interest, interacting_sequence)
    eps_vectors=get_epsilon_vectors(sequence_of_interest, interacting_sequence)
    # make a starting sequence equal to length of the starting sequence. 
    starting_sequence=create_starting_sequence(len(sequence_of_interest), exclude_aas=exclude_aas)
    # using the random sequence as a starting point, now optimize back to the epsilon value
    # we want
    new_sequence=optimize_to_epsilon_value(starting_sequence, interacting_sequence, epsilon, exclude_aas=exclude_aas)
    return new_sequence



def create_attractive_or_repulsive_seq(objective_seq_length, interacting_sequence, 
    attractive_or_repulsive, exclude_aas=[]):
    '''
    Function that creates a sequence that is attractive to the interacting sequence.
    
    Parameters
    -----------
    objective_seq_length : int
        The length of the sequence that you want to generate.

    interacting_sequence : str
        The sequence that the sequence of interest is interacting with.

    attractive_or_repulsive : str
        Whether the sequence should be attractive or repulsive to the interacting sequence.

    exclude_aas : list
        A list of amino acids to exclude from the generated sequence. 
        Default = []
    '''

    if attractive_or_repulsive not in ['attractive', 'repulsive']:
        raise goose_exceptions.GooseInputError("attractive_or_repulsive must be either 'attractive' or 'repulsive'.")

    # pairwise interaction values for each amino acid and U (RNA) against all
    # amino acids and U. Uses mPiPi-GGv1.
    epsilon_aas={'A': {'A': 0.1729, 'C': 0.1479, 'D': 0.1122, 'E': 0.0991, 'F': -0.1814, 'G': 0.1355, 'H': -0.1472, 'I': 0.1184, 'K': 0.1819, 'L': 0.1302, 'M': 0.1617, 'N': 0.0295, 'P': 0.0992, 'Q': 0.0159, 'R': -0.1, 'S': 0.1578, 'T': 0.1822, 'V': 0.1273, 'W': -0.3609, 'Y': -0.2128}, 'C': {'A': 0.1479, 'C': 0.1244, 'D': 0.0871, 'E': 0.0731, 'F': -0.2188, 'G': 0.1127, 'H': -0.1831, 'I': 0.1859, 'K': 0.1586, 'L': 0.1751, 'M': 0.1461, 'N': 0.0009, 'P': 0.0673, 'Q': -0.0136, 'R': -0.1343, 'S': 0.1352, 'T': 0.16, 'V': 0.1827, 'W': -0.4051, 'Y': -0.2514}, 'D': {'A': 0.1122, 'C': 0.0871, 'D': 2.574, 'E': 2.4772, 'F': 0.0495, 'G': 0.0783, 'H': -1.1931, 'I': 0.1469, 'K': -2.5449, 'L': 0.136, 'M': 0.1068, 'N': -0.0382, 'P': 0.0172, 'Q': -0.0538, 'R': -2.4973, 'S': 0.0989, 'T': 0.1225, 'V': 0.1447, 'W': -0.0261, 'Y': 0.0362}, 'E': {'A': 0.0991, 'C': 0.0731, 'D': 2.4772, 'E': 2.3835, 'F': 0.0398, 'G': 0.0643, 'H': -1.1415, 'I': 0.1344, 'K': -2.4518, 'L': 0.1232, 'M': 0.0931, 'N': -0.056, 'P': -0.0018, 'Q': -0.0721, 'R': -2.3995, 'S': 0.0854, 'T': 0.1096, 'V': 0.1323, 'W': -0.0378, 'Y': 0.0262}, 'F': {'A': -0.1814, 'C': -0.2188, 'D': 0.0495, 'E': 0.0398, 'F': -0.6114, 'G': -0.2051, 'H': -0.566, 'I': -0.1716, 'K': 0.0014, 'L': -0.1832, 'M': -0.2144, 'N': -0.3574, 'P': -0.393, 'Q': -0.3819, 'R': -0.2339, 'S': -0.199, 'T': -0.1842, 'V': -0.1662, 'W': -0.8217, 'Y': -0.6489}, 'G': {'A': 0.1355, 'C': 0.1127, 'D': 0.0783, 'E': 0.0643, 'F': -0.2051, 'G': -0.0747, 'H': -0.1708, 'I': 0.1667, 'K': 0.1411, 'L': 0.1568, 'M': 0.1302, 'N': -0.0007, 'P': 0.0525, 'Q': -0.0154, 'R': -0.1286, 'S': -0.0354, 'T': 0.1446, 'V': 0.1647, 'W': -0.3788, 'Y': -0.2356}, 'H': {'A': -0.1472, 'C': -0.1831, 'D': -1.1931, 'E': -1.1415, 'F': -0.566, 'G': -0.1708, 'H': 0.8857, 'I': -0.1362, 'K': 1.3994, 'L': -0.1475, 'M': -0.178, 'N': -0.3183, 'P': -0.3451, 'Q': -0.3419, 'R': 1.289, 'S': -0.1641, 'T': -0.1492, 'V': -0.1313, 'W': -0.7715, 'Y': -0.6027}, 'I': {'A': 0.1184, 'C': 0.1859, 'D': 0.1469, 'E': 0.1344, 'F': -0.1716, 'G': 0.1667, 'H': -0.1362, 'I': 0.0416, 'K': 0.2281, 'L': 0.0644, 'M': 0.1023, 'N': 0.056, 'P': 0.149, 'Q': 0.0432, 'R': -0.0806, 'S': 0.1954, 'T': 0.2247, 'V': 0.0536, 'W': -0.3651, 'Y': -0.2054}, 'K': {'A': 0.1819, 'C': 0.1586, 'D': -2.5449, 'E': -2.4518, 'F': 0.0014, 'G': 0.1411, 'H': 1.3994, 'I': 0.2281, 'K': 2.2789, 'L': 0.2165, 'M': 0.1853, 'N': 0.0263, 'P': 0.1123, 'Q': 0.0126, 'R': 2.0644, 'S': 0.1688, 'T': 0.1977, 'V': 0.2233, 'W': 0.017, 'Y': 0.0275}, 'L': {'A': 0.1302, 'C': 0.1751, 'D': 0.136, 'E': 0.1232, 'F': -0.1832, 'G': 0.1568, 'H': -0.1475, 'I': 0.0644, 'K': 0.2165, 'L': 0.0794, 'M': 0.1153, 'N': 0.045, 'P': 0.1346, 'Q': 0.0319, 'R': -0.0924, 'S': 0.1849, 'T': 0.2138, 'V': 0.0761, 'W': -0.3771, 'Y': -0.2171}, 'M': {'A': 0.1617, 'C': 0.1461, 'D': 0.1068, 'E': 0.0931, 'F': -0.2144, 'G': 0.1302, 'H': -0.178, 'I': 0.1023, 'K': 0.1853, 'L': 0.1153, 'M': 0.1501, 'N': 0.0155, 'P': 0.0957, 'Q': 0.0016, 'R': -0.1241, 'S': 0.1566, 'T': 0.1844, 'V': 0.1118, 'W': -0.4093, 'Y': -0.2485}, 'N': {'A': 0.0295, 'C': 0.0009, 'D': -0.0382, 'E': -0.056, 'F': -0.3574, 'G': -0.0007, 'H': -0.3183, 'I': 0.056, 'K': 0.0263, 'L': 0.045, 'M': 0.0155, 'N': -0.127, 'P': -0.0983, 'Q': -0.1452, 'R': -0.2736, 'S': 0.0151, 'T': 0.0356, 'V': 0.0562, 'W': -0.5511, 'Y': -0.3917}, 'P': {'A': 0.0992, 'C': 0.0673, 'D': 0.0172, 'E': -0.0018, 'F': -0.393, 'G': 0.0525, 'H': -0.3451, 'I': 0.149, 'K': 0.1123, 'L': 0.1346, 'M': 0.0957, 'N': -0.0983, 'P': 0.0539, 'Q': -0.1179, 'R': -0.2801, 'S': 0.0821, 'T': 0.1148, 'V': 0.1451, 'W': -0.6428, 'Y': -0.4368}, 'Q': {'A': 0.0159, 'C': -0.0136, 'D': -0.0538, 'E': -0.0721, 'F': -0.3819, 'G': -0.0154, 'H': -0.3419, 'I': 0.0432, 'K': 0.0126, 'L': 0.0319, 'M': 0.0016, 'N': -0.1452, 'P': -0.1179, 'Q': -0.1638, 'R': -0.2956, 'S': 0.001, 'T': 0.0222, 'V': 0.0435, 'W': -0.5805, 'Y': -0.417}, 'R': {'A': -0.1, 'C': -0.1343, 'D': -2.4973, 'E': -2.3995, 'F': -0.2339, 'G': -0.1286, 'H': 1.289, 'I': -0.0806, 'K': 2.0644, 'L': -0.0924, 'M': -0.1241, 'N': -0.2736, 'P': -0.2801, 'Q': -0.2956, 'R': 2.08, 'S': -0.1167, 'T': -0.098, 'V': -0.0777, 'W': -0.3863, 'Y': -0.2948}, 'S': {'A': 0.1578, 'C': 0.1352, 'D': 0.0989, 'E': 0.0854, 'F': -0.199, 'G': -0.0354, 'H': -0.1641, 'I': 0.1954, 'K': 0.1688, 'L': 0.1849, 'M': 0.1566, 'N': 0.0151, 'P': 0.0821, 'Q': 0.001, 'R': -0.1167, 'S': 0.1455, 'T': 0.1698, 'V': 0.1921, 'W': -0.3809, 'Y': -0.2308}, 'T': {'A': 0.1822, 'C': 0.16, 'D': 0.1225, 'E': 0.1096, 'F': -0.1842, 'G': 0.1446, 'H': -0.1492, 'I': 0.2247, 'K': 0.1977, 'L': 0.2138, 'M': 0.1844, 'N': 0.0356, 'P': 0.1148, 'Q': 0.0222, 'R': -0.098, 'S': 0.1698, 'T': 0.1965, 'V': 0.2204, 'W': -0.371, 'Y': -0.2168}, 'V': {'A': 0.1273, 'C': 0.1827, 'D': 0.1447, 'E': 0.1323, 'F': -0.1662, 'G': 0.1647, 'H': -0.1313, 'I': 0.0536, 'K': 0.2233, 'L': 0.0761, 'M': 0.1118, 'N': 0.0562, 'P': 0.1451, 'Q': 0.0435, 'R': -0.0777, 'S': 0.1921, 'T': 0.2204, 'V': 0.0726, 'W': -0.3555, 'Y': -0.1992}, 'W': {'A': -0.3609, 'C': -0.4051, 'D': -0.0261, 'E': -0.0378, 'F': -0.8217, 'G': -0.3788, 'H': -0.7715, 'I': -0.3651, 'K': 0.017, 'L': -0.3771, 'M': -0.4093, 'N': -0.5511, 'P': -0.6428, 'Q': -0.5805, 'R': -0.3863, 'S': -0.3809, 'T': -0.371, 'V': -0.3555, 'W': -1.0436, 'Y': -0.8617}, 'Y': {'A': -0.2128, 'C': -0.2514, 'D': 0.0362, 'E': 0.0262, 'F': -0.6489, 'G': -0.2356, 'H': -0.6027, 'I': -0.2054, 'K': 0.0275, 'L': -0.2171, 'M': -0.2485, 'N': -0.3917, 'P': -0.4368, 'Q': -0.417, 'R': -0.2948, 'S': -0.2308, 'T': -0.2168, 'V': -0.1992, 'W': -0.8617, 'Y': -0.687}}

    # list of all amino acids
    all_aas=['A', 'C', 'D', 'E', 'F', 'G', 'H', 'I', 'K', 'L', 'M', 'N', 'P', 'Q', 'R', 'S', 'T', 'V', 'W', 'Y']
    
    # check number of amino acids excluded
    if len(set(exclude_aas))==20:
        raise goose_exceptions.GooseInputError('Cannot exclude all aminon acids.')


    # string to hold new sequence we are generating
    new_sequence=''

    # loop through the length of the sequence we want to generate
    for i in range(0, objective_seq_length):
        # hold possible aas
        possible_aas={}

        # if seq over 8 amino acids, start constrainin all_aas
        if len(new_sequence) > 8:
            all_aas=return_constrained_aa_list(new_sequence)

        # loop through all aas
        for aa in all_aas:
            total_interaction=0
            # loop through all aas in the interacting sequence
            for aa2 in interacting_sequence:
                # get the epsilon value between the aa and the interacting sequence
                interaction_value=total_interaction+epsilon_aas[aa][aa2]
            possible_aas[aa]=interaction_value
        # choose all amino acids with + or - value depending on if we want attractive or repulsive
        if attractive_or_repulsive=='attractive':
            use_these_aas={k:v for k,v in possible_aas.items() if v<0}
        else:
            use_these_aas={k:v for k,v in possible_aas.items() if v>0}

        # if we don't have any amino acids that are attractive or repulsive, choose the most negative or positive
        if use_these_aas=={}:
            if attractive_or_repulsive=='attractive':
                # sort by most negative to positive
                sorted_aas=sorted(possible_aas.items(), key=lambda x: x[1])
            else:
                # sort by most negative to positive
                sorted_aas=sorted(possible_aas.items(), key=lambda x: x[1], reverse=True)
            # add the best aa we can
            new_sequence+=sorted_aas[0][0]
        else:
            if attractive_or_repulsive=='attractive':
                sorted_aas=sorted(use_these_aas.items(), key=lambda x: x[1])
            else:
                sorted_aas=sorted(use_these_aas.items(), key=lambda x: x[1], reverse=True)
            if len(sorted_aas) > 8:
                sorted_aas=sorted_aas[:8]
            new_sequence+=random.choice(sorted_aas)[0]
    return new_sequence

