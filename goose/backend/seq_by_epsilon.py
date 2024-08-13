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
'''

from collections import Counter
import itertools
import random

import matplotlib.pyplot as plt
import numpy as np

from finches import epsilon_calculation
from finches.forcefields.mpipi import Mpipi_model
from finches.forcefields.calvados import calvados_model
from goose.backend import lists
from goose.backend.lists import precomputed_epsilon
from goose import goose_exceptions


# -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
#        -=-=-=-=-=-= Code for epsilon related stuffs =-=-=-=-=-=-=
# -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=

def load_IMC_object(model='Mpipi_GGv1'):
    '''
    Function to load a FINCHES IMC object.
    Lets us only load it once and then can use it iteratively. 

    Parameters
    -----------
    model : string
        The specific model parameters we are using
        default = Mpipi_GGv1
        options are 'Mpipi_GGv1', 'CALVADOS2'

    '''

    # check the current implementations of forcefields.
    if model not in ['Mpipi_GGv1', 'CALVADOS2']:
        raise goose_exceptions.GooseInputError('Only Mpipi_GGv1 and CALVADOS2 forcefields have been implemented.')
    
    # initialize forcefield parameters
    if model in ['Mpipi_GGv1']:
        ff_model = Mpipi_model(model)
    else:
        ff_model = calvados_model(model)
    
    # make IMC Object
    IMC_object = epsilon_calculation.InteractionMatrixConstructor(ff_model)
    return IMC_object

def get_interaction_vectors(seq1, seq2, IMC_object, approach='mean'):
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

    IMC_object : Object
        A loaded FINCHES IMC object.

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


def get_epsilon_vectors(seq1, seq2, IMC_object):
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

    IMC_object : Object
        A loaded FINCHES IMC object.

    Returns
    -------
    list of np.arrays
        Returns a list of np.arrays. First item in the list is the
        interaction vector for seq1, the second is for seq2. 
    '''
    # use IMC_object to calculate heterotypic matrix
    interaction_matrix=IMC_object.calculate_epsilon_vectors(seq1,seq2)

    return [interaction_matrix[0], interaction_matrix[1]]


def get_epsilon_value(seq1, seq2, IMC_object):
    '''
    Function to get the epsilon value using the Mpipi_GGv1 model. 

    Parameters
    -----------
    seq1 : string
        the amino acid sequence for seq1 as a string

    seq2 : string
        the amino acid sequence for seq2 as a string

    IMC_object : Object
        A loaded FINCHES IMC object.

    Returns
    -------
    list of np.arrays
        Returns a list of np.arrays. First item in the list is the
        interaction vector for seq1, the second is for seq2. 
    '''    
    # use IMC_object to calculate heterotypic matrix
    epslon_value=IMC_object.calculate_epsilon_value(seq1,seq2)

    return epslon_value


# -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
#    -=-=-=-=-=-=-=-=- Code for modifiying sequences =-=-=-=-=-=-=-=-=-
# -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=



def optimize_to_epsilon_value(starting_sequence, interacting_sequence, objective_epsilon,
    allowed_error=0.1, optimization_iterations=None, exclude_aas=[], model='Mpipi_GGv1',
    preloaded_IMC_object=None, return_best_sequence=False, maximal_optimization=False):
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

    exclude_aas : list
        A list of amino acids to exclude from being added to the sequence during optimization.

    model : str
        The specific model parameters we are using
        default = Mpipi_GGv1
        options are 'Mpipi_GGv1', 'CALVADOS2'

    preloaded_IMC_object : object
        If we are using this optimizer and we already loaded the FINCHES model, we can
        skip that here to speed things up. 
        Default is None. 

    return_best_sequence : bool
        whether to just return the best sequence we could get to
        the objective epsilon value. 
        Default is False. 

    maximal_optimization : bool
        Whether to optimize to the maximal extent possible. 
        Reduces the sequence space explored, so default is False

    Returns
    -------
    str
        The optimized sequence that has the closest epsilon value to the objective epsilon value.
    '''

    # check the current implementations of forcefields.
    if model not in lists.implimented_finches_models:
        raise goose_exceptions.GooseInputError(f'Only {lists.implimented_finches_models} forcefields have been implemented.')

    # pairwise interaction values for each amino acid and U (RNA) against all
    # amino acids and U. Uses mPiPi-GGv1.
    epsilon_aas=precomputed_epsilon[model]
    
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

    if preloaded_IMC_object==None:
        # load IMC object
        loaded_model = load_IMC_object(model)
    else:
        loaded_model=preloaded_IMC_object

    # calculate the epsilon value of the starting sequence
    eps_vectors=get_epsilon_vectors(starting_sequence, interacting_sequence, loaded_model)
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
            eps_vectors=get_epsilon_vectors(current_sequence, interacting_sequence, loaded_model)
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

        # see if we are trying to maximize the optimziation
        if maximal_optimization==True:
            new_aa=possible_aas[0]
        else:
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
            new_aa=random.choice(possible_aas[:randomness])


        # change the amino acid
        current_sequence=current_sequence[:aa_to_change]+new_aa+current_sequence[aa_to_change+1:]
        # recalculate the epsilon value
        eps_vectors=get_epsilon_vectors(current_sequence, interacting_sequence, loaded_model)
        current_epsilon=np.sum(eps_vectors)
        # calculate the new difference
        current_epsilon_diff=objective_epsilon-current_epsilon
        # if the difference is less than the allowed error, break
        if abs(current_epsilon_diff) < allowed_error:
            return current_sequence
    if return_best_sequence==False:
        raise goose_exceptions.GooseFail('Failed to optimize to epsilon value.')
    else:
        return current_sequence



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
    
    if input_sequence=='':
        finaas=[]
        for aa in all_aas:
            if aa not in exclude_aas:
                finaas.append(aa)
        return finaas

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
    if 'C' in all_aas:
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
    exclude_aas=[], model='Mpipi_GGv1'):
    '''
    Function to create a sequence with approximately similar 
    interaction vectors to a sequence of interest. This function
    has some tolerance in that it will randomly choose from the
    best 2 amino acids at each position, so the total epsilon value
    will differ. However, the patterna cross the sequence should be quite close. 

    Parameters
    -----------
    sequence_of_interest : str
        The sequence that you want your generated sequence to have an epsilon value relative to.

    interacting_sequence : str
        If you want to make a sequence that has an epsilon equal
        to that between your sequence_of_interest and some interacting sequence,
        set this equal to the interacting sequence. If None, then it will be the same as sequence of interest (
        so will make something with the same homotypic interactions).

        ex. If you want to make a FUS variant that interacts with hnRNPA1 in a specific way, 
        set sequence_of_interest=FUS and interacting_sequence=hnRNPA1.

        Default = None

    exclude_aas : list
        A list of amino acids to exclude from the generated sequence. 
        Default = []

    model : str
        The specific model parameters we are using
        default = Mpipi_GGv1
        options are 'Mpipi_GGv1', 'CALVADOS2'

    '''
    # check the current implementations of forcefields.
    if model not in lists.implimented_finches_models:
        raise goose_exceptions.GooseInputError(f'Only {lists.implimented_finches_models} forcefields have been implemented.')

    # pairwise interaction values for each amino acid and U (RNA) against all
    # amino acids and U. Uses mPiPi-GGv1.
    epsilon_aas=precomputed_epsilon[model]
    
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

    # load imc object
    loaded_model = load_IMC_object(model)

    # calculate eps vectors
    eps_vectors=get_epsilon_vectors(sequence_of_interest, interacting_sequence, loaded_model)
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
    exclude_aas=[], model='Mpipi_GGv1'):
    '''
    Function that creates a sequence with the same total matching epsilon. Doesn't
    Care about linear space as much as previous function. 

    Parameters
    -----------
    sequence_of_interest : str
        The sequence that you want your generated sequence to have an epsilon value relative to.

    interacting_sequence : str
        If you want to make a sequence that has an epsilon equal
        to that between your sequence_of_interest and some interacting sequence,
        set this equal to the interacting sequence. If None, then it will be the same as sequence of interest (
        so will make something with the same homotypic interactions).
        
        ex. If you want to make a FUS variant that interacts with hnRNPA1 in a specific way, 
        set sequence_of_interest=FUS and interacting_sequence=hnRNPA1.

        Default = None

    exclude_aas : list
        A list of amino acids to exclude from the generated sequence. 
        Default = []        

    model : str
        The specific model parameters we are using
        default = Mpipi_GGv1
        options are 'Mpipi_GGv1', 'CALVADOS2'

    '''
    # check the current implementations of forcefields.
    if model not in lists.implimented_finches_models:
        raise goose_exceptions.GooseInputError(f'Only {lists.implimented_finches_models} forcefields have been implemented.')

    # pairwise interaction values for each amino acid and U (RNA) against all
    # amino acids and U. Uses mPiPi-GGv1.
    epsilon_aas=precomputed_epsilon[model]
    
    # list of all amino acids
    all_aas=['A', 'C', 'D', 'E', 'F', 'G', 'H', 'I', 'K', 'L', 'M', 'N', 'P', 'Q', 'R', 'S', 'T', 'V', 'W', 'Y']
    
    if len(exclude_aas)==20:
        raise goose_exceptions.GooseInputError('Cannot exclude all aminon acids.')

    # remove excluded amino acids
    for aa in exclude_aas:
        all_aas.remove(aa)

    # make sequence interactor self if None is specified.
    if interacting_sequence==None:
        interacting_sequence=sequence_of_interest
    
    # load imc object
    loaded_model = load_IMC_object(model)

    # get epsilon value between sequence and interacting_sequence
    eps_vectors=get_epsilon_vectors(sequence_of_interest, interacting_sequence, loaded_model)
    epsilon = np.sum(eps_vectors)
    # make a starting sequence equal to length of the starting sequence. 
    starting_sequence=create_starting_sequence(len(sequence_of_interest), exclude_aas=exclude_aas)
    # using the random sequence as a starting point, now optimize back to the epsilon value
    # we want
    new_sequence=optimize_to_epsilon_value(starting_sequence, interacting_sequence, 
        epsilon, exclude_aas=exclude_aas, model=model, preloaded_IMC_object=loaded_model)
    return new_sequence



def create_attractive_or_repulsive_seq(objective_seq_length, interacting_sequence, 
    attractive_or_repulsive, exclude_aas=[], model='Mpipi_GGv1'):
    '''
    Function that creates a sequence that is attractive
    or repulsive to the interacting sequence. You can also
    specify the length of the sequence. 
    
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

    model : str
        The specific model parameters we are using
        default = Mpipi_GGv1
        options are 'Mpipi_GGv1', 'CALVADOS2'

    '''
    # check the current implementations of forcefields.
    if model not in lists.implimented_finches_models:
        raise goose_exceptions.GooseInputError(f'Only {lists.implimented_finches_models} forcefields have been implemented.')

    # load the model
    loaded_model = load_IMC_object(model)

    if attractive_or_repulsive not in ['attractive', 'repulsive']:
        raise goose_exceptions.GooseInputError("attractive_or_repulsive must be either 'attractive' or 'repulsive'.")

    # pairwise interaction values for each amino acid and U (RNA) against all
    # amino acids and U. Uses mPiPi-GGv1.
    epsilon_aas=precomputed_epsilon[model]

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
            all_aas=return_constrained_aa_list(new_sequence, exclude_aas=exclude_aas)
        else:
            all_aas=return_constrained_aa_list('', exclude_aas=exclude_aas)



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

        # if we have fewer than 3 amino acids...
        if len(use_these_aas) < 3:
            if attractive_or_repulsive=='attractive':
                # sort by most negative to positive
                sorted_aas=sorted(possible_aas.items(), key=lambda x: x[1])
            else:
                # sort by most negative to positive
                sorted_aas=sorted(possible_aas.items(), key=lambda x: x[1], reverse=True)
            # add from best 5
            new_sequence+=random.choice(sorted_aas[:3])[0]
        else:
            if attractive_or_repulsive=='attractive':
                sorted_aas=sorted(use_these_aas.items(), key=lambda x: x[1])
            else:
                sorted_aas=sorted(use_these_aas.items(), key=lambda x: x[1], reverse=True)
            if len(sorted_aas) > 8:
                sorted_aas=sorted_aas[:8]
            new_sequence+=random.choice(sorted_aas)[0]

    # check to see if sequence is attractive or repulsive.
    final_eps = get_epsilon_value(new_sequence, interacting_sequence, loaded_model)
    if final_eps < 0 and attractive_or_repulsive=='repulsive':
        # try optimizing the sequence to become positive
        attempted_sequence=optimize_to_epsilon_value(new_sequence, interacting_sequence, 10,
            exclude_aas=exclude_aas, model=model, preloaded_IMC_object=loaded_model, return_best_sequence=True)
    elif final_eps > 0 and attractive_or_repulsive=='attractive':
        # try optimizing the sequence to become positive
        attempted_sequence=optimize_to_epsilon_value(new_sequence, interacting_sequence, -10,
            exclude_aas=exclude_aas, model=model, preloaded_IMC_object=loaded_model, return_best_sequence=True)
    else:
        # if we wanted it to be attractive and it was, 
        # or if we wanted it to be repulsive and it was, 
        # return it.
        return new_sequence

    # now check the attempted sequence
    attempted_eps = get_epsilon_value(attempted_sequence, interacting_sequence, loaded_model)
    if attempted_eps > 0 and attractive_or_repulsive=='repulsive':
        return attempted_sequence
    elif attempted_eps < 0 and attractive_or_repulsive=='attractive':
        return attempted_sequence
    else:
        # if we couldn't optimize it to be more attractive or repulsive,
        # raise goose_exceptions.GooseFail
        raise goose_exceptions.GooseFail(f'Failed to create {attractive_or_repulsive} sequence.')
            

def increase_epsilon(sequence_of_interest, interacting_sequence=None,
                    num_iterations=None, exclude_aas=[], model='Mpipi_GGv1',
                    maximal_optimization=False):
    '''
    Function to increase the epsilon value of a sequence relative to another sequence. 

    Parameters
    -----------
    sequence_of_interest : str
        The sequence that you want to increase or decrease the epsilon value of.

    interacting_sequence : str
        The sequence that the sequence of interest is interacting with.
        If None, this will just do homotypic interactions (will do 
        sequence_of_interest with itself)

    num_iterations : int
        The number of iterations to run the optimization. 
        If none, default is length of sequence_of_interest/10

    exclude_aas : list
        A list of amino acids to exclude from the generated sequence. 
        Default = []

    model : str
        The specific model parameters we are using
        default = Mpipi_GGv1
        options are 'Mpipi_GGv1', 'CALVADOS2'

    maximal_optimization : bool
        Whether to optimize to the maximal extent possible. 
        Reduces the sequence space explored, so default is False

    '''
    # check the current implementations of forcefields.
    if model not in lists.implimented_finches_models:
        raise goose_exceptions.GooseInputError(f'Only {lists.implimented_finches_models} forcefields have been implemented.')

    # load the model
    loaded_model = load_IMC_object(model)

    # if interacting sequence is None, set it to the sequence of interest
    if num_iterations is None:
        num_iterations=int(len(sequence_of_interest)/10)
    if num_iterations==0:
        num_iterations=1

    return optimize_to_epsilon_value(sequence_of_interest, interacting_sequence, 1000,
        allowed_error=0.1, optimization_iterations=num_iterations, exclude_aas=exclude_aas,
        model=model, preloaded_IMC_object=loaded_model, return_best_sequence=True,
        maximal_optimization=maximal_optimization)

def decrease_epsilon(sequence_of_interest, interacting_sequence=None,
                    num_iterations=None, exclude_aas=[], model='Mpipi_GGv1',
                    maximal_optimization=False):
    '''
    Function to decrease the epsilon value of a sequence relative to another sequence. 

    Parameters
    -----------
    sequence_of_interest : str
        The sequence that you want to increase or decrease the epsilon value of.

    interacting_sequence : str
        The sequence that the sequence of interest is interacting with.
        If None, this will just do homotypic interactions (will do 
        sequence_of_interest with itself)

    num_iterations : int
        The number of iterations to run the optimization. 
        If none, default is length of sequence_of_interest/10

    exclude_aas : list
        A list of amino acids to exclude from the generated sequence. 
        Default = []

    model : str
        The specific model parameters we are using
        default = Mpipi_GGv1
        options are 'Mpipi_GGv1', 'CALVADOS2'

    maximal_optimization : bool
        Whether to optimize to the maximal extent possible. 
        Reduces the sequence space explored, so default is False

    '''
    # check the current implementations of forcefields.
    if model not in lists.implimented_finches_models:
        raise goose_exceptions.GooseInputError(f'Only {lists.implimented_finches_models} forcefields have been implemented.')

    # load the model
    loaded_model = load_IMC_object(model)

    # if interacting sequence is None, set it to the sequence of interest
    if num_iterations is None:
        num_iterations=int(len(sequence_of_interest)/10)
    if num_iterations==0:
        num_iterations=1

    return optimize_to_epsilon_value(sequence_of_interest, interacting_sequence, -1000,
        allowed_error=0.1, optimization_iterations=num_iterations, exclude_aas=exclude_aas,
        model=model, preloaded_IMC_object=loaded_model, return_best_sequence=True,
        maximal_optimization=maximal_optimization)

