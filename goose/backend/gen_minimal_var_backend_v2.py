'''
objective of this sequence is to change the fewest amino acids
possible to get the desired sequence. 
'''

from goose.backend.protein import Protein as pr
from goose.backend.amino_acids import AminoAcid as AA
from goose.backend.goose_tools import check_valid_kwargs
import random


def identify_sequence_change(input_sequence, **kwargs):
    '''
    Basically looks at the difference in the input sequence and specified
    parameters to identify what to change.

    Parameters
    -----------
    input_sequence: string
        The sequence to be changed

    **kwargs: dictionary:
        the properties specified by the user

    '''
    valid_kwargs=['FCR', 'NCPR', 'hydropathy', 'kappa']
    check_valid_kwargs(kwargs, valid_kwargs)
    # first calculate the values for the input sequence
    initial_parameters = pr(input_sequence).properties

    # now figure out differences between input sequence and the objective parameters
    sequence_changes={}
    for param in kwargs.keys():
        if initial_parameters[param] != kwargs[param]:
            sequence_changes[param]=kwargs[param]-initial_parameters[param]
    return sequence_changes


def identify_best_hydro_change_diff_class(input_sequence, hydro_change):
    '''
    Function to identify the best residue to change for altering the
    hydropathy of a sequence when you cannot keep the residue the same
    within class. Basically prioritizes to minimze sequence chemistry
    changes but does allow chemistry to change. 
    '''
    # first checks hydrophobics and polar. This is a big change but less
    # bad than changing aromatics or charged, etc. This could change as far as order in the future.
    # order of residues to check 
    aromatics = ['F', 'W', 'Y']
    polar = ['Q', 'N', 'S', 'T']
    aliphatics = ['I', 'V', 'L', 'A', 'M']
    positive = ['K', 'R']
    negative = ['D', 'E']
    special_cases = ['C', 'P', 'H','G']


    best_aas = []
    # iterate through input sequence and choose best aa(s) to change.
    for ind, res in enumerate(input_sequence):
        if hydro_change > 0:
            if res in polar:
                best_aas.append(ind)
                change_to = aliphatics
        else:
            if res in aliphatics:
                best_aas.append(ind)
                change_to = polar

    if best_aas == []:
        for ind, res in enumerate(input_sequence):
            if res in special_cases:
                best_aas.append(ind)
        change_to = polar
        change_to.extend(aliphatics)
    if best_aas == []:
        for ind, res in enumerate(input_sequence):
            if hydro_change > 0:
                if res in negative or res in positive:
                    best_aas.append(ind)
            else:
                if res in aromatics:
                    best_aas.append(ind)

    # if for some reason we haven't found anything easy to change, just go for anything.
    if best_aas == []:
        best_aas = [aa for aa in range(0, len(input_sequence))]
        change_to = ['S', 'T', 'N', 'Q', 'G', 'A', 'L', 'V', 'I', 'M', 'P', 'C', 'Y', 'W', 'F', 'H', 'D', 'E', 'K', 'R']

    # keep track of best change
    best_change=10000000

    # first try to change a residue within class.
    # just try to find the biggest within class change possible.
    for ind in best_aas:
        res = input_sequence[ind]
        cur_hydro = AA.hydro(res)
        if len(change_to)>1:
            for possible_res in change_to:
                if possible_res != res:
                    hydro_diff = AA.hydro(possible_res) - cur_hydro
                    if abs(hydro_change-hydro_diff) < best_change:
                        best_change = abs(hydro_change-hydro_diff)
                        possible_in_class_changes=[[ind, possible_res]]
                    if abs(hydro_change-hydro_diff) == best_change:
                        possible_in_class_changes.append([ind, possible_res])
    if len(possible_in_class_changes)==1:
        return possible_in_class_changes
    elif len(possible_in_class_changes)>1:
        return possible_in_class_changes[random.randint(0, len(possible_in_class_changes)-1)]
    else:
        # can't do an in class change. Now need to prioritize a bit...
        # Going to choose class to change in specific order to see if I can get that first.
        raise Exception('Unable to change residue within class.')    
    


def identify_best_within_class_hydro_change(input_sequence, hydro_change):
    '''
    Function to identify the best hydropathy change for min var. Different than the 
    within class hydro optimization function in that this is meant to find a single
    residue and change it whereas the hydropathy optimization allows for changing as
    many residues as possible. Falls back to identify_best_hydro_change_diff_class function.

    Parameters
    ----------
    input_sequence: string
        The sequence to be changed

    hydro_change: float
        The desired change in hydropathy

    returns
    -------
    res_index : int
        returns the index for the residue to change
    '''
    
    # residues that you can change to make a difference in hydropathy within class
    aa_class_dict = {'aromatic' : ['F', 'W', 'Y'], 'polar' : ['Q', 'N', 'S', 'T'], 'positive' : ['K', 'R'], 'negative' : ['D', 'E'], 'hydrophobic' : ['I', 'V', 'L', 'A', 'M'], 'G':['G'], 'H':['H'], 'P':['P'], 'C':['C']}

    # keep track of best change
    best_change=10000000
    # keep track of possible in track changes
    possible_in_class_changes=[]

    # first try to change a residue within class.
    # just try to find the biggest within class change possible.
    for ind, res in enumerate(input_sequence):
        cur_hydro = AA.hydro(res)
        cur_class = AA(res).AA_class
        within_class_changes=aa_class_dict[cur_class]
        if len(within_class_changes)>1:
            if within_class_changes!=['D', 'E']:
                for possible_res in within_class_changes:
                    if possible_res != res:
                        hydro_diff = AA.hydro(possible_res) - cur_hydro
                        if abs(hydro_change-hydro_diff) < best_change:
                            best_change = abs(hydro_change-hydro_diff)
                            possible_in_class_changes=[[ind, possible_res]]
                        if abs(hydro_change-hydro_diff) == best_change:
                            possible_in_class_changes.append([ind, possible_res])

    if len(possible_in_class_changes)==1:
        return possible_in_class_changes
    elif len(possible_in_class_changes)>1:
        return possible_in_class_changes[random.randint(0, len(possible_in_class_changes)-1)]
    else:
        # can't do an in class change. Now need to prioritize a bit...
        # Going to choose class to change in specific order to see if I can get that first.
        raise Exception('Unable to change residue within class.')


def identify_best_FCR_change(input_sequence, FCR_change, charge_change=None):
    '''
    Function to identify the best residues to change in a sequence
    when adjusting FCR. 

    Basically uses a simple 'change hierarchy' in that it tries to
    minimize chemical changes in the output sequence. For example, increasing
    FCR at the expense of aromatics would be a huge change, so we try to avoid this. 

    Parameters
    ----------
    input_sequence: string
        The sequence to be changed

    FCR_change: string
        The desired change in FCR
        Options: 'increase' or 'decrease'

    charge_change : string
        the desired charged residue to change. Let's you choose
        negative or positive. Default is random. 


    Returns
    -------
    res_index : int
        returns the index for the residue to change
    '''
    # TO DO: Add in some checks here.

    if FCR_change == 'decrease':
        if charge_change == None:
            change_residues=['D', 'E', 'K']
        elif charge_change == 'positive':
            change_residues = ['K']
        elif charge_Change == 'negative':
            change_residues= ['D', 'E']
        else:
            raise Exception('Can only specify charge_change as None, negative, or positive.')
        # find target residues
        target_residues=[]
        arg_pos=[]
        for ind, res in enumerate(input_sequence):
            if res in change_residues:
                target_residues.append(ind)
            if res == 'R':
                arg_pos.append(ind)

        # shuffle the list
        random.shuffle(target_residues)
        # add on R last
        if charge_change == None or charge_change=='positive':
            target_residues.extend(arg_pos)

    elif FCR_change == 'increase':
        # order of residues to change to charged based on how much it would fundamentally change the seq.
        # this list could be changed TBH, but this is my best guess for now.
        to_charged = ['S', 'T', 'N', 'Q', 'G', 'A', 'L', 'V', 'I', 'M', 'P', 'C', 'Y', 'W', 'F', 'H']
        target_residues=[]
        for target_res in to_charged:
            if target_residues==[]:
                for ind, res in enumerate(input_sequence):
                    if res == target_res:
                        target_residues.append(ind)
        random.shuffle(target_residues)

    else:
        raise Exception('Parameter FCR_change can only be decrease or increase.')
        
    # return the first element in the target residues list
    if target_residues != []:
        return target_residues[0]
    else:
        raise Exception('Was not able to find a residue to change using the identify_best_FCR_change function.')



def identify_best_NCPR_change(input_sequence, NCPR_change):
    '''
    Function to identify the best position for changing NCPR.
    Basically tries to just change lysine before arg.

    Parameters
    ----------
    input_sequence: string
        The sequence to be changed

    NCPR_change: string
        The desired change in NCPR. Specify increase or decrease.
        options: 'increase' or 'decrease'

    Returns
    -------
    res_index : int
        returns the index for the residue to change
    '''
    # TO DO: Add in some checks here.
    if NCPR_change == 'increase':
        targets = ['D', 'E']
    elif NCPR_change == 'decrease':
        targets = ['K']
    else:
        raise Exception('Parameter NCPR_change can only be decrease or increase.')

    # find target residues
    target_residues=[]
    for ind, res in enumerate(input_sequence):
        if res in targets:
            target_residues.append(ind)
    # shuffle the list
    random.shuffle(target_residues)
    if NCPR_change == 'decrease':
        for ind, res in enumerate(input_sequence):
            if res == 'R':
                target_residues.append(ind)
    # return the first element in the target residues list
    if target_residues != []:
        return target_residues[0]
    else:
        raise Exception('Was not able to find a residue to change using the identify_best_NCPR_change function.')





"""
# dumb function to test the function above, will delete later.
def change_hydro(input_sequence, objective_hydro):
    cur_hydro=3
    objective_hydro=4
    total_change=(objective_hydro-cur_hydro)*len(input_sequence)
    for i in range(0, 1000):
        curchange = identify_best_hydro_change(input_sequence, total_change)
        newseq = f'{input_sequence[:curchange[0]]}{curchange[1]}{input_sequence[curchange[0]+1:]}'
        input_sequence = newseq
        cur_hydro = pr(input_sequence).hydropathy
        total_change=(objective_hydro-cur_hydro)*len(input_sequence)
        print(cur_hydro)
        print(input_sequence)
        print()

"""


          







