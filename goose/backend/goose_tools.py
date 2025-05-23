'''
various essential tools to help GOOSE take flight
'''
import random
import csv

from goose.goose_exceptions import GooseError, GooseInputError
from goose.backend import parameters
from goose.backend.sequence_generation_backend import calculate_max_charge


def check_valid_kwargs(kwargs_dict, valid_keywords):
    """
    Function which takes a kwargs dict and a list of allowed keywords and
    ensures no non-allowed keywords have been passed. This should be run as
    the first function whereever kwargs** is invoked to ensure no silent errors
    are introduced (where invalid keywords are passed but then ignored).
    
    Parameters
    ------------------
    kwargs_dict : dict
        Dictionary of keyword arguments which maps keywords (keys) to passed
        values (values).

    valid_keywords : list
        A list of keywords this function is allowed to accept. This does not
        throw an error if a keyword is missing, but does throw an error if 
        an unexpected keyword is passed.

    Returns
    -----------------
    None
        This function does not return anything nor alter the state of any
        variables.
    

    Raises
    -----------------
    GooseInputError
        If an invalid keyword is identified, this function raises a 
        GooseInputError exception.

    """

    
    for kw in kwargs_dict:
        if kw not in valid_keywords:
            raise GooseInputError(f'Error: Invalid keyword [{kw}] was passed')
    

def length_check(input_val, length_max=parameters.MAXIMUM_LENGTH, length_min=parameters.MINIMUM_LENGTH):
    '''
    function to make sure length parameters don't go out of bounds
    '''
    if type(input_val) == str:
        length = len(input_val)
    elif type(input_val) == int:
        length = input_val
    else:
        raise GooseError('length check function in goose tools got an input value not string or int.')

    if length > length_max or length < length_min:
        error_message = f'\nLength of {length} is not within the allowed range for length of between {parameters.MINIMUM_LENGTH} and {parameters.MAXIMUM_LENGTH}'
        raise GooseInputError(error_message)


def write_csv(input_list, output_file, properties):
    """
    Function that writes the scores in an input dictionary out to a standardized CVS file format.

    Parameters
    -----------
    input_list : List
        List of dictionaries

    output_file : str
        Location and filename for the output file. Assumes .csv is provided.

    properties : Bool
        Whether or not to save the sequence properties. 

        If True
            save an autogenerated sequence name, length, the FCR, NCPR, hydropathy,
            sigma, delta values, fractions of the sequence followed by the sequence.

        If false
            just saves the sequence

    Returns
    --------
    None
        No return value, but writes a .csv file to disk

    """

    header_names=[]

    if properties == True:
        for header, values in input_list[0].items():
            header_names.append(header)
    
    else:
        # if properties was false, input_list is not a list of dicts
        # but rather a list of strings (the seqeunces). Need to make
        # into a list of dicts and give each sequence a name!
        amino_acids = ['A', 'C', 'D', 'E', 'F', 'G', 'H', 'I', 'K', 'L', 'M', 'N', 'P', 'Q', 'R', 'S', 'T', 'V', 'W', 'Y']
        final_list = []
        header_names = ['Name', 'sequence']
        for sequence in input_list:
            temp_dict = {}
            seq_name = '>'        
            for i in range(0, 5):
                seq_name += amino_acids[random.randint(0, len(amino_acids)-1)]
                seq_name += str(random.randint(0, 9))
            temp_dict['Name'] = seq_name
            temp_dict['sequence'] = sequence
            final_list.append(temp_dict)

    # try and open the file and throw exception if anything goes wrong    
    try:    
        with open(output_file, 'w', newline='') as csvfile:
            writer=csv.DictWriter(csvfile, fieldnames=header_names)
            writer.writeheader()
            if properties == True:
                for i in input_list:
                    writer.writerow(i)

            else:
                for i in final_list:
                    writer.writerow(i)

    except Exception:
            raise GooseError('Unable to write to file destination %s' % (output_file))                  







def remove_None(**kwargs):
    '''
    function to remove none values from kwargs to make the CLI work
    make a list to keep track of kwargs to remove because can't remove
    during iteration as it changes the dictionary length
    '''

    remove_from_kwargs = []
    for i in kwargs.keys():
        if (kwargs[i]) == None:
           remove_from_kwargs.append(i)

    # for everything identified as None, remove it form kwargs
    for i in remove_from_kwargs:
        del kwargs[i]

    return kwargs


def check_props_parameters(**kwargs):
    '''
    function that makes sure the values for a given seq are not 
    outside possible bounds. Only for generating
    sequences by specifying properties.
    '''

    # check disorder disorder cutoff value bounds
    if 'cutoff' in list(kwargs.keys()):
        curval = kwargs['cutoff']
        if curval != None:
            if curval > parameters.MAXIMUM_DISORDER:
                error_message = f'The disorder cutoff value {curval} is greater than the max allowed value of {parameters.MAXIMUM_DISORDER}'
                print(error_message)
                raise GooseInputError(error_message)
            if curval < parameters.MINIMUM_DISORDER:
                error_message = f'The disorder cutoff {curval} is less than the minimum allowed value of {parameters.MINIMUM_DISORDER}'
                raise GooseInputError(error_message)

    # check FCR bounds
    if kwargs['FCR'] != None:
        curval = kwargs['FCR']
        if curval > parameters.MAXIMUM_FCR:
            error_message = f'The FCR value {curval} is greater than the allowed value of {parameters.MAXIMUM_FCR}'
            raise GooseInputError(error_message)
        if curval < parameters.MINIMUM_FCR:
            error_message = f'The FCR value {curval} is less than the allowed value of {parameters.MINIMUM_FCR}'
            raise GooseInputError(error_message)

    # check NCPR bounds
    if kwargs['NCPR'] != None:
        curval = kwargs['NCPR']
        if curval > parameters.MAXIMUM_NCPR:
            error_message = f'The NCPR value {curval} is greater than the allowed value of {parameters.MAXIMUM_NCPR}'
            raise GooseInputError(error_message)
        if curval < parameters.MINIMUM_NCPR:
            error_message = f'The NCPR value {curval} is less than the allowed value of {parameters.MINIMUM_NCPR}'
            raise GooseInputError(error_message)

    # check hydropathy bounds
    if kwargs['hydropathy'] != None:
        curval = kwargs['hydropathy']
        if curval > parameters.MAXIMUM_HYDRO:
            error_message = f'The hydropathy value {curval} is greater than the allowed value of {parameters.MAXIMUM_HYDRO}'
            raise GooseInputError(error_message)
        if curval < parameters.MINIMUM_HYDRO:
            error_message = f'The hydropathy value {curval} is less than the allowed value of {parameters.MINIMUM_HYDRO}'
            raise GooseInputError(error_message)


    # check FCR and NPCR error
    if kwargs['NCPR'] != None and kwargs['FCR'] != None:
        if abs(kwargs['NCPR']) > kwargs['FCR']:
            valncpr = kwargs['NCPR']
            valfcr = kwargs['FCR']
            error_message = f'The NCPR value of {valncpr} is not possible given the FCR value of {valfcr}'
            raise GooseInputError(error_message)

    # check kappa
    if kwargs['kappa'] != None:
        val_kappa = kwargs['kappa']
        if kwargs['kappa'] > 1 or kwargs['kappa'] < 0:
            error_message = f'The kappa value of {val_kappa} is not possible. Values are between 0 and 1'
            raise GooseInputError(error_message)        

    # check kappa based on FCR stuff
    if kwargs['kappa'] != None:
        if kwargs['FCR'] == 0:
            raise GooseInputError('Cannot specifiy kappa values when FCR is 0.')
        if kwargs['FCR'] != None:
            if kwargs['NCPR'] != None:
                if kwargs['NCPR'] == kwargs['FCR']:
                    raise GooseInputError('Cannot have NCPR = FCR for kappa to be a value other than -1.')


    # now check for charge AND hydro
    if kwargs['NCPR'] != None and kwargs['hydropathy'] != None:
        hydro_and_charge = True
    elif kwargs['FCR'] != None and kwargs['hydropathy'] != None:
        hydro_and_charge = True
    else:
        hydro_and_charge = False

    # get the charge value to use as the val for
    # checking hydro + charge errors
    if hydro_and_charge == True:
        if kwargs['FCR'] == None:
            check_charge_val = abs(kwargs['NCPR'])
        else:
            check_charge_val = kwargs['FCR']

        # now calculate maximum charge value
        max_possible_charge_value = calculate_max_charge(kwargs['hydropathy'])
        # now see if charge value used is within bounds
        if check_charge_val > max_possible_charge_value:
            curhydro = kwargs['hydropathy']
            error_message = f'The specified hydropathy value if {curhydro} is not compatible with the specified charge value.'
            raise GooseInputError(error_message)



def check_and_correct_props_kwargs(**kwargs):
    '''
    Checks the kwargs for creating a sequence by
    properties.
    
    Parameters
    -----------
    **kwargs : Dict
        The input dict to evaluate for completion and errors

    Returns
    -------
    Dict
        A dictionary holding all necessary
        values for sequence generation.
    '''
    # first remove any 'None' arguments passed from the command line
    # from the kwargs dict.
    kwargs = remove_None(**kwargs)

    # make sure the six essential keywords have been initialized to their default values if they were not provided
    essential_kwargs = {'cutoff': parameters.DISORDER_THRESHOLD,'FCR': None, 'NCPR': None, 'hydropathy': None, 'kappa': None, 
                        'attempts':parameters.DEFAULT_ATTEMPTS, 'exclude':None, 'strict_disorder': False,
                        'metapredict_version': parameters.METAPREDICT_DEFAULT_VERSION, 'return_all_sequences': False, 
                        'use_weighted_probabilities': False, 'custom_probabilities':None}
    
    for kw in essential_kwargs:
        if kw not in kwargs:
            kwargs[kw] = essential_kwargs[kw]
    
    # return the corrected dict
    return kwargs


def check_and_correct_fracs_kwargs(**kwargs):
    '''
    Checks the kwargs for creating a sequence by
    fractions.
    
    Parameters
    -----------
    **kwargs : Dict
        The input dict to evaluate for completion and errors

    Returns
    -------
    Dict
        A dictionary holding all necessary (and corrected) 
        values for sequence generation.
    '''
    # first remove any 'None' arguments passed from the command line
    # from the kwargs dict.
    kwargs = remove_None(**kwargs)
    

    # make sure the six essential keywords have been initialized to their default values if they were not provided
    essential_kwargs = {'cutoff': parameters.DISORDER_THRESHOLD, 
                        'strict_disorder': False, 
                        'attempts': parameters.DEFAULT_ATTEMPTS, 
                        'max_aa_fractions': {},
                        'metapredict_version': parameters.METAPREDICT_DEFAULT_VERSION,
                        'return_all_sequences': False,
                        'use_weighted_probabilities': False,
                        'custom_probabilities': None,
                        'exclude': None}

    for kw in essential_kwargs:
        if kw not in kwargs:
            kwargs[kw] = essential_kwargs[kw]

    # return the corrected dict
    return kwargs



def check_fracs_parameters(**kwargs):
    '''
    function that makes sure the values for a given seq are not 
    outside possible bounds. Only for generating
    sequences by specifying fractions of amino acids.

    Parameters
    -------------------
    TBD

    Returns
    -------------------
    None
        No return variable

    Raises 
    --------------------
    GooseInputError  
        Returns a Goose input error in a passed fracion is invalid.

    '''

    # check disorder disorder cutoff value bounds (max and min are 1 and 0)
    if 'cutoff' in kwargs:
        curval = kwargs['cutoff']
        if curval != None:
            if curval > parameters.MAXIMUM_DISORDER:
                error_message = f'The disorder cutoff value {curval} is greater than the max allowed value of {parameters.MAXIMUM_DISORDER}'
                print(error_message)
                raise GooseInputError(error_message)
            if curval < parameters.MINIMUM_DISORDER:
                error_message = f'The disorder cutoff {curval} is less than the minimum allowed value of {parameters.MINIMUM_DISORDER}'
                raise GooseInputError(error_message)

    amino_acids = ['A', 'C', 'D', 'E', 'F', 'G', 'H', 'I', 'K', 'L', 'M', 'N', 'P', 'Q', 'R', 'S', 'T', 'V', 'W', 'Y']

    # cycle over each keyword (kw) in the kwargs dictionary
    for kw in kwargs:

        # if the keyword is an amino acid...
        if kw in amino_acids:
            
            # get current fraction value
            current_value = kwargs[kw]

            # if a an override max value was provided
            if kw in kwargs['max_aa_fractions']:
                max_value = kwargs['max_aa_fractions'][kw]

            # otherwise use the default
            else:                
                max_value = parameters.MAX_FRACTION_DICT[kw]
            # make sure that current value not greater than max value
            
            if current_value > max_value:
                error_message = f'Current value of {current_value} is greater than max possible fraction value for amino acid {kw}. Max value is {max_value}.'
                raise GooseInputError(error_message)

    # make sure these bounds make sense...
    for aa in kwargs['max_aa_fractions']:
        if kwargs['max_aa_fractions'][aa] < 0 or kwargs['max_aa_fractions'][aa] > 1:
            raise GooseInputError(f"Max aa fraction requested for {aa} is invalid ({kwargs['max_aa_fractions'][aa]}). Must be between 0 and 1")
            

    # finally - check the requested fractions don't sum to more than 1 and all values are valid
    summed_fraction = 0
    for aa in amino_acids:
        
        if aa in kwargs:
            if kwargs[aa] < 0 or kwargs[aa] > 1:
                raise GooseInputError(f'Fractional content requested for {aa} is invalid ({kwargs[aa]}). Must be between 0 and 1')
            
            summed_fraction = summed_fraction + kwargs[aa]

    # round to 10th fractional value to avoid dumb errors
    if round(summed_fraction,10) > 1:
        raise GooseInputError(f'Requested a sequence where the sum of the fractional components is greater than 1. This will not work!')
        



def gen_random_name():
    '''
    generates a random name for a protein.
    '''
    random_name = '>'
    amino_acids = ['A', 'C', 'D', 'E', 'F', 'G', 'H', 'I', 'K', 'L', 'M', 'N', 'P', 'Q', 'R', 'S', 'T', 'V', 'W', 'Y']
    for i in range(0, 5):
        random_name += amino_acids[random.randint(0, len(amino_acids)-1)]
        random_name += str(random.randint(0, 9))
    return random_name    






