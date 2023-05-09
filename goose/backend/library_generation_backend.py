'''
backend code for library generation
'''
import numpy as np

# Uncomment later
from goose.backend import parameters
from goose import goose_exceptions
from goose.backend.sequence_generation_backend import calculate_max_charge
import itertools


def generate_library_by_parameter_ranges(length,
    FCR=None, NCPR=None, hydropathy=None, kappa=None, silent_failed_seqs=False):
    '''
    Function that returns a list of dictioanries where each
    dictionary has the values specified for a library of sequences.
    Bit of a long function because all of the paramter values are checked
    before adding the sequence to the library. This way the user can specify 
    ranges of, for example, hydropathy, that are not possible with certain charge
    values but are possible with other charge values. 


    parameters
    ----------
    length : Int 
        The length of the sequence

    FCR : list of float(s) or float
        The fraction of charged residues as a decimal, between 0 and 1.
        Input as a range where the first value is the lower bounds wanted
        for the library and the second value is the upper bounds wanted for
        the library. If a single number is given, then that value will be constant
        across the library.
        A third value can be specified, which will determine the
        interval between the lowest value and the highest value. 
        If this interval is not perfect, then GOOSE will take the 
        closest value to the desired maximum. 

    NCPR : list of float(s) or float
        The net charge of the sequence as a decimal, between -1 and 1
        Absolute value cannot be greater than FCR.
        Input as a range where the first value is the lower bounds wanted
        for the library and the second value is the upper bounds wanted for
        the library. If a single number is given, then that value will be constant
        across the library.
        A third value can be specified, which will determine the
        interval between the lowest value and the highest value. 
        If this interval is not perfect, then GOOSE will take the 
        closest value to the desired maximum. 

    hydropathy : list of float(s) or float
        The mean hydropathy of the sequence. Lower numbers are less
        hydrophobic. Between 0.6 and 6.1
        Input as a range where the first value is the lower bounds wanted
        for the library and the second value is the upper bounds wanted for
        the library. If a single number is given, then that value will be constant
        across the library.
        A third value can be specified, which will determine the
        interval between the lowest value and the highest value. 
        If this interval is not perfect, then GOOSE will take the 
        closest value to the desired maximum.         

    kappa : list of float(s) or float
        The charge asymmetry metric. Describes how oppositely charged
        residues are patterned across the sequeence. Between 0 and 1 
        where a higher value is a more asymmetrically distributed 
        positioning of oppositely charged residues.
        Input as a range where the first value is the lower bounds wanted
        for the library and the second value is the upper bounds wanted for
        the library. If a single number is given, then that value will be constant
        across the library.        
        A third value can be specified, which will determine the
        interval between the lowest value and the highest value. 
        If this interval is not perfect, then GOOSE will take the 
        closest value to the desired maximum. 

        Input as a range where the first value is the lower bounds wanted
        for the library and the second value is the upper bounds wanted for
        the library. If a single number is given, then that value will be constant
        across the library.

    silent_failed_seqs : bool
        Whether to silence any printed warnings of sequences that
        are not possible to generate due to incompatible charge/hydorpatyh values

    returns
    -------
    sequence_list : list of dictionaries
        Returns a list of dicts where each dictionary in the list 
        can be input into the sequence generators downstream. 
        The dict holds the parameter specifications for each sequence in the library.
    '''
    # empty dictionary to hold parameter values
    parameter_dictionary={}

    # go through and get the lists of values for all specified parameters
    # main reason for this is to check the values for input errors. 
    params_possible = [FCR, NCPR, hydropathy, kappa]
    param_names = ['FCR', 'NCPR', 'hydropathy', 'kappa']
    for param_ind in range(0, len(params_possible)):
        param = params_possible[param_ind]
        param_name = param_names[param_ind]
        if param != None:
            if type(param)==float or type(param)==int:
                parameter_dictionary[param_name]=[param]
            else:
                if type(param) == list:
                    if len(param) == 1:
                        parameter_dictionary[param_name]= param
                    elif len(param)==2:
                        parameter_dictionary[param_name] = param
                    elif len(param)==3:
                        # make sure the specified interval value is possible.
                        if param[2]> max(param):
                            message = f'\n\tThe third number specified in the list for {param_name} of {param[2]} is\n\t greater than the range specified of {param[0]} and {param[1]}.\n\t Specify a third value in the list between the parameter\n\t value range of {param[0]} and {param[1]}\n'
                            raise Exception(message)
                        if param[2]<=0:
                            message = f'\n\tThe third number specified in the list for {param_name} of {param[2]} is\n\t less than 0. You must specify a value greater than 0 between {param[0]} and {param[1]}.\n'
                        # now that that's taken care of, make sure all values in list are int or float
                        for val in param:
                            # make sure values in the parameters are integers or floats.
                            if type(val) != int:
                                if type(val) != float:
                                    message = f'\n\tThe input type for {param_name} for value {val}\n\tis not an integer or decimal value!'
                                    raise Exception(message)
                        # finally can get on with it.
                        interval = param[2]
                        values_list = list(np.arange(param[0], param[1], interval))
                        # add max value
                        possible_max_value=param[1]+interval
                        if possible_max_value > param[1]:
                            values_list.append(param[1])
                        else:
                            values_list.append(possible_max_value)
                        # add the values to the parameter dictionary
                        parameter_dictionary[param_name] = values_list
                        # get the minimum number of sequences

                    else:
                        # The range / range+interval list cannot be over 3 values.
                        message = f'\n\tThe number of inputs for {param_name} cannot be over 3!\n\tPlease specify 3 or fewer numbers for each parameter.\n'
                        raise Exception(message)
                else:
                    # if something other than float, int, or list is input, ...
                    message = f'The input for parameter {param_name} has an invalid input data type of {type(param)}!'
                    raise Exception(message)
        else:
            parameter_dictionary[param_name]=[None]
    # now have the parameter_dictionary populated with values specified by
    # the user. Now look at values specified and add possible values to the final list.
    list_of_seqs = []
    # list to hold sequences GOOSE could not make due to hydropathy value
    # incompatiblity
    failed_seqs=[]
    for cur_fcr_val in parameter_dictionary['FCR']:
        fcr_ok=True
        if cur_fcr_val != None:
            if cur_fcr_val < 0 or cur_fcr_val > 1:
                fcr_ok=False
        for cur_ncpr_val in parameter_dictionary['NCPR']:
            ncpr_ok=True
            if cur_ncpr_val != None:
                if cur_ncpr_val < -1 or cur_ncpr_val > 1:
                    ncpr_ok=False
                if cur_fcr_val != None:
                    if abs(cur_ncpr_val)>cur_fcr_val:
                        ncpr_ok = False
            for cur_kappa_val in parameter_dictionary['kappa']:
                kappa_ok=True
                if cur_kappa_val != None:
                    if cur_kappa_val < 0 or cur_kappa_val > 1:
                        kappa_ok = False
                    if cur_fcr_val != None:
                        if cur_fcr_val == 0:
                            kappa_ok = False
                    if cur_ncpr_val != None:
                        if abs(cur_ncpr_val) ==1:
                            kappa_ok = False    
                        if cur_fcr_val != None:
                            if cur_fcr_val == abs(cur_ncpr_val):
                                kappa_ok = False
                for cur_hydro_val in parameter_dictionary['hydropathy']:
                    hydro_ok = True
                    if cur_hydro_val != None:
                        if cur_hydro_val > parameters.MAXIMUM_HYDRO or cur_hydro_val < parameters.MINIMUM_HYDRO:
                            hydro_ok=False
                        max_charge = calculate_max_charge(cur_hydro_val)
                        if cur_fcr_val != None:
                            if cur_fcr_val > max_charge:
                                hydro_ok=False
                                failed_seqs.append({'FCR':cur_fcr_val, 'NCPR':cur_ncpr_val, 'kappa':cur_kappa_val, 'hydropathy': cur_hydro_val})
                        if cur_ncpr_val != None:
                            if abs(cur_ncpr_val)>max_charge:
                                hydro_ok=False
                                failed_seqs.append({'FCR':cur_fcr_val, 'NCPR':cur_ncpr_val, 'kappa':cur_kappa_val, 'hydropathy': cur_hydro_val})
                    # make sure none of the paramters had any issues...
                    if fcr_ok==True:
                        if ncpr_ok==True:
                            if kappa_ok==True:
                                if hydro_ok==True:
                                    list_of_seqs.append({'FCR':cur_fcr_val, 'NCPR':cur_ncpr_val, 'kappa':cur_kappa_val, 'hydropathy': cur_hydro_val})
    if failed_seqs != []:
        if silent_failed_seqs==False:
            print('The following sequences were not made because the hydropathy and charge values were not possible:')
            for seq in failed_seqs:
                print(seq)
            print('End failed seqs.\n\n\n')

    return list_of_seqs


# to get the combos of frac values from the library fractions function
def dict_combinations(input_dict):
    """
    Takes in a dict and returns a list of dicts with each possible value from that dict.
    """
    keys = list(input_dict.keys())
    values = list(input_dict.values())
    combinations = list(itertools.product(*values))
    output_vals = []
    for combo in combinations:
        temp_dict={}
        for vals in range(0, len(combo)):
            temp_dict[keys[vals]]=combo[vals]
        output_vals.append(temp_dict)
    return output_vals

def generate_library_by_fraction_ranges(length, **kwargs):
    '''
    Function that returns a list of dictioanries where each
    dictionary has the values specified for a library of sequences.
    Bit of a long function because all of the paramter values are checked
    before adding the sequence to the library. This way the user can specify 
    ranges of, for example, aliphatics that add up to an impossible value. 

    Parameters
    ------------
    length : int
        the length of the sequence to generate

    <each of the 20 amino acids> : float
        Specify the fraction of the sequence that should be made up of one or more
        of the 20 natural amino acids (e.g. A=0.2, Y=0.05) etc.

    max_aa_fractions : dict 
        Dictionary which, if provided, allows the user to over-ride the 
        fraction of a sequence which can be made up of any given amino
        acid. The passed dictionary should contain key/value pairs, where
        keys are one of the twenty standard amino acids and values is a
        float between 0 and 1. If amino acids are missing then the default
        thresholds set by GOOSE are used.

    Returns 
    --------
    length, final_seq_specifications : list
        returns a list where the first element is the lenght of the sequence
        and the second element is the     

    '''
    # empty dictionary to hold parameter values
    parameter_dictionary={}

    # go through and get the lists of values for all specified parameters
    # main reason for this is to check the values for input errors. 
    # set all to None by default before checking
    init_dict = {}
    param_names=['A', 'C', 'D', 'E', 'F', 'G', 'H', 'I', 'K', 'L', 'M', 'N', 'P', 'Q', 'R', 'S', 'T', 'V', 'W', 'Y']
    for val in param_names:
        if val not in kwargs.keys():
            init_dict[val]=None
        else:
            init_dict[val]=kwargs[val]

    for param_ind in init_dict:
        param = init_dict[param_ind]
        param_name=(str(param_ind))
        if param != None:
            if type(param)==float or type(param)==int:
                parameter_dictionary[param_name]=[param]
            else:
                if type(param) == list:
                    if len(param) == 1:
                        parameter_dictionary[param_name]= param
                    elif len(param)==2:
                        parameter_dictionary[param_name] = param
                    elif len(param)==3:
                        # make sure the specified interval value is possible.
                        if param[2]> max(param):
                            message = f'\n\tThe third number specified in the list for {param_name} of {param[2]} is\n\t greater than the range specified of {param[0]} and {param[1]}.\n\t Specify a third value in the list between the parameter\n\t value range of {param[0]} and {param[1]}\n'
                            raise Exception(message)
                        if param[2]<=0:
                            message = f'\n\tThe third number specified in the list for {param_name} of {param[2]} is\n\t less than 0. You must specify a value greater than 0 between {param[0]} and {param[1]}.\n'
                        # now that that's taken care of, make sure all values in list are int or float
                        for val in param:
                            # make sure values in the parameters are integers or floats.
                            if type(val) != int:
                                if type(val) != float:
                                    message = f'\n\tThe input type for {param_name} for value {val}\n\tis not an integer or decimal value!'
                                    raise Exception(message)
                        # finally can get on with it.
                        interval = param[2]
                        values_list = list(np.arange(param[0], param[1], interval))
                        # add max value
                        possible_max_value=param[1]+interval
                        if possible_max_value > param[1]:
                            values_list.append(param[1])
                        else:
                            values_list.append(possible_max_value)
                        # add the values to the parameter dictionary
                        parameter_dictionary[param_name] = values_list
                        # get the minimum number of sequences

                    else:
                        # The range / range+interval list cannot be over 3 values.
                        message = f'\n\tThe number of inputs for {param_name} cannot be over 3!\n\tPlease specify 3 or fewer numbers for each parameter.\n'
                        raise Exception(message)
                else:
                    # if something other than float, int, or list is input, ...
                    message = f'The input for parameter {param_name} has an invalid input data type of {type(param)}!'
                    raise Exception(message)

    # now have the parameter_dictionary populated with values specified by
    # the user. Now get all possible combinations.
    all_combos = dict_combinations(parameter_dictionary)

    # now need to verify that everything in range of allowable fractions.
    # first see if any max fractions are overridden:
    if 'max_aa_fractions' in kwargs.keys():
        max_fracs = {}
        for val in parameters.MAX_FRACTION_DICT.keys():
            if val in kwargs['max_aa_fractions'].keys():
                max_fracs[val]=kwargs['max_aa_fractions'][val]
            else:
                max_fracs[val]=parameters.MAX_FRACTION_DICT[val]
    else:
        max_fracs = parameters.MAX_FRACTION_DICT
    
    # now go through each value in all_combos and make sure total value of specified fracs < 1 and 
    # that each amino acid doesn't have a greater evalue than specified max value.
    final_seq_specifications=[]
    for com in all_combos:
        # first make sure that total not over 1. 
        if sum(com.values())<= 1:
            # now go through each aa
            # first set add to true
            add_to_list=True
            for aa in com.keys():
                if com[aa] > max_fracs[aa]:
                    add_to_list=False
            if add_to_list==True:
                final_seq_specifications.append(com)

    return length, final_seq_specifications









