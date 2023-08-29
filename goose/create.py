##
## create.py
## 
## create.py contains all the user-facing function associated with metapredict. If a new function is added it should be included
## here and added to the __all__ list
## 

##Handles the primary functions

# if any new functions are added to create.py, you need to add them here.
__all__ =  ['seq_fractions', 'sequence', 'minimal_var', 'new_seq_constant_class_var', 'new_var', 'constant_class_var', 'hydro_class_var', 'constant_residue_var', 'shuffle_var', 'kappa_var', 'asymmetry_var', 'fcr_class_var', 'ncpr_class_var', 'all_props_class_var', 'alpha_helix', 'beta_strand', 'beta_sheet', 'seq_property_library']

import os
import sys
import random

# note - we import packages below with a leading _ which means they are ignored in the import

#for sequence generation
from goose.backend.sequence_generation import generate_disordered_seq_by_fractions as _generate_disordered_seq_by_fractions
from goose.backend.sequence_generation import generate_disordered_seq_by_props as _generate_disordered_seq_by_props

# goose tools for checking and fixing parameters
from goose.backend.goose_tools import check_and_correct_props_kwargs as _check_and_correct_props_kwargs
from goose.backend.goose_tools import check_props_parameters as _check_props_parameters
from goose.backend.goose_tools import check_and_correct_fracs_kwargs as _check_and_correct_fracs_kwargs
from goose.backend.goose_tools import check_fracs_parameters as _check_fracs_parameters
from goose.backend.goose_tools import length_check as _length_check
from goose.backend.goose_tools import gen_random_name as _gen_random_name
from goose.backend.goose_tools import check_valid_kwargs as _check_valid_kwargs

# variant generation stuff
from goose.backend.variant_generation import gen_kappa_variant as _gen_kappa_variant
from goose.backend.variant_generation import gen_shuffle_variant as _gen_shuffle_variant
from goose.backend.variant_generation import gen_constant_residue_variant as _gen_constant_residue_variant
from goose.backend.variant_generation import gen_hydropathy_class_variant as _gen_hydropathy_class_variant
from goose.backend.variant_generation import gen_new_variant as _gen_new_variant
from goose.backend.variant_generation import gen_constant_class_variant as _gen_constant_class_variant
from goose.backend.variant_generation import gen_new_var_constant_class as _gen_new_var_constant_class
from goose.backend.variant_generation import gen_asymmetry_variant as _gen_asymmetry_variant
from goose.backend.variant_generation import gen_fcr_class_variant as _gen_fcr_class_variant
from goose.backend.variant_generation import gen_ncpr_class_variant as _gen_ncpr_class_variant
from goose.backend.variant_generation import gen_all_props_class_variant as _gen_all_props_class_variant
from goose.backend.variant_generation import gen_targeted_shuffle_variant as _gen_targeted_shuffle_variant
from goose.backend.variant_generation import gen_excluded_shuffle_variant as _gen_excluded_shuffle_variant

# for minimal variant, need to... get rid of this old code.
from goose.backend.gen_minimal_variant_backend import gen_minimal_sequence_variant as _gen_minimal_sequence_variant

# library creation
from goose.backend.library_generation_backend import generate_library_by_parameter_ranges as _generate_library_by_parameter_ranges
from goose.backend.library_generation_backend import generate_library_by_fraction_ranges as _generate_library_by_fraction_ranges

# for folded structure generation
from goose.backend.folded_region_generation import gen_helix as _gen_helix
from goose.backend.folded_region_generation import gen_beta_strand as _gen_beta_strand
from goose.backend.folded_region_generation import gen_beta_sheet as _gen_beta_sheet

# FOR WHEN THINGS GO WRONG
from goose import goose_exceptions

# because every good package needs some parameters
from goose.backend import parameters

#-\|/-\|/-\|/-\|/-\|/-\|/-\|/-\|/-\|/-\|/-\|/-\|/-\|/-\|/-\|/-\|/-\|/-\|/-\|/-\|/-\|/-
#-\|/-\|/-\|/-\|/-\|/-\|/-\|/-\|/-\|/             \|/-\|/-\|/-\|/-\|/-\|/-\|/-\|/-\|/-
#-\|/-\|/-\|/-\|/-\|/-\|/-\|/-\|/-\|/  Create     \|/-\|/-\|/-\|/-\|/-\|/-\|/-\|/-\|/-
#-\|/-\|/-\|/-\|/-\|/-\|/-\|/-\|/-\|/  sequence   \|/-\|/-\|/-\|/-\|/-\|/-\|/-\|/-\|/-
#-\|/-\|/-\|/-\|/-\|/-\|/-\|/-\|/-\|/  By         \|/-\|/-\|/-\|/-\|/-\|/-\|/-\|/-\|/-
#-\|/-\|/-\|/-\|/-\|/-\|/-\|/-\|/-\|/  Specifying \|/-\|/-\|/-\|/-\|/-\|/-\|/-\|/-\|/-
#-\|/-\|/-\|/-\|/-\|/-\|/-\|/-\|/-\|/  Properties \|/-\|/-\|/-\|/-\|/-\|/-\|/-\|/-\|/-
#-\|/-\|/-\|/-\|/-\|/-\|/-\|/-\|/-\|/             \|/-\|/-\|/-\|/-\|/-\|/-\|/-\|/-\|/-
#-\|/-\|/-\|/-\|/-\|/-\|/-\|/-\|/-\|/-\|/-\|/-\|/-\|/-\|/-\|/-\|/-\|/-\|/-\|/-\|/-\|/-

def sequence(length, **kwargs):
    """
    Stand alone function that takes care of creating sequences with specified
    properties. Takes in various sequence parameters as arguments. 

    NB:

    1. You can specify NCPR and FCR simultaneously
    2. You can specify NCPR, FCR, and hydropathy simultaneously
    3. If you specify sigma you cannot specify NCPR, FCR, or hydropathy

    Parameters
    ------------
    length : int
        Defines the length of desired disordered sequence

    FCR : float
        Defines the requested fraction of charged residues (between 0 and 1)

    NCPR : float 
        Defines the net charge per residue of generated sequences (between -1 and 1)

    sigma : float
        Defines the sigma value of generated sequences(between 0 and 1). Sigma reports
        on the charge asymmetry (between 0 and 1).

    hydropathy : float 
        Defines the mean hydropathy of generated sequence (between 0 and 6.1).

    kappa : float
        specify the kappa value of generated sequence. Kappa reports on the charge
        patterning (between 0 and 1 if there are both positive and negative residues).

    attempts : int
        Specify the number of times to make the sequence. Default is 20. Greater numbers
        of attempts increase the odds that a sequence will be generated but will increase
        the duration of attempting to make the sequence. 

    cutoff : float
        The disorder cutoff threshold. This ensures a sequence has a mean disorder above
        this cutoff

    attempts : int
        Number of attempts to try
            

    Returns
    -----------
    generated_seq : String
        Returns a string that is the amino acid sequence

    """
    # check length. Note this checks min/max length as well as
    # casts a string length to an int
    _length_check(length)

    # check we passed in acceptable keyword arguments. At this stage, if a keyword
    # was passed that is not found in the list passed to _check_valid_kwargs then
    # an exception is raised. 
    _check_valid_kwargs(kwargs, ['FCR','NCPR','sigma', 'hydropathy', 'kappa', 'cutoff', 'attempts'])
    

    # First correct kwargs. Do this first because
    # the next function that looks over kwargs values
    # can only take in corrected kwargs.
    kwargs = _check_and_correct_props_kwargs(**kwargs)

    # now make sure that the input vals are within appropriate bounds
    _check_props_parameters(**kwargs)

    # make the sequence
    try:
        generated_seq = _generate_disordered_seq_by_props(length, FCR=kwargs['FCR'], NCPR=kwargs['NCPR'], hydropathy=kwargs['hydropathy'],
            sigma = kwargs['sigma'], attempts = kwargs['attempts'], allowed_hydro_error = parameters.HYDRO_ERROR, disorder_threshold = kwargs['cutoff'])
    except:
        raise goose_exceptions.GooseFail('Unable to generate sequence. Please try again with different parameters or a lower cutoff value.')

    # this is a bit hacky for now, but it works.
    if kwargs['kappa'] != None:
        try:
            generated_seq = _gen_kappa_variant(generated_seq, kappa=kwargs['kappa'], allowed_kappa_error = parameters.MAXIMUM_KAPPA_ERROR, disorder_threshold=kwargs['cutoff'], strict_disorder=False)
        except:
            raise goose_exceptions.GooseFail('Unable to get kappa of generated sequence correct')
    # return the seq
    return generated_seq


#-\|/-\|/-\|/-\|/-\|/-\|/-\|/-\|/-\|/-\|/-\|/-\|/-\|/-\|/-\|/-\|/-\|/-\|/-\|/-\|/-\|/-
#-\|/-\|/-\|/-\|/-\|/-\|/-\|/-\|/-\|/             \|/-\|/-\|/-\|/-\|/-\|/-\|/-\|/-\|/-
#-\|/-\|/-\|/-\|/-\|/-\|/-\|/-\|/-\|/  Create     \|/-\|/-\|/-\|/-\|/-\|/-\|/-\|/-\|/-
#-\|/-\|/-\|/-\|/-\|/-\|/-\|/-\|/-\|/  sequence   \|/-\|/-\|/-\|/-\|/-\|/-\|/-\|/-\|/-
#-\|/-\|/-\|/-\|/-\|/-\|/-\|/-\|/-\|/  By         \|/-\|/-\|/-\|/-\|/-\|/-\|/-\|/-\|/-
#-\|/-\|/-\|/-\|/-\|/-\|/-\|/-\|/-\|/  Specifying \|/-\|/-\|/-\|/-\|/-\|/-\|/-\|/-\|/-
#-\|/-\|/-\|/-\|/-\|/-\|/-\|/-\|/-\|/  Fractions  \|/-\|/-\|/-\|/-\|/-\|/-\|/-\|/-\|/-
#-\|/-\|/-\|/-\|/-\|/-\|/-\|/-\|/-\|/             \|/-\|/-\|/-\|/-\|/-\|/-\|/-\|/-\|/-
#-\|/-\|/-\|/-\|/-\|/-\|/-\|/-\|/-\|/-\|/-\|/-\|/-\|/-\|/-\|/-\|/-\|/-\|/-\|/-\|/-\|/-


def seq_fractions(length, **kwargs):
    """
    Stand-alone function that takes care of creating sequences with specified
    fractions of amino acids. 

    Parameters
    ------------
    length : int
        length of the desired disordered sequence

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

    cutoff : float
        The disorder cutoff threshold. Default used is 0.6
    
    attempts : int
        Number of attempts that will be made to generate a given sequence. Default
        is 100.
    
    strict_disorder : bool
        if set to true, will not count a sequence as disordered even if a single amino
        acid falls below the cutoff value. Default = False.
             
    Returns
    -----------
    generated_seq : string
        Returns a string that is the amino acid sequence

    """

    # check length. Note this checks min/max length as well as
    # casts a string length to an int
    _length_check(length)

    # check we passed in acceptable keyword arguments. At this stage, if a keyword
    # was passed that is not found in the list passed to _check_valid_kwargs then
    # an exception is raised. 
    _check_valid_kwargs(kwargs, ['cutoff', 'attempts',  'strict_disorder',  'max_aa_fractions', 'A','C','D','E','F','G','H','I','K','L','M','N','P','Q','R','S','T','V','W','Y'])

    # First correct kwargs. Do this first because
    # the next function that looks over kwargs values
    # can only take in corrected kwargs.
    kwargs = _check_and_correct_fracs_kwargs(**kwargs)

    # now make sure that the input vals are within appropriate bounds
    _check_fracs_parameters(**kwargs)

    generated_seq = _generate_disordered_seq_by_fractions(length, **kwargs)

    # return the seq
    return generated_seq


'''
/-/-/-/-/-/-/-/-/-/-/-/-/-/-/-/-/-/-/
/-/-/-/-/-/-/-/-/-/-/-/-/-/-/-/-/-/-/
          VARIANT GENERATORS
/-/-/-/-/-/-/-/-/-/-/-/-/-/-/-/-/-/-/
/-/-/-/-/-/-/-/-/-/-/-/-/-/-/-/-/-/-/
'''



def minimal_var(input_sequence, hydropathy = '', fcr = '', 
    ncpr = '', scd='', cutoff=parameters.DISORDER_THRESHOLD, strict=False):
    '''
    User facing function for generating the minimal sequence variant. This variant
    tries to make a sequence as similar to the input sequence as possible all while
    minimizing the number of amino acids changed. 
    
    Parameters
    ------------
    input_sequence : str
        the sequence to make a variant of

    hydropathy : float
        the mean hydropathy of the sequence. If not specified does not change.

    fcr : float
        the fraction of charged residues in the sequence. If not specified does not change.

    ncpr : float
        the net charge per residue of the sequence. If not specified does not change.

    scd : float
        the charge asymmetry of the sequence. If not specified does not change.

    cutoff : float
        the disorder cutoff threshold. Closer to 1 has a higher chance of being disordered.
        Must be betwee 0 and 1. 

    strict : bool
        whether to use a strict disorder calculation. By default, variants are allowed
        to have regions below the disorder threshold *for regions where the input sequence
        is also below the threshold*. 

    Returns
    -----------
    final_sequence : str
        the final sequence variant
    '''

    # make sure that the input sequence is all caps
    input_sequence = input_sequence.upper()

    # check length
    _length_check(input_sequence)

    if cutoff > 1 or cutoff < 0:
        raise goose_exceptions.GooseInputError('cutoff value must be between 0 and 1 for disorder threshold')    
    try:
        final_sequence = _gen_minimal_sequence_variant(input_sequence, mean_hydro = hydropathy, fraction = fcr, 
        net_charge = ncpr, charge_asymmetry=scd, cutoff=cutoff, strict=strict)
    except:
        raise goose_exceptions.GooseFail('Sorry! GOOSE was unable to generate the sequence. Please try again or try with different parameters or a lower cutoff value.')
    return final_sequence



def new_seq_constant_class_var(sequence, attempts=5,
    cutoff = parameters.DISORDER_THRESHOLD, strict=False):
    '''
    A function to generate a variant where the sequence composition is new but
    the numbers of each residue from each class is the same. The overall properties
    of the generated sequence will also be constant.

    Parameters
    ------------
    sequence : str
        the sequence to make a variant of

    attempts : int
        Specify the number of times to make the sequence. Default is 5. Greater numbers
        of attempts increase the odds that a sequence will be generated but will increase
        the duration of attempting to make the sequence.         

    cutoff : float
        the disorder cutoff threshold. Closer to 1 has a higher chance of being disordered.

    strict : bool
        whether to use a strict disorder calculation. By default, variants are allowed

    Returns
    -----------
    final_sequence : str
        the final sequence variant

    '''
    # make sure that the input sequence is all caps
    sequence = sequence.upper()

    # check length
    _length_check(sequence)

    if cutoff > 1 or cutoff < 0:
        raise goose_exceptions.GooseInputError('cutoff value must be between 0 and 1 for disorder threshold')
    try:
        final_sequence = _gen_new_var_constant_class(sequence, attempts=attempts, disorder_threshold=cutoff, strict_disorder=strict)
    except:
        raise goose_exceptions.GooseFail('Sorry! GOOSE was unable to generate the sequence. Please try again or try with a lower cutoff value.')
    return final_sequence



def constant_class_var(sequence, attempts=5,
    cutoff = parameters.DISORDER_THRESHOLD, strict=False):
    '''
    function to generate a variant with the same properties as the 
    input variant as well as the same order of amino acids as
    far as class and the same number in each class
    
    Parameters
    ------------
    sequence : str
        the sequence to make a variant of

    attempts : int
        Specify the number of times to make the sequence. Default is 5. Greater numbers
        of attempts increase the odds that a sequence will be generated but will increase
        the duration of attempting to make the sequence.    

    cutoff : float
        the disorder cutoff threshold. Closer to 1 has a higher chance of being disordered.

    strict : bool
        whether to use a strict disorder calculation. By default, variants are allowed

    Returns
    -----------
    final_sequence : str
        the final sequence variant    
    '''
    # make sure that the input sequence is all caps
    sequence = sequence.upper()

    # check length
    _length_check(sequence)

    if cutoff > 1 or cutoff < 0:
        raise goose_exceptions.GooseInputError('cutoff value must be between 0 and 1 for disorder threshold')    
    try:
        final_sequence = _gen_constant_class_variant(sequence, attempts=attempts, disorder_threshold=cutoff, strict_disorder=strict)
    except:
        raise goose_exceptions.GooseFail('Sorry! GOOSE was unable to generate the sequence. Please try again or try with a lower cutoff value.')
    return final_sequence


def new_var(sequence, attempts=5,
    cutoff = parameters.DISORDER_THRESHOLD, strict=False):
    '''
    function to generate a variant that is completely different
    in sequence to the input but has all the same overall parameters.
    Does not account for specific classes of residues.
    
    Parameters
    ------------
    sequence : str
        the sequence to make a variant of

    attempts : int
        Specify the number of times to make the sequence. Default is 5. Greater numbers
        of attempts increase the odds that a sequence will be generated but will increase
        the duration of attempting to make the sequence.            

    cutoff : float
        the disorder cutoff threshold. Closer to 1 has a higher chance of being disordered.

    strict : bool
        whether to use a strict disorder calculation. By default, variants are allowed

    Returns
    -----------
    final_sequence : str
        the final sequence variant    
    '''
    # make sure that the input sequence is all caps
    sequence = sequence.upper()

    # check length
    _length_check(sequence)

    if cutoff > 1 or cutoff < 0:
        raise goose_exceptions.GooseInputError('cutoff value must be between 0 and 1 for disorder threshold')    
    try:
        final_sequence = _gen_new_variant(sequence, attempts=attempts, disorder_threshold=cutoff, strict_disorder=strict)
    except:
        raise goose_exceptions.GooseFail('Sorry! GOOSE was unable to generate the sequence. Please try again or try with a lower cutoff value.')
    return final_sequence


def hydro_class_var(sequence, hydropathy, hydro_error = parameters.HYDRO_ERROR,
    attempts=5, cutoff = parameters.DISORDER_THRESHOLD, strict=False):
    '''
    function to take in a sequence and make a variant that adjusts the
    hydropathy while keeping the position and number of amino acids the
    same by class of amino acid

    Parameters
    ------------
    sequence : str
        the sequence to make a variant of
    
    hydropathy : float
        the mean hydropathy of the sequence variant. 

    hydro_erorr : float
        the allowed error between the specified hydropathy and the
        hydropathy of the returned sequence

    attempts : int
        Specify the number of times to make the sequence. Default is 5. Greater numbers
        of attempts increase the odds that a sequence will be generated but will increase
        the duration of attempting to make the sequence.        

    cutoff : float
        the disorder cutoff threshold. Closer to 1 has a higher chance of being disordered.

    strict : bool
        whether to use a strict disorder calculation. By default, variants are allowed

    Returns
    -----------
    final_sequence : str
        the final sequence variant    

    '''
    # make sure that the input sequence is all caps
    sequence = sequence.upper()

    # check length
    _length_check(sequence)

    if cutoff > 1 or cutoff < 0:
        raise goose_exceptions.GooseInputError('cutoff value must be between 0 and 1 for disorder threshold')    
    try:
        final_sequence = _gen_hydropathy_class_variant(sequence, hydropathy=hydropathy, allowed_hydro_error = hydro_error, attempts=attempts, disorder_threshold=cutoff, strict_disorder=strict)
    except:
        raise goose_exceptions.GooseFail('Sorry! GOOSE was unable to generate the sequence. Please try again or try with a lower cutoff value.')
    return final_sequence




def constant_residue_var(sequence, constant=[], attempts=5, 
    cutoff = parameters.DISORDER_THRESHOLD, strict=False):
    '''
    function that will generate a new sequence variant
    where specific residues are held constant. The 
    variant will have the same aggregate properties
    as the original sequence.

    Parameters
    ------------
    sequence : str
        the sequence to make a variant of

    constant : list
        A list of residues in the sequence to hold constant

    attempts : int
        Specify the number of times to make the sequence. Default is 5. Greater numbers
        of attempts increase the odds that a sequence will be generated but will increase
        the duration of attempting to make the sequence.   

    cutoff : float
        the disorder cutoff threshold. Closer to 1 has a higher chance of being disordered.

    strict : bool
        whether to use a strict disorder calculation. By default, variants are allowed

    Returns
    -----------
    final_sequence : str
        the final sequence variant 

    '''
    # make sure that the input sequence is all caps
    sequence = sequence.upper()

    # check length
    _length_check(sequence)

    if cutoff > 1 or cutoff < 0:
        raise goose_exceptions.GooseInputError('cutoff value must be between 0 and 1 for disorder threshold')
    try:
        final_sequence = _gen_constant_residue_variant(sequence, constant_residues=constant, attempts=attempts, disorder_threshold=cutoff, strict_disorder=strict)
    except:
        raise goose_exceptions.GooseFail('Sorry! GOOSE was unable to generate the sequence. Please try again or try with a lower cutoff value or with different constant residues.')
    return final_sequence



def shuffle_var(sequence, shuffle=[], attempts=5,
    cutoff = parameters.DISORDER_THRESHOLD, strict=False):
    '''
    Function that will shuffle specific regions of an IDR.
    Multiple regions can be specified simultaneously.

    Parameters
    ------------
    sequence : str
        the sequence to make a variant of

    shuffle : list
        A list of regions to shuffle

    attempts : int
        Specify the number of times to make the sequence. Default is 5. Greater numbers
        of attempts increase the odds that a sequence will be generated but will increase
        the duration of attempting to make the sequence.  

    cutoff : float
        the disorder cutoff threshold. Closer to 1 has a higher chance of being disordered.

    strict : bool
        whether to use a strict disorder calculation. By default, variants are allowed

    Returns
    -----------
    final_sequence : str
        the final sequence variant 
    '''
    # make sure that the input sequence is all caps
    sequence = sequence.upper()

    # check length
    _length_check(sequence)

    if cutoff > 1 or cutoff < 0:
        raise goose_exceptions.GooseInputError('cutoff value must be between 0 and 1 for disorder threshold')    

    all_vals = []
    if shuffle != []:
        for val in shuffle:
            if type(val) == list:
                for subval in val:
                    all_vals.append(subval)
            else:
                all_vals.append(val)
    
    for val in all_vals:
        if val < 1:
            raise goose_exceptions.GooseInputError('Cannot have a value to shuffle below 1')
        if val > len(sequence):
            raise goose_exceptions.GooseInputError('Cannot specify to shuffle a region greater than the length of your sequence')

    curvals = []
    if type(shuffle[0])==list:
        if len(shuffle) >= 2:
            for sublist in shuffle:
                for i in range(sublist[0], sublist[1]+1):
                    if i in curvals:
                        raise goose_exceptions.GooseInputError('Cannot have overlapping regions to shuffle.')
                    else:
                        curvals.append(i)

    try:
        final_sequence = _gen_shuffle_variant(sequence, shuffle_regions=shuffle, attempts=attempts, disorder_threshold=cutoff, strict_disorder=strict)
    except:
        raise goose_exceptions.GooseFail('Sorry! GOOSE was unable to generate the sequence. Please try again or try with a lower cutoff value or with different constant residues.')
    return final_sequence



def kappa_var(sequence, kappa, kappa_error=parameters.MAXIMUM_KAPPA_ERROR, 
    attempts=10, cutoff=parameters.DISORDER_THRESHOLD, strict=False):
    '''
    Function to generate a sequence with a user-defined
    kappa value. Requires kappa calculation using 
    SPARROW. Kappa is a function of charge asymmetry, larger
    kappa values have more asymmetrically distributed
    charged residues.


    Parameters
    ------------
    sequence : str
        the sequence to make a variant of

    kappa : float
        The desired kappa value

    kappa_error : float
        The allowed error between the desired kappa value and the kappa value of the
        returned sequence

    attempts : int
        Specify the number of times to make the sequence. Default is 10. Greater numbers
        of attempts increase the odds that a sequence will be generated but will increase
        the duration of attempting to make the sequence.  

    cutoff : float
        the disorder cutoff threshold. Closer to 1 has a higher chance of being disordered.

    strict : bool
        whether to use a strict disorder calculation. By default, variants are allowed

    Returns
    -----------
    final_sequence : str
        the final sequence variant 
    '''
    # make sure that the input sequence is all caps
    sequence = sequence.upper()

    # check length
    _length_check(sequence)

    if cutoff > 1 or cutoff < 0:
        raise goose_exceptions.GooseInputError('cutoff value must be between 0 and 1 for disorder threshold')    
    if len(sequence) < 6:
        raise GooseInputError('Cannot have sequence with a length less than 6')
    if kappa > 1 or kappa < 0:
        raise GooseInputError('Kappa values must be between 0 and 1')

    try:
        final_sequence = _gen_kappa_variant(sequence, kappa=kappa, allowed_kappa_error = kappa_error, attempts=attempts, disorder_threshold=cutoff, strict_disorder=strict)
    except:
        raise goose_exceptions.GooseFail('Sorry! GOOSE was unable to generate the sequence. Please try again or try with a lower cutoff value or with different kappa value')
    return final_sequence


def asymmetry_var(sequence, increase_decrease, aa_class, number_changes=None,
    attempts=10, cutoff = parameters.DISORDER_THRESHOLD, strict=False):
    '''
    user facing function for generating a variant that has 
    altered asymmetry of some class of residues or user specificed
    list of residues.

    Parameters
    ------------
    sequence : str
        the sequence to make a variant of

    increase_decrease : str
        whether to increase or decrease the asymmetry of a specific 
        residue or residue class. Set to 'increase' or 'decrease'.

    aa_class : str or list
        the class of amino acids to alter asymmetry of. Can specify a custom
        class by passing in a list of residues to target or you can just type in
        the name of a canonical class as a string. Possible classes include:
            'negative' - negative residues
            'positive' - positive residues    
            'proline' - prolines
            'aromatic' - W Y F
            'aliphatic' - I V L A M
            'polar' - Q N S T        

    number_changes : bool
        the number of amino acids to change in the sequence to change asymmetry. 
        By default, will attempt to change asymmetry a number of times equal to
        1/5th the sequence length. Higher numbers increase the amount that the sequence
        becomes more or less symmetrical for the specified residue.

    attempts : int
        Specify the number of times to make the sequence. Default is 10. Greater numbers
        of attempts increase the odds that a sequence will be generated but will increase
        the duration of attempting to make the sequence.          

    cutoff : float
        the disorder cutoff threshold. Closer to 1 has a higher chance of being disordered.

    strict : bool
        whether to use a strict disorder calculation. By default, variants are allowed

    Returns
    -----------
    final_sequence : str
        the final sequence variant     
    '''
    # make sure that the input sequence is all caps
    sequence = sequence.upper()

    # check length
    _length_check(sequence)

    if cutoff > 1 or cutoff < 0:
        raise goose_exceptions.GooseInputError('cutoff value must be between 0 and 1 for disorder threshold')    
    if len(sequence) < 6:
        raise GooseInputError('Cannot have sequence with a length less than 6')

    try:
        final_sequence = _gen_asymmetry_variant(sequence, increase_decrease, aa_class, num_change=number_changes, attempts=attempts, disorder_threshold=cutoff, strict_disorder=strict)
    except:
        raise goose_exceptions.GooseFail('Sorry! GOOSE was unable to generate the sequence. Please try again or try with a lower cutoff value or with different specified amino acids')
    return final_sequence    


def fcr_class_var(sequence, fcr, attempts=10, cutoff=parameters.DISORDER_THRESHOLD, strict=False):
    '''
    user facing funcitonality to generate variants where
    the fcr is changed and the changes to classes of amino
    acids is minimized
    
    Parameters
    ------------
    sequence : str
        the sequence to make a variant of

    fcr : float
        The desired fcr value

    attempts : int
        Specify the number of times to make the sequence. Default is 10. Greater numbers
        of attempts increase the odds that a sequence will be generated but will increase
        the duration of attempting to make the sequence.  

    cutoff : float
        the disorder cutoff threshold. Closer to 1 has a higher chance of being disordered.

    strict : bool
        whether to use a strict disorder calculation. By default, variants are allowed

    Returns
    -----------
    final_sequence : str
        the final sequence variant     
    '''
    # make sure that the input sequence is all caps
    sequence = sequence.upper()

    # check length
    _length_check(sequence)

    if cutoff > 1 or cutoff < 0:
        raise goose_exceptions.GooseInputError('cutoff value must be between 0 and 1 for disorder threshold')    
    if len(sequence) < 6:
        raise GooseInputError('Cannot have sequence with a length less than 6')

    if fcr > 1 or fcr < 0:
        raise goose_exceptions.GooseInputError('fcr values must be between 0 and 1.')

    try:
        final_sequence = _gen_fcr_class_variant(sequence, fcr=fcr, attempts=attempts, disorder_threshold=cutoff, strict_disorder=strict)
    except:
        raise goose_exceptions.GooseFail('Sorry! GOOSE was unable to generate the sequence. Please try again or try with a different fcr value or a different cutoff value.')
    return final_sequence       
    


def ncpr_class_var(sequence, ncpr, attempts=10, cutoff=parameters.DISORDER_THRESHOLD, strict=False):
    '''
    user facing funcitonality to generate variants where
    the ncpr is changed and the changes to classes of amino
    acids is minimized

    Parameters
    ------------
    sequence : str
        the sequence to make a variant of

    ncpr : float
        The desired ncpr value

    attempts : int
        Specify the number of times to make the sequence. Default is 10. Greater numbers
        of attempts increase the odds that a sequence will be generated but will increase
        the duration of attempting to make the sequence.          

    cutoff : float
        the disorder cutoff threshold. Closer to 1 has a higher chance of being disordered.

    strict : bool
        whether to use a strict disorder calculation. By default, variants are allowed

    Returns
    -----------
    final_sequence : str
        the final sequence variant      
    '''
    # make sure that the input sequence is all caps
    sequence = sequence.upper()

    # check length
    _length_check(sequence)

    if cutoff > 1 or cutoff < 0:
        raise goose_exceptions.GooseInputError('cutoff value must be between 0 and 1 for disorder threshold')    
    if len(sequence) < 6:
        raise GooseInputError('Cannot have sequence with a length less than 6')

    if ncpr > 1 or ncpr < -1:
        raise goose_exceptions.GooseInputError('fcr values must be between -1 and 1.')

    try:
        final_sequence = _gen_ncpr_class_variant(sequence, ncpr=ncpr, attempts=attempts, disorder_threshold=cutoff, strict_disorder=strict, constant_fcr=False)
    except:
        raise goose_exceptions.GooseFail('Sorry! GOOSE was unable to generate the sequence. Please try again or try with a different ncpr value or a different cutoff value.')
    return final_sequence   


def all_props_class_var(sequence, hydropathy=None, fcr=None, ncpr=None, kappa=None,
    attempts=10, cutoff=parameters.DISORDER_THRESHOLD, strict=False):
    '''
    user facing funcitonality to generate variants where
    the ncpr, fcr, and/or hydropathy is/are changed and 
    the changes to classes of amino acids is minimized
    
    Parameters
    ------------
    sequence : str
        the sequence to make a variant of

    hydropathy : float
        The desired hydropathy value

    fcr : float
        The desired fcr value

    ncpr : float
        The desired ncpr value

    kappa : float
        The desired kappa value

    attempts : int
        Specify the number of times to make the sequence. Default is 10. Greater numbers
        of attempts increase the odds that a sequence will be generated but will increase
        the duration of attempting to make the sequence.  

    cutoff : float
        the disorder cutoff threshold. Closer to 1 has a higher chance of being disordered.

    strict : bool
        whether to use a strict disorder calculation. By default, variants are allowed

    Returns
    -----------
    final_sequence : str
        the final sequence variant       
    '''
    # make sure that the input sequence is all caps
    sequence = sequence.upper()

    # check length
    _length_check(sequence)

    if cutoff > 1 or cutoff < 0:
        raise goose_exceptions.GooseInputError('cutoff value must be between 0 and 1 for disorder threshold')    
    if len(sequence) < 6:
        raise GooseInputError('Cannot have sequence with a length less than 6')

    if ncpr != None:
        if ncpr > 1 or ncpr < -1:
            raise goose_exceptions.GooseInputError('fcr values must be between -1 and 1.')

    if hydropathy != None:
        if hydropathy < 0 or hydropathy > 6.4:
            raise goose_exceptions.GooseInputError('hydropathy values must be between 0 and 6.4')

    if fcr != None:
        if fcr > 1 or fcr < 0:
            raise goose_exceptions.GooseInputError('fcr values must be between 0 and 1.')

    if kappa != None:
        if kappa > 1 or kappa < 0:
            raise goose_exceptions.GooseInputError('kappa values must be between 0 and 1')


    try:
        final_sequence = _gen_all_props_class_variant(sequence, hydropathy=hydropathy, fcr=fcr,
        ncpr=ncpr, kappa=kappa, attempts=attempts, disorder_threshold=cutoff, strict_disorder=strict)
    except:
        raise goose_exceptions.GooseFail('Sorry! GOOSE was unable to generate the sequence. Please try again or try with a different input values or a different cutoff value.')
    return final_sequence   


def targeted_shuffle_variant(sequence, target_aas, attempts=10,
    cutoff=parameters.DISORDER_THRESHOLD, strict=False):
    '''
    user facing funcitonality to generate variants where
    you can specify residues or classes of residues to shuffle.

    parameters
    ----------
    sequence : str
        the amino acid sequence as a string

    target_aas : str or list
        a list of amino acids to shuffle
        or a class of amino acids to shuffle
        Possible target classes:
            charged : DEKR
            polar : QNST
            aromatic : FYW
            aliphatic : IVLAM
            negative: DE
            positive : KR

    attempts : int
        Specify the number of times to make the sequence. Default is 10. Greater numbers
        of attempts increase the odds that a sequence will be generated but will increase
        the duration of attempting to make the sequence.  

    cutoff : float
        the cutoff value for disorder between 0 and 1. 
        Higher values have a higher likelihood of being disordered. 

    strict : bool
        whether to enforce a strict disorder cutoff value. For variants,
        Lower than cutoff values are allowed for areas of the input sequence
        that have below cutoff values. Nothing in the returned sequence
        will have a higher diosrder value ata ny residue location than the value at that
        position for the input sequence. 
    
    Returns
    -------
    returns the amino acid sequence as a string. 
    '''
    # make sure that the input sequence is all caps
    sequence = sequence.upper()

    # check length
    _length_check(sequence)

    if cutoff > 1 or cutoff < 0:
        raise goose_exceptions.GooseInputError('cutoff value must be between 0 and 1 for disorder threshold')    
    if len(sequence) < 6:
        raise GooseInputError('Cannot have sequence with a length less than 6')

    # make sequence.
    try:
        final_sequence = _gen_targeted_shuffle_variant(sequence, target_aas, attempts=attempts, disorder_threshold=cutoff, strict_disorder=strict)
    except:
        raise goose_exceptions.GooseFail('Sorry! GOOSE was unable to generate the sequence. Please try again or try with a different input values or a different cutoff value.')
    return final_sequence

def excluded_shuffle_variant(sequence, exclude_aas, attempts=10,
    cutoff=parameters.DISORDER_THRESHOLD, strict=False):
    '''
    user facing funcitonality to generate variants where
    you can specify residues or classes of residues to NOT shuffle

    parameters
    ----------
    sequence : str
        the amino acid sequence as a string

    exclude_aas : str or list
        a list of amino acids to not shuffle
        or a class of amino acids to NOT shuffle
        Possible target classes:
            charged : DEKR
            polar : QNST
            aromatic : FYW
            aliphatic : IVLAM
            negative: DE
            positive : KR

    attempts : int
        Specify the number of times to make the sequence. Default is 10. Greater numbers
        of attempts increase the odds that a sequence will be generated but will increase
        the duration of attempting to make the sequence.  

    cutoff : float
        the cutoff value for disorder between 0 and 1. 
        Higher values have a higher likelihood of being disordered. 

    strict : bool
        whether to enforce a strict disorder cutoff value. For variants,
        Lower than cutoff values are allowed for areas of the input sequence
        that have below cutoff values. Nothing in the returned sequence
        will have a higher diosrder value ata ny residue location than the value at that
        position for the input sequence. 
    
    Returns
    -------
    returns the amino acid sequence as a string. 

    '''
    # make sure that the input sequence is all caps
    sequence = sequence.upper()

    # check length
    _length_check(sequence)

    if cutoff > 1 or cutoff < 0:
        raise goose_exceptions.GooseInputError('cutoff value must be between 0 and 1 for disorder threshold')    
    if len(sequence) < 6:
        raise GooseInputError('Cannot have sequence with a length less than 6')

    # make sequence.
    try:
        final_sequence = _gen_excluded_shuffle_variant(sequence, exclude_aas, attempts=attempts, disorder_threshold=cutoff, strict_disorder=strict)
    except:
        raise goose_exceptions.GooseFail('Sorry! GOOSE was unable to generate the sequence. Please try again or try with a different input values or a different cutoff value.')
    return final_sequence


'''
/-/-/-/-/-/-/-/-/-/-/-/-/-/-/-/-/-/-/
/-/-/-/-/-/-/-/-/-/-/-/-/-/-/-/-/-/-/
 PREDICTED FOLDED REGION GENERATORS
/-/-/-/-/-/-/-/-/-/-/-/-/-/-/-/-/-/-/
/-/-/-/-/-/-/-/-/-/-/-/-/-/-/-/-/-/-/
'''

def alpha_helix(length, attempts=500):
    '''
    user facing function for generating a sequence predicted to be an alpha
    helix based on DSSP scores.

    Parameters
    ----------
    length : int
        the length of the desired alpha helix

    attempts : int
        Specify the number of times to make the sequence. Default is 500. Greater numbers
        of attempts increase the odds that a sequence will be generated but will increase
        the duration of attempting to make the sequence.  

    returns
    -------
    final_seq : str

    '''
    if length > 150:
        raise goose_exceptions.GooseInputError('Unable to make alpha helix with length greater than 150.')    
    elif length < 8:
        raise goose_exceptions.GooseInputError('Unable to make alpha helix with length less than 8.')    
    else:
        try:
            final_seq = _gen_helix(length, max_iters=attempts)
        except:
            raise goose_exceptions.GooseFail('Sorry! Goose was unable to make that helix. Try again or try a different length.')
        return final_seq


def beta_strand(length, attempts=5000):
    '''
    user facing function for generating a sequence predicted to be 
    a beta strand based on DSSP scores.

    Parameters
    ----------
    length : int
        the length of the desired beta strand
    
    attempts : int
        Specify the number of times to make the sequence. Default is 5000. Greater numbers
        of attempts increase the odds that a sequence will be generated but will increase
        the duration of attempting to make the sequence.  
    
    returns
    -------
    final_seq : str    
    '''
    if length > 34:
        raise goose_exceptions.GooseInputError('Unable to make beta strands with length greater than 34.')
    elif length < 5:
        raise goose_exceptions.GooseInputError('Unable to make beta strands with length less than 5.')
    else:
        try:
            final_seq = _gen_beta_strand(length, max_iters=attempts)
        except:
            raise goose_exceptions.GooseFail('Sorry! Goose was unable to make that strand. Try again or try a different length.')
        return final_seq



def beta_sheet(length):
    '''
    user facing function for generating a sequence predicted to be 
    a beta sheet based on DSSP scores. uses coils to connect strands.

    Parameters
    ----------
    length : int
        the length of the desired beta_sheet helix

    returns
    -------
    final_seq : str    
    '''
    if length < 18:
        raise goose_exceptions.GooseInputError('cannot generate beta sheet less than 18 amino acids.')
    elif length > 400:
        raise goose_exceptions.GooseInputError('cannot generate beta sheet greater than 400 amino acids.')
    try:
        final_seq = _gen_beta_sheet(length)
    except:
        raise goose_exceptions.GooseFail('Sorry! Goose was unable to make that beta sheet. Try again or try a different length.')
    return final_seq



'''
/-/-/-/-/-/-/-/-/-/-/-/-/-/-/-/-/-/-/
/-/-/-/-/-/-/-/-/-/-/-/-/-/-/-/-/-/-/
        LIBRARY GENERATION
/-/-/-/-/-/-/-/-/-/-/-/-/-/-/-/-/-/-/
/-/-/-/-/-/-/-/-/-/-/-/-/-/-/-/-/-/-/
'''


def seq_property_library(length,
    FCR=None, NCPR=None, hydropathy=None, kappa=None, 
    cutoff=parameters.DISORDER_THRESHOLD, silent_failed_seqs=False, beta_mode=False, random_name=False):
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

    cutoff : float
        The cutoff value for disorder as a float. Higher values
        lead to a more 'strict' cutoff for what is considered to be disordered.

    silent_failed_seqs : bool
        Whether to silence any printed warnings of sequences that
        are not possible to generate due to incompatible charge/hydorpatyh values

    beta_mode : bool
        For testing. Set to True to get some printouts as seqs are being generated.

    random_name : bool
        whether to use a randomly generated string as the protein name 

    returns
    -------
    sequence_list : list of dictionaries
        Returns a list of dicts where each dictionary in the list 
        can be input into the sequence generators downstream. 
        The dict holds the parameter specifications for each sequence in the library.
    '''
    # make the list of sequences to generate
    sequence_list = _generate_library_by_parameter_ranges(length=length,
    FCR=FCR, NCPR=NCPR, hydropathy=hydropathy, kappa=kappa, silent_failed_seqs=silent_failed_seqs)
    # dict to return
    seq_dict={}
    # make a dict to hold the library
    for seq_specified in sequence_list:
        if beta_mode==True:
            print(f'{seq_specified}')
        # make the sequence
        seq = sequence(length, **seq_specified)
        # make a sequence name
        seq_name='>'
        name_values = ['FCR', 'NCPR', 'kappa', 'hydropathy']
        for vals in range(0,len(name_values)):
            param_val = name_values[vals]
            curval = seq_specified[param_val]
            if curval != None:
                if seq_name=='>':
                    seq_name+=f'{param_val}_{round(curval, 4)}'
                else:
                    seq_name+=f'_{param_val}_{round(curval, 4)}'
        # if user wants a randodm name, make it happen.
        if random_name==True:
            seq_name+=f'_{_gen_random_name()[1:]}'

        # add to the dict
        seq_dict[seq_name]=seq
    # return the dict
    return seq_dict

def seq_fractions_library(length, random_name=False, warn_user=True, robust_warning=False, **kwargs):
    '''
    User facing function for generating a library of sequences that have varying specified fractions
    of aminon acids. This function will only make fractions that are possible (ie. less than 1.0 specified
    and all values less than the respective value for each seq as far as the maximum possible). By default will
    alert the user to any 'failed' generate sequences...


    Parameters
    ------------
    length : int
        the length of the sequence to generate

    random_name : bool
        whether to make a random name for each sequence.

    warn_user : bool
        whether to warn the user if the fractions of a sequence do not match 
        what the user input. This is typically due to length problems (bascially,
        you can specify F=0.95 but will never get 0.95 for F if your sequence is 10
        amino acids long because you can't have 1/2 a F)

    robust_warning : bool
        Whether to print out a robust message for each incorrect sequence highlighting what was wrong

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
    seq_dict : dict
        returns a dictionary of the sequences.   
    '''

    # make the list of sequences to generate
    sequence_list = _generate_library_by_fraction_ranges(length=length, **kwargs)
    # need dict for naming...
    possible_vals_name=['A', 'C', 'D', 'E', 'F', 'G', 'H', 'I', 'K', 'L', 'M', 'N', 'P', 'Q', 'R', 'S', 'T', 'V', 'W', 'Y']
    name_vals = []
    for specified_arg in kwargs.keys():
        if specified_arg in possible_vals_name:
            name_vals.append(specified_arg)
    # keep track of sequences with erros
    sequence_warnings=0
    # dict to return
    seq_dict={}
    # list of additional kwargs
    additional_kwargs=['cutoff', 'attempts',  'strict_disorder',  'max_aa_fractions']
    for seq_specified in sequence_list[1]:
        # add in any additional kwargs to input into seq_fractions function
        for cur_kwarg in additional_kwargs:
            if cur_kwarg in kwargs.keys():
                seq_specified[cur_kwarg]=kwargs[cur_kwarg]
        sequence = seq_fractions(length, **seq_specified)
        # make a sequence name
        seq_name='>'
        sequence_correct=True
        error_warning=f'For:\n{sequence}\n'
        if name_vals != []:
            for vals in name_vals:
                actual_val = sequence.count(vals)/len(sequence)
                param_val = seq_specified[vals]
                # see if actual val equals param val
                if round(float(actual_val), 8) != round(float(param_val), 8):
                    sequence_correct=False
                    error_warning += f'Objective fraction for {vals}: {round(float(param_val), 8)}, Actual value: {round(float(actual_val), 8)}\n'
                if random_name==False:                
                    if seq_name=='>':
                        seq_name+=f'{vals}_{round(actual_val, 4)}'
                    else:
                        seq_name+=f'_{vals}_{round(actual_val, 4)}'
                # if user wants a randodm name, make it happen.
                else:
                    seq_name+=f'_{_gen_random_name()[1:]}'
            if sequence_correct==False:
                sequence_warnings+=1
                if robust_warning==True:
                    print(error_warning)
                    print()

        # make sure seq_name not already in dict
        if seq_name not in seq_dict.keys():
            # add to seq dict
            seq_dict[seq_name]=sequence
        else:
            # this shouldn't ever be a problem, but just to be safe add in a random 4 digit ID.
            add_id = '_ID_'
            for i in range(0, 8):
                add_id += f'{random.randint(0, 9)}'
            seq_name+=add_id
            seq_dict[seq_name]=sequence
    # warn user if needed and if wanted by user
    if warn_user == True:
        if sequence_warnings > 0:
            print(f'\nWARNING!\n GOOSE detected {sequence_warnings} sequences out of {len(seq_dict)} sequences where the generated sequence does not have the exact value as the input.\n This is often due to a mismatch in the specified length and possible fractions.\n')
    return seq_dict




