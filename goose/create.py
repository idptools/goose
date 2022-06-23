'''
user-facing functionality
'''


#__all__ =  ['sequence', 'seq_by_fractions', 'sequence_library', 'seq_fraction_library']

import os
import sys

from goose.backend.sequence_generation import generate_disordered_seq_by_fractions as _generate_disordered_seq_by_fractions
from goose.backend.sequence_generation import generate_disordered_seq_by_props as _generate_disordered_seq_by_props
from goose.backend.goose_tools import check_and_correct_props_kwargs as _check_and_correct_props_kwargs
from goose.backend.goose_tools import check_props_parameters as _check_props_parameters
from goose.backend.goose_tools import check_and_correct_fracs_kwargs as _check_and_correct_fracs_kwargs
from goose.backend.goose_tools import check_fracs_parameters as _check_fracs_parameters



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

    Parameters
    ------------
    length : Int
        length of desired disordered sequence
    **kwargs : parameter and parameter value
        Takes in possible parameters for GOOSE. Possible parameters/parameter combinations are:
            FCR : Float
                specify the fraction of charged residues (between 0 and 1)
            NCPR : Float 
                specify the net charge per residue of generated sequences (between -1 and 1)
            sigma : Float
                specify the sigma value of generated sequences(between 0 and 1)
            hydro : float 
                specify the mean hydropathy of generated sequences
            cutoff : Float
                the disorder cutoff threshold
            **Note** can specify NCPR and FCR simultaneously
                     can specify NCPR and FCR and hydro simultaneously


    Returns
    -----------
    generated_seq : String
        Returns a string that is the amino acid sequence

    """

    # First correct kwargs. Do this first because
    # the next function that looks over kwargs values
    # can only take in corrected kwargs.
    kwargs = _check_and_correct_props_kwargs(**kwargs)

    # now make sure that the input vals are within appropriate bounds
    _check_props_parameters(**kwargs)

    # make sure length is within bounds
    if length > parameters.MAXIMUM_LENGTH:
        error_message = f'length of {length} is greater than maximum allowed value of {parameters.MAXIMUM_LENGTH}'
        raise goose_exceptions.GooseInputError(error_message)
    if length < parameters.MINIMUM_LENGTH:
        error_message = f'length of {length} is less than maximum allowed value of {parameters.MINIMUM_LENGTH}'
        raise goose_exceptions.GooseInputError(error_message)

    # make the sequence
    generated_seq = _generate_disordered_seq_by_props(length, FCR=kwargs['FCR'], NCPR=kwargs['NCPR'], hydropathy=kwargs['hydropathy'],
        sigma = kwargs['sigma'], attempts = 1, allowed_hydro_error = parameters.HYDRO_ERROR, disorder_threshold = kwargs['cutoff'])

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


def sequence_fractions(length, **kwargs):
    """
    Stand alone function that takes care of creating sequences with specified
    fractions of amino acids. 

    Parameters
    ------------
    length : Int
        length of desired disordered sequence
    **kwargs : parameter and parameter value
        Takes in amino acids followed by the fraction value as a decimal.


    Returns
    -----------
    generated_seq : String
        Returns a string that is the amino acid sequence

    """

    # First correct kwargs. Do this first because
    # the next function that looks over kwargs values
    # can only take in corrected kwargs.
    kwargs = _check_and_correct_fracs_kwargs(**kwargs)

    # now make sure that the input vals are within appropriate bounds
    _check_fracs_parameters(**kwargs)

    # make sure length is within bounds
    # length is the only thing not checked by my check / check and correct functions.
    if length > parameters.MAXIMUM_LENGTH:
        error_message = f'length of {length} is greater than maximum allowed value of {parameters.MAXIMUM_LENGTH}'
        raise goose_exceptions.GooseInputError(error_message)
    if length < parameters.MINIMUM_LENGTH:
        error_message = f'length of {length} is less than maximum allowed value of {parameters.MINIMUM_LENGTH}'
        raise goose_exceptions.GooseInputError(error_message)

    generated_seq = _generate_disordered_seq_by_fractions(length, **kwargs)

    # return the seq
    return generated_seq


"""
now need to add variant creation to user facing functionality.


def gen_minimal_sequence_variant(input_sequence, mean_hydro = '',
 fraction = '', net_charge = '', charge_asymmetry='', cutoff=parameters.DISORDER_THRESHOLD, strict=False):
    '''
    tries to change the paramters you want 
    while minimally changing the starting sequence
    '''


def gen_new_var_constant_class(sequence, attempts=5, 
    disorder_threshold = parameters.DISORDER_THRESHOLD, strict_disorder=False):
    '''
    A function to generate a variant where the sequence composition is new but
    the numbers of each residue from each class is the same. The overall properties
    of the generated sequence will also be constant.

def gen_constant_class_variant(sequence, attempts=5, 
    disorder_threshold = parameters.DISORDER_THRESHOLD, strict_disorder=False):
    '''
    function to generate a variant with the same properties as the 
    input variant as well as the same order of amino acids as
    far as class and the same number in each class

def gen_new_variant(sequence, attempts=5, 
    disorder_threshold = parameters.DISORDER_THRESHOLD, strict_disorder=False):
    '''
    function to generate a variant that is completely different
    in sequence to the input but has all the same overall parameters.
    Does not account for specific classes of residues.

def gen_hydropathy_class_variant(sequence, hydropathy, allowed_hydro_error = parameters.HYDRO_ERROR,
    attempts=5, disorder_threshold = parameters.DISORDER_THRESHOLD, strict_disorder=False):
    '''
    function to take in a sequence and make a variant that adjusts the
    hydropathy while keeping the position and nuimber of amino acids the
    same by class of amino acid


def gen_constant_residue_variant(sequence, constant_residues = [],
    attempts=5, disorder_threshold = parameters.DISORDER_THRESHOLD, strict_disorder=False):
    '''
    function that will generate a new sequence variant
    where specific residues are held constant. The 
    variant will have the same aggregate properties
    as the original sequence.


def gen_shuffle_variant(sequence, shuffle_regions = [], use_index=False,
    attempts=5, disorder_threshold = parameters.DISORDER_THRESHOLD, strict_disorder=False):
    '''
    Function that will shuffle specific regions of an IDR.
    Multiple regions can be specified simultaneously.


def gen_kappa_variant(sequence, kappa, allowed_kappa_error = parameters.MAXIMUM_KAPPA_ERROR,
    attempts=5, disorder_threshold = parameters.DISORDER_THRESHOLD, strict_disorder=False):
    '''
    Function to generate a sequence with a user-defined
    kappa value. Requires kappa calculation using 
    SPARROW. Kappa is a function of charge asymmetry, larger
    kappa values have more asymmetrically distributed
    charged residues.

Once done with that, update the docs!


"""