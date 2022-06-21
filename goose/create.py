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

'''
MIGHT BE WORTH TESTING...

CHANGE ALL GENERATORS TO MAKE A BUNCH OF SEQEUNCES WITH NECECSSARY PARAMETERS THAN TEST FOR DISORDER
'''



'''
# testing code below

next steps:
implement variant generation


2 types of variants:
class variants
    accounts for classes of amino acids.

non-class variants
    doesn't care about the class of amino acids.


Also want:
minimal change variant
charge asymmetry variants



****
For charge asymmetery variants, can basically just use the function
from OG GOOSE and let the user specify how many times to run
it where increased number of times increasinlgy changes the asymetry
****


need to import:



from goose.backend.protein import Protein

test = sequence(100, FCR=0.1, NCPR=0.1)
print(Protein.all_properties_library(test))

'''
