'''
Uses the variant_generation_functions to make a sequence. Important thing here is that
it also checks for disorder. These things are separated to keep things faster and more modular.
'''

from goose.backend_variant_generation.helper_functions import check_variant_disorder_vectorized
from goose.backend_property_optimization.optimize_dimensions import predict_rg, predict_re
from goose.backend_variant_generation import variant_generation_functions as vgf


# code for variant generation with disorder checking
def gen_constant_class_variant(input_sequence,
                               num_attempts=5,
                               strict_disorder=False,
                               disorder_cutoff=0.5,
                               metapredict_version=3):
    """
    Generate a constant class variant of the input sequence.
    Parameters
    ----------
    input_sequence : str
        The input sequence to generate a variant from.
    num_attempts : int
        Number of attempts to generate a variant.
        default is 5
    strict_disorder : bool
        If True, all residues must be above the disorder cutoff to be considered disordered.
        default is False
    disorder_cutoff : float
        Cutoff for disorder. Above this value is considered disordered.
        default is 0.5
    metapredict_version : int
        Version of MetaPredict to use (1, 2, or 3)
        default is 3

    Returns
    -------
    str
        A generated variant sequence that meets the disorder criteria.
    """
    for _ in range(num_attempts):
        # Generate a variant sequence
        variant_sequence = vgf.gen_constant_class_variant(input_sequence)
        if variant_sequence == None:
            continue
        
        # Check if the variant sequence meets the disorder criteria
        disordered_sequence = check_variant_disorder_vectorized(
            original_sequence=input_sequence,
            sequences=variant_sequence,
            strict_disorder=strict_disorder,
            disorder_cutoff=disorder_cutoff,
            metapredict_version=metapredict_version,
            return_best_sequence=False
        )
        if disordered_sequence is None:
            continue
    
    # return the disordered sequence. This could be None if no valid sequence was found
    # but we will handle that exception in the frontend to avoid an annoying
    # string of exceptions. 
    return disordered_sequence

