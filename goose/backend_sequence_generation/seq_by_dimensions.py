import random
import numpy as np

from sparrow.predictors import batch_predict
from sparrow.protein import Protein

from goose.backend import parameters
from goose.backend.parameters import re_error, rg_error, rg_re_attempt_num
from goose.backend import lists
from goose.backend.lists import disordered_list, disordered_list_reduced_charge
from goose.backend_property_optimization.optimize_dimensions import optimize_seq_dims


def create_seq_by_dims(seq_length, objective_dim, rg_or_re='rg',
                       allowed_error=None, num_attempts=rg_re_attempt_num,
                       reduce_pos_charged=True, exclude_aas=None,
                       variants_per_iteration=32, mutation_fraction=0.0125):
    """
    Create a sequence of a specified length that meets a target radius of gyration (Rg)
    or end-to-end distance (Re) with a given error tolerance.

    Parameters
    ----------
    seq_length : int
        Length of the sequence to generate
    objective_dim : float
        Target value for the specified dimensional property
    rg_or_re : {'rg', 're'}, default='rg'
        Target property: 'rg' for radius of gyration, 're' for end-to-end distance
    allowed_error : float or None, default=None
        Maximum allowed deviation from target. If None, uses backend parameter defaults
    num_attempts : int, default=rg_re_attempt_num
        Maximum number of optimization iterations
    reduce_pos_charged : bool, default=True
        Whether to reduce positively charged residues (K, R) in the sequence.
        Based on in vivo data suggesting positive charges may not drive expansion
        as predicted by current models.
    exclude_aas : list of str, optional
        Amino acids to exclude from the optimization process
    variants_per_iteration : int, default=32    
        Number of sequence variants to generate per optimization iteration
    mutation_fraction : float, default=0.0125
        Fraction of sequence length to mutate per iteration (minimum 1 residue)
        
    Returns 
    -------
    str
        Generated sequence that meets the target dimensional property
        
    Raises
    ------
    RuntimeError
        If optimization fails to find a suitable sequence
    ValueError
        If input parameters are invalid or constraints cannot be satisfied
    """

    # Set default error tolerance
    if allowed_error is None:
        allowed_error = parameters.rg_error if rg_or_re == 'rg' else parameters.re_error
    
    # Configure amino acid pools based on constraints
    if reduce_pos_charged:
        available_aas = list(disordered_list_reduced_charge)
    else:
        available_aas = list(disordered_list)
    
    # Handle amino acid exclusions
    if exclude_aas is not None:
        exclude_set = set(exclude_aas)
        available_aas = [aa for aa in available_aas if aa not in exclude_set]
        
        if not available_aas:
            raise ValueError("Cannot exclude all available amino acids")
    
    # Generate initial random sequence
    initial_sequence = ''.join(random.choices(available_aas, k=seq_length))
    
    # Use the optimize_seq_dims function to refine the sequence
    optimized_sequence = optimize_seq_dims(
        input_sequence=initial_sequence,
        rg_or_re=rg_or_re,
        objective_dim=objective_dim,
        allowed_error=allowed_error,
        num_attempts=num_attempts,
        reduce_pos_charged=reduce_pos_charged,
        exclude_aas=exclude_aas,
        variants_per_iteration=variants_per_iteration,
        mutation_fraction=mutation_fraction
    )
    return optimized_sequence
        
    