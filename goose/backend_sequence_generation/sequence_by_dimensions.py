import random

from goose.backend import parameters
from goose.backend.lists import disordered_list, disordered_list_reduced_charge
from goose.backend_property_optimization.optimize_dimensions import optimize_seq_dims
from goose.backend_property_optimization.helper_functions import generate_random_sequence
from goose.data import aa_list_probabilities as aa_probs


def create_seq_by_dims(seq_length, objective_dim, rg_or_re='rg',
                       allowed_error=parameters.MAXIMUM_RG_RE_ERROR, 
                       num_attempts_dimensions=parameters.RG_RE_ATTEMPT_NUMBER,
                       reduce_pos_charged=True, exclude_aas=None,
                       variants_per_iteration=64, mutation_fraction=0.0125,
                       num_iterations=10, starting_probabilities=None):
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
    allowed_error : float or None, default=parameters.MAXIMUM_RG_RE_ERROR
        Allowed deviation from the target dimension. If None, uses backend parameter defaults.
        For 'rg', default is parameters.rg_error; for 're', default is parameters.re_error.
    num_attempts_dimensions : int, default=parameters.RG_RE_ATTEMPT_NUMBER
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
    num_iterations : int, default=10
        Number of iterations to attempt to make the sequence starting from random. 
    starting_probabilities : dict, optional
        Initial amino acid probabilities to use for sequence generation.
        If None, uses default probabilities from aa_list_probabilities.
    Returns
    -------
    str or None
        Generated sequence that meets the target dimensional property and is disordered.
        If no suitable sequence is found, returns None.

    """

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
    
    for _ in range(num_iterations):
        # Generate initial random sequence
        initial_sequence = generate_random_sequence(seq_length,
                                                    probabilities=starting_probabilities)
        
        # Use the optimize_seq_dims function to refine the sequence
        optimized_sequence = optimize_seq_dims(
            input_sequence=initial_sequence,
            rg_or_re=rg_or_re,
            objective_dim=objective_dim,
            allowed_error=allowed_error,
            num_attempts=num_attempts_dimensions,
            reduce_pos_charged=reduce_pos_charged,
            exclude_aas=exclude_aas,
            variants_per_iteration=variants_per_iteration,
            mutation_fraction=mutation_fraction,
            return_all_sequences=True
        )

        # if no sequence was found, continue to the next attempt
        if optimized_sequence is not None:
            return optimized_sequence

    # If no sequence meets the criteria after all attempts, return None
    return None
        
    