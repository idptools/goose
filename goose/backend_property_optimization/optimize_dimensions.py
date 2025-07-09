import random

from sparrow.predictors import batch_predict
from sparrow.protein import Protein

from goose.backend import parameters
from goose.backend.parameters import rg_re_attempt_num
from goose.backend.lists import disordered_list, disordered_list_reduced_charge

def predict_rg(sequence): return Protein(sequence).predictor.radius_of_gyration(use_scaled=True)
def predict_re(sequence): return Protein(sequence).predictor.end_to_end_distance(use_scaled=True)
def batch_predict_rg(sequences,show_progress_bar=False, return_seq2prediction=False,
                        batch_size=32): 
    seqobs=[Protein(seq) for seq in sequences]
    rgs=batch_predict.batch_predict(seqobs, 'scaled_rg',
        show_progress_bar=show_progress_bar, batch_size=batch_size,
        return_seq2prediction=return_seq2prediction)
    return rgs

def batch_predict_re(sequences, show_progress_bar=False, return_seq2prediction=False,
                        batch_size=32):
    seqobs=[Protein(seq) for seq in sequences]
    re=batch_predict.batch_predict(seqobs, 'scaled_re',
        show_progress_bar=show_progress_bar, batch_size=batch_size,
        return_seq2prediction=return_seq2prediction)
    return re


def optimize_seq_dims(input_sequence, rg_or_re, objective_dim, allowed_error=None,
                      num_attempts=rg_re_attempt_num, reduce_pos_charged=False, exclude_aas=None,
                      variants_per_iteration=32, mutation_fraction=0.0125,
                      return_all_sequences=False):
    """
    Optimize a sequence to achieve a target radius of gyration (Rg) or end-to-end distance (Re).

    This function uses an iterative evolutionary approach to modify an input sequence
    to achieve the desired dimensional property while maintaining amino acid composition
    constraints.

    Parameters
    ----------
    input_sequence : str
        The input amino acid sequence to optimize
    rg_or_re : {'rg', 're'}
        Target property: 'rg' for radius of gyration, 're' for end-to-end distance
    objective_dim : float
        Target value for the specified dimensional property
    allowed_error : float or 'default_error'
        Maximum allowed deviation from target. If 'default_error', uses
        backend parameter defaults (re_error or rg_error)
    num_attempts : int
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
    return_all_sequences : bool, default=False
        If True, returns all generated sequences that meet the target criteria.
        If False, returns only the best sequence found.

    Returns
    -------
    str
        Optimized sequence with the target dimensional property

    Raises
    ------
    ValueError
        If input parameters are invalid or optimization constraints cannot be satisfied
    RuntimeError
        If optimization fails to converge within the specified number of attempts

    Notes
    -----
    The optimization uses bias dictionaries to favor amino acids that promote
    either sequence collapse or expansion based on the current vs target dimensions.
    """
    # Set default error tolerance
    if allowed_error == None:
        allowed_error = parameters.rg_error if rg_or_re == 'rg' else parameters.re_error
    
    # Configure amino acid pools and bias dictionaries
    if reduce_pos_charged:
        available_aas = list(disordered_list_reduced_charge)
        bias_dict = {
            'collapse': ['W', 'Y', 'G', 'F', 'Q', 'N'], 
            'expand': ['D', 'E', 'P', 'S', 'T']
        }
    else:
        available_aas = list(disordered_list)
        bias_dict = {
            'collapse': ['W', 'Y', 'G', 'F', 'Q', 'N'], 
            'expand': ['D', 'E', 'K', 'R', 'P', 'S', 'T']
        }
    
    # Handle amino acid exclusions
    if exclude_aas is not None:
        exclude_set = set(exclude_aas)
        
        # Remove excluded amino acids from available pool
        available_aas = [aa for aa in available_aas if aa not in exclude_set]
        
        # Remove excluded amino acids from bias dictionaries
        for bias_type in bias_dict:
            bias_dict[bias_type] = [aa for aa in bias_dict[bias_type] if aa not in exclude_set]
    
    # Choose prediction function
    predict_func = batch_predict_rg if rg_or_re == 'rg' else batch_predict_re
    
    # Get initial dimensional property
    current_sequence = input_sequence
    current_dim = predict_func([current_sequence], return_seq2prediction=True)[current_sequence]
    
    # Check if sequence already meets criteria
    if abs(current_dim - objective_dim) <= allowed_error:
        return current_sequence
    
    # Optimization variables
    best_sequence = current_sequence
    best_error = abs(current_dim - objective_dim)
    mutation_size = max(1, int(len(input_sequence) * mutation_fraction))
    
    # list for all sequences
    all_sequences = []

    # Optimization loop
    for attempt in range(num_attempts):
        # Determine which amino acids to favor based on current vs target
        if current_dim > objective_dim:
            favored_aas = bias_dict['collapse']
        else:
            favored_aas = bias_dict['expand']
        
        # Generate sequence variants
        variants = []
        current_list = list(current_sequence)
        
        for _ in range(variants_per_iteration):
            variant = current_list.copy()
            
            # Randomly select positions to mutate
            mutation_positions = random.sample(range(len(variant)), mutation_size)
            
            # Replace selected positions with favored amino acids
            for pos in mutation_positions:
                variant[pos] = random.choice(favored_aas)
            
            variants.append(''.join(variant))
        
        # Predict dimensional properties for all variants
        predictions = predict_func(variants, return_seq2prediction=True,
                                   batch_size=len(variants))
        
        # Find the best variant
        for sequence, predicted_dim in predictions.items():
            error = abs(predicted_dim - objective_dim)
            
            # Check if we've reached the target
            if error <= allowed_error:
                if return_all_sequences:
                    all_sequences.append(sequence)
                    return all_sequences
                else:
                    return sequence
            
            # Update best candidate
            if error < best_error:
                best_error = error
                best_sequence = sequence
                current_sequence = sequence
                current_dim = predicted_dim
                if return_all_sequences:
                    all_sequences.append(sequence)

    return best_sequence if not return_all_sequences else all_sequences

