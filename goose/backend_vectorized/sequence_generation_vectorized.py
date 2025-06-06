'''
Version of sequence_generation.py but using the new vectorized functions.
The strategy here is to generate many sequences rapidly and then filter for ones we want.
'''
import random
import numpy as np
import metapredict as meta
from sparrow.protein import Protein
from goose import goose_exceptions
from goose.backend_vectorized.optimize_kappa_vectorized import optimize_kappa_vectorized
from goose.backend_vectorized.optimize_hydropathy_vectorized import optimize_hydropathy_vectorized
from goose.backend_vectorized.seq_by_probability_vectorized import SequenceGenerator
from goose.backend_vectorized.seq_by_fractions_vectorized import FractionBasedSequenceGenerator



def check_disorder_vectorized(sequences, 
                              strict_disorder=False,
                              disorder_cutoff=0.5,
                              max_consecutive_ordered=3,
                              max_total_ordered=0.05,
                              metapredict_version=3):
    """
    Check disorder of sequences using vectorized operations.
    
    Parameters
    ----------
    sequences : list, str, dict
        list, str, or dict of sequences to check for disorder
    strict_disorder : bool
        If True, all residues must be above the disorder cutoff
        default is False
    disorder_cutoff : float
        Cutoff for disorder. Above this value is considered disordered.
    max_consecutive_ordered : int
        Maximum number of consecutive residues allowed to be below the disorder cutoff
    max_total_ordered : float or int
        If float (0-1): Maximum fraction of residues allowed to be below the cutoff
        If int (>1): Maximum absolute number of residues allowed to be below the cutoff
    metapredict_version : int
        Version of MetaPredict to use (1, 2, or 3)
        default is 3
    
    Returns
    -------
    list
        List of sequences that are disordered
    """
    # if sequences are a numpy array, convert to list
    if isinstance(sequences, np.ndarray):
        sequences = sequences.tolist()
    # convert to list if str
    if isinstance(sequences, str):
        sequences = [sequences]
    # check if sequences is a dict
    if isinstance(sequences, dict):
        sequences = list(sequences.values())
    
    # predict disorder for sequences
    disorder_predictions = meta.predict_disorder(sequences, version=metapredict_version)
    
    # separate out sequences and predictions into their own lists.
    seqs = [seq[0] for seq in disorder_predictions]
    disorder_scores = np.array([pred[1] for pred in disorder_predictions])

    # if strict_disorder is True, check if all residues are above the cutoff
    if strict_disorder:
        disorder = np.all(disorder_scores > disorder_cutoff, axis=1)

    else:    
        N, seq_length = disorder_scores.shape
        
        # Create binary mask where True means the residue is ordered (below cutoff)
        ordered_mask = disorder_scores < disorder_cutoff
        
        # Check total number of ordered residues for each sequence
        total_ordered = np.sum(ordered_mask, axis=1)
        
        # Determine max allowed ordered residues
        if isinstance(max_total_ordered, float) and 0 <= max_total_ordered <= 1:
            # If it's a fraction, calculate absolute number
            max_ordered_allowed = int(max_total_ordered * seq_length)
        else:
            # Use as absolute number
            max_ordered_allowed = int(max_total_ordered)
        
        # Check if total ordered residues are within limit
        total_condition = total_ordered <= max_ordered_allowed
        
        # Initialize array for consecutive check
        consecutive_condition = np.ones(N, dtype=bool)
        
        # Use efficient sliding window approach for checking consecutive residues
        if max_consecutive_ordered < seq_length:
            max_consecutive_to_check = max_consecutive_ordered + 1
            
            # Check windows of size max_consecutive_ordered+1
            # If any window is all ordered, then we have too many consecutive ordered residues
            for i in range(N):
                # For each possible window of size max_consecutive_to_check
                for j in range(seq_length - max_consecutive_to_check + 1):
                    if np.all(ordered_mask[i, j:j+max_consecutive_to_check]):
                        consecutive_condition[i] = False
                        break
        
        # Combine conditions
        disorder = total_condition & consecutive_condition

    # return sequences that are disordered
    return [seqs[i] for i in range(len(seqs)) if disorder[i]]


# make function to generate sequences with specific properties
def generate_seq_by_props(length, 
                           fcr=None, 
                           ncpr=None, 
                           hydropathy=None, 
                           kappa=None, 
                           exclude_residues=None,
                           num_attempts=100,
                           strict_disorder=False,
                           hydropathy_tolerance = 0.05,
                           kappa_tolerance=0.03,
                           disorder_cutoff=0.5,
                           max_consecutive_ordered=3,
                           max_total_ordered=0.05,
                           metapredict_version=3,
                           return_all_sequences=False,
                           use_weighted_probabilities=True,
                           chosen_probabilities=None,
                           batch_size=1000):
    """
    Generate sequences with specific properties and check for disorder.
    
    Parameters
    ----------
    length : int
        Length of sequences to generate
    fcr : float
        FCR value to generate sequences with
    ncpr : float
        NCPR value to generate sequences with
    hydropathy : float
        Hydropathy value to generate sequences with
    kappa : float
        Kappa value to generate sequences with
    exclude_residues : list
        List of residues to exclude from the sequence
    num_attempts : int
        Number of attempts to make a sequence
    strict_disorder : bool
        If True, all residues must be above the disorder cutoff
        default is False
    hydropathy_tolerance : float
        Tolerance for hydropathy optimization
        default is 0.05
    kappa_tolerance : float
        Tolerance for kappa optimization
        default is 0.03
    disorder_cutoff : float
        Cutoff for disorder. Above this value is considered disordered.
    max_consecutive_ordered : int
        Maximum number of consecutive residues allowed to be below the disorder cutoff
    max_total_ordered : float or int
        If float (0-1): Maximum fraction of residues allowed to be below the cutoff
        If int (>1): Maximum absolute number of residues allowed to be below the cutoff
    metapredict_version : int
        Version of MetaPredict to use (1, 2, or 3)
        default is 3
    return_all_sequences : bool
        If True, return all sequences generated, not just the ones that pass the disorder check
        default is False
    use_weighted_probabilities : bool
        If True, use weighted probabilities for sequence generation. 
        If False, uniform probabilities are used. This is very slow. 
        default is True
    chosen_probabilities : dict
        Dictionary of probabilities for each amino acid. 
        Only used if use_weighted_probabilities is True.
    batch_size : int
        Number of sequences to generate in each batch
        default is 1000
    Returns
    -------
    list
        List of sequences that meet the criteria
    """

    # initialize the sequence generators
    seq_gen = SequenceGenerator( 
                 use_weighted_probabilities=use_weighted_probabilities)
    
    # iterate over the number of attempts
    for _ in range(num_attempts):
        # generate starter sequences with specified properties
        seqs = seq_gen.generate_sequences_vectorized(length, fcr=fcr, ncpr=ncpr, hydropathy=hydropathy, num_sequences=batch_size,
                                                     specific_probabilities=chosen_probabilities,
                                                     exclude_residues=exclude_residues)

        # if hydropathy is not None, we need to optimize hydropathy. This will only return seqs within tolerance.
        if hydropathy is not None:
            # set preserve charge to false as default, can change to true if we need to. 
            preserve_charge=False
            if fcr is not None:
                preserve_charge=True
            if ncpr is not None:
                preserve_charge=True
            seqs = optimize_hydropathy_vectorized(seqs, hydropathy,
                                                preserve_charged=preserve_charge, 
                                                tolerance=hydropathy_tolerance,
                                                return_when_num_hit=250,
                                                only_return_within_tolernace=True,
                                                exclude_residues=exclude_residues)
        # make sure we still have sequences
        if len(seqs) == 0:
            continue

        # cap number needed for kappa. 
        if len(seqs) < 25:
            num_for_kappa = len(seqs)
        else:
            num_for_kappa=25

        if len(seqs) > 250:
            seqs=seqs[:250]


        # see if we need to optimize kappa
        if kappa is not None:
            # get first successful kappa seq. This is needed because kappa is computationally intensive. 
            seqs = optimize_kappa_vectorized(seqs, kappa, return_when_num_hit=25, only_return_within_tolerance=True,
                                             tolerance=kappa_tolerance)

        # make sure we still have sequences
        if len(seqs) == 0:
            continue

        # finally, check for disorder.
        seqs = check_disorder_vectorized(seqs, strict_disorder=strict_disorder,
                                        disorder_cutoff=disorder_cutoff,
                                        max_consecutive_ordered=max_consecutive_ordered,
                                        max_total_ordered=max_total_ordered,
                                        metapredict_version=metapredict_version)
        
        if len(seqs) == 0:
            continue

        # if return_all_sequences, return seqs
        if return_all_sequences:
            return seqs
        else:
            # return a single sequence
            return random.choice(seqs)

    # if we get here, we didn't find any sequences that met the criteria
    raise goose_exceptions.GooseFail('Failed to generate sequence!')



# similar function for sequence fractions
def generate_seq_by_fractions(length,
                                fractions,
                                randomize_unspecified=False,
                                remaining_probabilities=None,
                                num_attempts=100,
                                strict_disorder=False,
                                disorder_cutoff=0.5,
                                max_consecutive_ordered=3,
                                max_total_ordered=0.05,
                                metapredict_version=3,
                                return_all_sequences=False,
                                batch_size=None):
    """"
    Generate sequences with specific fractions of amino acids and check for disorder."

    Parameters
    ----------
    length : int
        Length of sequences to generate
    fractions : dict
        Dictionary of fractions for each amino acid. 
        e.g. {'A': 0.1, 'C': 0.2, 'D': 0.3, 'E': 0.4}
    randomize_unspecified : bool
        If True, randomize the unspecified residues
        default is False
    remaining_probabilities : dict
        Dictionary of probabilities for the unspecified residues. 
        These are simply probabilities, not exact fractions. 
        Only used if randomize_unspecified is False.
    num_attempts : int
        Number of attempts to make a sequence
    strict_disorder : bool
        If True, all residues must be above the disorder cutoff
        default is False
    disorder_cutoff : float
        Cutoff for disorder. Above this value is considered disordered.
    max_consecutive_ordered : int
        Maximum number of consecutive residues allowed to be below the disorder cutoff
    max_total_ordered : float or int
        If float (0-1): Maximum fraction of residues allowed to be below the cutoff
        If int (>1): Maximum absolute number of residues allowed to be below the cutoff
    metapredict_version : int
        Version of MetaPredict to use (1, 2, or 3)
        default is 3
    return_all_sequences : bool
        If True, return all sequences generated, not just the ones that pass the disorder check
        default is False
    batch_size : int
        Number of sequences to generate in each batch
        default is None
    """
    # if batch size not specified, set it based on metapredict version
    if batch_size is None:
        if metapredict_version==3:
            batch_size=256
        else:
            batch_size=32

    # initialize the sequence generators
    seq_gen = FractionBasedSequenceGenerator(length=length,
                                             fractions=fractions,
                                                randomize_unspecified=randomize_unspecified,
                                                default_remaining_probabilities=remaining_probabilities)

    # iterate over the number of attempts
    for _ in range(num_attempts):
        # generate starter sequences with specified properties
        seqs = seq_gen.generate_sequences(num_sequences=batch_size)

        # finally, check for disorder.
        seqs = check_disorder_vectorized(seqs, strict_disorder=strict_disorder,
                                        disorder_cutoff=disorder_cutoff,
                                        max_consecutive_ordered=max_consecutive_ordered,
                                        max_total_ordered=max_total_ordered,
                                        metapredict_version=metapredict_version)
        
        if len(seqs) == 0:
            continue

        # if return_all_sequences, return seqs
        if return_all_sequences:
            return seqs
        else:
            # return a single sequence
            return random.choice(seqs)
    
    # if we get here, we didn't find any sequences that met the criteria
    raise goose_exceptions.GooseFail('Failed to generate sequence!')
