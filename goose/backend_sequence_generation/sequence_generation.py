'''
Version of sequence_generation.py but using the new vectorized functions.
The strategy here is to generate many sequences rapidly and then filter for ones we want.
'''
import random
import numpy as np
import metapredict as meta
from numpy.lib.stride_tricks import sliding_window_view
from goose import goose_exceptions
from goose.backend import parameters
from goose.backend_property_calculation.calculate_properties_batch import matrices_to_sequences
from goose.backend_property_calculation.calculate_kappa import kappa as calculate_kappa
from goose.backend_property_optimization.optimize_kappa import optimize_kappa
from goose.backend_property_optimization.optimize_hydropathy import optimize_hydropathy
from goose.backend_sequence_generation.sequence_by_probability import SequenceGenerator
from goose.backend_sequence_generation.sequence_by_fractions import FractionBasedSequenceGenerator
from goose.backend_sequence_generation.sequence_by_class import create_sequence_by_class
from goose.backend_sequence_generation.sequence_by_dimensions import create_seq_by_dims
from goose.backend_property_optimization.optimize_disorder import optimize_disorder



def check_disorder(sequences, 
                              strict_disorder=False,
                              disorder_cutoff=parameters.DISORDER_THRESHOLD,
                              max_consecutive_ordered=parameters.ALLOWED_CONSECUTIVE_ORDERED,
                              max_total_ordered=parameters.ALLOWED_TOTAL_ORDERED_FRACTION,
                              metapredict_version=parameters.METAPREDICT_DEFAULT_VERSION,
                              return_best_sequence=False):
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
    return_best_sequence : bool
        If True, return the sequence with the best disorder score
        default is False
    
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
        
        # Vectorized consecutive residue check
        if max_consecutive_ordered < seq_length:
            max_consecutive_to_check = max_consecutive_ordered + 1
            
            # Create sliding windows for all sequences at once using broadcasting
            # This creates a 3D array: (N, num_windows, window_size)
            num_windows = seq_length - max_consecutive_to_check + 1
            
            if num_windows > 0:
                # Create sliding windows for all sequences at once
                windowed_mask = sliding_window_view(ordered_mask, 
                                                  window_shape=max_consecutive_to_check, 
                                                  axis=1)
                
                # Check if any window in each sequence has all True values
                # windowed_mask shape: (N, num_windows, window_size)
                # np.all(axis=2) checks if all values in each window are True
                # np.any(axis=1) checks if any window in each sequence violates the condition
                has_consecutive_violation = np.any(np.all(windowed_mask, axis=2), axis=1)
                consecutive_condition = ~has_consecutive_violation
        
        # Combine conditions
        disorder = total_condition & consecutive_condition

    # return sequences that are disordered
    disordered_seqs=[seqs[i] for i in range(len(seqs)) if disorder[i]]

    # if no disorderd sequences and return_best_sequence is True,
    # return the sequence with the highest mean disorder score
    if return_best_sequence and len(disordered_seqs) == 0:
        # return the sequence with the highest mean disorder score
        best_index = np.argmax(np.mean(disorder_scores, axis=1))
        return seqs[best_index]
    
    # if disordered_seqs is empty, return None
    if len(disordered_seqs) == 0:
        return None
    
    # return disordered_seqs
    return disordered_seqs


# make function to generate sequences with specific properties
def by_properties(length, 
                           fcr=None, 
                           ncpr=None, 
                           hydropathy=None, 
                           kappa=None, 
                           exclude_residues=None,
                           num_attempts=5000,
                           strict_disorder=False,
                           hydropathy_tolerance = parameters.MAXIMUM_HYDRO_ERROR,
                           kappa_tolerance=parameters.MAXIMUM_KAPPA_ERROR,
                           disorder_cutoff=parameters.DISORDER_THRESHOLD,
                           max_consecutive_ordered=parameters.ALLOWED_CONSECUTIVE_ORDERED,
                           max_total_ordered=parameters.ALLOWED_TOTAL_ORDERED_FRACTION,
                           metapredict_version=parameters.METAPREDICT_DEFAULT_VERSION,
                           return_all_sequences=False,
                           use_weighted_probabilities=False,
                           chosen_probabilities=None,
                           batch_size=None,
                           check_sequence_disorder=True):
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
        set by parameters.HYDRO_ERROR
    kappa_tolerance : float
        Tolerance for kappa optimization
        default is 0.03
        set by parameters.MAXIMUM_KAPPA_ERROR
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
    check_sequence_disorder : bool
        If True, check the generated sequences for disorder
        default is True

    Returns
    -------
    list
        List of sequences that meet the criteria
    """
    # set batch size based on what is specified.
    if batch_size is None:
        # dynamic batch size. Basically slowly increases until it reaches a maximum size.
        # idea here is to increase the odds of a batch having a disordered sequence. However,
        # larger batches take longer, so if we can use a smaller batch size, we will.
        use_dynamic_batching = True
        dynamic_batch_sizes = [(10, 2, 1), (10, 2, 1), (10, 2, 1), (10, 2, 1),  (10, 2, 1), 
                             (20, 4, 1), (20, 6, 1), (20, 4, 1), (20, 4, 1),
                             (40, 8, 2), (40, 8, 2),(40, 8, 2),
                             (50, 10, 3),(50, 10, 3),
                             (100, 20, 6), 
                             (200, 40, 12),
                             (300, 60, 24), 
                             (400, 80, 48)]
    else:
        use_dynamic_batching = False
        if batch_size < 1:
            raise goose_exceptions.GooseFail('Batch size must be at least 1!')
        # if batch size is specified, use it.
        if batch_size < 2:
            required_hydro_batch_size = 1
            required_kappa_batch_size = 1
        else:
            required_hydro_batch_size = int(batch_size/10)
            if required_hydro_batch_size < 2:
                required_hydro_batch_size = 2
            required_kappa_batch_size = int(required_hydro_batch_size/10)
            if required_kappa_batch_size < 1:
                required_kappa_batch_size = 1

    # initialize the sequence generators
    seq_gen = SequenceGenerator( 
                 use_weighted_probabilities=use_weighted_probabilities)
    

    # iterate over the number of attempts
    batch_to_use = 0
    for _ in range(num_attempts):
        if use_dynamic_batching:
            # set batch size based on the number of attempts
            if batch_to_use < len(dynamic_batch_sizes):
                batch_size, required_hydro_batch_size, required_kappa_batch_size = dynamic_batch_sizes[batch_to_use]
            else:
                batch_size, required_hydro_batch_size, required_kappa_batch_size = dynamic_batch_sizes[-1]
            batch_to_use += 1
        # generate starter sequences with specified properties
        seqs = seq_gen.generate_sequences_vectorized(length, fcr=fcr, ncpr=ncpr, 
                                                     hydropathy=hydropathy, num_sequences=batch_size,
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
            seqs = optimize_hydropathy(seqs, hydropathy,
                                                preserve_charged=preserve_charge, 
                                                hydropathy_tolerance=hydropathy_tolerance,
                                                return_when_num_hit=required_hydro_batch_size,
                                                only_return_within_tolernace=True,
                                                exclude_residues=exclude_residues)
        # make sure we still have sequences
        if len(seqs) == 0:
            continue

        # see if we need to optimize kappa
        if kappa is not None and kappa != -1:
            # set charge placement variable 
            preserve_charge_placement = True
            # make sure kappa is not -1
            kappa_values = calculate_kappa(seqs)
            # filter out sequences with kappa -1
            valid_indices = np.where(kappa_values != -1)[0]
            # get valid sequences
            seqs = seqs[valid_indices]

            # iterate over sequences that have made it this far. 
            # get first successful kappa seq. This is needed because kappa is computationally intensive. 
            seqs = optimize_kappa(seqs, kappa, 
                                            return_when_num_hit=required_kappa_batch_size, 
                                            only_return_within_tolerance=True,
                                            kappa_tolerance=kappa_tolerance)
            
            # make sure we still have sequences
            if len(seqs) == 0:
                continue

        else:
            # if we are not optimizing kappa, we need to convert sequences back from matrices. 
            seqs = matrices_to_sequences(seqs)
            # set charge preservation for disorder optimization
            preserve_charge_placement = False

        # if we are checking disorder, do that. 
        if check_sequence_disorder:
            seqs = check_disorder(seqs, strict_disorder=strict_disorder,
                                        disorder_cutoff=disorder_cutoff,
                                        max_consecutive_ordered=max_consecutive_ordered,
                                        max_total_ordered=max_total_ordered,
                                        metapredict_version=metapredict_version,
                                        return_best_sequence=True)
            
            # if we still have a string, try to optimize it
            if isinstance(seqs, str):
                # only use optimization when we have failed to make a sequence consistently. 
                if batch_size >= 100:
                    seqs=optimize_disorder(seqs, disorder_cutoff=disorder_cutoff,
                            max_iterations=50, 
                            preserve_charge_placement=preserve_charge_placement,
                            metapredict_version=metapredict_version)
                    # check if the optimized sequence meets the disorder cutoffa
                    cur_disorder = meta.predict_disorder(seqs, version=metapredict_version)
                    if np.min(cur_disorder) >= disorder_cutoff:
                        return seqs  
                else:
                    continue
            else:
                if seqs == []:
                    continue
                else:
                    # if return_all_sequences, return seqs
                    if return_all_sequences:
                        return seqs
                    else:
                        # return a single sequence
                        return random.choice(seqs)
        else:
            # return all the sequences. 
            return seqs

    # if we get here, we didn't find any sequences that met the criteria
    return None



# similar function for sequence fractions
def by_fractions(length,
                                fractions,
                                randomize_unspecified=False,
                                remaining_probabilities=None,
                                num_attempts=100,
                                strict_disorder=False,
                                disorder_cutoff=parameters.DISORDER_THRESHOLD,
                                max_consecutive_ordered=parameters.ALLOWED_CONSECUTIVE_ORDERED,
                                max_total_ordered=parameters.ALLOWED_TOTAL_ORDERED_FRACTION,
                                metapredict_version=parameters.METAPREDICT_DEFAULT_VERSION,
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
        # use dynamic batching
        use_dynamic_batching = True
        dynamic_batch_sizes = [1, 1, 1, 1, 1,
                               2, 2, 2, 2,
                               4, 4, 4, 8, 8, 
                               16, 32, 64, 128, 256]
    else:
        use_dynamic_batching = False
        if batch_size < 1:
            raise goose_exceptions.GooseFail('Batch size must be at least 1!')
        
    if remaining_probabilities is not None:
        # make sure randomize_unspecified is false
        randomize_unspecified = False

    # initialize the sequence generators
    seq_gen = FractionBasedSequenceGenerator(length=length,
                                             fractions=fractions,
                                                randomize_unspecified=randomize_unspecified,
                                                default_remaining_probabilities=remaining_probabilities)

    cur_batch_num=0
    # iterate over the number of attempts
    for _ in range(num_attempts):
        if use_dynamic_batching:
            # set batch size based on the number of attempts
            cur_batch_size = dynamic_batch_sizes[cur_batch_num]
            if cur_batch_num < len(dynamic_batch_sizes) - 1:
                cur_batch_num += 1
        else:
            cur_batch_size = batch_size
    
        # generate starter sequences with specified properties
        seqs = seq_gen.generate_sequences(num_sequences=cur_batch_size)


        # finally, check for disorder.
        seqs = check_disorder(seqs, strict_disorder=strict_disorder,
                                        disorder_cutoff=disorder_cutoff,
                                        max_consecutive_ordered=max_consecutive_ordered,
                                        max_total_ordered=max_total_ordered,
                                        metapredict_version=metapredict_version,
                                        return_best_sequence=True)

        # the check_disordered_vectorized only returns a string if the sequence doesn't make
        # the cutoff for bieng disordered. Try to optimize it. 
        # only attempt to optimize if we have a batch size larger than 128.
        if isinstance(seqs, str):
            if cur_batch_size > 128:
                # means we got back a single sequence as a string and therefore have 
                # a sequence not disordered. Try to optimize. 
                # first try making 100 shuffles of the sequence we got back and checking disorder for that.
                # this is a bit of a hack, but it works.
                shuff_seqs = [seqs]
                seqs=list(seqs) # convert to list for consistency
                # shuffle the sequence 256 times
                for _ in range(256):
                    # shuffle the sequence
                    random.shuffle(seqs)
                    shuff_seqs.append(''.join(seqs))

                # check disorder
                seqs = check_disorder(shuff_seqs, strict_disorder=strict_disorder,
                                            disorder_cutoff=disorder_cutoff,
                                            max_consecutive_ordered=max_consecutive_ordered,
                                            max_total_ordered=max_total_ordered,
                                            metapredict_version=metapredict_version,
                                            return_best_sequence=True)
                
                # if we still have a string, try to optimize it
                if isinstance(seqs, str):
                    seqs=optimize_disorder(seqs, disorder_cutoff=disorder_cutoff,
                            max_iterations=500, preserve_charge_placement=False,
                            metapredict_version=metapredict_version)
                    # check if the optimized sequence meets the disorder cutoffa
                    cur_disorder = meta.predict_disorder(seqs, version=metapredict_version)
                    if np.min(cur_disorder) >= disorder_cutoff:
                        return seqs 
                    else:
                        # if we still have a string, we failed to make a sequence that meets the disorder cutoff
                        continue
            else:
                continue
        
        # if we have a list of sequences, check if they are empty
        if seqs == []:
            continue

        if seqs is None:
            # if we got back None, we failed to make a sequence that meets the disorder cutoff
            continue
        
        # if return_all_sequences, return seqs
        if return_all_sequences:
            return seqs
        else:
            # return a single sequence
            if isinstance(seqs, list):
                if len(seqs) == 0:
                    continue
            # if we have a list of sequences, return a random one
            if isinstance(seqs, list):
                # if we have a string, return it as a list
                seqs = random.choice(seqs)
            # return a random sequence from the list
            return seqs
    
    # if we get here, we didn't find any sequences that met the criteria
    return None


    
def by_class(length: int,
                            aromatic_fraction: float = None,
                            aliphatic_fraction: float = None,
                            polar_fraction: float = None,
                            positive_fraction: float = None,
                            negative_fraction: float = None,
                            glycine_fraction: float = None,
                            proline_fraction: float = None,
                            cysteine_fraction: float = None,
                            histidine_fraction: float = None,
                            num_attempts=100, strict_disorder=False,
                            disorder_cutoff=parameters.DISORDER_THRESHOLD,
                            metapredict_version=parameters.METAPREDICT_DEFAULT_VERSION,
                            max_consecutive_ordered=parameters.ALLOWED_CONSECUTIVE_ORDERED,
                            max_total_ordered=parameters.ALLOWED_TOTAL_ORDERED_FRACTION,
                            remaining_probabilities=None):
    """
    Generate a sequence of a specified length with specific amino acid class fractions.
    Non-specified classes will be randomly filled in. 
    
    Parameters
    ----------
    seq_length : int
        Length of the sequence to generate
    aromatic_fraction : float, default=None
        Fraction of aromatic amino acids (F, W, Y)
    aliphatic_fraction : float, default=None
        Fraction of aliphatic amino acids (A, I, L, V)
    polar_fraction : float, default=None
        Fraction of polar amino acids (N, Q, S, T)
    positive_fraction : float, default=None
        Fraction of positively charged amino acids (K, R)
    negative_fraction : float, default=None
        Fraction of negatively charged amino acids (D, E)
    glycine_fraction : float, default=None  
        Fraction of glycine (G)
    proline_fraction : float, default=None
        Fraction of proline (P)
    cysteine_fraction : float, default=None
        Fraction of cysteine (C)
    histidine_fraction : float, default=None
        Fraction of histidine (H)
    num_attempts : int, default=100
        Number of attempts to make the sequence
    strict_disorder : bool, default=False
        If True, applies strict disorder checks using MetaPredict
    disorder_cutoff : float, default=parameters.DISORDER_THRESHOLD
        Cutoff for disorder. Above this value is considered disordered.
    metapredict_version : int, default=parameters.METAPREDICT_DEFAULT_VERSION
        Version of MetaPredict to use (1, 2, or 3),
    max_consecutive_ordered : int, default=3
        Maximum number of consecutive residues allowed to be below the disorder cutoff
    max_total_ordered : float or int, default=0.05
        If float (0-1): Maximum fraction of residues allowed to be below the cutoff
        If int (>1): Maximum absolute number of residues allowed to be below the cutoff
    remaining_probabilities : dict, optional
        Dictionary of probabilities for the unspecified amino acid classes.
        If not provided, the unspecified classes will be filled randomly.

    Returns
    -------
    str or None
        Generated sequence that meets the target amino acid class fractions and is disordered.
        If no suitable sequence is found, returns None.
    """
    for i in range(num_attempts):
        # Call the create_sequence_by_class function to
        # generate the sequence with specified class fractions
        seq = create_sequence_by_class(length,
                            aromatic_fraction=aromatic_fraction,
                            aliphatic_fraction=aliphatic_fraction,
                            polar_fraction=polar_fraction,
                            positive_fraction=positive_fraction,
                            negative_fraction=negative_fraction,
                            glycine_fraction=glycine_fraction,
                            proline_fraction=proline_fraction,
                            cysteine_fraction=cysteine_fraction,
                            histidine_fraction=histidine_fraction,
                            num_sequences=32,
                            remaining_probabilities=remaining_probabilities)
        
        # check disorder
        disordered_seqs = check_disorder(seq, strict_disorder=strict_disorder,
                                                        disorder_cutoff=disorder_cutoff,
                                                        max_consecutive_ordered=max_consecutive_ordered,
                                                        max_total_ordered=max_total_ordered,
                                                        metapredict_version=metapredict_version)
        # if we have a disordered sequence, return it
        if disordered_seqs is not None:
            return disordered_seqs[0]
    
    # if we get here, we didn't find any sequences that met the criteria
    return None



def by_dimensions(seq_length, objective_dim, rg_or_re='rg',
                       allowed_error=parameters.MAXIMUM_RG_RE_ERROR,
                       reduce_pos_charged=True, exclude_aas=None,
                       variants_per_iteration=64, mutation_fraction=0.0125,
                       num_attempts=100, strict_disorder=False,
                       disorder_cutoff=parameters.DISORDER_THRESHOLD,
                       metapredict_version=parameters.METAPREDICT_DEFAULT_VERSION,
                       max_consecutive_ordered=parameters.ALLOWED_CONSECUTIVE_ORDERED,
                       max_total_ordered=parameters.ALLOWED_TOTAL_ORDERED_FRACTION,
                       num_attempts_dimensions=parameters.RG_RE_ATTEMPT_NUMBER,):
    """
    Create a sequence of a specified length that meets a target radius of gyration (Rg)
    or end-to-end distance (Re) with a given error tolerance.

    Parameters
    ----------
    seq_length : int
        Length of the sequence to generate
    objective_dim : float
        Target value for the specified dimensional property
    rg_or_re : {'rg', 're'}
        Target property: 'rg' for radius of gyration, 're' for end-to-end distance
    allowed_error : float or 'default_error'
        Maximum allowed deviation from target. If 'default_error', uses
        backend parameter defaults (re_error or rg_error)
    num_attempts_dimensions : int, default=rg_re_attempt_num
        Maximum number of optimization iterations
    reduce_pos_charged : bool, default=True
        Whether to reduce positively charged residues (K, R) in the sequence.
        Based on in vivo data suggesting positive charges may not drive expansion   
        as predicted by current models.
    exclude_aas : list of str, optional
        Amino acids to exclude from the optimization process
    variants_per_iteration : int, default=64
        Number of sequence variants to generate per optimization iteration
    mutation_fraction : float, default=0.0125
        Fraction of sequence length to mutate per iteration (minimum 1 residue)
    num_attempts : int, default=10
        Number of iterations to attempt to make the sequence starting from random. 
    strict_disorder : bool, default=False
        If True, applies strict disorder checks using MetaPredict
    disorder_cutoff : float, default=parameters.DISORDER_THRESHOLD
        Cutoff for disorder. Above this value is considered disordered.
    metapredict_version : int, default=parameters.METAPREDICT_DEFAULT_VERSION
        Version of MetaPredict to use (1, 2, or 3)
    max_consecutive_ordered : int, default=3
        Maximum number of consecutive residues allowed to be below the disorder cutoff
    max_total_ordered : float or int, default=0.05
        If float (0-1): Maximum fraction of residues allowed to be below the cutoff
        If int (>1): Maximum absolute number of residues allowed to be below the cutoff

    Returns
    -------
    str or None
        Generated sequence that meets the target dimensional property and is disordered.
        If no suitable sequence is found, returns None.
    """
    for _ in range(num_attempts):
        # Call the create_seq_by_dims function to generate the sequence
        seq = create_seq_by_dims(seq_length, objective_dim, rg_or_re=rg_or_re,
                                allowed_error=allowed_error,
                                num_attempts_dimensions=num_attempts_dimensions,
                                reduce_pos_charged=reduce_pos_charged,
                                exclude_aas=exclude_aas,
                                variants_per_iteration=variants_per_iteration,
                                mutation_fraction=mutation_fraction)
        
        # note: the create_seq_dims function retruns a list or None. Therefore, we don't need
        # to wrap seq in a list when predicting disorder.
        # if seq is None, we didn't find a sequence that meets the criteria
        # if seq is not None, check disorder
        if seq is not None:
            # check disorder of the sequence
            disordered_seqs = check_disorder(seq, strict_disorder=strict_disorder,
                                                        disorder_cutoff=disorder_cutoff,
                                                        max_consecutive_ordered=max_consecutive_ordered,
                                                        max_total_ordered=max_total_ordered,
                                                        metapredict_version=metapredict_version)
            # if we have a disordered sequence, return it
            if disordered_seqs is not None:
                return disordered_seqs[0]