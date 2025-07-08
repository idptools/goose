'''
Uses the variant_generation_functions to make a sequence. Important thing here is that
it also checks for disorder. These things are separated to keep things faster and more modular.
'''
from sparrow.protein import Protein
from goose.backend_sequence_generation.sequence_generation_vectorized import generate_seq_by_props
from goose.backend_variant_generation.helper_functions import check_variant_disorder_vectorized
from goose.backend_property_optimization.optimize_dimensions import predict_rg, predict_re
from goose.backend_variant_generation import variant_generation_functions as vgf
from goose.backend import parameters


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
        else:
            # if we have a disordered sequence, we can break out of the loop
            break
    
    # return the disordered sequence. This could be None if no valid sequence was found
    # but we will handle that exception in the frontend to avoid an annoying
    # string of exceptions. 
    if isinstance(disordered_sequence, list):
        # if the disordered_sequence is a list, we will return the first one
        disordered_sequence = disordered_sequence[0]
    return disordered_sequence



# code for variant generation with disorder checking
def gen_minimal_variant(input_sequence,
                            target_hydropathy: float = None,
                            target_kappa: float = None,
                            target_FCR: float = None,
                            target_NCPR: float = None,
                            num_attempts=5,
                            hydropathy_tolerance: float = parameters.HYDRO_ERROR,
                            kappa_tolerance: float = parameters.MAXIMUM_KAPPA_ERROR,
                            strict_disorder=False,
                            disorder_cutoff=0.5,
                            metapredict_version=3):
    """
    Generate a constant class variant of the input sequence.
    Parameters
    ----------
    input_sequence : str
        The input sequence to generate a variant from.
    target_hydropathy : float
        Target mean hydropathy value to achieve.
        If None, hydropathy optimization is not performed.
    target_kappa : float
        Target kappa value to achieve.
        If None, kappa optimization is not performed.
    target_FCR : float
        Target FCR value to achieve.
        If None, FCR optimization is not performed.
    target_NCPR : float
        Target NCPR value to achieve.
        If None, NCPR optimization is not performed.
    num_attempts : int
        Number of attempts to generate a variant.
        default is 5
    hydropathy_tolerance : float
        Acceptable difference between achieved and target hydropathy.
        default is parameters.HYDRO_ERROR
    kappa_tolerance : float
        Acceptable difference between achieved and target kappa.
        default is parameters.MAXIMUM_KAPPA_ERROR
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
        variant_sequence = vgf.generate_minimal_variant(input_sequence,
                                                        target_hydropathy=target_hydropathy,
                                                        target_kappa=target_kappa,
                                                        target_FCR=target_FCR,
                                                        target_NCPR=target_NCPR,
                                                        hydropathy_tolerance=hydropathy_tolerance,
                                                        kappa_tolerance=kappa_tolerance)
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
        else:
            # if we have a disordered sequence, we can break out of the loop
            break
    
    # make sure we return a sequence
    if isinstance(disordered_sequence, list):
        # if the disordered_sequence is a list, we will return the first one
        disordered_sequence = disordered_sequence[0]
    
    return disordered_sequence

def gen_region_shuffle_variant(input_sequence: str,
                               shuffle_regions: list,
                               num_attempts: int = 5,
                                strict_disorder=False,
                                disorder_cutoff=0.5,
                                metapredict_version=3                               
                               ) -> str:
    """    Generate a variant sequence by shuffling specified regions of the input sequence.
    
    Parameters
    ----------
    input_sequence : str
        The input sequence to generate a variant from.
    shuffle_regions : list
        List of tuples specifying regions to shuffle. Each tuple should contain
        (start_index, end_index) for the region to shuffle.
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
        A generated variant sequence with specified regions shuffled.
    """
    for _ in range(num_attempts):
        # generate multiple variants per round because these are easy to generate
        shuffle_vars=[]
        for sub_loop in range(10):
            shuffle_vars.append(vgf.generate_region_shuffle_variant(input_sequence, shuffle_regions))

        # Check if the variant sequence meets the disorder criteria
        disordered_sequence = check_variant_disorder_vectorized(
            original_sequence=input_sequence,
            sequences=shuffle_vars,
            strict_disorder=strict_disorder,
            disorder_cutoff=disorder_cutoff,
            metapredict_version=metapredict_version,
            return_best_sequence=False
        )
        if disordered_sequence is None:
            continue
        else:
            # if we have a disordered sequence, we can break out of the loop
            break
    
    # make sure we return a sequence
    if isinstance(disordered_sequence, list):
        # if the disordered_sequence is a list, we will return the first one
        disordered_sequence = disordered_sequence[0]

    return disordered_sequence

def gen_excluded_shuffle_variant(input_sequence: str,
                                 excluded_residues: list,
                                 num_attempts: int = 5,
                                 strict_disorder=False,
                                 disorder_cutoff=0.5,
                                 metapredict_version=3) -> str:
    """
    Generate a variant sequence by shuffling residues not specified in excluded_residues.

    Parameters
    ----------
    input_sequence : str
        The input sequence to generate a variant from.
    excluded_residues : list
        List of residues to exclude from shuffling.
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
        A generated variant sequence with specified residues excluded from shuffling.
    """
    for _ in range(num_attempts):
        # Generate a variant sequence
        shuff_seqs=[]
        for sub_loop in range(10):
            # generate multiple variants per round because these are easy to generate
            shuff_seqs.append(vgf.generate_excluded_shuffle_variant(input_sequence, 
                                                                    excluded_residues=excluded_residues))
        
        # check disroder
        disordered_sequence = check_variant_disorder_vectorized(
            original_sequence=input_sequence,
            sequences=shuff_seqs,
            strict_disorder=strict_disorder,
            disorder_cutoff=disorder_cutoff,
            metapredict_version=metapredict_version,
            return_best_sequence=False
        )
        if disordered_sequence is None:
            continue
        else:
            # if we have a disordered sequence, we can break out of the loop
            break
    
    # make sure we return a sequence
    if isinstance(disordered_sequence, list):
        # if the disordered_sequence is a list, we will return the first one
        disordered_sequence = disordered_sequence[0]

    return disordered_sequence


def gen_targeted_shuffle_variant(input_sequence: str,
                                 target_residues: list,
                                 num_attempts: int = 5,
                                 strict_disorder=False,
                                 disorder_cutoff=0.5,
                                 metapredict_version=3) -> str:
    '''
    Generate a variant sequence by shuffling residues within a specified target residues.

    Parameters
    ----------
    input_sequence : str
        The input sequence to generate a variant from.
    target_residues : list
        List of residues to shuffle within the sequence.
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
        A generated variant sequence with specified residues shuffled.  
    '''
    for _ in range(num_attempts):
        # Generate a variant sequence
        shuff_seqs=[]
        for sub_loop in range(10):
            # generate multiple variants per round because these are easy to generate
            shuff_seqs.append(vgf.generate_targeted_shuffle_variant(input_sequence, target_residues))
        
        # check disroder
        disordered_sequence = check_variant_disorder_vectorized(
            original_sequence=input_sequence,
            sequences=shuff_seqs,
            strict_disorder=strict_disorder,
            disorder_cutoff=disorder_cutoff,
            metapredict_version=metapredict_version,
            return_best_sequence=False
        )
        if disordered_sequence is None:
            continue
        else:
            # if we have a disordered sequence, we can break out of the loop
            break
    
    # make sure we return a sequence
    if isinstance(disordered_sequence, list):
        # if the disordered_sequence is a list, we will return the first one
        disordered_sequence = disordered_sequence[0]

    return disordered_sequence    

def gen_new_var_constant_class(input_sequence: str,
                               num_attempts: int = 50,
                               kappa_tolerance: float = parameters.MAXIMUM_KAPPA_ERROR,
                               hydropathy_tolerance: float = parameters.HYDRO_ERROR,
                               strict_disorder=False,
                               disorder_cutoff=0.5,
                               metapredict_version=3) -> str:
    """
    Generate a new variant sequence by constant class mutation.

    Parameters
    ----------
    input_sequence : str
        The input sequence to generate a variant from.
    num_attempts : int
        Number of attempts to generate a variant.
        default is 5
    kappa_tolerance : float
        Acceptable difference between achieved and target kappa.
        default is parameters.MAXIMUM_KAPPA_ERROR
    hydropathy_tolerance : float
        Acceptable difference between achieved and target hydropathy.
        default is parameters.HYDRO_ERROR
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
        variant_sequence = vgf.generate_new_seq_constant_class_variant(input_sequence,
                                                          kappa_tolerance=kappa_tolerance,
                                                          hydropathy_tolerance=hydropathy_tolerance)
        
        # if we fail at the variant generation step, try again. 
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
        # if we fail at the disorder check, try again
        if disordered_sequence is None:
            continue
        else:
            # if we have a disordered sequence, we can break out of the loop
            break   

    # make sure we return a sequence
    if isinstance(disordered_sequence, list):
        # if the disordered_sequence is a list, we will return the first one
        disordered_sequence = disordered_sequence[0]

    return disordered_sequence

def gen_new_variant(input_sequence: str,
                    num_attempts: int = 5,
                    hydropathy_tolerance: float = parameters.HYDRO_ERROR,
                    kappa_tolerance: float = parameters.MAXIMUM_KAPPA_ERROR,
                    strict_disorder=False,
                    disorder_cutoff=0.5,
                    metapredict_version=3,
                    exclude_residues: list = None
                    ) -> str:
    """
    Generate a new variant where the sequence is constrained by hyrdropathy, 
    NCPR, FCR, and kappa values.
    
    Parameters
    ----------
    input_sequence : str
        The input sequence to generate a variant from.
    num_attempts : int
        Number of attempts to generate a variant.
        default is 5
    hydropathy_tolerance : float
        Acceptable difference between achieved and target hydropathy.
        default is parameters.HYDRO_ERROR
    kappa_tolerance : float
        Acceptable difference between achieved and target kappa.
        default is parameters.MAXIMUM_KAPPA_ERROR
    strict_disorder : bool
        If True, all residues must be above the disorder cutoff to be considered disordered.
        default is False
    disorder_cutoff : float
        Cutoff for disorder. Above this value is considered disordered.
        default is 0.5
    metapredict_version : int
        Version of MetaPredict to use (1, 2, or 3)
        default is 3
    exclude_residues : list
        List of residues to exclude from the variant generation.
        If None, no residues are excluded.
        default is None

    Returns
    -------
    str
        A generated variant sequence that meets the disorder criteria.
    """

    # no backend function for this because it just uses the generate_seq_by_props functionality. 
    # get starting properties
    protein = Protein(input_sequence)
    original_hydropathy = protein.hydrophobicity
    original_kappa = protein.kappa
    original_FCR = protein.FCR
    original_NCPR = protein.NCPR
    length = len(input_sequence)

    for i in range(num_attempts):
        seq = generate_seq_by_props(
                            length, 
                            fcr=original_FCR, 
                            ncpr=original_NCPR, 
                            hydropathy=original_hydropathy, 
                            kappa=original_kappa, 
                            exclude_residues=exclude_residues,
                            num_attempts=5000,
                            hydropathy_tolerance = hydropathy_tolerance,
                            kappa_tolerance=kappa_tolerance,
                            metapredict_version=metapredict_version,
                            return_all_sequences=True,
                            check_disorder=False,
                            batch_size=200)
        if seq is None:
            continue
        
        # Check if the variant sequence meets the disorder criteria
        disordered_sequence = check_variant_disorder_vectorized(
            original_sequence=input_sequence,
            sequences=seq,
            strict_disorder=strict_disorder,
            disorder_cutoff=disorder_cutoff,
            metapredict_version=metapredict_version,
            return_best_sequence=False
        )
        if disordered_sequence is None:
            continue
        else:
            # if we have a disordered sequence, we can break out of the loop
            break
    
    # make sure we return a sequence
    if isinstance(disordered_sequence, list):
        # if the disordered_sequence is a list, we will return the first one
        disordered_sequence = disordered_sequence[0]

    return disordered_sequence

def gen_constant_residue_variant(input_sequence: str,
                                 constant_residues: list,
                                 hydropathy_tolerance: float = parameters.HYDRO_ERROR,
                                 kappa_tolerance: float = parameters.MAXIMUM_KAPPA_ERROR,
                                 num_attempts: int = 50,
                                 strict_disorder=False,
                                 disorder_cutoff=0.5,
                                 metapredict_version=3) -> str:
    """
    Generate a variant sequence by keeping specified residues constant.
    Hydropathy, kappa, ncpr and fcr are also held constant. 

    Parameters
    ----------
    input_sequence : str
        The input sequence to generate a variant from.
    constant_residues : list
        List of residues to keep constant in the variant sequence.
    hydropathy_tolerance : float
        Acceptable difference between achieved and target hydropathy.
        default is parameters.HYDRO_ERROR
    kappa_tolerance : float
        Acceptable difference between achieved and target kappa.
        default is parameters.MAXIMUM_KAPPA_ERROR
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
        variant_sequence = vgf.generate_constant_residue_variant(
            input_sequence,
            constant_residues=constant_residues,
            hydropathy_tolerance=hydropathy_tolerance,
            kappa_tolerance=kappa_tolerance
        )
        
        # if we fail at the variant generation step, try again. 
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
        # if we fail at the disorder check, try again
        if disordered_sequence is None:
            continue
        else:
            # if we have a disordered sequence, we can break out of the loop
            break
    
    # make sure we return a sequence
    if isinstance(disordered_sequence, list):
        # if the disordered_sequence is a list, we will return the first one
        disordered_sequence = disordered_sequence[0]

    return disordered_sequence

def gen_asymmetry_variant(input_sequence: str,
                          target_residues: list,
                          increase_or_decrease: str = 'increase',
                          num_attempts: int = 50,
                          strict_disorder=False,
                          disorder_cutoff=0.5,
                          metapredict_version=3) -> str:
    """
    Generate a variant sequence by introducing asymmetry in specified residues.

    Parameters
    ----------
    input_sequence : str
        The input sequence to generate a variant from.
    target_residues : list
        List of residues to introduce asymmetry in.
    increase_or_decrease : str
        Whether to increase or decrease the asymmetry.
        default is 'increase'
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
        A generated variant sequence with specified residues made asymmetric.
    """
    for _ in range(num_attempts):
        # Generate a variant sequence
        variant_sequence = vgf.generate_asymmetry_variant(
            input_sequence,
            target_residues=target_residues,
            increase_or_decrease=increase_or_decrease
        )

        # if we fail at the variant generation step, try again. 
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
        # if we fail at the disorder check, try again
        if disordered_sequence is None:
            continue
        else:
            # if we have a disordered sequence, we can break out of the loop
            break
    
    # make sure we return a sequence
    if isinstance(disordered_sequence, list):
        # if the disordered_sequence is a list, we will return the first one
        disordered_sequence = disordered_sequence[0]

    return disordered_sequence


def gen_hydropathy_class_variant(input_sequence,
                                 target_hydropathy: float,
                                 num_attempts: int = 5,
                                 hydropathy_tolerance: float = parameters.HYDRO_ERROR,
                                 strict_disorder=False,
                                 disorder_cutoff=0.5,
                                 metapredict_version=3):
    """
    Generate a hydropathy class variant of the input sequence.

    Parameters
    ----------
    input_sequence : str
        The input sequence to generate a variant from.
    target_hydropathy : float
        Target mean hydropathy value to achieve.
    num_attempts : int
        Number of attempts to generate a variant.
        default is 5
    hydropathy_tolerance : float
        Acceptable difference between achieved and target hydropathy.
        default is parameters.HYDRO_ERROR
    strict_disorder : bool
        If True, all residues must be above the disorder cutoff to be considered disordered.
        default is False
    disorder_cutoff : float
        Cutoff for disorder. Above  this value is considered disordered.
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
        variant_sequence = vgf.generate_hydro_class_variant(
            input_sequence,
            target_hydropathy=target_hydropathy,
            tolerance=hydropathy_tolerance)
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
        else:
            # if we have a disordered sequence, we can break out of the loop
            break
    
    # make sure we return a sequence
    if isinstance(disordered_sequence, list):
        # if the disordered_sequence is a list, we will return the first one
        disordered_sequence = disordered_sequence[0]

    return disordered_sequence


def gen_fcr_class_variant(input_sequence: str,
                          target_FCR: float,
                          hydropathy_tolerance: float = parameters.HYDRO_ERROR,
                          kappa_tolerance: float = parameters.MAXIMUM_KAPPA_ERROR,
                          num_attempts: int = 5,
                          strict_disorder=False,
                          disorder_cutoff=0.5,
                          metapredict_version=3) -> str:
    """
    Generate a FCR class variant of the input sequence.

    Parameters
    ----------
    input_sequence : str
        The input sequence to generate a variant from.
    target_FCR : float
        Target FCR value to achieve.
    hydropathy_tolerance : float
        Acceptable difference between achieved and target hydropathy.
        default is parameters.HYDRO_ERROR
    kappa_tolerance : float
        Acceptable difference between achieved and target kappa.
        default is parameters.MAXIMUM_KAPPA_ERROR
    num_attempts : int
        Number of attempts to generate a variant.
        default is 5
    strict_disorder : bool
        If True, all residues must be   above the disorder cutoff to be considered disordered.  
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
        variant_sequence = vgf.generate_fcr_class_variant(
            input_sequence,
            target_FCR=target_FCR,
            hydropathy_tolerance=hydropathy_tolerance,
            kappa_tolerance=kappa_tolerance)
        
        # if we fail at the variant generation step, try again. 
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
        # if we fail at the disorder check, try again
        if disordered_sequence is None:
            continue
        else:
            # if we have a disordered sequence, we can break out of the loop
            break
    
    # make sure we return a sequence
    if isinstance(disordered_sequence, list):
        # if the disordered_sequence is a list, we will return the first one
        disordered_sequence = disordered_sequence[0]

    return disordered_sequence

def gen_ncpr_class_variant(input_sequence: str,
                          target_NCPR: float,
                          hydropathy_tolerance: float = parameters.HYDRO_ERROR,
                          kappa_tolerance: float = parameters.MAXIMUM_KAPPA_ERROR,
                          num_attempts: int = 50,
                          strict_disorder=False,
                          disorder_cutoff=0.5,
                          metapredict_version=3) -> str:
    """
    Generate a NCPR class variant of the input sequence.

    Parameters
    ----------
    input_sequence : str
        The input sequence to generate a variant from.
    target_NCPR : float
        Target NCPR value to achieve.
    hydropathy_tolerance : float
        Acceptable difference between achieved and target hydropathy.
        default is parameters.HYDRO_ERROR
    kappa_tolerance : float
        Acceptable difference between achieved and target kappa.
        default is parameters.MAXIMUM_KAPPA_ERROR
    num_attempts : int
        Number of attempts to generate a variant.
        default is 5
    strict_disorder : bool
        If True, all residues must be above the disorder cutoff to be considered disordered.
        default is  False
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
        variant_sequence = vgf.generate_ncpr_class_variant(
            input_sequence,
            target_NCPR=target_NCPR,
            hydropathy_tolerance=hydropathy_tolerance,
            kappa_tolerance=kappa_tolerance)
        
        # if we fail at the variant generation step, try again. 
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
        # if we fail at the disorder check, try again
        if disordered_sequence is None:
            continue
        else:
            # if we have a disordered sequence, we can break out of the loop
            break
    
    # make sure we return a sequence
    if isinstance(disordered_sequence, list):
        # if the disordered_sequence is a list, we will return the first one
        disordered_sequence = disordered_sequence[0]
    
    return disordered_sequence


def gen_kappa_variant(input_sequence: str,
                      target_kappa: float,
                      kappa_tolerance: float = parameters.MAXIMUM_KAPPA_ERROR,
                      num_attempts: int = 5,
                        strict_disorder=False,
                        disorder_cutoff=0.5,
                        metapredict_version=3) -> str:
    """
    Generate a kappa class variant of the input sequence.
    Parameters
    ----------
    input_sequence : str
        The input sequence to generate a variant from.
    target_kappa : float
        Target kappa value to achieve.
    kappa_tolerance : float
        Acceptable difference between achieved and target kappa.
        default is parameters.MAXIMUM_KAPPA_ERROR
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
        A generated variant sequence that meets the kappa criteria. 
    """
    for _ in range(num_attempts):
        # Generate a variant sequence
        variant_sequence = vgf.generate_kappa_variant(
            input_sequence,
            target_kappa=target_kappa,
            kappa_tolerance=kappa_tolerance)
        
        # if we fail at the variant generation step, try again. 
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
        # if we fail at the disorder check, try again
        if disordered_sequence is None:
            continue
        else:
            # if we have a disordered sequence, we can break out of the loop
            break
            
    # make sure we return a sequence
    if isinstance(disordered_sequence, list):
        # if the disordered_sequence is a list, we will return the first one
        disordered_sequence = disordered_sequence[0]

    return disordered_sequence

def gen_all_props_class_variant(input_sequence: str,
                                 target_FCR: float,
                                 target_NCPR: float,
                                 target_kappa: float,
                                 target_hydropathy: float,
                                 hydropathy_tolerance: float = parameters.HYDRO_ERROR,
                                 kappa_tolerance: float = parameters.MAXIMUM_KAPPA_ERROR,
                                 num_attempts: int = 5,
                                 strict_disorder=False,
                                 disorder_cutoff=0.5,
                                 metapredict_version=3) -> str:
    """
    Generate a variant sequence by constraining FCR, NCPR, kappa, and hydropathy.
    
    Parameters
    ----------
    input_sequence : str
        The input sequence to generate a variant from.
    target_FCR : float
        Target FCR value to achieve.
    target_NCPR : float
        Target NCPR value to achieve.
    target_kappa : float
        Target kappa value to achieve.
    target_hydropathy : float
        Target hydropathy value to achieve.
    hydropathy_tolerance : float
        Acceptable difference between achieved and target hydropathy.
        default is parameters.HYDRO_ERROR   
    kappa_tolerance : float
        Acceptable difference between achieved and target kappa.
        default is parameters.MAXIMUM_KAPPA_ERROR
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
        variant_sequence = vgf.generate_all_props_class_var(
            input_sequence,
            target_FCR=target_FCR,
            target_NCPR=target_NCPR,
            target_kappa=target_kappa,
            target_hydropathy=target_hydropathy,
            hydropathy_tolerance=hydropathy_tolerance,
            kappa_tolerance=kappa_tolerance
        )
        
        # if we fail at the variant generation step, try again. 
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
        # if we fail at the disorder check, try again
        if disordered_sequence is None:
            continue
        else:
            # if we have a disordered sequence, we can break out of the loop
            break
    
    # make sure we return a sequence
    if isinstance(disordered_sequence, list):
        # if the disordered_sequence is a list, we will return the first one
        disordered_sequence = disordered_sequence[0]
        
    return disordered_sequence