##
## create.py
## 
## create.py contains all the user-facing functions associated with GOOSE. 
## This module provides high-level interfaces for generating sequences and variants
## with specific properties. If a new function is added it should be included here 
## and added to the __all__ list.
## 

# If any new functions are added to create.py, you need to add them here.
__all__ = ['sequence', 'seq_by_fractions', 'seq_by_classes', 'seq_by_re', 'seq_by_rg', 'variant', 'seq_fractions']

# Import required modules and functions
from goose import goose_exceptions
from goose.backend_sequence_generation import sequence_generation
from goose.backend_variant_generation.variant_generator import VariantGenerator
from goose.backend import goose_tools
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
    Generate a disordered sequence with specified physicochemical properties.

    This is the main function for creating intrinsically disordered sequences
    with specific characteristics. You can specify multiple properties 
    simultaneously to create sequences with desired combinations of NCPR, 
    FCR, kappa, and hydropathy values.

    Parameters
    ----------
    length : int
        Length of the desired disordered sequence. Must be between the minimum
        and maximum allowed lengths as defined in the parameters module.

    FCR : float, optional
        Fraction of charged residues (between 0 and 1). This includes both
        positively and negatively charged residues.

    NCPR : float, optional
        Net charge per residue (between -1 and 1). Positive values indicate
        net positive charge, negative values indicate net negative charge.

    hydropathy : float, optional
        Mean hydropathy of the sequence (between 0 and 6.1). Higher values
        indicate more hydrophobic sequences.

    kappa : float, optional
        Kappa value describing charge patterning (between 0 and 1). Values
        closer to 1 indicate more even charge distribution.

    attempts : int, optional
        Number of attempts to generate the sequence. Default is 20. Higher
        values increase success probability but take longer.

    disorder_cutoff : float, optional
        Disorder threshold for sequence validation. Sequences must have
        disorder scores above this threshold. Default from parameters module.

    exclude : list, optional
        List of amino acid residues to exclude from sequence generation.
        Cannot exclude charged residues if FCR is specified.

    use_weighted_probabilities : bool, optional
        Whether to use weighted amino acid probabilities. This can increase
        generation success but may reduce sequence diversity. Default is False.

    strict_disorder : bool, optional
        Whether to use strict disorder checking. If True, all residues must
        be above the disorder threshold. Default is False.

    return_all_sequences : bool, optional
        Whether to return all generated sequences. If False, returns only
        the first successful sequence. Default is False.

    custom_probabilities : dict, or string optional
        Custom amino acid probabilities for sequence generation. Keys should
        be single-letter amino acid codes, values should be probabilities.
        String options include the specified organisms in idr_probabilities.py
        These are:
        'mouse', 'fly', 'neurospora', 'yeast', 'arabidopsis', 'e_coli', 'worm', 
        'zebrafish', 'frog', 'dictyostelium', 'human', 'unbiased', 'all'

    metapredict_version : int, optional
        Version of MetaPredict to use for disorder prediction. Default is 3.

    max_consecutive_ordered : int, optional
        Maximum number of consecutive ordered residues allowed. Default from
        parameters module.

    max_total_ordered : float, optional
        Maximum fraction of ordered residues allowed. Default from parameters
        module.

    batch_size : int, optional
        Number of sequences to generate in each batch. Default from parameters
        module.

    hydropathy_tolerance : float, optional
        Tolerance for hydropathy matching. Default from parameters module.

    kappa_tolerance : float, optional
        Tolerance for kappa matching. Default from parameters module.

    Returns
    -------
    str or list
        Generated amino acid sequence as a string if return_all_sequences is
        False, or list of sequences if return_all_sequences is True.

    Raises
    ------
    GooseInputError
        If invalid parameters are provided.
    GooseFail
        If sequence generation fails after all attempts.

    Examples
    --------
    >>> # Generate a 100-residue sequence with specific properties
    >>> seq = sequence(100, FCR=0.3, NCPR=0.1, hydropathy=3.0)
    >>> 
    >>> # Generate sequence excluding certain residues
    >>> seq = sequence(50, exclude=['C', 'M'])
    """
    # Validate sequence length and convert to int if needed
    goose_tools.length_check(length)

    # handle custom_probabilities
    if 'custom_probabilities' in kwargs:
        if kwargs['custom_probabilities'] is not None:
            kwargs['custom_probabilities'] = goose_tools.handle_custom_probabilities(kwargs['custom_probabilities'])

    # Handle legacy parameter name: convert 'hydrophobicity' to 'hydropathy'
    if 'hydrophobicity' in kwargs:
        kwargs['hydropathy'] = kwargs['hydrophobicity']
        # Remove the legacy parameter name
        del kwargs['hydrophobicity']

    # handle legacy parameter name: convert 'cutoff' to 'disorder_cutoff'
    if 'cutoff' in kwargs:
        # Handle legacy parameter name: convert 'cutoff' to 'disorder_cutoff'
        kwargs['disorder_cutoff'] = kwargs['cutoff']
        # Remove the legacy parameter name
        del kwargs['cutoff']

    # Validate that only acceptable keyword arguments were passed
    # Raises exception if unknown parameters are provided
    goose_tools.check_valid_kwargs(kwargs, ['FCR','NCPR', 'hydropathy', 'kappa', 'disorder_cutoff',
                                  'attempts', 'exclude', 'use_weighted_probabilities', 
                                  'strict_disorder', 'return_all_sequences', 'custom_probabilities', 
                                  'metapredict_version', 'max_consecutive_ordered',
                                  'max_total_ordered', 'batch_size',
                                  'hydropathy_tolerance', 'kappa_tolerance'])
    
    # Validate and correct parameter values, setting defaults where needed
    # This must be done before parameter validation to ensure correct types
    kwargs = goose_tools.check_and_correct_props_kwargs(**kwargs)

    # Validate that parameter values are within acceptable ranges
    goose_tools.check_props_parameters(**kwargs)

    # Validate common parameters shared across multiple functions
    # Ensures proper types and ranges for basic sequence generation parameters
    goose_tools.check_basic_parameters(num_attempts=kwargs['attempts'],
                            strict_disorder=kwargs['strict_disorder'],
                            disorder_cutoff=kwargs['disorder_cutoff'],
                            metapredict_version=kwargs['metapredict_version'],
                            return_all_sequences=kwargs['return_all_sequences'],
                            use_weighted_probabilities=kwargs['use_weighted_probabilities'],
                            max_consecutive_ordered=kwargs['max_consecutive_ordered'],
                            max_total_ordered=kwargs['max_total_ordered'],
                            batch_size=kwargs['batch_size'],
                            custom_probabilities=kwargs['custom_probabilities'],
                            exclude=kwargs['exclude'])

    # Generate the sequence using the backend sequence generation engine
    generated_seq = sequence_generation.by_properties(length, fcr=kwargs['FCR'], 
                                                        ncpr=kwargs['NCPR'], 
                                                        hydropathy=kwargs['hydropathy'], 
                                                        kappa=kwargs['kappa'], 
                                                        exclude_residues=kwargs['exclude'], 
                                                        num_attempts=kwargs['attempts'],
                                                        strict_disorder=kwargs['strict_disorder'], 
                                                        disorder_cutoff=kwargs['disorder_cutoff'],
                                                        metapredict_version=kwargs['metapredict_version'], 
                                                        return_all_sequences=kwargs['return_all_sequences'],
                                                        use_weighted_probabilities=kwargs['use_weighted_probabilities'], 
                                                        chosen_probabilities=kwargs['custom_probabilities'],
                                                        max_consecutive_ordered= kwargs['max_consecutive_ordered'],
                                                        max_total_ordered=kwargs['max_total_ordered'],
                                                        batch_size=kwargs['batch_size'],
                                                        hydropathy_tolerance=kwargs['hydropathy_tolerance'],
                                                        kappa_tolerance=kwargs['kappa_tolerance'])
        
    # Check if sequence generation failed and raise appropriate error
    if generated_seq is None:
        raise goose_exceptions.GooseFail('Unable to generate sequence. Please try again with different parameters or a lower cutoff value.')

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


def seq_by_fractions(length, **kwargs):
    """
    Generate a disordered sequence with specified amino acid fractions.

    This function creates intrinsically disordered sequences where you can
    specify the exact fraction of each amino acid type. This provides fine-
    grained control over sequence composition.

    Parameters
    ----------
    length : int
        Length of the desired disordered sequence.

    A, C, D, E, F, G, H, I, K, L, M, N, P, Q, R, S, T, V, W, Y : float, optional
        Fraction of the sequence that should be made up of each specific
        amino acid (e.g., A=0.2, Y=0.05). Values should be between 0 and 1.
        The sum of all specified fractions should not exceed 1.

    max_aa_fractions : dict, optional
        Dictionary to override the maximum allowed fraction of any amino acid.
        Keys should be single-letter amino acid codes, values should be floats
        between 0 and 1. If not specified, default GOOSE thresholds are used.

    disorder_cutoff : float, optional
        Disorder threshold for sequence validation. Default is 0.6.

    attempts : int, optional
        Number of attempts to generate the sequence. Default is 100.

    strict_disorder : bool, optional
        Whether to use strict disorder checking. If True, all residues must
        be above the disorder threshold. Default is False.

    remaining_probabilities : dict, or string optional
        Custom amino acid probabilities for sequence generation. Keys should
        be single-letter amino acid codes, values should be probabilities.
        String options include the specified organisms in idr_probabilities.py
        These are:
        'mouse', 'fly', 'neurospora', 'yeast', 'arabidopsis', 'e_coli', 'worm', 
        'zebrafish', 'frog', 'dictyostelium', 'human', 'unbiased', 'all'

    return_all_sequences : bool, optional
        Whether to return all generated sequences. Default is False.

    metapredict_version : int, optional
        Version of MetaPredict to use for disorder prediction. Default is 3.

    max_consecutive_ordered : int, optional
        Maximum number of consecutive ordered residues allowed.

    max_total_ordered : float, optional
        Maximum fraction of ordered residues allowed.

    batch_size : int, optional
        Number of sequences to generate in each batch.

    Returns
    -------
    str or list
        Generated amino acid sequence as a string, or list of sequences if
        return_all_sequences is True.

    Raises
    ------
    GooseInputError
        If invalid parameters are provided.
    GooseFail
        If sequence generation fails after all attempts.

    Examples
    --------
    >>> # Generate sequence with 30% alanine and 10% glycine
    >>> seq = seq_by_fractions(100, A=0.3, G=0.1)
    >>> 
    >>> # Generate sequence with custom max fractions
    >>> seq = seq_by_fractions(50, A=0.4, max_aa_fractions={'A': 0.5})
    """

    # Validate sequence length and convert to int if needed
    goose_tools.length_check(length)

    # handle remaining_probabilities
    if 'remaining_probabilities' in kwargs:
        # If provided, handle custom probabilities
        if kwargs['remaining_probabilities'] is not None:
            kwargs['remaining_probabilities'] = goose_tools.handle_custom_probabilities(kwargs['remaining_probabilities'])

    # handle legacy parameter name: convert 'cutoff' to 'disorder_cutoff'
    if 'cutoff' in kwargs:
        # Handle legacy parameter name: convert 'cutoff' to 'disorder_cutoff'
        kwargs['disorder_cutoff'] = kwargs['cutoff']
        # Remove the legacy parameter name
        del kwargs['cutoff']

    # Validate that only acceptable keyword arguments were passed
    # Raises exception if unknown parameters are provided
    goose_tools.check_valid_kwargs(kwargs, ['remaining_probabilities', 'attempts', 
                                  'strict_disorder', 'disorder_cutoff', 'max_consecutive_ordered',
                                  'max_total_ordered','max_aa_fractions', 
                                 'A','C','D','E','F','G','H','I','K','L','M','N',
                                 'P','Q','R','S','T','V','W','Y', 
                                 'return_all_sequences', 'metapredict_version',
                                  'batch_size'])


    # Validate and correct parameter values, setting defaults where needed
    # This must be done before parameter validation to ensure correct types
    kwargs = goose_tools.check_and_correct_fracs_kwargs(**kwargs)

    # Extract amino acid fractions from kwargs
    # Build dictionary of explicitly specified amino acid fractions
    fractions = {}
    for f in kwargs:
        if f in ['A','C','D','E','F','G','H','I','K','L','M','N','P','Q','R','S','T','V','W','Y']:
            fractions[f] = kwargs[f]

    # Validate that fraction parameter values are within acceptable ranges
    goose_tools.check_fracs_parameters(**kwargs)
    
    # Validate common parameters shared across multiple functions
    # Ensures proper types and ranges for basic sequence generation parameters
    goose_tools.check_basic_parameters(num_attempts=kwargs['attempts'],
                            strict_disorder=kwargs['strict_disorder'],
                            disorder_cutoff=kwargs['disorder_cutoff'],
                            metapredict_version=kwargs['metapredict_version'],
                            return_all_sequences=kwargs['return_all_sequences'],
                            max_consecutive_ordered=kwargs['max_consecutive_ordered'],
                            max_total_ordered=kwargs['max_total_ordered'],
                            batch_size=kwargs['batch_size'],
                            custom_probabilities=kwargs['remaining_probabilities'])

    
    # Generate the sequence using the backend fraction-based generation engine
    generated_seq = sequence_generation.by_fractions(length,fractions=fractions,
                                                     remaining_probabilities=kwargs['remaining_probabilities'],
                                                     num_attempts=kwargs['attempts'],
                                                     strict_disorder=kwargs['strict_disorder'],
                                                     disorder_cutoff=kwargs['disorder_cutoff'],
                                                     max_consecutive_ordered=kwargs['max_consecutive_ordered'],
                                                     max_total_ordered=kwargs['max_total_ordered'],
                                                     metapredict_version=kwargs['metapredict_version'],
                                                     return_all_sequences=kwargs['return_all_sequences'],
                                                     batch_size = kwargs['batch_size'])

    # Return the generated sequence (no need to check for None as backend handles this)
    return generated_seq

def seq_fractions(length, **kwargs):
    """
    Generate a disordered sequence with specified amino acid fractions.
    
    This function is a backwards compatibility wrapper around seq_by_fractions.
    Please use seq_by_fractions for new code.
    
    Parameters
    ----------
    length : int
        Length of the desired disordered sequence.
    **kwargs : dict
        All keyword arguments are passed directly to seq_by_fractions.
        See seq_by_fractions documentation for full parameter details.
    
    Returns
    -------
    str or list
        Generated amino acid sequence(s) - see seq_by_fractions for details.
    
    See Also
    --------
    seq_by_fractions : The main function for generating sequences by fractions.
    
    Examples
    --------
    >>> # Generate sequence with 30% alanine and 10% glycine
    >>> seq = seq_fractions(100, A=0.3, G=0.1)
    """
    return seq_by_fractions(length, **kwargs)


#-\|/-\|/-\|/-\|/-\|/-\|/-\|/-\|/-\|/-\|/-\|/-\|/-\|/-\|/-\|/-\|/-\|/-\|/-\|/-\|/-\|/-
#-\|/-\|/-\|/-\|/-\|/-\|/-\|/-\|/-\|/             \|/-\|/-\|/-\|/-\|/-\|/-\|/-\|/-\|/-
#-\|/-\|/-\|/-\|/-\|/-\|/-\|/-\|/-\|/  Create     \|/-\|/-\|/-\|/-\|/-\|/-\|/-\|/-\|/-
#-\|/-\|/-\|/-\|/-\|/-\|/-\|/-\|/-\|/  sequence   \|/-\|/-\|/-\|/-\|/-\|/-\|/-\|/-\|/-
#-\|/-\|/-\|/-\|/-\|/-\|/-\|/-\|/-\|/  By         \|/-\|/-\|/-\|/-\|/-\|/-\|/-\|/-\|/-
#-\|/-\|/-\|/-\|/-\|/-\|/-\|/-\|/-\|/  classes    \|/-\|/-\|/-\|/-\|/-\|/-\|/-\|/-\|/-
#-\|/-\|/-\|/-\|/-\|/-\|/-\|/-\|/-\|/             \|/-\|/-\|/-\|/-\|/-\|/-\|/-\|/-\|/-
#-\|/-\|/-\|/-\|/-\|/-\|/-\|/-\|/-\|/-\|/-\|/-\|/-\|/-\|/-\|/-\|/-\|/-\|/-\|/-\|/-\|/-

def seq_by_classes(length: int,
                    aromatic: float = 0.0,
                    aliphatic: float = 0.0,
                    polar: float = 0.0,
                    positive: float = 0.0,
                    negative: float = 0.0,
                    glycine: float = 0.0,
                    proline: float = 0.0,
                    cysteine: float = 0.0,
                    histidine: float = 0.0,
                    num_attempts=10, strict_disorder=False,
                    disorder_cutoff=parameters.DISORDER_THRESHOLD,
                    cutoff = None, # legacy parameter name
                    metapredict_version=parameters.METAPREDICT_DEFAULT_VERSION,
                    max_consecutive_ordered=parameters.ALLOWED_CONSECUTIVE_ORDERED,
                    max_total_ordered=parameters.ALLOWED_TOTAL_ORDERED_FRACTION,
                    remaining_probabilities=None,
                    max_class_fractions=None):
    """
    Generate a disordered sequence with specified amino acid class fractions.

    This function creates intrinsically disordered sequences where you can
    specify the fraction of different amino acid classes (aromatic, aliphatic,
    polar, charged, etc.) rather than individual amino acids. This provides
    a higher-level approach to sequence composition control.

    Parameters
    ----------
    length : int
        Length of the desired disordered sequence.
    aromatic : float, optional
        Fraction of aromatic amino acids (F, W, Y) in the sequence 
        (between 0 and 1). Default is 0.0.
    aliphatic : float, optional
        Fraction of aliphatic amino acids (A, I, L, V) in the sequence 
        (between 0 and 1). Default is 0.0.
    polar : float, optional
        Fraction of polar amino acids (N, Q, S, T) in the sequence 
        (between 0 and 1). Default is 0.0.
    positive : float, optional
        Fraction of positively charged amino acids (K, R) in the sequence 
        (between 0 and 1). Default is 0.0.
    negative : float, optional
        Fraction of negatively charged amino acids (D, E) in the sequence 
        (between 0 and 1). Default is 0.0.
    glycine : float, optional
        Fraction of glycine (G) in the sequence (between 0 and 1). Default is 0.0.
    proline : float, optional
        Fraction of proline (P) in the sequence (between 0 and 1). Default is 0.0.
    cysteine : float, optional
        Fraction of cysteine (C) in the sequence (between 0 and 1). Default is 0.0.
    histidine : float, optional
        Fraction of histidine (H) in the sequence (between 0 and 1). Default is 0.0.
    num_attempts : int, optional
        Number of attempts to generate the sequence. Default is 10.
    strict_disorder : bool, optional
        Whether to use strict disorder checking. If True, all residues must
        be above the disorder threshold. Default is False.
    disorder_cutoff : float, optional
        Disorder threshold for sequence validation. Default from parameters module.
    metapredict_version : int, optional
        Version of MetaPredict to use for disorder prediction. Default is 3.
    max_consecutive_ordered : int, optional
        Maximum number of consecutive ordered residues allowed. Default from
        parameters module.
    max_total_ordered : float, optional
        Maximum fraction of ordered residues allowed. Default from parameters
        module.
    remaining_probabilities : dict, or string optional
        Custom amino acid probabilities for sequence generation. Keys should
        be single-letter amino acid codes, values should be probabilities.
        String options include the specified organisms in idr_probabilities.py
        These are:
        'mouse', 'fly', 'neurospora', 'yeast', 'arabidopsis', 'e_coli', 'worm', 
        'zebrafish', 'frog', 'dictyostelium', 'human', 'unbiased', 'all'
    cutoff : float, optional
        Legacy parameter name for disorder cutoff. If provided, it will override
        the default disorder_cutoff value.
    max_class_fractions : dict, optional
        Dictionary to override the maximum allowed fraction of any amino acid class.
        Keys should be class names ('aromatic', 'aliphatic', 'polar', 'positive',
        'negative', 'glycine', 'proline', 'cysteine', 'histidine'), values should
        be floats between 0 and 1. If not specified, default GOOSE thresholds are used.

    Returns
    -------
    str
        Generated amino acid sequence as a string.

    Raises
    ------
    GooseInputError
        If invalid parameters are provided.
    GooseFail
        If sequence generation fails after all attempts.

    Examples
    --------
    >>> # Generate sequence with 20% aromatic and 10% positive residues
    >>> seq = seq_by_classes(100, aromatic=0.2, positive=0.1)
    >>> 
    >>> # Generate sequence with multiple class constraints
    >>> seq = seq_by_classes(75, aromatic=0.15, polar=0.25, glycine=0.1)
    """
    # Validate sequence length and convert to int if needed
    goose_tools.length_check(length)

    if max_class_fractions is not None:
        # overwrite values in the default max_class_fractions
        original_max_class_fractions = parameters.MAX_CLASS_FRACTIONS.copy()
        for key in max_class_fractions:
            if key in original_max_class_fractions:
                original_max_class_fractions[key] = max_class_fractions[key]
        max_class_fractions = original_max_class_fractions
    else:
        # Use default max_class_fractions if not provided
        max_class_fractions = parameters.MAX_CLASS_FRACTIONS.copy()

    # Validate that the specified class fractions are within acceptable bounds
    goose_tools.check_class_values(max_class_fractions, aromatic=aromatic, aliphatic=aliphatic, polar=polar,
                        positive=positive, negative=negative, glycine=glycine,
                        proline=proline, cysteine=cysteine, histidine=histidine)
    

    # handle legacy parameter name: convert 'cutoff' to 'disorder_cutoff'
    if cutoff is not None:
        # Handle legacy parameter name: convert 'cutoff' to 'disorder_cutoff'
        disorder_cutoff = cutoff

    # handle remaining_probabilities
    if remaining_probabilities is not None:
        remaining_probabilities = goose_tools.handle_custom_probabilities(remaining_probabilities)

    # Validate common parameters shared across multiple functions
    # Ensures proper types and ranges for basic sequence generation parameters
    goose_tools.check_basic_parameters(num_attempts=num_attempts,
                        strict_disorder=strict_disorder,
                        disorder_cutoff=disorder_cutoff,
                        metapredict_version=metapredict_version,
                        max_consecutive_ordered=max_consecutive_ordered,
                        max_total_ordered=max_total_ordered,
                        custom_probabilities=remaining_probabilities)

    # Generate the sequence using the backend class-based generation engine
    generated_seq = sequence_generation.by_class(length, 
        aromatic_fraction=aromatic,
        aliphatic_fraction=aliphatic,
        polar_fraction=polar,
        positive_fraction=positive,
        negative_fraction=negative,
        glycine_fraction=glycine,
        proline_fraction=proline,
        cysteine_fraction=cysteine,
        histidine_fraction=histidine,
        num_attempts=num_attempts, strict_disorder=strict_disorder,
        disorder_cutoff=disorder_cutoff, metapredict_version=metapredict_version,
        max_consecutive_ordered=max_consecutive_ordered,
        max_total_ordered=max_total_ordered,
        remaining_probabilities=remaining_probabilities)

    # Check if sequence generation failed and raise appropriate error
    if generated_seq is None:
        raise goose_exceptions.GooseFail('Unable to generate sequence. Please try again with different parameters or a lower cutoff value.')

    return generated_seq

#-\|/-\|/-\|/-\|/-\|/-\|/-\|/-\|/-\|/-\|/-\|/-\|/-\|/-\|/-\|/-\|/-\|/-\|/-\|/-\|/-\|/-
#-\|/-\|/-\|/-\|/-\|/-\|/-\|/-\|/-\|/             \|/-\|/-\|/-\|/-\|/-\|/-\|/-\|/-\|/-
#-\|/-\|/-\|/-\|/-\|/-\|/-\|/-\|/-\|/  Create     \|/-\|/-\|/-\|/-\|/-\|/-\|/-\|/-\|/-
#-\|/-\|/-\|/-\|/-\|/-\|/-\|/-\|/-\|/  sequence   \|/-\|/-\|/-\|/-\|/-\|/-\|/-\|/-\|/-
#-\|/-\|/-\|/-\|/-\|/-\|/-\|/-\|/-\|/  By         \|/-\|/-\|/-\|/-\|/-\|/-\|/-\|/-\|/-
#-\|/-\|/-\|/-\|/-\|/-\|/-\|/-\|/-\|/  dimensions \|/-\|/-\|/-\|/-\|/-\|/-\|/-\|/-\|/-
#-\|/-\|/-\|/-\|/-\|/-\|/-\|/-\|/-\|/             \|/-\|/-\|/-\|/-\|/-\|/-\|/-\|/-\|/-
#-\|/-\|/-\|/-\|/-\|/-\|/-\|/-\|/-\|/-\|/-\|/-\|/-\|/-\|/-\|/-\|/-\|/-\|/-\|/-\|/-\|/-

def seq_by_re(length, objective_re, allowed_error=parameters.MAXIMUM_RG_RE_ERROR, 
           attempts=20, disorder_cutoff=parameters.DISORDER_THRESHOLD, 
           strict_disorder=False, reduce_pos_charged=False, exclude_aas=None,
            metapredict_version=parameters.METAPREDICT_DEFAULT_VERSION,
            max_consecutive_ordered=parameters.ALLOWED_CONSECUTIVE_ORDERED,
            max_total_ordered=parameters.ALLOWED_TOTAL_ORDERED_FRACTION,
            cutoff=None):
    """
    Generate a disordered sequence with a specified end-to-end distance (Re).

    This function creates intrinsically disordered sequences with a target
    end-to-end distance (Re) in Angstroms. The end-to-end distance is the
    average distance between the N and C termini of the sequence.

    Parameters
    ----------
    length : int
        Length of the sequence to generate.
    objective_re : float
        Target end-to-end distance in Angstroms.
    allowed_error : float, optional
        Allowed error between the target and actual Re value. Default from
        parameters module.
    attempts : int, optional
        Number of attempts to generate the sequence. Default is 20.
    disorder_cutoff : float, optional
        Disorder threshold for sequence validation. Default from parameters module.
    strict_disorder : bool, optional
        Whether to use strict disorder checking. If True, all residues must
        be above the disorder threshold. Default is False.
    reduce_pos_charged : bool, optional
        Whether to reduce positively charged amino acids in the sequence.
        Default is False. In vivo data suggests positively charged residues
        may not drive sequence expansion as much as predicted by the model.
    exclude_aas : list, optional
        List of amino acids to exclude from the sequence. Default is None.
    metapredict_version : int, optional
        Version of MetaPredict to use for disorder prediction. Default is 3.
    max_consecutive_ordered : int, optional
        Maximum number of consecutive ordered residues allowed. Default from
        parameters module.
    max_total_ordered : float, optional
        Maximum fraction of ordered residues allowed. Default from parameters
        module.
    cutoff : float, optional
        Legacy parameter name for disorder cutoff. If provided, it will override
        the default disorder_cutoff value.

    Returns
    -------
    str
        Generated amino acid sequence as a string.

    Raises
    ------
    GooseInputError
        If the objective_re is outside the possible range for the given length,
        or if other invalid parameters are provided.
    GooseFail
        If sequence generation fails after all attempts.

    Examples
    --------
    >>> # Generate a 100-residue sequence with Re = 50 Å
    >>> seq = seq_by_re(100, 50.0)
    >>> 
    >>> # Generate with custom error tolerance
    >>> seq = seq_by_re(75, 40.0, allowed_error=2.0)
    """
    # Validate sequence length and convert to int if needed
    goose_tools.length_check(length)

    # Validate that the objective Re is within the possible range for this length
    if objective_re < parameters.get_min_re(length) or objective_re > parameters.get_max_re(length):
        min_possible_value = parameters.get_min_re(length)
        max_possible_value = parameters.get_max_re(length)
        raise goose_exceptions.GooseInputError(f'Cannot generate sequence, for length {length}, min Re = {min_possible_value}, max Re = {max_possible_value}.')

    # handle legacy parameter name: convert 'cutoff' to 'disorder_cutoff'
    if cutoff is not None:
        # Handle legacy parameter name: convert 'cutoff' to 'disorder_cutoff'
        disorder_cutoff = cutoff

    # Validate that allowed error is positive
    if allowed_error < 0:
        raise goose_exceptions.GooseInputError('Allowed error must be a positive number.')
    
    # Validate that reduce_pos_charged is a boolean
    if not isinstance(reduce_pos_charged, bool):
        raise goose_exceptions.GooseInputError('reduce_pos_charged must be a boolean value.')

    # Validate common parameters shared across multiple functions
    # Ensures proper types and ranges for basic sequence generation parameters
    goose_tools.check_basic_parameters(num_attempts=attempts,
                            strict_disorder=strict_disorder,
                            disorder_cutoff=disorder_cutoff,
                            metapredict_version=metapredict_version,
                            max_consecutive_ordered=max_consecutive_ordered,
                            max_total_ordered=max_total_ordered,
                            exclude=exclude_aas)
    
    # Generate the sequence using the backend dimensional constraint engine
    sequence = sequence_generation.by_dimensions(
        length, objective_re, rg_or_re='re',  
        allowed_error=allowed_error, 
        reduce_pos_charged=reduce_pos_charged,
        exclude_aas=exclude_aas,
        num_attempts=attempts, 
        strict_disorder=strict_disorder,
        disorder_cutoff=disorder_cutoff,
        metapredict_version=metapredict_version,
        max_consecutive_ordered=max_consecutive_ordered,
        max_total_ordered=max_total_ordered
       )
    return sequence


def seq_by_rg(length, objective_rg, allowed_error=parameters.MAXIMUM_RG_RE_ERROR, 
           attempts=20, disorder_cutoff=parameters.DISORDER_THRESHOLD, 
           strict_disorder=False, reduce_pos_charged=False, exclude_aas=None,
            metapredict_version=parameters.METAPREDICT_DEFAULT_VERSION,
            max_consecutive_ordered=parameters.ALLOWED_CONSECUTIVE_ORDERED,
            max_total_ordered=parameters.ALLOWED_TOTAL_ORDERED_FRACTION,
            cutoff=None):
    """
    Generate a disordered sequence with a specified radius of gyration (Rg).

    This function creates intrinsically disordered sequences with a target
    radius of gyration (Rg) in Angstroms. The radius of gyration is a measure
    of the compactness of the sequence's ensemble of conformations.

    Parameters
    ----------
    length : int
        Length of the sequence to generate.
    objective_rg : float
        Target radius of gyration in Angstroms.
    allowed_error : float, optional
        Allowed error between the target and actual Rg value. Default from
        parameters module.
    attempts : int, optional
        Number of attempts to generate the sequence. Default is 20.
    disorder_cutoff : float, optional
        Disorder threshold for sequence validation. Default from parameters module.
    strict_disorder : bool, optional
        Whether to use strict disorder checking. If True, all residues must
        be above the disorder threshold. Default is False.
    reduce_pos_charged : bool, optional
        Whether to reduce positively charged amino acids in the sequence.
        Default is False. In vivo data suggests positively charged residues
        may not drive sequence expansion as much as predicted by the model.
    exclude_aas : list, optional
        List of amino acids to exclude from the sequence. Default is None.
    metapredict_version : int, optional
        Version of MetaPredict to use for disorder prediction. Default is 3.
    max_consecutive_ordered : int, optional
        Maximum number of consecutive ordered residues allowed. Default from
        parameters module.
    max_total_ordered : float, optional
        Maximum fraction of ordered residues allowed. Default from parameters
        module.
    cutoff : float, optional
        Legacy parameter name for disorder cutoff. If provided, it will override
        the default disorder_cutoff value.

    Returns
    -------
    str
        Generated amino acid sequence as a string.

    Raises
    ------
    GooseInputError
        If the objective_rg is outside the possible range for the given length,
        or if other invalid parameters are provided.
    GooseFail
        If sequence generation fails after all attempts.

    Examples
    --------
    >>> # Generate a 100-residue sequence with Rg = 25 Å
    >>> seq = seq_by_rg(100, 25.0)
    >>> 
    >>> # Generate with reduced positive charges
    >>> seq = seq_by_rg(75, 20.0, reduce_pos_charged=True)
    """
    # Validate sequence length and convert to int if needed
    goose_tools.length_check(length)

    # Validate that the objective Rg is within the possible range for this length
    if objective_rg < parameters.get_min_rg(length) or objective_rg > parameters.get_max_rg(length):
        min_possible_value = parameters.get_min_rg(length)
        max_possible_value = parameters.get_max_rg(length)
        raise goose_exceptions.GooseInputError(f'Cannot generate sequence, for length {length}, min Rg = {min_possible_value}, max Rg = {max_possible_value}.')

    # handle legacy parameter name: convert 'cutoff' to 'disorder_cutoff'
    if cutoff is not None:
        # Handle legacy parameter name: convert 'cutoff' to 'disorder_cutoff'
        disorder_cutoff = cutoff

    # Validate that allowed error is positive
    if allowed_error < 0:
        raise goose_exceptions.GooseInputError('Allowed error must be a positive number.')
    
    # Validate that reduce_pos_charged is a boolean
    if not isinstance(reduce_pos_charged, bool):
        raise goose_exceptions.GooseInputError('reduce_pos_charged must be a boolean value.')

    # Validate common parameters shared across multiple functions
    # Ensures proper types and ranges for basic sequence generation parameters
    goose_tools.check_basic_parameters(num_attempts=attempts,
                            strict_disorder=strict_disorder,
                            disorder_cutoff=disorder_cutoff,
                            metapredict_version=metapredict_version,
                            max_consecutive_ordered=max_consecutive_ordered,
                            max_total_ordered=max_total_ordered,
                            exclude=exclude_aas)
    
    # Generate the sequence using the backend dimensional constraint engine
    sequence = sequence_generation.by_dimensions(
        length, objective_rg, rg_or_re='rg',  
        allowed_error=allowed_error, 
        reduce_pos_charged=reduce_pos_charged,
        exclude_aas=exclude_aas,
        num_attempts=attempts, 
        strict_disorder=strict_disorder,
        disorder_cutoff=disorder_cutoff,
        metapredict_version=metapredict_version,
        max_consecutive_ordered=max_consecutive_ordered,
        max_total_ordered=max_total_ordered
       )
    return sequence


#===============================================================================
#===============================================================================
#                             VARIANT GENERATORS
#===============================================================================
#===============================================================================

def variant(sequence,
            variant_type,
            **kwargs):
    """
    Generate variants of an input sequence using various transformation methods.

    This function provides a unified interface for creating sequence variants
    using different algorithms. It supports shuffling, repositioning, and
    property-based modifications of amino acid sequences while maintaining
    disorder characteristics.

    Parameters
    ----------
    sequence : str
        The amino acid sequence to generate variants from. Must be a non-empty
        string containing valid amino acid codes.
    
    variant_type : str
        The type of variant to generate. Available options:
        
        Shuffling methods:
        - 'shuffle_specific_regions': Shuffle only specified regions
        - 'shuffle_except_specific_regions': Shuffle all except specified regions
        - 'shuffle_specific_residues': Shuffle only specific residue types
        - 'shuffle_except_specific_residues': Shuffle all except specific residue types
        - 'weighted_shuffle_specific_residues': Weighted shuffle of specific residues
        - 'targeted_reposition_specific_residues': Reposition specific residues
        
        Asymmetry and property methods:
        - 'change_residue_asymmetry': Change residue asymmetry patterns
        - 'constant_properties': Generate variant with constant properties
        - 'constant_residues_and_properties': Keep specified residues and properties constant
        - 'constant_properties_and_class': Generate variant with constant properties and class
        - 'constant_properties_and_class_by_order': Generate variant with constant properties and class by order
        
        Property modification methods:
        - 'change_hydropathy_constant_class': Change hydropathy while keeping class constant
        - 'change_fcr_minimize_class_changes': Change FCR while minimizing class changes
        - 'change_ncpr_constant_class': Change NCPR while keeping class constant
        - 'change_kappa': Change kappa value
        - 'change_properties_minimize_differences': Change properties while minimizing differences
        - 'change_any_properties': Change any combination of properties
        - 'change_dimensions': Change sequence dimensions (Rg/Re)
    
    **kwargs : dict
        Additional parameters specific to the variant type. Common parameters include:
        
        General parameters:
        - num_attempts (int): Number of attempts to generate variant (default: 100)
        - strict_disorder (bool): Whether to use strict disorder checking (default: False)
        - disorder_cutoff (float): Disorder cutoff threshold (default: from parameters)
        - metapredict_version (int): MetaPredict version to use (default: 3)
        - hydropathy_tolerance (float): Hydropathy tolerance (default: from parameters)
        - kappa_tolerance (float): Kappa tolerance (default: from parameters)
        
        Variant-specific parameters:
        - shuffle_regions (list): Regions to shuffle (tuple pairs of start/end positions)
        - excluded_regions (list): Regions to exclude from shuffling
        - target_residues (list): Specific residues to target
        - excluded_residues (list): Specific residues to exclude
        - shuffle_weight (float): Weight for shuffling operations
        - num_changes (int): Number of changes to make
        - increase_or_decrease (str): Direction of change ('increase' or 'decrease')
        - exclude_residues (list): Residues to exclude from modifications
        - constant_residues (list): Residues to keep constant
        - target_hydropathy (float): Target hydropathy value
        - target_FCR (float): Target FCR value
        - target_NCPR (float): Target NCPR value
        - target_kappa (float): Target kappa value
        - rg_or_re (str): Whether to optimize 'rg' or 're'
        - num_dim_attempts (int): Number of dimensional optimization attempts
        - allowed_error (float): Allowed error for dimensional constraints
        - reduce_pos_charged (bool): Whether to reduce positive charges
        - exclude_aas (list): Amino acids to exclude from generation

    Returns
    -------
    str
        Generated variant sequence as a string.

    Raises
    ------
    GooseInputError
        If invalid parameters are provided, including:
        - Empty or invalid sequence
        - Invalid variant_type
        - Missing required parameters for the specified variant type
        - Invalid parameter values
    GooseFail
        If variant generation fails after all attempts.

    Examples
    --------
    >>> # Shuffle specific regions of a sequence
    >>> original = "MSEDKQRTYHLNVAIGPKWF"
    >>> variant = variant(original, 'shuffle_specific_regions', 
    ...                  shuffle_regions=[(0, 5), (10, 15)])
    >>> 
    >>> # Change hydropathy while keeping amino acid classes constant
    >>> variant = variant(original, 'change_hydropathy_constant_class',
    ...                  target_hydropathy=3.5)
    >>> 
    >>> # Generate variant with constant properties but different sequence
    >>> variant = variant(original, 'constant_properties', num_attempts=50)
    
    Notes
    -----
    The function uses the VariantGenerator class from the backend to perform
    the actual sequence modifications. Each variant type has specific parameter
    requirements - consult the documentation for detailed parameter descriptions.
    
    Region specifications use 0-based indexing where (start, end) includes
    positions from start to end-1, following Python slice conventions.
    """
    # Validate that the input sequence is a non-empty string
    if not isinstance(sequence, str):
        raise goose_exceptions.GooseInputError('Sequence must be a string')
    
    if len(sequence) == 0:
        raise goose_exceptions.GooseInputError('Sequence cannot be empty')
    
    # Define all valid variant types for validation
    valid_variant_types = {
        'shuffle_specific_regions',
        'shuffle_except_specific_regions', 
        'shuffle_specific_residues',
        'shuffle_except_specific_residues',
        'weighted_shuffle_specific_residues',
        'targeted_reposition_specific_residues',
        'change_residue_asymmetry',
        'constant_properties',
        'constant_residues_and_properties',
        'constant_properties_and_class',
        'constant_properties_and_class_by_order',
        'change_hydropathy_constant_class',
        'change_fcr_minimize_class_changes',
        'change_ncpr_constant_class',
        'change_kappa',
        'change_properties_minimize_differences',
        'change_any_properties',
        'change_dimensions'
    }
    
    # Validate that the variant type is supported
    if variant_type not in valid_variant_types:
        raise goose_exceptions.GooseInputError(f'Invalid variant_type: {variant_type}. Must be one of: {", ".join(sorted(valid_variant_types))}')
    
    # Extract common parameters with default values
    # These parameters are used across multiple variant generation methods
    common_params = {
        'num_attempts': kwargs.get('num_attempts', 100),
        'strict_disorder': kwargs.get('strict_disorder', False),
        'disorder_cutoff': kwargs.get('disorder_cutoff', parameters.DISORDER_THRESHOLD),
        'metapredict_version': kwargs.get('metapredict_version', parameters.METAPREDICT_DEFAULT_VERSION),
        'hydropathy_tolerance': kwargs.get('hydropathy_tolerance', parameters.MAXIMUM_HYDRO_ERROR),
        'kappa_tolerance': kwargs.get('kappa_tolerance', parameters.MAXIMUM_KAPPA_ERROR)
    }
    
    # Create VariantGenerator instance with common parameters
    generator = VariantGenerator(**common_params)
    
    # Define method dispatch mapping: variant_type -> (method_name, required_params, optional_params)
    # This maps user-facing variant types to backend method names and their parameter requirements
    method_dispatch = {
        'shuffle_specific_regions': ('shuffle_specific_regions', ['shuffle_regions']),
        'shuffle_except_specific_regions': ('shuffle_except_specific_regions', ['excluded_regions']),
        'shuffle_specific_residues': ('shuffle_specific_residues', ['target_residues']),
        'shuffle_except_specific_residues': ('shuffle_except_specific_residues', ['excluded_residues']),
        'weighted_shuffle_specific_residues': ('weighted_shuffle_specific_residues', ['target_residues', 'shuffle_weight']),
        'targeted_reposition_specific_residues': ('targeted_reposition_specific_residues', ['target_residues']),
        'change_residue_asymmetry': ('change_residue_asymmetry', ['target_residues'], ['num_changes', 'increase_or_decrease']),
        'constant_properties': ('constant_properties', [], ['exclude_residues']),
        'constant_residues_and_properties': ('constant_residues_and_properties', ['constant_residues']),
        'constant_properties_and_class': ('constant_properties_and_class', []),
        'constant_properties_and_class_by_order': ('constant_properties_and_class_by_order', []),
        'change_hydropathy_constant_class': ('change_hydropathy_constant_class', ['target_hydropathy']),
        'change_fcr_minimize_class_changes': ('change_fcr_minimize_class_changes', ['target_FCR']),
        'change_ncpr_constant_class': ('change_ncpr_constant_class', ['target_NCPR']),
        'change_kappa': ('change_kappa', ['target_kappa']),
        'change_properties_minimize_differences': ('change_properties_minimze_differences', [], ['target_hydropathy', 'target_kappa', 'target_FCR', 'target_NCPR']),
        'change_any_properties': ('change_any_properties', ['target_FCR', 'target_NCPR', 'target_kappa', 'target_hydropathy']),
        'change_dimensions': ('change_dimensions', ['increase_or_decrease', 'rg_or_re'], ['num_dim_attempts', 'allowed_error', 'reduce_pos_charged', 'exclude_aas'])
    }
    
    # Get method information for the requested variant type
    method_info = method_dispatch[variant_type]
    method_name = method_info[0]
    required_params = method_info[1]
    optional_params = method_info[2] if len(method_info) > 2 else []
    
    # Validate that all required parameters are provided
    for param in required_params:
        if param not in kwargs:
            raise goose_exceptions.GooseInputError(f'Missing required parameter for {variant_type}: {param}')
    
    # Prepare arguments for the backend method call
    method_args = {'input_sequence': sequence}
    
    # Add required parameters to method arguments
    for param in required_params:
        method_args[param] = kwargs[param]
    
    # Add optional parameters if they were provided
    for param in optional_params:
        if param in kwargs:
            method_args[param] = kwargs[param]
    
    # Call the appropriate backend method and handle potential errors
    try:
        method = getattr(generator, method_name)
        result = method(**method_args)
        
        # Check if variant generation failed
        if result is None:
            raise goose_exceptions.GooseFail(f'Failed to generate variant of type {variant_type}. Try adjusting parameters or increasing num_attempts.')
        
        return result
        
    except AttributeError:
        # This should not happen if method_dispatch is correct, but catch it just in case
        raise goose_exceptions.GooseInputError(f'Method {method_name} not found in VariantGenerator')
    except Exception as e:
        # Catch any other errors from the backend and re-raise as GooseFail
        raise goose_exceptions.GooseFail(f'Error generating variant: {str(e)}')