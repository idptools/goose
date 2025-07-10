##
## create.py
## 
## create.py contains all the user-facing function associated with GOOSE. 
## If a new function is added it should be included here and added to the __all__ list
## 

# if any new functions are added to create.py, you need to add them here.
__all__ =  ['sequence', 'seq_by_fractions', 'seq_by_classes', 'seq_by_re', 'seq_by_rg', 'variant']

# imports
from goose.backend_sequence_generation import sequence_generation
from goose.backend_variant_generation.variant_generator import VariantGenerator
from goose.backend import goose_tools
from goose import goose_exceptions
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

    NB:

    1. You can specify NCPR, FCR, kappa, and hydropathy simultaneously.
    2. You don't need to specify anything other than length. 

    Parameters
    ------------
    length : int
        Defines the length of desired disordered sequence. (Required)

    FCR : float
        Defines the requested fraction of charged residues (between 0 and 1). (Optional)

    NCPR : float 
        Defines the net charge per residue of generated sequences (between -1 and 1). (Optional)

    hydropathy : float 
        Defines the mean hydropathy of generated sequence (between 0 and 6.1). (Optional)

    kappa : float
        specify the kappa value of generated sequence. Kappa reports on the charge
        patterning (between 0 and 1 if there are both positive and negative residues). (Optional)

    attempts : int
        Specify the number of times to make the sequence. Default is 20. Greater numbers
        of attempts increase the odds that a sequence will be generated but will increase
        the duration of attempting to make the sequence. (Optional)

    disorder_cutoff : float
        The disorder cutoff threshold. This ensures a sequence has a mean disorder above
        this cutoff. (Optional)

    exclude : list
        A list of residues to exclude from sequence generation. 
        Cannot exclude charged residues if FCR is specified. 

    use_weighted_probabilities : bool
        whether to use weighted probabilities for sequence generation. This can
        increase the likelihood of successfully generating a sequence but does 
        decrease sequence diversity. Default is False. (Optional)

    strict_disorder : bool
        Whether to use a strict disorder calculation. By default, sequences are allowed
        to have short regions below the disorder threshold.

    return_all_sequences : bool
        Whether to return all sequences generated. Default is False. (Optional)
    
    custom_probabilities : dict
        A dictionary of custom probabilities to use for sequence generation. 
        Default is None. (Optional)

    metapredict_version : int
        The version of metapredict to use for sequence generation. Default is 3. (Optional)

    Returns
    -----------
    generated_seq : string or list
        Returns a string that is the amino acid sequence if return_all_sequences is False.
        If return_all_sequences is True, returns a list of all generated sequences.

    """
    # check length. Note this checks min/max length as well as
    # casts a string length to an int
    goose_tools.length_check(length)

    # if hydrophobicity in kwargs, change to hydropathy
    if 'hydrophobicity' in kwargs:
        kwargs['hydropathy'] = kwargs['hydrophobicity']
        # delete the old hydrophobicity key
        del kwargs['hydrophobicity']

    # check we passed in acceptable keyword arguments. At this stage, if a keyword
    # was passed that is not found in the list passed to goose_tools.check_valid_kwargs then
    # an exception is raised. 
    goose_tools.check_valid_kwargs(kwargs, ['FCR','NCPR', 'hydropathy', 'kappa', 'disorder_cutoff',
                                  'attempts', 'exclude', 'use_weighted_probabilities', 
                                  'strict_disorder', 'return_all_sequences', 'custom_probabilities', 
                                  'metapredict_version', 'max_consecutive_ordered',
                                  'max_total_ordered', 'batch_size',
                                  'hydropathy_tolerance', 'kappa_tolerance'])
    
    # First correct kwargs. Do this first because
    # the next function that looks over kwargs values
    # can only take in corrected kwargs.
    kwargs = goose_tools.check_and_correct_props_kwargs(**kwargs)

    # now make sure that the input vals are within appropriate bounds
    goose_tools.check_props_parameters(**kwargs)

    # finally, make sure all the random kwargs shared between various things are 
    # the correct type and within the correct bounds. Raise an exception if not.
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

    # make the sequence
    
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
    Stand-alone function that takes care of creating sequences with specified
    fractions of amino acids. 

    Parameters
    ------------
    length : int
        Length of the desired disordered sequence.

    <each of the 20 amino acids> : float
        Specify the fraction of the sequence that should be made up of one or more
        of the 20 natural amino acids (e.g. A=0.2, Y=0.05) etc.

    max_aa_fractions : dict 
        Dictionary which, if provided, allows the user to over-ride the 
        fraction of a sequence which can be made up of any given amino
        acid. The passed dictionary should contain key/value pairs, where
        keys are one of the twenty standard amino acids and values is a
        float between 0 and 1. If amino acids are missing then the default
        thresholds set by GOOSE are used.  (Optional)

    disorder_cutoff : float
        The disorder cutoff threshold. Default used is 0.6.  (Optional)
    
    attempts : int
        Number of attempts that will be made to generate a given sequence. Default
        is 100.  (Optional)
    
    strict_disorder : bool
        If set to true, will not count a sequence as disordered even if a single amino
        acid falls below the cutoff value. Default = False.  (Optional)
             
    Returns
    -----------
    generated_seq : string
        Returns a string that is the amino acid sequence

    """

    # check length. Note this checks min/max length as well as
    # casts a string length to an int
    goose_tools.length_check(length)

    # check we passed in acceptable keyword arguments. At this stage, if a keyword
    # was passed that is not found in the list passed to goose_tools.check_valid_kwargs then
    # an exception is raised. 
    goose_tools.check_valid_kwargs(kwargs, ['remaining_probabilities', 'attempts', 
                                  'strict_disorder', 'disorder_cutoff', 'max_consecutive_ordered',
                                  'max_total_ordered','max_aa_fractions', 
                                 'A','C','D','E','F','G','H','I','K','L','M','N',
                                 'P','Q','R','S','T','V','W','Y', 
                                 'return_all_sequences', 'metapredict_version',
                                  'batch_size'])

    

    # First correct kwargs. Do this first because
    # the next function that looks over kwargs values
    # can only take in corrected kwargs.
    kwargs = goose_tools.check_and_correct_fracs_kwargs(**kwargs)

    fractions={}
    for f in kwargs:
        if f in ['A','C','D','E','F','G','H','I','K','L','M','N','P','Q','R','S','T','V','W','Y']:
            fractions[f] = kwargs[f]

    # now make sure that the input vals are within appropriate bounds
    goose_tools.check_fracs_parameters(**kwargs)
    
    # finally, make sure all the random kwargs shared between various things are 
    # the correct type and within the correct bounds. Raise an exception if not.
    goose_tools.check_basic_parameters(num_attempts=kwargs['attempts'],
                            strict_disorder=kwargs['strict_disorder'],
                            disorder_cutoff=kwargs['disorder_cutoff'],
                            metapredict_version=kwargs['metapredict_version'],
                            return_all_sequences=kwargs['return_all_sequences'],
                            max_consecutive_ordered=kwargs['max_consecutive_ordered'],
                            max_total_ordered=kwargs['max_total_ordered'],
                            batch_size=kwargs['batch_size'],
                            custom_probabilities=kwargs['remaining_probabilities'])

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

    # return the seq
    return generated_seq

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
                    metapredict_version=parameters.METAPREDICT_DEFAULT_VERSION,
                    max_consecutive_ordered=parameters.ALLOWED_CONSECUTIVE_ORDERED,
                    max_total_ordered=parameters.ALLOWED_TOTAL_ORDERED_FRACTION,
                    remaining_probabilities=None):
    """
    Stand-alone function that takes care of creating sequences with specified
    classes of amino acids.

    Parameters
    ------------
    length : int
        Length of the desired disordered sequence. (Required)
    aromatic : float
        Fraction of aromatic amino acids in the sequence (between 0 and 1). (Optional)
    aliphatic : float
        Fraction of aliphatic amino acids in the sequence (between 0 and 1). (Optional)
    polar : float
        Fraction of polar amino acids in the sequence (between 0 and 1). (Optional)
    positive : float
        Fraction of positively charged amino acids in the sequence (between 0 and 1). (Optional)
    negative : float
        Fraction of negatively charged amino acids in the sequence (between 0 and 1). (Optional)
    glycine : float
        Fraction of glycine in the sequence (between 0 and 1). (Optional)
    proline : float
        Fraction of proline in the  sequence (between 0 and 1). (Optional)  
    cysteine : float
        Fraction of cysteine in the sequence (between 0 and 1). (Optional)
    histidine : float
        Fraction of histidine in the sequence (between 0 and 1). (Optional)
    num_attempts : int
        Number of attempts to generate a sequence. Default is 10. (Optional)
    strict_disorder : bool
        If set to true, will not count a sequence as disordered even if a single amino
        acid falls below the disorder cutoff value. Default = False. (Optional)
    disorder_cutoff : float
        The disorder cutoff threshold. Default used is 0.6. (Optional)
    metapredict_version : int
        The version of metapredict to use for sequence generation. Default is 3. (Optional)
    max_consecutive_ordered : int
        The maximum number of consecutive ordered residues allowed in the sequence.
        Default is 5. (Optional)
    max_total_ordered : float
        The maximum fraction of ordered residues allowed in the sequence.
        Default is 0.2. (Optional)
    remaining_probabilities : dict
        A dictionary of custom probabilities to use for sequence generation.

    Returns
    -----------
    generated_seq : string
        Returns a string that is the amino acid sequence.
    """
    # check length. Note this checks min/max length as well as
    # casts a string length to an int
    goose_tools.length_check(length)

    # make sure the values specified for each class are within the bounds.
    goose_tools.check_class_values(aromatic=aromatic, aliphatic=aliphatic, polar=polar,
                        positive=positive, negative=negative, glycine=glycine,
                        proline=proline, cysteine=cysteine, histidine=histidine)
    
    # check basic params
    goose_tools.check_basic_parameters(num_attempts=num_attempts,
                        strict_disorder=strict_disorder,
                        disorder_cutoff=disorder_cutoff,
                        metapredict_version=metapredict_version,
                        max_consecutive_ordered=max_consecutive_ordered,
                        max_total_ordered=max_total_ordered,
                        custom_probabilities=remaining_probabilities)

    # make the sequence
    generated_seq = sequence_generation.by_classes(
        length, aromatic=aromatic,
        aliphatic=aliphatic,
        polar=polar,
        positive=positive,
        negative=negative,
        glycine=glycine,
        proline=proline,
        cysteine=cysteine,
        histidine=histidine,
        num_attempts=num_attempts, strict_disorder=strict_disorder,
        disorder_cutoff=disorder_cutoff, metapredict_version=metapredict_version,
        max_consecutive_ordered=max_consecutive_ordered,
        max_total_ordered=max_total_ordered,
        remaining_probabilities=remaining_probabilities)

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
            max_total_ordered=parameters.ALLOWED_TOTAL_ORDERED_FRACTION):
    '''
    Parameters
    ----------
    length: int
        Length of sequence to generate
    objective_re: float
        the wanted Re in Å
    allowed_error: float
        Allowed error between the specified radius of gyration and the actual radius of gyration
        default is from the backend.parameters module, re_error or rg_error
    attempts : Int
        The number of times to attempt to build a sequence before throwing
        in the towel
    disorder_cutoff : Float
        The value for a residue to be considered disordered.
    strict_disorder : Bool
        Whether to have a strict cutoff for disorder where if any single residue
        falls below cutoff, the sequence is not considered disordered.
        Set to False by default allowing single residues to occassionally drop below
        the disorder theshold provided it is minimal. See check_disorder for more
        details.
    individual_rg_re_attempts : int
        Number of attempts to make the objective rg or re. 
        Does not account for disorder
    reduce_pos_charged : bool
        Whether to reduce positively charged amino acids in the sequence. 
        default is True
        Reason for this is that in vivo data suggests that positively charged
        residues may not drive sequence expansion as much as was predicted by
        the model used here for predicted rg / re. Therefore, when set to True,
        this function will largely avoid (+) charged residues. 
    exclude_aas : list
        list of amino acids to exclude from the sequence. 
        default is None

    Returns
    -------
    sequence : str
        The generated sequence.
    '''
    # check length. Note this checks min/max length as well as
    # casts a string length to an int
    goose_tools.length_check(length)

    # check that the rg or re range is possible. 
    if objective_re < parameters.get_min_re(length) or objective_re > parameters.get_max_re(length):
        min_possible_value=parameters.get_min_re(length)
        max_possible_value=parameters.get_max_re(length)
        raise goose_exceptions.GooseInputError(f'Cannot generate sequence, for length {length}, min Re = {min_possible_value}, max Re = {max_possible_value}.')

    if allowed_error < 0:
        raise goose_exceptions.GooseInputError('Allowed error must be a positive number.')
    
    if not isinstance(reduce_pos_charged, bool):
        raise goose_exceptions.GooseInputError('reduce_pos_charged must be a boolean value.')

    # check basic params
    goose_tools.check_basic_parameters(num_attempts=attempts,
                            strict_disorder=strict_disorder,
                            disorder_cutoff=disorder_cutoff,
                            metapredict_version=metapredict_version,
                            max_consecutive_ordered=max_consecutive_ordered,
                            max_total_ordered=max_total_ordered,
                            exclude=exclude_aas)
    
    # try to make the sequence.
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
            max_total_ordered=parameters.ALLOWED_TOTAL_ORDERED_FRACTION):
    '''
    Parameters
    ----------
    length: int
        Length of sequence to generate
    objective_rg: float
        the wanted Re in Å
    allowed_error: float
        Allowed error between the specified radius of gyration and the actual radius of gyration
        default is from the backend.parameters module, re_error or rg_error
    attempts : Int
        The number of times to attempt to build a sequence before throwing
        in the towel
    disorder_cutoff : Float
        The value for a residue to be considered disordered.
    strict_disorder : Bool
        Whether to have a strict cutoff for disorder where if any single residue
        falls below cutoff, the sequence is not considered disordered.
        Set to False by default allowing single residues to occassionally drop below
        the disorder theshold provided it is minimal. See check_disorder for more
        details.
    individual_rg_rg_attempts : int
        Number of attempts to make the objective rg or re. 
        Does not account for disorder
    reduce_pos_charged : bool
        Whether to reduce positively charged amino acids in the sequence. 
        default is True
        Reason for this is that in vivo data suggests that positively charged
        residues may not drive sequence expansion as much as was predicted by
        the model used here for predicted rg / re. Therefore, when set to True,
        this function will largely avoid (+) charged residues. 
    exclude_aas : list
        list of amino acids to exclude from the sequence. 
        default is None

    Returns
    -------
    sequence : str
        The generated sequence.
    '''
    # check length. Note this checks min/max length as well as
    # casts a string length to an int
    goose_tools.length_check(length)

    # check that the rg or re range is possible. 
    if objective_rg < parameters.get_min_rg(length) or objective_rg > parameters.get_max_rg(length):
        min_possible_value=parameters.get_min_rg(length)
        max_possible_value=parameters.get_max_rg(length)
        raise goose_exceptions.GooseInputError(f'Cannot generate sequence, for length {length}, min Rg = {min_possible_value}, max Rg = {max_possible_value}.')

    if allowed_error < 0:
        raise goose_exceptions.GooseInputError('Allowed error must be a positive number.')
    
    if not isinstance(reduce_pos_charged, bool):
        raise goose_exceptions.GooseInputError('reduce_pos_charged must be a boolean value.')

    # check basic params
    goose_tools.check_basic_parameters(num_attempts=attempts,
                            strict_disorder=strict_disorder,
                            disorder_cutoff=disorder_cutoff,
                            metapredict_version=metapredict_version,
                            max_consecutive_ordered=max_consecutive_ordered,
                            max_total_ordered=max_total_ordered,
                            exclude=exclude_aas)
    
    # try to make the sequence.
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



'''
/-/-/-/-/-/-/-/-/-/-/-/-/-/-/-/-/-/-/
/-/-/-/-/-/-/-/-/-/-/-/-/-/-/-/-/-/-/
          VARIANT GENERATORS
/-/-/-/-/-/-/-/-/-/-/-/-/-/-/-/-/-/-/
/-/-/-/-/-/-/-/-/-/-/-/-/-/-/-/-/-/-/
'''

def variant(sequence,
            variant_type,
            **kwargs):
    """
    Stand-alone function that takes care of creating variants of a sequence.
    
    Parameters
    ------------
    sequence : str
        The amino acid sequence to generate variants from. (Required)
    
    variant_type : str
        The type of variant to generate. Available options:
        - 'shuffle_specific_regions': Shuffle specified regions
        - 'shuffle_except_specific_regions': Shuffle all except specified regions
        - 'shuffle_specific_residues': Shuffle specific residues
        - 'shuffle_except_specific_residues': Shuffle all except specific residues
        - 'weighted_shuffle_specific_residues': Weighted shuffle of specific residues
        - 'targeted_reposition_specific_residues': Reposition specific residues
        - 'change_residue_asymmetry': Change residue asymmetry
        - 'constant_properties': Generate variant with constant properties
        - 'constant_residues_and_properties': Keep specified residues and properties constant
        - 'constant_properties_and_class': Generate variant with constant properties and class
        - 'constant_properties_and_class_by_order': Generate variant with constant properties and class by order
        - 'change_hydropathy_constant_class': Change hydropathy while keeping class constant
        - 'change_fcr_minimize_class_changes': Change FCR while minimizing class changes
        - 'change_ncpr_constant_class': Change NCPR while keeping class constant
        - 'change_kappa': Change kappa value
        - 'change_properties_minimize_differences': Change properties while minimizing differences
        - 'change_any_properties': Change any combination of properties
        - 'change_dimensions': Change sequence dimensions (Rg/Re)
        (Required)
    
    **kwargs : dict
        Additional parameters specific to the variant type:
        
        Common parameters:
        - num_attempts (int): Number of attempts to generate variant (default: 50)
        - strict_disorder (bool): Whether to use strict disorder checking (default: False)
        - disorder_cutoff (float): Disorder cutoff threshold (default: 0.5)
        - metapredict_version (int): MetaPredict version to use (default: 3)
        - hydropathy_tolerance (float): Hydropathy tolerance (default: from parameters)
        - kappa_tolerance (float): Kappa tolerance (default: from parameters)
        
        Variant-specific parameters:
        - shuffle_regions (list): Regions to shuffle (for shuffle_specific_regions)
        - excluded_regions (list): Regions to exclude (for shuffle_except_specific_regions)
        - target_residues (list): Residues to target (for various methods)
        - excluded_residues (list): Residues to exclude (for shuffle_except_specific_residues)
        - shuffle_weight (float): Weight for shuffling (for weighted_shuffle_specific_residues)
        - num_changes (int): Number of changes to make (for change_residue_asymmetry)
        - increase_or_decrease (str): Direction of change (for change_residue_asymmetry, change_dimensions)
        - exclude_residues (list): Residues to exclude (for constant_properties)
        - constant_residues (list): Residues to keep constant (for constant_residues_and_properties)
        - target_hydropathy (float): Target hydropathy value
        - target_FCR (float): Target FCR value
        - target_NCPR (float): Target NCPR value
        - target_kappa (float): Target kappa value
        - rg_or_re (str): Whether to optimize 'rg' or 're' (for change_dimensions)
        - num_dim_attempts (int): Number of dimension attempts (for change_dimensions)
        - allowed_error (float): Allowed error for dimensions (for change_dimensions)
        - reduce_pos_charged (bool): Whether to reduce positive charges (for change_dimensions)
        - exclude_aas (list): Amino acids to exclude (for change_dimensions)
        
    Returns
    -------
    str or None
        The generated variant sequence, or None if generation failed.
        
    Raises
    ------
    GooseInputError
        If invalid parameters are provided.
    GooseFail
        If variant generation fails.
    """
    
    # Validate sequence
    if not isinstance(sequence, str):
        raise goose_exceptions.GooseInputError('Sequence must be a string')
    
    if len(sequence) == 0:
        raise goose_exceptions.GooseInputError('Sequence cannot be empty')
    
    # Validate variant type
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
    
    if variant_type not in valid_variant_types:
        raise goose_exceptions.GooseInputError(f'Invalid variant_type: {variant_type}. Must be one of: {", ".join(sorted(valid_variant_types))}')
    
    # Extract common parameters with defaults
    common_params = {
        'num_attempts': kwargs.get('num_attempts', 100),
        'strict_disorder': kwargs.get('strict_disorder', False),
        'disorder_cutoff': kwargs.get('disorder_cutoff', parameters.DISORDER_THRESHOLD),
        'metapredict_version': kwargs.get('metapredict_version', parameters.METAPREDICT_DEFAULT_VERSION),
        'hydropathy_tolerance': kwargs.get('hydropathy_tolerance', parameters.MAXIMUM_HYDRO_ERROR),
        'kappa_tolerance': kwargs.get('kappa_tolerance', parameters.MAXIMUM_KAPPA_ERROR)
    }
    
    # Create VariantGenerator instance
    generator = VariantGenerator(**common_params)
    
    # Method dispatch dictionary
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
    
    # Get method info
    method_info = method_dispatch[variant_type]
    method_name = method_info[0]
    required_params = method_info[1]
    optional_params = method_info[2] if len(method_info) > 2 else []
    
    # Validate required parameters
    for param in required_params:
        if param not in kwargs:
            raise goose_exceptions.GooseInputError(f'Missing required parameter for {variant_type}: {param}')
    
    # Prepare method arguments
    method_args = {'input_sequence': sequence}
    
    # Add required parameters
    for param in required_params:
        method_args[param] = kwargs[param]
    
    # Add optional parameters if provided
    for param in optional_params:
        if param in kwargs:
            method_args[param] = kwargs[param]
    
    # Call the appropriate method
    try:
        method = getattr(generator, method_name)
        result = method(**method_args)
        
        if result is None:
            raise goose_exceptions.GooseFail(f'Failed to generate variant of type {variant_type}. Try adjusting parameters or increasing num_attempts.')
        
        return result
        
    except AttributeError:
        raise goose_exceptions.GooseInputError(f'Method {method_name} not found in VariantGenerator')
    except Exception as e:
        raise goose_exceptions.GooseFail(f'Error generating variant: {str(e)}')