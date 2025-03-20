##
## create.py
## 
## create.py contains all the user-facing function associated with metapredict. If a new function is added it should be included
## here and added to the __all__ list
## 

##Handles the primary functions

# if any new functions are added to create.py, you need to add them here.
__all__ =  ['seq_fractions', 'sequence', 'seq_re', 'seq_rg', 'minimal_var', 'new_seq_constant_class_var', 'constant_properties_var', 'constant_class_var', 'hydro_class_var', 'constant_residue_var', 'region_shuffle_var', 'kappa_var', 'asymmetry_var', 'fcr_class_var', 'ncpr_class_var', 'all_props_class_var', 're_var', 'rg_var', 'alpha_helix', 'beta_strand', 'beta_sheet', 'seq_property_library', 'excluded_shuffle_var', 'targeted_shuffle_var', 'targeted_reposition_var',  'weighted_shuffle_var']

import os
import sys
import random

# note - we import packages below with a leading _ which means they are ignored in the import

#for sequence generation
from goose.backend.sequence_generation import generate_disordered_seq_by_fractions as _generate_disordered_seq_by_fractions
from goose.backend.sequence_generation import generate_disordered_seq_by_props as _generate_disordered_seq_by_props
from goose.backend.sequence_generation import generate_disordered_seq_by_dimensions as _generate_disordered_seq_by_dimensions
from goose.backend.sequence_generation_backend import calculate_max_charge as _calculate_max_charge

# goose tools for checking and fixing parameters
from goose.backend.goose_tools import check_and_correct_props_kwargs as _check_and_correct_props_kwargs
from goose.backend.goose_tools import check_props_parameters as _check_props_parameters
from goose.backend.goose_tools import check_and_correct_fracs_kwargs as _check_and_correct_fracs_kwargs
from goose.backend.goose_tools import check_fracs_parameters as _check_fracs_parameters
from goose.backend.goose_tools import length_check as _length_check
from goose.backend.goose_tools import gen_random_name as _gen_random_name
from goose.backend.goose_tools import check_valid_kwargs as _check_valid_kwargs

# variant generation stuff
from goose.backend.variant_generation import gen_kappa_variant as _gen_kappa_variant
from goose.backend.variant_generation import gen_region_shuffle_variant as _gen_region_shuffle_variant
from goose.backend.variant_generation import gen_constant_residue_variant as _gen_constant_residue_variant
from goose.backend.variant_generation import gen_hydropathy_class_variant as _gen_hydropathy_class_variant
from goose.backend.variant_generation import gen_new_variant as _gen_new_variant
from goose.backend.variant_generation import gen_constant_class_variant as _gen_constant_class_variant
from goose.backend.variant_generation import gen_new_var_constant_class as _gen_new_var_constant_class
from goose.backend.variant_generation import gen_asymmetry_variant as _gen_asymmetry_variant
from goose.backend.variant_generation import gen_fcr_class_variant as _gen_fcr_class_variant
from goose.backend.variant_generation import gen_ncpr_class_variant as _gen_ncpr_class_variant
from goose.backend.variant_generation import gen_all_props_class_variant as _gen_all_props_class_variant
from goose.backend.variant_generation import gen_targeted_shuffle_variant as _gen_targeted_shuffle_variant
from goose.backend.variant_generation import gen_targeted_reposition_variant as _gen_targeted_reposition_variant
from goose.backend.variant_generation import gen_weighted_shuffle_variant as _gen_weighted_shuffle_variant
from goose.backend.variant_generation import gen_excluded_shuffle_variant as _gen_excluded_shuffle_variant
from goose.backend.variant_generation import gen_dimensions_variant as _gen_dimensions_variant

# for minimal variant, need to... get rid of this old code.
from goose.backend.gen_minimal_variant_backend import gen_minimal_sequence_variant as _gen_minimal_sequence_variant

# library creation
from goose.backend.library_generation_backend import generate_library_by_parameter_ranges as _generate_library_by_parameter_ranges
from goose.backend.library_generation_backend import generate_library_by_fraction_ranges as _generate_library_by_fraction_ranges

# for folded structure generation
from goose.backend.folded_region_generation import gen_helix as _gen_helix
from goose.backend.folded_region_generation import gen_beta_strand as _gen_beta_strand
from goose.backend.folded_region_generation import gen_beta_sheet as _gen_beta_sheet

# FOR WHEN THINGS GO WRONG
from goose import goose_exceptions

# because every good package needs some parameters
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

    cutoff : float
        The disorder cutoff threshold. This ensures a sequence has a mean disorder above
        this cutoff. (Optional)

    exclude : list
        A list of residues to exclude from sequence generation. 
        Cannot exclude charged residues if FCR is specified. 
            

    Returns
    -----------
    generated_seq : string
        Returns a string that is the amino acid sequence.

    """
    # check length. Note this checks min/max length as well as
    # casts a string length to an int
    _length_check(length)

    # if just specifying kappa, add in some charged residues. Otherwise can't adjust kappa.
    '''
    This section has expanded more or less into a hack to keep things possible as far as sequence
    space. When only a few parameters are specified, there are inherent ranges on some other parameters.
    There is a better way to do this, but there's also this way. I'm doing this way because it will be faster. 
    I'll probably rewrite this sooner or later because it's pretty inefficiently written at the moment. 
    '''

    # increase attempts if hydropathy over 5.9
    if 'hydropathy' in kwargs:
        if kwargs['hydropathy']!=None:
            if kwargs['hydropathy']>=5.9:
                if 'attempts' not in kwargs:
                    kwargs['attempts']=200
                else:
                    kwargs['attempts']=kwargs['attempts']+200

    # verify that charged residues not in exclude if FCR or NCPR specified.
    if 'exclude' in kwargs:
        if kwargs['exclude'] != None:
            if 'FCR' in kwargs or 'NCPR' in kwargs:
                for val in kwargs['exclude']:
                    if val in ['R', 'K', 'D', 'E']:
                        raise goose_exceptions.GooseInputError('Cannot exclude charged residues if FCR or NCPR specified.')
            if len(kwargs['exclude']) > 10:
                raise goose_exceptions.GooseInputError('Cannot exclude more than 10 residues.')

    if 'kappa' in kwargs:
        if kwargs['kappa'] != None:
            if 'FCR' in kwargs and 'NCPR' in kwargs:
                if kwargs['FCR'] != None and kwargs['NCPR']!= None:
                    if kwargs['FCR']==kwargs['NCPR']:
                        raise goose_exceptions.GooseInputError('Cannot specify FCR and NCPR to be the same value and specify kappa. Kappa requires the presence of oppositely charged residues to be specified.')
            if 'FCR' in kwargs:
                if kwargs['FCR']==0:
                    raise goose_exceptions.GooseInputError('When specifying kappa, FCR must be greater than 0. FCR must be a high enough value to result in at least 2 charged residues to be in the sequence so oppositely charged residue spacing (kappa) can be specified.')


    # check we passed in acceptable keyword arguments. At this stage, if a keyword
    # was passed that is not found in the list passed to _check_valid_kwargs then
    # an exception is raised. 
    _check_valid_kwargs(kwargs, ['FCR','NCPR','sigma', 'hydropathy', 'kappa', 'cutoff', 'attempts', 'exclude', 'no_constraints'])
    

    # First correct kwargs. Do this first because
    # the next function that looks over kwargs values
    # can only take in corrected kwargs.
    kwargs = _check_and_correct_props_kwargs(**kwargs)


    # now make sure that the input vals are within appropriate bounds
    _check_props_parameters(**kwargs)

    # add exclude to kwargs if not in. 
    if 'exclude' not in kwargs:
        kwargs['exclude']=[]

    # placeholder for sequence we haven't managed to make
    generated_seq=None
    # make the sequence
    if kwargs['kappa'] == None:
        for sub_attempt in range(0, 20):
            try:
                generated_seq = _generate_disordered_seq_by_props(length, FCR=kwargs['FCR'], NCPR=kwargs['NCPR'], hydropathy=kwargs['hydropathy'],
                    sigma = kwargs['sigma'], attempts = kwargs['attempts'], allowed_hydro_error = parameters.HYDRO_ERROR,
                    disorder_threshold = kwargs['cutoff'], exclude=kwargs['exclude'])
            except:
                continue
            # return the seq
            if generated_seq != None:
                return generated_seq
    
    else:        
        for sub_attempt in range(0, 30):
            # if specifying kappa and FCR or NCPR are not specified, this code
            # helps use increase our odds of making a sequence successfully.
            # You can disable these constraints using 'no_constraints=True.'
            if 'no_constraints' not in kwargs:
                remove_constraints=False
            else:
                if kwargs['no_constraints']==False:
                    remove_constraints=False
                else:
                    remove_constraints=True
            # if we have not set remove_constraints to True...
            if remove_constraints==False:
                # first let's see if we can modify FCR. 
                if kwargs['FCR']!=None:
                    fcr_val = kwargs['FCR']
                else:
                    # this is hydropathy dependent, so check that. 
                    if kwargs['hydropathy']==None:
                        # set max number of charged resiudes to be the length.
                        max_FCR = 1
                    else:
                        # set max to be dependent on hydropathy. 
                        max_FCR = min([1, _calculate_max_charge(kwargs['hydropathy'])])
                        if max_FCR-0.04 < 0.2:
                            max_FCR=0.21
                        else:
                            max_FCR = max_FCR-0.04
                    # check for NCPR. This will help us decide a minimum FCR val.
                    if kwargs['NCPR']==None:
                        min_FCR = 0.2
                    else:
                        absncpr=abs(kwargs['NCPR'])
                        if absncpr>=0.75:
                            min_FCR=1
                        elif absncpr>=0.5 and absncpr<0.75:
                            min_FCR = min([absncpr+0.35,0.99,1])
                        elif absncpr >= 0.25 and absncpr < 0.5:
                            min_FCR = min([absncpr+0.25, 1])
                        else:
                            min_FCR = min([absncpr+0.15, 1])
                    # make sure min_FCR not greater than max FCR. Then set FCR val. 
                    if min_FCR >= max_FCR:
                        fcr_val = max_FCR
                    else:
                        fcr_val = random.uniform(min_FCR, max_FCR)
                    
                # now let's see if we can modify NCPR. 
                if kwargs['NCPR']!=None:
                    ncpr_val = kwargs['NCPR']
                else:
                    # we can modulate the NCPR. 
                    if length <= 15:
                        ncpr_val=0.0
                    else:
                        min_ncpr = -fcr_val+(2/length)+0.05
                        if min_ncpr < -0.3:
                            min_ncpr=-0.3
                        max_ncpr = fcr_val-(2/length)-0.05
                        if max_ncpr > 0.3:
                            max_ncpr=0.3
                        if min_ncpr >= max_ncpr:
                            ncpr_val = max_ncpr
                        else:
                            ncpr_val = random.uniform(min_ncpr, max_ncpr)
                # modulate the NCPR value if we have a specified FCR but not NCPR. 
                if kwargs['FCR']!=None and kwargs['NCPR']==None:
                    FCR_res = length*fcr_val
                    ncpr_res = abs(ncpr_val)*length
                    if round(FCR_res-ncpr_res)%2 != 0:
                        if ncpr_res < FCR_res:
                            ncpr_res+=1
                        else:
                            ncpr_res-=1
                    if ncpr_val < 0:
                        ncpr_val = -ncpr_res/length
                    else:
                        ncpr_val = ncpr_res/length
            else:
                # make sure we have fcr_val, ncpr_val specified.
                if 'FCR' in kwargs:
                    fcr_val = kwargs['FCR']
                else:
                    fcr_val = None
                if 'NCPR' in kwargs:
                    ncpr_val = kwargs['NCPR']
                else:
                    ncpr_val = None


            generated_seq=None
            try:
                generated_seq = _generate_disordered_seq_by_props(length, FCR=fcr_val, NCPR=ncpr_val, hydropathy=kwargs['hydropathy'],
                    sigma = kwargs['sigma'], attempts = kwargs['attempts'], allowed_hydro_error = parameters.HYDRO_ERROR,
                    disorder_threshold = kwargs['cutoff'], exclude=kwargs['exclude'])
            except:
                continue
            if generated_seq != None:
                try:
                    generated_seq = _gen_kappa_variant(generated_seq, kappa=kwargs['kappa'], 
                        allowed_kappa_error = parameters.MAXIMUM_KAPPA_ERROR, 
                        disorder_threshold=kwargs['cutoff'], strict_disorder=False)
                except:
                    continue
            if generated_seq != None:
                return generated_seq


    raise goose_exceptions.GooseFail('Unable to generate sequence. Please try again with different parameters or a lower cutoff value.')


#-\|/-\|/-\|/-\|/-\|/-\|/-\|/-\|/-\|/-\|/-\|/-\|/-\|/-\|/-\|/-\|/-\|/-\|/-\|/-\|/-\|/-
#-\|/-\|/-\|/-\|/-\|/-\|/-\|/-\|/-\|/             \|/-\|/-\|/-\|/-\|/-\|/-\|/-\|/-\|/-
#-\|/-\|/-\|/-\|/-\|/-\|/-\|/-\|/-\|/  Create     \|/-\|/-\|/-\|/-\|/-\|/-\|/-\|/-\|/-
#-\|/-\|/-\|/-\|/-\|/-\|/-\|/-\|/-\|/  sequence   \|/-\|/-\|/-\|/-\|/-\|/-\|/-\|/-\|/-
#-\|/-\|/-\|/-\|/-\|/-\|/-\|/-\|/-\|/  By         \|/-\|/-\|/-\|/-\|/-\|/-\|/-\|/-\|/-
#-\|/-\|/-\|/-\|/-\|/-\|/-\|/-\|/-\|/  Specifying \|/-\|/-\|/-\|/-\|/-\|/-\|/-\|/-\|/-
#-\|/-\|/-\|/-\|/-\|/-\|/-\|/-\|/-\|/  Fractions  \|/-\|/-\|/-\|/-\|/-\|/-\|/-\|/-\|/-
#-\|/-\|/-\|/-\|/-\|/-\|/-\|/-\|/-\|/             \|/-\|/-\|/-\|/-\|/-\|/-\|/-\|/-\|/-
#-\|/-\|/-\|/-\|/-\|/-\|/-\|/-\|/-\|/-\|/-\|/-\|/-\|/-\|/-\|/-\|/-\|/-\|/-\|/-\|/-\|/-


def seq_fractions(length, **kwargs):
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

    cutoff : float
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
    _length_check(length)

    # check we passed in acceptable keyword arguments. At this stage, if a keyword
    # was passed that is not found in the list passed to _check_valid_kwargs then
    # an exception is raised. 
    _check_valid_kwargs(kwargs, ['cutoff', 'attempts',  'strict_disorder',  'max_aa_fractions', 'A','C','D','E','F','G','H','I','K','L','M','N','P','Q','R','S','T','V','W','Y'])

    # First correct kwargs. Do this first because
    # the next function that looks over kwargs values
    # can only take in corrected kwargs.
    kwargs = _check_and_correct_fracs_kwargs(**kwargs)

    # now make sure that the input vals are within appropriate bounds
    _check_fracs_parameters(**kwargs)

    generated_seq = _generate_disordered_seq_by_fractions(length, **kwargs)

    # return the seq
    return generated_seq


#-\|/-\|/-\|/-\|/-\|/-\|/-\|/-\|/-\|/-\|/-\|/-\|/-\|/-\|/-\|/-\|/-\|/-\|/-\|/-\|/-\|/-
#-\|/-\|/-\|/-\|/-\|/-\|/-\|/-\|/-\|/             \|/-\|/-\|/-\|/-\|/-\|/-\|/-\|/-\|/-
#-\|/-\|/-\|/-\|/-\|/-\|/-\|/-\|/-\|/  Create     \|/-\|/-\|/-\|/-\|/-\|/-\|/-\|/-\|/-
#-\|/-\|/-\|/-\|/-\|/-\|/-\|/-\|/-\|/  sequence   \|/-\|/-\|/-\|/-\|/-\|/-\|/-\|/-\|/-
#-\|/-\|/-\|/-\|/-\|/-\|/-\|/-\|/-\|/  By         \|/-\|/-\|/-\|/-\|/-\|/-\|/-\|/-\|/-
#-\|/-\|/-\|/-\|/-\|/-\|/-\|/-\|/-\|/  dimensions \|/-\|/-\|/-\|/-\|/-\|/-\|/-\|/-\|/-
#-\|/-\|/-\|/-\|/-\|/-\|/-\|/-\|/-\|/             \|/-\|/-\|/-\|/-\|/-\|/-\|/-\|/-\|/-
#-\|/-\|/-\|/-\|/-\|/-\|/-\|/-\|/-\|/-\|/-\|/-\|/-\|/-\|/-\|/-\|/-\|/-\|/-\|/-\|/-\|/-

def seq_re(length, objective_re, allowed_error=parameters.re_error, attempts=20, 
    cutoff=parameters.DISORDER_THRESHOLD, strict_disorder=False,
    individual_rg_re_attempts=parameters.rg_re_attempt_num,
    reduce_pos_charged=True, exclude_aas=None):
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
    cutoff : Float
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
    _length_check(length)

    # check that the rg or re range is possible. 
    if objective_re < parameters.get_min_re(length) or objective_re > parameters.get_max_re(length):
        min_possible_value=parameters.get_min_re(length)
        max_possible_value=parameters.get_max_re(length)
        raise goose_exceptions.GooseInputError(f'Cannot generate sequence, for length {length}, min Re = {min_possible_value}, max Re = {max_possible_value}.')

    # try to make the sequence.
    sequence=_generate_disordered_seq_by_dimensions(
        length, 're', objective_re, attempts=attempts, allowed_error=allowed_error,
        disorder_threshold=cutoff, strict_disorder=strict_disorder,
        individual_rg_re_attempts=individual_rg_re_attempts,
        reduce_pos_charged=reduce_pos_charged, exclude_aas=exclude_aas)
    return sequence



def seq_rg(length, objective_rg, allowed_error=parameters.rg_error, attempts=20, 
    cutoff=parameters.DISORDER_THRESHOLD, strict_disorder=False,
    individual_rg_re_attempts=parameters.rg_re_attempt_num,
    reduce_pos_charged=True, exclude_aas=None):
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
    cutoff : Float
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
    _length_check(length)

    # check that the rg or re range is possible. 
    if objective_rg < parameters.get_min_rg(length) or objective_rg > parameters.get_max_rg(length):
        min_possible_value=parameters.get_min_rg(length)
        max_possible_value=parameters.get_max_rg(length)
        raise goose_exceptions.GooseInputError(f'Cannot generate sequence, for length {length}, min Rg = {min_possible_value}, max Rg = {max_possible_value}.')

    # try to make the sequence.
    sequence=_generate_disordered_seq_by_dimensions(
        length, 'rg', objective_rg, attempts=attempts, allowed_error=allowed_error,
        disorder_threshold=cutoff, strict_disorder=strict_disorder,
        individual_rg_re_attempts=individual_rg_re_attempts,
        reduce_pos_charged=reduce_pos_charged, exclude_aas=None)
    return sequence


'''
/-/-/-/-/-/-/-/-/-/-/-/-/-/-/-/-/-/-/
/-/-/-/-/-/-/-/-/-/-/-/-/-/-/-/-/-/-/
          VARIANT GENERATORS
/-/-/-/-/-/-/-/-/-/-/-/-/-/-/-/-/-/-/
/-/-/-/-/-/-/-/-/-/-/-/-/-/-/-/-/-/-/
'''



def minimal_var(input_sequence, hydropathy = '', FCR = '', 
    NCPR = '', SCD='', cutoff=parameters.DISORDER_THRESHOLD, strict=False):
    '''
    User facing function for generating the minimal sequence variant. This variant
    tries to make a sequence as similar to the input sequence as possible all while
    minimizing the number of amino acids changed. 
    
    Parameters
    ------------
    input_sequence : str
        The sequence to make a variant of.

    hydropathy : float
        The mean hydropathy of the sequence. If not specified does not change. (Optional)

    FCR : float
        The fraction of charged residues in the sequence. If not specified does not change. (Optional)

    NCPR : float
        The net charge per residue of the sequence. If not specified does not change. (Optional)

    SCD : float
        The charge asymmetry of the sequence. If not specified does not change. (Optional)

    cutoff : float
        The disorder cutoff threshold. Closer to 1 has a higher chance of being disordered.
        Must be betwee 0 and 1. (Optional)

    strict : bool
        Whether to use a strict disorder calculation. By default, variants are allowed
        to have regions below the disorder threshold *for regions where the input sequence
        is also below the threshold*. (Optional)

    Returns
    -----------
    final_sequence : str
        The final sequence variant.
    '''

    # make sure that the input sequence is all caps
    input_sequence = input_sequence.upper()

    # check length
    _length_check(input_sequence)

    if cutoff > 1 or cutoff < 0:
        raise goose_exceptions.GooseInputError('cutoff value must be between 0 and 1 for disorder threshold')    
    try:
        final_sequence = _gen_minimal_sequence_variant(input_sequence, mean_hydro = hydropathy, fraction = FCR, 
        net_charge = NCPR, charge_asymmetry=SCD, cutoff=cutoff, strict=strict)
    except:
        raise goose_exceptions.GooseFail('Sorry! GOOSE was unable to generate the sequence. Please try again or try with different parameters or a lower cutoff value.')
    return final_sequence



def new_seq_constant_class_var(sequence, attempts=5,
    cutoff = parameters.DISORDER_THRESHOLD, strict=False):
    '''
    A function to generate a variant where the sequence composition is new but
    the numbers of each residue from each class is the same. The overall properties
    of the generated sequence will also be constant.

    Parameters
    ------------
    sequence : str
        The sequence to make a variant of.

    attempts : int
        Specify the number of times to make the sequence. Default is 5. Greater numbers
        of attempts increase the odds that a sequence will be generated but will increase
        the duration of attempting to make the sequence. (Optional)

    cutoff : float
        the disorder cutoff threshold. Closer to 1 has a higher chance of being disordered. (Optional)

    strict : bool
        Whether to use a strict disorder calculation. By default, variants are allowed
        to have regions below the disorder threshold *for regions where the input sequence
        is also below the threshold*. (Optional)

    Returns
    -----------
    final_sequence : str
        The final sequence variant

    '''
    # make sure that the input sequence is all caps
    sequence = sequence.upper()

    # check length
    _length_check(sequence)

    if cutoff > 1 or cutoff < 0:
        raise goose_exceptions.GooseInputError('cutoff value must be between 0 and 1 for disorder threshold')
    try:
        final_sequence = _gen_new_var_constant_class(sequence, attempts=attempts, disorder_threshold=cutoff, strict_disorder=strict)
    except:
        raise goose_exceptions.GooseFail('Sorry! GOOSE was unable to generate the sequence. Please try again or try with a lower cutoff value.')
    return final_sequence



def constant_class_var(sequence, attempts=5,
    cutoff = parameters.DISORDER_THRESHOLD, strict=False):
    '''
    Function to generate a variant with the same properties as the 
    input variant as well as the same order of amino acids as
    far as class and the same number in each class.
    
    Parameters
    ------------
    sequence : str
        The sequence to make a variant of

    attempts : int
        Specify the number of times to make the sequence. Default is 5. Greater numbers
        of attempts increase the odds that a sequence will be generated but will increase
        the duration of attempting to make the sequence. (Optional)  

    cutoff : float
        The disorder cutoff threshold. Closer to 1 has a higher chance of being disordered. (Optional)

    strict : bool
        Whether to use a strict disorder calculation. By default, variants are allowed
        to have regions below the disorder threshold *for regions where the input sequence
        is also below the threshold*. (Optional)

    Returns
    -----------
    final_sequence : str
        The final sequence variant.
    '''
    # make sure that the input sequence is all caps
    sequence = sequence.upper()

    # check length
    _length_check(sequence)

    if cutoff > 1 or cutoff < 0:
        raise goose_exceptions.GooseInputError('cutoff value must be between 0 and 1 for disorder threshold')    
    try:
        final_sequence = _gen_constant_class_variant(sequence, attempts=attempts, disorder_threshold=cutoff, strict_disorder=strict)
    except:
        raise goose_exceptions.GooseFail('Sorry! GOOSE was unable to generate the sequence. Please try again or try with a lower cutoff value.')
    return final_sequence


def constant_properties_var(sequence, attempts=5,
    cutoff = parameters.DISORDER_THRESHOLD, strict=False):
    '''
    Function to generate a variant that is completely different
    in sequence to the input but has all the same overall parameters.
    Does not account for specific classes of residues.
    
    Parameters
    ------------
    sequence : str
        The sequence to make a variant of.

    attempts : int
        Specify the number of times to make the sequence. Default is 5. Greater numbers
        of attempts increase the odds that a sequence will be generated but will increase
        the duration of attempting to make the sequence. (Optional)

    cutoff : float
        The disorder cutoff threshold. Closer to 1 has a higher chance of being disordered. (Optional)

    strict : bool
        Whether to use a strict disorder calculation. By default, variants are allowed
        to have regions below the disorder threshold *for regions where the input sequence
        is also below the threshold*. (Optional)

    Returns
    -----------
    final_sequence : str
        The final sequence Variant.
    '''
    # make sure that the input sequence is all caps
    sequence = sequence.upper()

    # check length
    _length_check(sequence)

    if cutoff > 1 or cutoff < 0:
        raise goose_exceptions.GooseInputError('cutoff value must be between 0 and 1 for disorder threshold')    
    try:
        final_sequence = _gen_new_variant(sequence, attempts=attempts, disorder_threshold=cutoff, strict_disorder=strict)
    except:
        raise goose_exceptions.GooseFail('Sorry! GOOSE was unable to generate the sequence. Please try again or try with a lower cutoff value.')
    return final_sequence


def hydro_class_var(sequence, hydropathy, hydro_error = parameters.HYDRO_ERROR,
    attempts=5, cutoff = parameters.DISORDER_THRESHOLD, strict=False):
    '''
    Function to take in a sequence and make a variant that adjusts the
    hydropathy while keeping the position and number of amino acids the
    same by class of amino acid.

    Parameters
    ------------
    sequence : str
        The sequence to make a variant of.
    
    hydropathy : float
        The mean hydropathy of the sequence variant. 

    hydro_erorr : float
        The allowed error between the specified hydropathy and the
        hydropathy of the returned sequence. (Optional)

    attempts : int
        Specify the number of times to make the sequence. Default is 5. Greater numbers
        of attempts increase the odds that a sequence will be generated but will increase
        the duration of attempting to make the sequence. (Optional)   

    cutoff : float
        The disorder cutoff threshold. Closer to 1 has a higher chance of being disordered. (Optional)

    strict : bool
        Whether to use a strict disorder calculation. By default, variants are allowed
        to have regions below the disorder threshold *for regions where the input sequence
        is also below the threshold*. (Optional)

    Returns
    -----------
    final_sequence : str
        The final sequence variant.

    '''
    # make sure that the input sequence is all caps
    sequence = sequence.upper()

    # check length
    _length_check(sequence)

    if cutoff > 1 or cutoff < 0:
        raise goose_exceptions.GooseInputError('cutoff value must be between 0 and 1 for disorder threshold')    
    try:
        final_sequence = _gen_hydropathy_class_variant(sequence, hydropathy=hydropathy, allowed_hydro_error = hydro_error, attempts=attempts, disorder_threshold=cutoff, strict_disorder=strict)
    except:
        raise goose_exceptions.GooseFail('Sorry! GOOSE was unable to generate the sequence. Please try again or try with a lower cutoff value.')
    return final_sequence




def constant_residue_var(sequence, constant=[], attempts=5, 
    cutoff = parameters.DISORDER_THRESHOLD, strict=False):
    '''
    function that will generate a new sequence variant
    where specific residues are held constant. The 
    variant will have the same aggregate properties
    as the original sequence.

    Parameters
    ------------
    sequence : str
        The sequence to make a variant of.

    constant : list
        A list of residues in the sequence to hold constant.

    attempts : int
        Specify the number of times to make the sequence. Default is 5. Greater numbers
        of attempts increase the odds that a sequence will be generated but will increase
        the duration of attempting to make the sequence. (Optional)

    cutoff : float
        the disorder cutoff threshold. Closer to 1 has a higher chance of being disordered. (Optional)

    strict : bool
        Whether to use a strict disorder calculation. By default, variants are allowed
        to have regions below the disorder threshold *for regions where the input sequence
        is also below the threshold*. (Optional)

    Returns
    -----------
    final_sequence : str
        the final sequence variant 

    '''
    # make sure that the input sequence is all caps
    sequence = sequence.upper()

    # check length
    _length_check(sequence)

    if cutoff > 1 or cutoff < 0:
        raise goose_exceptions.GooseInputError('cutoff value must be between 0 and 1 for disorder threshold')
    try:
        final_sequence = _gen_constant_residue_variant(sequence, constant_residues=constant, attempts=attempts, disorder_threshold=cutoff, strict_disorder=strict)
    except:
        raise goose_exceptions.GooseFail('Sorry! GOOSE was unable to generate the sequence. Please try again or try with a lower cutoff value or with different constant residues.')
    return final_sequence



def region_shuffle_var(sequence, shuffle=[], attempts=5,
    cutoff = parameters.DISORDER_THRESHOLD, strict=False):
    '''
    Function that will shuffle specific regions of an IDR.
    Multiple regions can be specified simultaneously.

    Parameters
    ------------
    sequence : str
        The sequence to make a variant of.

    shuffle : list
        A list of regions to shuffle.

    attempts : int
        Specify the number of times to make the sequence. Default is 5. Greater numbers
        of attempts increase the odds that a sequence will be generated but will increase
        the duration of attempting to make the sequence. (Optional) 

    cutoff : float
        The disorder cutoff threshold. Closer to 1 has a higher chance of being disordered. (Optional)

    strict : bool
        Whether to use a strict disorder calculation. By default, variants are allowed
        to have regions below the disorder threshold *for regions where the input sequence
        is also below the threshold*. (Optional)

    Returns
    -----------
    final_sequence : str
        The final sequence variant.
    '''

    # make sure that the input sequence is all caps
    sequence = sequence.upper()

    # check length
    _length_check(sequence)

    if cutoff > 1 or cutoff < 0:
        raise goose_exceptions.GooseInputError('cutoff value must be between 0 and 1 for disorder threshold')    

    all_vals = []
    if shuffle != []:
        for val in shuffle:
            if type(val) == list:
                for subval in val:
                    all_vals.append(subval)
            else:
                all_vals.append(val)
    
    for val in all_vals:
        if val < 1:
            raise goose_exceptions.GooseInputError('Cannot have a value to shuffle below 1')
        if val > len(sequence):
            raise goose_exceptions.GooseInputError('Cannot specify to shuffle a region greater than the length of your sequence')

    curvals = []
    if type(shuffle[0])==list:
        if len(shuffle) >= 2:
            for sublist in shuffle:
                for i in range(sublist[0], sublist[1]+1):
                    if i in curvals:
                        raise goose_exceptions.GooseInputError('Cannot have overlapping regions to shuffle.')
                    else:
                        curvals.append(i)

    try:
        final_sequence = _gen_region_shuffle_variant(sequence, shuffle_regions=shuffle, attempts=attempts, disorder_threshold=cutoff, strict_disorder=strict)
    except:
        raise goose_exceptions.GooseFail('Sorry! GOOSE was unable to generate the sequence. Please try again or try with a lower cutoff value or with different constant residues.')
    return final_sequence



def kappa_var(sequence, kappa, kappa_error=parameters.MAXIMUM_KAPPA_ERROR, 
    attempts=10, cutoff=parameters.DISORDER_THRESHOLD, strict=False):
    '''
    Function to generate a sequence with a user-defined
    kappa value. Requires kappa calculation using 
    SPARROW. Kappa is a function of charge asymmetry, larger
    kappa values have more asymmetrically distributed
    charged residues.

    Parameters
    ------------
    sequence : str
        The sequence to make a variant of.

    kappa : float
        The desired kappa value. 

    kappa_error : float
        The allowed error between the desired kappa value and the kappa value of the
        returned sequence. (Optional)

    attempts : int
        Specify the number of times to make the sequence. Default is 10. Greater numbers
        of attempts increase the odds that a sequence will be generated but will increase
        the duration of attempting to make the sequence. (Optional)

    cutoff : float
        the disorder cutoff threshold. Closer to 1 has a higher chance of being disordered. (Optional)

    strict : bool
        Whether to use a strict disorder calculation. By default, variants are allowed
        to have regions below the disorder threshold *for regions where the input sequence
        is also below the threshold*. (Optional)

    Returns
    -----------
    final_sequence : str
        The final sequence variant.
    '''
    # make sure that the input sequence is all caps
    sequence = sequence.upper()

    # check length
    _length_check(sequence)

    if cutoff > 1 or cutoff < 0:
        raise goose_exceptions.GooseInputError('cutoff value must be between 0 and 1 for disorder threshold')    
    if len(sequence) < 6:
        raise goose_exceptions.GooseInputError('Cannot have sequence with a length less than 6')
    if kappa > 1 or kappa < 0:
        raise goose_exceptions.GooseInputError('Kappa values must be between 0 and 1')

    try:
        final_sequence = _gen_kappa_variant(sequence, kappa=kappa, allowed_kappa_error = kappa_error, attempts=attempts, disorder_threshold=cutoff, strict_disorder=strict)
    except:
        raise goose_exceptions.GooseFail('Sorry! GOOSE was unable to generate the sequence. Please try again or try with a lower cutoff value or with different kappa value')
    return final_sequence


def asymmetry_var(sequence, increase_decrease, aa_class, number_changes=None,
    attempts=10, cutoff = parameters.DISORDER_THRESHOLD, strict=False):
    '''
    User facing function for generating a variant that has 
    altered asymmetry of some class of residues or user specificed
    list of residues.

    Parameters
    ------------
    sequence : str
        The sequence to make a variant of.

    increase_decrease : str
        Whether to increase or decrease the asymmetry of a specific 
        residue or residue class. Set to 'increase' or 'decrease'.

    aa_class : str or list
        the class of amino acids to alter asymmetry of. Can specify a custom
        class by passing in a list of residues to target or you can just type in
        the name of a canonical class as a string. Possible classes include:
        'negative' - negative residues
        'positive' - positive residues    
        'proline' - prolines
        'aromatic' - W Y F
        'aliphatic' - I V L A M
        'polar' - Q N S T        

    number_changes : bool
        The number of amino acids to change in the sequence to change asymmetry. 
        By default, will attempt to change asymmetry a number of times equal to
        1/5th the sequence length. Higher numbers increase the amount that the sequence
        becomes more or less symmetrical for the specified residue.

    attempts : int
        Specify the number of times to make the sequence. Default is 10. Greater numbers
        of attempts increase the odds that a sequence will be generated but will increase
        the duration of attempting to make the sequence.          

    cutoff : float
        The disorder cutoff threshold. Closer to 1 has a higher chance of being disordered.

    strict : bool
        Whether to use a strict disorder calculation. By default, variants are allowed
        to have regions below the disorder threshold *for regions where the input sequence
        is also below the threshold*. (Optional)

    Returns
    -----------
    final_sequence : str
        the final sequence variant     
    '''
    # make sure that the input sequence is all caps
    sequence = sequence.upper()

    # check length
    _length_check(sequence)

    if cutoff > 1 or cutoff < 0:
        raise goose_exceptions.GooseInputError('cutoff value must be between 0 and 1 for disorder threshold')    
    if len(sequence) < 6:
        raise goose_exceptions.GooseInputError('Cannot have sequence with a length less than 6')

    try:
        final_sequence = _gen_asymmetry_variant(sequence, increase_decrease, aa_class, num_change=number_changes, attempts=attempts, disorder_threshold=cutoff, strict_disorder=strict)
    except:
        raise goose_exceptions.GooseFail('Sorry! GOOSE was unable to generate the sequence. Please try again or try with a lower cutoff value or with different specified amino acids')
    return final_sequence    


def fcr_class_var(sequence, FCR, attempts=10, cutoff=parameters.DISORDER_THRESHOLD, strict=False):
    '''
    User facing funcitonality to generate variants where
    the fcr is changed and the changes to classes of amino
    acids is minimized.
    
    Parameters
    ------------
    sequence : str
        The sequence to make a variant of.

    FCR : float
        The desired fcr value.

    attempts : int
        Specify the number of times to make the sequence. Default is 10. Greater numbers
        of attempts increase the odds that a sequence will be generated but will increase
        the duration of attempting to make the sequence. (Optional)

    cutoff : float
        the disorder cutoff threshold. Closer to 1 has a higher chance of being disordered. (Optional)

    strict : bool
        Whether to use a strict disorder calculation. By default, variants are allowed
        to have regions below the disorder threshold *for regions where the input sequence
        is also below the threshold*. (Optional)

    Returns
    -----------
    final_sequence : str
        The final sequence variant.
    '''

    # make sure that the input sequence is all caps
    sequence = sequence.upper()

    # check length
    _length_check(sequence)

    if cutoff > 1 or cutoff < 0:
        raise goose_exceptions.GooseInputError('cutoff value must be between 0 and 1 for disorder threshold')    
    if len(sequence) < 6:
        raise goose_exceptions.GooseInputError('Cannot have sequence with a length less than 6')

    if FCR > 1 or FCR < 0:
        raise goose_exceptions.GooseInputError('fcr values must be between 0 and 1.')

    try:
        final_sequence = _gen_fcr_class_variant(sequence, fcr=FCR, attempts=attempts, disorder_threshold=cutoff, strict_disorder=strict)
    except:
        raise goose_exceptions.GooseFail('Sorry! GOOSE was unable to generate the sequence. Please try again or try with a different fcr value or a different cutoff value.')
    return final_sequence       
    


def ncpr_class_var(sequence, NCPR, attempts=10, cutoff=parameters.DISORDER_THRESHOLD, strict=False):
    '''
    User facing funcitonality to generate variants where
    the ncpr is changed and the changes to classes of amino
    acids is minimized.

    Parameters
    ------------
    sequence : str
        The sequence to make a variant of.

    NCPR : float
        The desired ncpr value.

    attempts : int
        Specify the number of times to make the sequence. Default is 10. Greater numbers
        of attempts increase the odds that a sequence will be generated but will increase
        the duration of attempting to make the sequence. (Optional)

    cutoff : float
        the disorder cutoff threshold. Closer to 1 has a higher chance of being disordered. (Optional)

    strict : bool
        Whether to use a strict disorder calculation. By default, variants are allowed
        to have regions below the disorder threshold *for regions where the input sequence
        is also below the threshold*. (Optional)

    Returns
    -----------
    final_sequence : str
        The final sequence variant.
    '''
    # make sure that the input sequence is all caps
    sequence = sequence.upper()

    # check length
    _length_check(sequence)

    if cutoff > 1 or cutoff < 0:
        raise goose_exceptions.GooseInputError('cutoff value must be between 0 and 1 for disorder threshold')    
    if len(sequence) < 6:
        raise goose_exceptions.GooseInputError('Cannot have sequence with a length less than 6')

    if NCPR > 1 or NCPR < -1:
        raise goose_exceptions.GooseInputError('NCPR values must be between -1 and 1.')

    try:
        final_sequence = _gen_ncpr_class_variant(sequence, ncpr=NCPR, attempts=attempts, disorder_threshold=cutoff, strict_disorder=strict, constant_fcr=False)
    except:
        raise goose_exceptions.GooseFail('Sorry! GOOSE was unable to generate the sequence. Please try again or try with a different ncpr value or a different cutoff value.')
    return final_sequence   


def all_props_class_var(sequence, hydropathy=None, FCR=None, NCPR=None, kappa=None,
    attempts=10, cutoff=parameters.DISORDER_THRESHOLD, strict=False):
    '''
    User facing funcitonality to generate variants where
    the ncpr, fcr, and/or hydropathy is/are changed and 
    the changes to classes of amino acids is minimized.
    
    Parameters
    ------------
    sequence : str
        The sequence to make a variant of.

    hydropathy : float
        The desired hydropathy value. (Optional)

    FCR : float
        The desired FCR value. (Optional)

    NCPR : float
        The desired NCPR value. (Optional)

    kappa : float
        The desired kappa value. (Optional)

    attempts : int
        Specify the number of times to make the sequence. Default is 10. Greater numbers
        of attempts increase the odds that a sequence will be generated but will increase
        the duration of attempting to make the sequence. (Optional)

    cutoff : float
        The disorder cutoff threshold. Closer to 1 has a higher chance of being disordered.

    strict : bool
        Whether to use a strict disorder calculation. By default, variants are allowed
        to have regions below the disorder threshold *for regions where the input sequence
        is also below the threshold*. (Optional)

    Returns
    -----------
    final_sequence : str
        The final sequence variant.
    '''
    # make sure that the input sequence is all caps
    sequence = sequence.upper()

    # check length
    _length_check(sequence)

    if cutoff > 1 or cutoff < 0:
        raise goose_exceptions.GooseInputError('cutoff value must be between 0 and 1 for disorder threshold')    
    if len(sequence) < 6:
        raise goose_exceptions.GooseInputError('Cannot have sequence with a length less than 6')

    if NCPR != None:
        if NCPR > 1 or NCPR < -1:
            raise goose_exceptions.GooseInputError('NCPR values must be between -1 and 1.')

    if hydropathy != None:
        if hydropathy < 0 or hydropathy > 6.4:
            raise goose_exceptions.GooseInputError('hydropathy values must be between 0 and 6.4')

    if FCR != None:
        if FCR > 1 or FCR < 0:
            raise goose_exceptions.GooseInputError('FCR values must be between 0 and 1.')

    if kappa != None:
        if kappa > 1 or kappa < 0:
            raise goose_exceptions.GooseInputError('kappa values must be between 0 and 1')


    try:
        final_sequence = _gen_all_props_class_variant(sequence, hydropathy=hydropathy, fcr=FCR,
        ncpr=NCPR, kappa=kappa, attempts=attempts, disorder_threshold=cutoff, strict_disorder=strict)
    except:
        raise goose_exceptions.GooseFail('Sorry! GOOSE was unable to generate the sequence. Please try again or try with a different input values or a different cutoff value.')
    return final_sequence   


def targeted_shuffle_var(sequence, target_aas, attempts=10,
    cutoff=parameters.DISORDER_THRESHOLD, strict=False):
    '''
    User facing funcitonality to generate variants where
    you can specify residues or classes of residues to shuffle.

    parameters
    ----------
    sequence : str
        The amino acid sequence as a string.

    target_aas : str or list
        A list of amino acids to shuffle
        or a class of amino acids to shuffle
        Possible target classes:
        charged : DEKR
        polar : QNST
        aromatic : FYW
        aliphatic : IVLAM
        negative: DE
        positive : KR

        If a list is specified, format should be (for example): 
        target_aas['K', 'W', 'R', 'Y']


    attempts : int
        Specify the number of times to make the sequence. Default is 10. Greater numbers
        of attempts increase the odds that a sequence will be generated but will increase
        the duration of attempting to make the sequence. (Optional)

    cutoff : float
        the cutoff value for disorder between 0 and 1. 
        Higher values have a higher likelihood of being disordered. (Optional)

    strict : bool
        Whether to use a strict disorder calculation. By default, variants are allowed
        to have regions below the disorder threshold *for regions where the input sequence
        is also below the threshold*. (Optional)
    
    Returns
    -------
    Returns the amino acid sequence as a string. 
    '''

    # make sure that the input sequence is all caps
    sequence = sequence.upper()

    # check length
    _length_check(sequence)

    if cutoff > 1 or cutoff < 0:
        raise goose_exceptions.GooseInputError('cutoff value must be between 0 and 1 for disorder threshold')    
    if len(sequence) < 6:
        raise goose_exceptions.GooseInputError('Cannot have sequence with a length less than 6')

    # make sequence.
    try:
        final_sequence = _gen_targeted_shuffle_variant(sequence, target_aas, attempts=attempts, disorder_threshold=cutoff, strict_disorder=strict)
    except:
        raise goose_exceptions.GooseFail('Sorry! GOOSE was unable to generate the sequence. Please try again or try with a different input values or a different cutoff value.')
    return final_sequence

def weighted_shuffle_var(sequence, target_aas, shuffle_weight=1.0, attempts=10,
    cutoff=parameters.DISORDER_THRESHOLD, strict=False):
    '''
    User facing funcitonality to generate variants where
    you can specify residues or classes of residues to shuffle 
    along with a weight to dictate how severely to shuffle
    the sequence.
    parameters
    ----------
    sequence : str
        The amino acid sequence as a string.
    target_aas : str or list
        A list of amino acids to shuffle
        or a class of amino acids to shuffle
        Possible target classes:
        charged : DEKR
        polar : QNST
        aromatic : FYW
        aliphatic : IVLAM
        negative: DE
        positive : KR
        If a list is specified, format should be (for example): 
        target_aas['K', 'W', 'R', 'Y']
        
    shuffle_weight : float
        a weight between 0.0-1.0 representing the probability of 
        moving a residue during shuffling
    attempts : int
        Specify the number of times to make the sequence. Default is 10. Greater numbers
        of attempts increase the odds that a sequence will be generated but will increase
        the duration of attempting to make the sequence. (Optional)
    cutoff : float
        the cutoff value for disorder between 0 and 1. 
        Higher values have a higher likelihood of being disordered. (Optional)
    strict : bool
        Whether to use a strict disorder calculation. By default, variants are allowed
        to have regions below the disorder threshold *for regions where the input sequence
        is also below the threshold*. (Optional)
    
    Returns
    -------
    Returns the amino acid sequence as a string. 
    '''

    # make sure that the input sequence is all caps
    sequence = sequence.upper()

    # check length
    _length_check(sequence)

    if cutoff > 1 or cutoff < 0:
        raise goose_exceptions.GooseInputError('cutoff value must be between 0 and 1 for disorder threshold')    
    if len(sequence) < 6:
        raise goose_exceptions.GooseInputError('Cannot have sequence with a length less than 6')

    # make sequence.
    try:
        final_sequence = _gen_weighted_shuffle_variant(sequence, target_aas, shuffle_weight, attempts=attempts, disorder_threshold=cutoff, strict_disorder=strict)
    except:
        raise goose_exceptions.GooseFail('Sorry! GOOSE was unable to generate the sequence. Please try again or try with a different input values or a different cutoff value.')
    return final_sequence


def targeted_reposition_var(sequence, target_aas, attempts=10,
    cutoff=parameters.DISORDER_THRESHOLD, strict=False):
    '''
    User facing funcitonality to generate variants where
    you can specify residues or classes of residues to reposition.

    parameters
    ----------
    sequence : str
        The amino acid sequence as a string.

    target_aas : str or list
        A list of amino acids to reposition
        or a class of amino acids to reposition
        Possible target classes:
        charged : DEKR
        polar : QNST
        aromatic : FYW
        aliphatic : IVLAM
        negative: DE
        positive : KR

        If a list is specified, format should be (for example): 
        target_aas['K', 'W', 'R', 'Y']


    attempts : int
        Specify the number of times to make the sequence. Default is 10. Greater numbers
        of attempts increase the odds that a sequence will be generated but will increase
        the duration of attempting to make the sequence. (Optional)

    cutoff : float
        the cutoff value for disorder between 0 and 1. 
        Higher values have a higher likelihood of being disordered. (Optional)

    strict : bool
        Whether to use a strict disorder calculation. By default, variants are allowed
        to have regions below the disorder threshold *for regions where the input sequence
        is also below the threshold*. (Optional)
    
    Returns
    -------
    Returns the amino acid sequence as a string. 
    '''

    # make sure that the input sequence is all caps
    sequence = sequence.upper()

    # check length
    _length_check(sequence)

    if cutoff > 1 or cutoff < 0:
        raise goose_exceptions.GooseInputError('cutoff value must be between 0 and 1 for disorder threshold')    
    if len(sequence) < 6:
        raise goose_exceptions.GooseInputError('Cannot have sequence with a length less than 6')

    # make sequence.
    try:
        final_sequence = _gen_targeted_reposition_variant(sequence, target_aas, attempts=attempts, disorder_threshold=cutoff, strict_disorder=strict)
    except:
        raise goose_exceptions.GooseFail('Sorry! GOOSE was unable to generate the sequence. Please try again or try with a different input values or a different cutoff value.')
    return final_sequence

def excluded_shuffle_var(sequence, exclude_aas, attempts=10,
    cutoff=parameters.DISORDER_THRESHOLD, strict=False):
    '''
    User facing funcitonality to generate variants where
    you can specify residues or classes of residues to NOT shuffle.

    parameters
    ----------
    sequence : str
        The amino acid sequence as a string.

    exclude_aas : str or list
        A list of amino acids to not shuffle
        or a class of amino acids to NOT shuffle
        Possible target classes:
        charged : DEKR
        polar : QNST
        aromatic : FYW
        aliphatic : IVLAM
        negative: DE
        positive : KR

        If a list is specified, format should be (for example): 
        target_aas['K', 'W', 'R', 'Y']

    attempts : int
        Specify the number of times to make the sequence. Default is 10. Greater numbers
        of attempts increase the odds that a sequence will be generated but will increase
        the duration of attempting to make the sequence.  

    cutoff : float
        The cutoff value for disorder between 0 and 1. 
        Higher values have a higher likelihood of being disordered. 

    strict : bool
        Whether to use a strict disorder calculation. By default, variants are allowed
        to have regions below the disorder threshold *for regions where the input sequence
        is also below the threshold*. (Optional)
    
    Returns
    -------
    Returns the amino acid sequence as a string. 

    '''
    # make sure that the input sequence is all caps
    sequence = sequence.upper()

    # check length
    _length_check(sequence)

    if cutoff > 1 or cutoff < 0:
        raise goose_exceptions.GooseInputError('cutoff value must be between 0 and 1 for disorder threshold')    
    if len(sequence) < 6:
        raise goose_exceptions.GooseInputError('Cannot have sequence with a length less than 6')

    # make sequence.
    try:
        final_sequence = _gen_excluded_shuffle_variant(sequence, exclude_aas, attempts=attempts, disorder_threshold=cutoff, strict_disorder=strict)
    except:
        raise goose_exceptions.GooseFail('Sorry! GOOSE was unable to generate the sequence. Please try again or try with a different input values or a different cutoff value.')
    return final_sequence




def re_var(sequence, increase_or_decrease, return_all=False, return_all_interval=0.2,
    include_original=False, attempts=None, cutoff=parameters.DISORDER_THRESHOLD, 
    strict=False):
    '''
    User facing funcitonality to generate variants with increased or decreased Re. 

    parameters
    ----------
    sequence : str
        The amino acid sequence as a string.

    increase_or_decrease : str
        Whether to increase or decrease the Re value.   

    return_all : bool
        Whether to return all sequences that have an increased or decreased Re. 
        Default=False

    return_all_interval : float
        minimum difference between Re values for each sequence returned if
        return_all=True. 
        Default=0.2

    include_original : bool
        Whether to return the orignal sequence along with the variants.
        Default=False
        If set to True, a nested dictionary with the keys 'original' and 'variant'
        are returned with a corresponding dictionary as the value where the key value pairs
        in the nested dictionary are the sequence and the rg for that sequence.

    number_attempts : int
        The number of attempts to mkae the sequence. 
        Default is 50*sequence length. 

    cutoff : float
        the cutoff value for somethign to be considered disordered. This only
        really matters if you set 'strict' to 'False'. Otherwise, the disorder
        of your input sequence is what is important. 

    strict : bool
        Whether to use a strict disorder calculation. By default, variants are allowed
        to have regions below the disorder threshold *for regions where the input sequence
        is also below the threshold*. (Optional)
    
    Returns
    -------
    Dict or nested dict
        If 'include' is set to True, returns a nested dict with 'original' and 'variants' as the
        keys with dictionaries as values wher the dictionaries contain the sequence as the key and 
        the predicted re as the value. If 'include_original' is set to fault (Default), then the 
        returned dictionary just has sequences as keys and thier corresponding predicted Rg as values. 
    '''
    # make sure that the input sequence is all caps
    sequence = sequence.upper()

    # check length
    _length_check(sequence)

    if cutoff > 1 or cutoff < 0:
        raise goose_exceptions.GooseInputError('cutoff value must be between 0 and 1 for disorder threshold')    

    # make sequence.
    try:
        final_sequence = _gen_dimensions_variant(sequence, increase_or_decrease, 're', 
            return_all=return_all, return_all_interval=return_all_interval,
            disorder_threshold=cutoff, strict_disorder=strict,
            include_original=include_original,
            num_attempts=attempts)
    except:
        raise goose_exceptions.GooseFail('Sorry! GOOSE was unable to generate the sequence. Please try again or try with a different input values or a different cutoff value.')
    return final_sequence




def rg_var(sequence, increase_or_decrease, return_all=False, return_all_interval=0.2,
    include_original=False, attempts=None, cutoff=parameters.DISORDER_THRESHOLD, 
    strict=False):
    '''
    User facing funcitonality to generate variants with increased or decreased Rg. 

    parameters
    ----------
    sequence : str
        The amino acid sequence as a string.

    increase_or_decrease : str
        Whether to increase or decrease the Re value.   

    return_all : bool
        Whether to return all sequences that have an increased or decreased Re. 
        Default=False

    return_all_interval : float
        minimum difference between Re values for each sequence returned if
        return_all=True. 
        Default=0.2

    include_original : bool
        Whether to return the orignal sequence along with the variants.
        Default=False
        If set to True, a nested dictionary with the keys 'original' and 'variant'
        are returned with a corresponding dictionary as the value where the key value pairs
        in the nested dictionary are the sequence and the rg for that sequence.

    number_attempts : int
        The number of attempts to mkae the sequence. 
        Default is 50*sequence length. 

    cutoff : float
        the cutoff value for somethign to be considered disordered. This only
        really matters if you set 'strict' to 'False'. Otherwise, the disorder
        of your input sequence is what is important. 

    strict : bool
        Whether to use a strict disorder calculation. By default, variants are allowed
        to have regions below the disorder threshold *for regions where the input sequence
        is also below the threshold*. (Optional)
    
    Returns
    -------
    Dict or nested dict
        If 'include' is set to True, returns a nested dict with 'original' and 'variants' as the
        keys with dictionaries as values wher the dictionaries contain the sequence as the key and 
        the predicted re as the value. If 'include_original' is set to fault (Default), then the 
        returned dictionary just has sequences as keys and thier corresponding predicted Rg as values. 
    '''
    # make sure that the input sequence is all caps
    sequence = sequence.upper()

    # check length
    _length_check(sequence)

    if cutoff > 1 or cutoff < 0:
        raise goose_exceptions.GooseInputError('cutoff value must be between 0 and 1 for disorder threshold')    

    # make sequence.
    try:
        final_sequence = _gen_dimensions_variant(sequence, increase_or_decrease, 'rg', 
            return_all=return_all, return_all_interval=return_all_interval,
            disorder_threshold=cutoff, strict_disorder=strict,
            include_original=include_original,
            num_attempts=attempts)
    except:
        raise goose_exceptions.GooseFail('Sorry! GOOSE was unable to generate the sequence. Please try again or try with a different input values or a different cutoff value.')
    return final_sequence



'''
/-/-/-/-/-/-/-/-/-/-/-/-/-/-/-/-/-/-/
/-/-/-/-/-/-/-/-/-/-/-/-/-/-/-/-/-/-/
 PREDICTED FOLDED REGION GENERATORS
/-/-/-/-/-/-/-/-/-/-/-/-/-/-/-/-/-/-/
/-/-/-/-/-/-/-/-/-/-/-/-/-/-/-/-/-/-/
'''

def alpha_helix(length, attempts=500):
    '''
    User facing function for generating a sequence predicted to be an alpha
    helix based on DSSP scores. 

    Parameters
    ----------
    length : int
        The length of the desired alpha helix.
        Must be between 8 and 150.

    attempts : int
        Specify the number of times to make the sequence. Default is 500. Greater numbers
        of attempts increase the odds that a sequence will be generated but will increase
        the duration of attempting to make the sequence. (Optional) 

    Returns
    -------
    final_seq : str
        The final sequence generated.

    '''
    if length > 150:
        raise goose_exceptions.GooseInputError('Unable to make alpha helix with length greater than 150.')    
    elif length < 8:
        raise goose_exceptions.GooseInputError('Unable to make alpha helix with length less than 8.')    
    else:
        try:
            final_seq = _gen_helix(length, max_iters=attempts)
        except:
            raise goose_exceptions.GooseFail('Sorry! Goose was unable to make that helix. Try again or try a different length.')
        return final_seq


def beta_strand(length, attempts=5000):
    '''
    User facing function for generating a sequence predicted to be 
    a beta strand based on DSSP scores.

    Parameters
    ----------
    length : int
        The length of the desired beta strand.
        Must be between 5 and 34.
    
    attempts : int
        Specify the number of times to make the sequence. Default is 5000. Greater numbers
        of attempts increase the odds that a sequence will be generated but will increase
        the duration of attempting to make the sequence.  
    
    returns
    -------
    final_seq : str  
        Final generated sequence.  
    '''
    if length > 34:
        raise goose_exceptions.GooseInputError('Unable to make beta strands with length greater than 34.')
    elif length < 5:
        raise goose_exceptions.GooseInputError('Unable to make beta strands with length less than 5.')
    else:
        try:
            final_seq = _gen_beta_strand(length, max_iters=attempts)
        except:
            raise goose_exceptions.GooseFail('Sorry! Goose was unable to make that strand. Try again or try a different length.')
        return final_seq

def beta_sheet(length):
    '''
    User facing function for generating a sequence predicted to be 
    a beta sheet based on DSSP scores. Uses coils to connect strands.

    Parameters
    ----------
    length : int
        The length of the desired beta_sheet.
        Must be between 18 and 400 residues.

    returns
    -------
    final_seq : str  
        Final generated sequence.  
    '''

    if length < 18:
        raise goose_exceptions.GooseInputError('cannot generate beta sheet less than 18 amino acids.')
    elif length > 400:
        raise goose_exceptions.GooseInputError('cannot generate beta sheet greater than 400 amino acids.')
    try:
        final_seq = _gen_beta_sheet(length)
    except:
        raise goose_exceptions.GooseFail('Sorry! Goose was unable to make that beta sheet. Try again or try a different length.')
    return final_seq



'''
/-/-/-/-/-/-/-/-/-/-/-/-/-/-/-/-/-/-/
/-/-/-/-/-/-/-/-/-/-/-/-/-/-/-/-/-/-/
        LIBRARY GENERATION
/-/-/-/-/-/-/-/-/-/-/-/-/-/-/-/-/-/-/
/-/-/-/-/-/-/-/-/-/-/-/-/-/-/-/-/-/-/
'''


def seq_property_library(length,
    FCR=None, NCPR=None, hydropathy=None, kappa=None, 
    cutoff=parameters.DISORDER_THRESHOLD, silent_failed_seqs=False, beta_mode=False, random_name=False):
    '''
    Function that returns a list of dictioanries where each
    dictionary has the values specified for a library of sequences.
    Bit of a long function because all of the paramter values are checked
    before adding the sequence to the library. This way the user can specify 
    ranges of, for example, hydropathy, that are not possible with certain charge
    values but are possible with other charge values. 


    parameters
    ----------
    length : Int 
        The length of the sequence

    FCR : list of float(s) or float
        The fraction of charged residues as a decimal, between 0 and 1.
        Input as a range where the first value is the lower bounds wanted
        for the library and the second value is the upper bounds wanted for
        the library. If a single number is given, then that value will be constant
        across the library.
        A third value can be specified, which will determine the
        interval between the lowest value and the highest value. 
        If this interval is not perfect, then GOOSE will take the 
        closest value to the desired maximum. (Optional)

    NCPR : list of float(s) or float
        The net charge of the sequence as a decimal, between -1 and 1
        Absolute value cannot be greater than FCR.
        Input as a range where the first value is the lower bounds wanted
        for the library and the second value is the upper bounds wanted for
        the library. If a single number is given, then that value will be constant
        across the library.
        A third value can be specified, which will determine the
        interval between the lowest value and the highest value. 
        If this interval is not perfect, then GOOSE will take the 
        closest value to the desired maximum. (Optional)

    hydropathy : list of float(s) or float
        The mean hydropathy of the sequence. Lower numbers are less
        hydrophobic. Between 0.6 and 6.1
        Input as a range where the first value is the lower bounds wanted
        for the library and the second value is the upper bounds wanted for
        the library. If a single number is given, then that value will be constant
        across the library.
        A third value can be specified, which will determine the
        interval between the lowest value and the highest value. 
        If this interval is not perfect, then GOOSE will take the 
        closest value to the desired maximum. (Optional)     

    kappa : list of float(s) or float
        The charge asymmetry metric. Describes how oppositely charged
        residues are patterned across the sequeence. Between 0 and 1 
        where a higher value is a more asymmetrically distributed 
        positioning of oppositely charged residues.
        Input as a range where the first value is the lower bounds wanted
        for the library and the second value is the upper bounds wanted for
        the library. If a single number is given, then that value will be constant
        across the library.        
        A third value can be specified, which will determine the
        interval between the lowest value and the highest value. 
        If this interval is not perfect, then GOOSE will take the 
        closest value to the desired maximum. (Optional)

        Input as a range where the first value is the lower bounds wanted
        for the library and the second value is the upper bounds wanted for
        the library. If a single number is given, then that value will be constant
        across the library.

    cutoff : float
        The cutoff value for disorder as a float. Higher values
        lead to a more 'strict' cutoff for what is considered to be disordered. (Optional)

    silent_failed_seqs : bool
        Whether to silence any printed warnings of sequences that
        are not possible to generate due to incompatible charge/hydorpatyh values. (Optional)

    beta_mode : bool
        For testing. Set to True to get some printouts as seqs are being generated. (Optional)

    random_name : bool
        Whether to use a randomly generated string as the protein name. Helps avoid
        overwriting sequences. (Optional)

    returns
    -------
    sequence_list : list of dictionaries
        Returns a list of dicts where each dictionary in the list 
        can be input into the sequence generators downstream. 
        The dict holds the parameter specifications for each sequence in the library.
    '''

    # make the list of sequences to generate
    sequence_list = _generate_library_by_parameter_ranges(length=length,
    FCR=FCR, NCPR=NCPR, hydropathy=hydropathy, kappa=kappa, silent_failed_seqs=silent_failed_seqs)
    # dict to return
    seq_dict={}
    # make a dict to hold the library
    for seq_specified in sequence_list:
        if beta_mode==True:
            print(f'{seq_specified}')
        # make the sequence
        seq = sequence(length, **seq_specified)
        # make a sequence name
        seq_name='>'
        name_values = ['FCR', 'NCPR', 'kappa', 'hydropathy']
        for vals in range(0,len(name_values)):
            param_val = name_values[vals]
            curval = seq_specified[param_val]
            if curval != None:
                if seq_name=='>':
                    seq_name+=f'{param_val}_{round(curval, 4)}'
                else:
                    seq_name+=f'_{param_val}_{round(curval, 4)}'
        # if user wants a randodm name, make it happen.
        if random_name==True:
            seq_name+=f'_{_gen_random_name()[1:]}'

        # add to the dict
        seq_dict[seq_name]=seq
    # return the dict
    return seq_dict

def seq_fractions_library(length, random_name=False, warn_user=True, robust_warning=False, **kwargs):
    '''
    User facing function for generating a library of sequences that have varying specified fractions
    of aminon acids. This function will only make fractions that are possible (ie. less than 1.0 specified
    and all values less than the respective value for each seq as far as the maximum possible). By default will
    alert the user to any 'failed' generate sequences...

    Parameters
    ------------
    length : int
        The length of the sequence to generate.

    random_name : bool
        Whether to make a random name for each sequence. Helps avoid overwriting
        sequences. (Optional)

    warn_user : bool
        whether to warn the user if the fractions of a sequence do not match 
        what the user input. This is typically due to length problems (bascially,
        you can specify F=0.95 but will never get 0.95 for F if your sequence is 10
        amino acids long because you can't have 1/2 a F). (Optional)

    robust_warning : bool
        Whether to print out a robust message for each incorrect sequence highlighting what was wrong. (Optional)

    <each of the 20 amino acids> : float
        Specify the fraction of the sequence that should be made up of one or more
        of the 20 natural amino acids (e.g. A=0.2, Y=0.05) etc.

    max_aa_fractions : dict 
        Dictionary which, if provided, allows the user to over-ride the 
        fraction of a sequence which can be made up of any given amino
        acid. The passed dictionary should contain key/value pairs, where
        keys are one of the twenty standard amino acids and values is a
        float between 0 and 1. If amino acids are missing then the default
        thresholds set by GOOSE are used. (Optional.)

    Returns 
    --------
    seq_dict : dict
        returns a dictionary of the sequences.   
    '''

    # make the list of sequences to generate
    sequence_list = _generate_library_by_fraction_ranges(length=length, **kwargs)
    # need dict for naming...
    possible_vals_name=['A', 'C', 'D', 'E', 'F', 'G', 'H', 'I', 'K', 'L', 'M', 'N', 'P', 'Q', 'R', 'S', 'T', 'V', 'W', 'Y']
    name_vals = []
    for specified_arg in kwargs.keys():
        if specified_arg in possible_vals_name:
            name_vals.append(specified_arg)
    # keep track of sequences with erros
    sequence_warnings=0
    # dict to return
    seq_dict={}
    # list of additional kwargs
    additional_kwargs=['cutoff', 'attempts',  'strict_disorder',  'max_aa_fractions']
    for seq_specified in sequence_list[1]:
        # add in any additional kwargs to input into seq_fractions function
        for cur_kwarg in additional_kwargs:
            if cur_kwarg in kwargs.keys():
                seq_specified[cur_kwarg]=kwargs[cur_kwarg]
        sequence = seq_fractions(length, **seq_specified)
        # make a sequence name
        seq_name='>'
        sequence_correct=True
        error_warning=f'For:\n{sequence}\n'
        if name_vals != []:
            for vals in name_vals:
                actual_val = sequence.count(vals)/len(sequence)
                param_val = seq_specified[vals]
                # see if actual val equals param val
                if round(float(actual_val), 8) != round(float(param_val), 8):
                    sequence_correct=False
                    error_warning += f'Objective fraction for {vals}: {round(float(param_val), 8)}, Actual value: {round(float(actual_val), 8)}\n'
                if random_name==False:                
                    if seq_name=='>':
                        seq_name+=f'{vals}_{round(actual_val, 4)}'
                    else:
                        seq_name+=f'_{vals}_{round(actual_val, 4)}'
                # if user wants a randodm name, make it happen.
                else:
                    seq_name+=f'_{_gen_random_name()[1:]}'
            if sequence_correct==False:
                sequence_warnings+=1
                if robust_warning==True:
                    print(error_warning)
                    print()

        # make sure seq_name not already in dict
        if seq_name not in seq_dict.keys():
            # add to seq dict
            seq_dict[seq_name]=sequence
        else:
            # this shouldn't ever be a problem, but just to be safe add in a random 4 digit ID.
            add_id = '_ID_'
            for i in range(0, 8):
                add_id += f'{random.randint(0, 9)}'
            seq_name+=add_id
            seq_dict[seq_name]=sequence
    # warn user if needed and if wanted by user
    if warn_user == True:
        if sequence_warnings > 0:
            print(f'\nWARNING!\n GOOSE detected {sequence_warnings} sequences out of {len(seq_dict)} sequences where the generated sequence does not have the exact value as the input.\n This is often due to a mismatch in the specified length and possible fractions.\n')
    return seq_dict




