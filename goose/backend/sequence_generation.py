'''
functionality for generating sequences where
disorder is actually checked. The sequence_generation_backend.py
holds code simply for generating sequences that are used here. The
backend does not check disorder because it allows for more rapid generation
of input sequences that can then be checked here.
'''

import metapredict as meta

from goose.backend.sequence_generation_backend import create_seq_by_props, sigma_FCR_NCPR, create_seq_by_fracs
from goose.backend.seq_by_dimension_backend import build_seq_by_dimensions
from goose.backend import parameters
from goose.goose_exceptions import GooseFail

# can probably delete this line
#from goose.backend.protein import Protein


def check_disorder(sequence, disorder_threshold=parameters.DISORDER_THRESHOLD, strict=False):
    '''
    function to check whether a generated sequence is disordered.
    The function allows some 'ordered' residues provided that 
    'strict' is not set to True.

    Parameters
    ----------
    sequence : String
        The amino acid sequence of the protein to check as a string

    disorder_threshold : Float
        The disorder threshold value. Higher values are more 'strict'

    strict : Bool
        Whether to allow for any residues below the disorder_threshold.
        If set to True, then if any residue goes below threshold, the sequence
        will not be considered disordered. Default is False.

    Returns
    -------
    Bool
        Returns True if sequence is disordered and False if not disordered

    
    '''
    # calculate the sequence disorder using metapredict
    if len(sequence) > 20:
        sequence_disorder = meta.predict_disorder(sequence[2:len(sequence)-2])
    else:
        sequence_disorder=meta.predict_disorder(sequence)

    # see what the min disorder is. If min is over thresh, good to go. 
    min_disorder = min(sequence_disorder)
    if min_disorder >= disorder_threshold:
        return True
    
    # Otherwise, check if we are strict with the check disorder functionality. 
    else:
        # if strict set to True, nothing can be below threshold. Return False
        if strict==True:
            return False
        # if strict not set to True, can have some residues below thresh value. 
        else:
            # allow a little bit of 'order'
            cur_order_streak = 0
            # keep track of total number of ordered residues
            total_ordered_residues = 0
            # allow up to 5 % of residues to be 'ordered' provided they aren't consecutive
            allowed_order_residues = round(0.05*len(sequence))
            if allowed_order_residues < 1:
                allowed_order_residues = 1
            
            # calculate number of conseuctively allowed residues
            # going to allow up to the length of the sequence over 25 but 
            # not greater than 10
            consec_allowed = round(len(sequence)/25)
            if consec_allowed < 1:
                consec_allowed = 1
            if consec_allowed > 10:
                consec_allowed = 10

            # now going to go through the sequence
            for residue in sequence_disorder:
                if residue < disorder_threshold:
                    cur_order_streak += 1
                    total_ordered_residues += 1
                else:
                    cur_order_streak = 0
                # if our current order streak greater than the consecutive allowed or our
                # total number of residues above concescutive allowed, return False
                if cur_order_streak > consec_allowed or total_ordered_residues > allowed_order_residues:
                    return False
    # finally, return True. Only get here if we have consecutive ordered residues below
    # the above calculated max and total ordered below above calculated max. 
    return True


def generate_disordered_seq_by_props(length, FCR=None, NCPR=None, hydropathy=None, sigma=None, attempts=50, 
    allowed_hydro_error = parameters.HYDRO_ERROR, disorder_threshold = parameters.DISORDER_THRESHOLD, 
    strict_disorder=False, exclude = []):
    '''
    Function to actually generate a disordered sequence.
    General idea is to first generate the sequecne and see
    how the disorder ends up. If it's good, just return it. 
    If not, find a region of the sequence that is 
    ordered and replace it if possible. This is dependent
    on the sequence size.

    Parameters
    ----------
    length : Int
        The length of the wanted protein sequence as an integer value

    FCR : Float
        The fraction of charged residues wanted for the sequence as a 
        decimal value.

    NCPR : Float
        The wanted net charge of the sequence given as a decimal value

    hydropathy : Float 
        The wanted mean hydropathy value of the sequence. Uses ajusted
        Kyte-doolittle hydropathy scale that goes from 0 to 9 where
        0 is least hydrophobic and 9 is most hydrophobic.

    attempts : Int
        The number of times to attempt to build a sequence before throwing
        in the towel

    allowed_hydro_error : Float
        The allowed error for hydropathy between the value of hydropathy and
        the final hydropathy value of the generated sequence

    disorder_threshold : Float
        The value for a residue to be considered disordered.

    strict_disorder : Bool
        Whether to have a strict cutoff for disorder where if any single residue
        falls below disorder_threshold, the sequence is not considered disordered.
        Set to False by default allowing single residues to occassionally drop below
        the disorder theshold provided it is minimal. See check_disorder for more
        details.

    exclude : list
        A list of residues to exclude from sequence generation. 

    Returns
    -------
    final_seq : String
        Returns the final sequence that was specified as a string.
    '''
    
    # try the number of specified attempts to build the seq
    for attempt_num in range(0, attempts):
        attempted_seq=None
        # if sigma is specified, get the corresponding
        # FCR and NCPR values
        if sigma != None:
            FCR_NCPR_Dict = sigma_FCR_NCPR(length, sigma)
            FCR = FCR_NCPR_Dict['FCR']
            NCPR = FCR_NCPR_Dict['NCPR']

        # try to build the sequence
        try:
            attempted_seq = create_seq_by_props(length, FCR=FCR, NCPR=NCPR, hydropathy=hydropathy, attempts=200, 
            allowed_hydro_error = parameters.HYDRO_ERROR, exclude=exclude)

        # if attempt to build the sequence failed, continue back at the 
        # beginning of the loop
        except:
            continue

        # if the sequence is disordered, return it
        if attempted_seq != None:
            if check_disorder(attempted_seq, disorder_threshold=disorder_threshold, strict=strict_disorder):
                return attempted_seq

    # if no disordered sequence in number of attempts, raise GooseFail
    raise GooseFail('Unable to generate sequence! Try increasing attempts!')




def generate_disordered_seq_by_fractions(length, **kwargs):
    '''
    Function to actually generate a disordered sequence.
    General idea is to first generate the sequecne and see
    how the disorder ends up. If it's good, just return it. 
    If not, find a region of the sequence that is 
    ordered and replace it if possible. This is dependent
    on the sequence size.

    Parameters
    ----------
    length : int
        The length of the wanted protein sequence as an integer value

    <each of the 20 amino acids> : float
        Specify the fraction of the sequence that should be made up of one or more
        of the 20 natural amino acids (e.g. A=0.2, Y=0.05) etc.

    attempts : int
        The number of times to attempt to build a sequence before throwing
        in the towel. 

    disorder_threshold : Float
        The value for a residue to be considered disordered.

    strict_disorder : Bool
        Whether to have a strict cutoff for disorder where if any single residue
        falls below disorder_threshold, the sequence is not considered disordered.
        Set to False by default allowing single residues to occassionally drop below
        the disorder theshold provided it is minimal. See check_disorder for more
        details.

    max_aa_fractions : dict 
        Dictionary which, if provided, allows the user to over-ride the 
        fraction of a sequence which can be made up of any given amino
        acid. The passed dictionary should contain key/value pairs, where
        keys are one of the twenty standard amino acids and values is a
        float between 0 and 1. If amino acids are missing then the default
        thresholds set by GOOSE are used.

    Returns
    -------
    final_seq : String
        Returns the final sequence that was specified as a string.


    '''

    ## NOTE FOR FUTURE CODE CLEANUP
    ## These should have already been defined in the create.seq_fractions
    ## functions, so suggest we remove these initializations as they're
    ## redundant. Leaving for now. (~ash 2023-04-23)
    
    # check for necessary kwargs
    if 'attempts' not in list(kwargs.keys()):
        attempts = parameters.DEFAULT_ATTEMPTS
    else:
        attempts = kwargs['attempts']

    if 'cutoff' not in list(kwargs.keys()):
        cutoff = parameters.DISORDER_THRESHOLD
    else:
        cutoff = kwargs['cutoff']

    if 'strict_disorder' not in list(kwargs.keys()):
        strict_disorder = False
    else:
        strict_disorder = kwargs['strict_disorder']    

    # make input kwargs just amino acids
    input_kwargs = {}
    amino_acids = ['A', 'C', 'D', 'E', 'F', 'G', 'H', 'I', 'K', 'L', 'M', 'N', 'P', 'Q', 'R', 'S', 'T', 'V', 'W', 'Y']
    for kw in kwargs.keys():
        curkwval = kwargs[kw]
        if kw in amino_acids:
            input_kwargs[kw] = curkwval


    # try the number of specified attempts to build the seq
    for attempt_num in range(0, attempts):
        # set attempted seq to empty string
        attempted_seq=''

        # try to build the sequence
        try:
            attempted_seq = create_seq_by_fracs(length, max_aa_fractions=kwargs['max_aa_fractions'], **input_kwargs)

        # if attempt to build the sequence failed, continue back at the 
        # beginning of the loop
        except:
            continue

        # if the sequence is disordered, return it
        if attempted_seq != '':
            if check_disorder(attempted_seq, disorder_threshold=cutoff, strict=strict_disorder):
                return attempted_seq

    # new fallback that should expand the sequence space. 
    # only use if we have already failed using the optimized approach
    for attempt_num in range(0, attempts):
        # set attempted seq to empty string
        attempted_seq=''
        # try to build the sequence *Without choosing optimized sequences**
        try:
            attempted_seq = create_seq_by_fracs(length, max_aa_fractions=kwargs['max_aa_fractions'],choose_optimized_residue=False, **input_kwargs)

        # if attempt to build the sequence failed, continue back at the 
        # beginning of the loop
        except:
            continue

        # if the sequence is disordered, return it
        if attempted_seq != '':
            if check_disorder(attempted_seq, disorder_threshold=cutoff, strict=strict_disorder):
                return attempted_seq

    # if no disordered sequence in number of attempts, raise GooseFail
    raise GooseFail('Unable to generate sequence! Try increasing attempts!')


def generate_disordered_seq_by_dimensions(seq_length, rg_or_re, objective_dims, attempts=20, 
    allowed_error = 'default_error', disorder_threshold = parameters.DISORDER_THRESHOLD, 
    strict_disorder=False, individual_rg_re_attempts=parameters.rg_re_attempt_num):
    '''
    Parameters
    ----------
    seq_length: int
        Length of sequence to generate
    rg_or_re: str
        'rg' or 're' depending on whether you want to specify radius of gyration or end to end distance
    objective_dim : float
        objective rg or re value
    allowed_error: float
        Allowed error between the specified radius of gyration and the actual radius of gyration
        default is from the backend.parameters module, re_error or rg_error
    attempts : Int
        The number of times to attempt to build a sequence before throwing
        in the towel
    disorder_threshold : Float
        The value for a residue to be considered disordered.
    strict_disorder : Bool
        Whether to have a strict cutoff for disorder where if any single residue
        falls below disorder_threshold, the sequence is not considered disordered.
        Set to False by default allowing single residues to occassionally drop below
        the disorder theshold provided it is minimal. See check_disorder for more
        details.
    individual_rg_re_attempts : int
        Number of attempts to make the objective rg or re. 
        Does not account for disorder

    Returns
    -------
    attempted_seq : String
        Returns the final sequence that was specified as a string.
    '''
    # make sure rg or re is good. 
    if rg_or_re not in ['re', 'rg']:
        raise GooseError('rg_or_re must be "rg" or "re"')

    # try for number of attempts. 
    for attempt in range(0, attempts):
        # try to build the sequence
        try:
            attempted_seq = build_seq_by_dimensions(seq_length,
            rg_or_re = rg_or_re, objective_dim=objective_dims,
            allowed_error=allowed_error, num_attempts=individual_rg_re_attempts)

        # if attempt to build the sequence failed, continue back at the 
        # beginning of the loop
        except:
            continue

        # if the sequence is disordered, return it
        if check_disorder(attempted_seq, disorder_threshold=disorder_threshold, strict=strict_disorder):
            return attempted_seq

    # if no disordered sequence in number of attempts, raise GooseFail
    raise GooseFail('Unable to generate sequence! Try increasing attempts!')



