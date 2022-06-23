'''
functionality for generating sequences where
disorder is actually checked. The sequence_generation_backend.py
holds code simply for generating sequences that are used here. The
backend does not check disorder because it allows for more rapid generation
of input sequences that can then be checked here.
'''

import metapredict as meta

from goose.backend.sequence_generation_backend import create_seq_by_props, sigma_FCR_NCPR, create_seq_by_fracs
from goose.backend import parameters
from goose.goose_exceptions import GooseFail
from goose.backend.protein import Protein


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

    # see what the min disorder is
    min_disorder = min(sequence_disorder)
    if min_disorder >= disorder_threshold:
        return True

    # if the minimum disorder is very low, just return false.
    elif min_disorder < 0.25:
        return False
    # see if sequence is likely to be disorderd
    else:
        # if strict set to True, nothing can be below threshold. Return False
        if strict==True:
            return False
        # if strict not set to True, can have some 'order'
        else:
            # allow a little bit of 'order'
            cur_order_streak = 0
            # keep track of total number of disordered residues
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
                if cur_order_streak > consec_allowed or total_ordered_residues > allowed_order_residues:
                    return False
    return True





def generate_disordered_seq_by_props(length, FCR=None, NCPR=None, hydropathy=None, sigma=None, attempts=20, 
    allowed_hydro_error = parameters.HYDRO_ERROR, disorder_threshold = parameters.DISORDER_THRESHOLD, strict_disorder=False):
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

    Returns
    -------
    final_seq : String
        Returns the final sequence that was specified as a string.


    '''
    # try the number of specified attempts to build the seq
    for attempt_num in range(0, attempts):
        # if sigma is specified, get the corresponding
        # FCR and NCPR values
        if sigma != None:
            FCR_NCPR_Dict = sigma_FCR_NCPR(length, sigma)
            FCR = FCR_NCPR_Dict['FCR']
            NCPR = FCR_NCPR_Dict['NCPR']

        # try to build the sequence
        try:
            attempted_seq = create_seq_by_props(length, FCR=FCR, NCPR=NCPR, hydropathy=hydropathy, attempts=20, 
            allowed_hydro_error = parameters.HYDRO_ERROR)

        # if attempt to build the sequence failed, continue back at the 
        # beginning of the loop
        except:
            continue

        # if the sequence is disordered, return it
        if check_disorder(attempted_seq, disorder_threshold=disorder_threshold, strict=strict_disorder):
            return attempted_seq

    # if no disordered sequence in number of attempts, raise GooseFail
    raise GooseFail('Unable to generate sequence!')




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
    length : Int
        The length of the wanted protein sequence as an integer value

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

    **kwargs : Variable, float
        The desired amino acid as a variable (no need for quotations).     
        The fraction of amino acids as a float followed immediately by


    Returns
    -------
    final_seq : String
        Returns the final sequence that was specified as a string.


    '''

    # check for necessary kwargs
    if 'attempts' not in list(kwargs.keys()):
        attempts = 1
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

        # try to build the sequence
        try:
            attempted_seq = create_seq_by_fracs(length, **input_kwargs)

        # if attempt to build the sequence failed, continue back at the 
        # beginning of the loop
        except:
            continue

        # if the sequence is disordered, return it
        if check_disorder(attempted_seq, disorder_threshold=cutoff, strict=strict_disorder):
            return attempted_seq

    # if no disordered sequence in number of attempts, raise GooseFail
    raise GooseFail('Unable to generate sequence!')
