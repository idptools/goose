'''
Testing whether it's possible to generate sequences with 
specific radius of gyration or end to end distance...

Current bug with SPARROW - 
the predict_scaled_re appears to predict re... 
Will figure this out but doesn't impact functionality
for what we are trying to do here. 

'''

from matplotlib import pyplot as plt
import random
import numpy as np

from sparrow import Protein as pr
from sparrow.predictors import batch_predict

from goose.backend.sequence_generation_backend import gen_sequence
from goose.backend import parameters
from goose.backend.parameters import re_error, rg_error, rg_re_attempt_num, get_min_re, get_max_re
from goose.backend import lists
from goose.backend.lists import disordered_list, disordered_list_reduced_charge
from goose.goose_exceptions import GooseError, GooseException, GooseInputError, GooseFail
from goose.backend.variant_generation_backend import create_asymmetry_variant



def predict_rg(sequence): return pr(sequence).predictor.radius_of_gyration(use_scaled=True)
def predict_re(sequence): return pr(sequence).predictor.end_to_end_distance(use_scaled=True)
def batch_predict_rg(sequences,show_progress_bar=False, return_seq2prediction=False): 
    seqobs=[pr(seq) for seq in sequences]
    rgs=batch_predict.batch_predict(seqobs, 'scaled_rg',
        show_progress_bar=show_progress_bar, batch_size=32,
        return_seq2prediction=return_seq2prediction)
    return rgs

def batch_predict_re(sequences, show_progress_bar=False, return_seq2prediction=False):
    seqobs=[pr(seq) for seq in sequences]
    re=batch_predict.batch_predict(seqobs, 'scaled_re',
        show_progress_bar=show_progress_bar, batch_size=32,
        return_seq2prediction=return_seq2prediction)
    return re



def optimize_seq_dims(input_sequence, rg_or_re, objective_dim, allowed_error='default_error',
    num_attempts=rg_re_attempt_num, reduce_pos_charged=True, exclude_aas=None):
    '''
    Takes an input sequence and tries to optimize it such that it has the correct
    objective Re or Rg.

    Parameters
    -----------
    input_sequence: str
        the input sequence to optimzies
    rg_or_re: str
        'rg' or 're' depending on whether you want to specify radius of gyration or end to end distance
    objective_dim : float
        objective rg or re value
    allowed_error: float
        Allowed error between the specified radius of gyration and the actual radius of gyration
        default is from the backend.parameters module, re_error or rg_error
    num_attempts: int
        Number of attempts to try to generate a sequence with the specified radius of gyration
        default is from the backend.parameters module, rg_re_attempt_num
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
    --------
    str
        Sequence with specified radius of gyration or end to end distance
    '''
    # which list to use
    if reduce_pos_charged==True:
        dis_list_to_use=disordered_list_reduced_charge
        bias_dim_dict={'collapse':['W','Y', 'G', 'F', 'Q', 'N'], 'expand':['D', 'E', 'P', 'S', 'T']}
    else:
        dis_list_to_use=disordered_list
        bias_dim_dict={'collapse':['W','Y', 'G', 'F', 'Q', 'N'], 'expand':['D', 'E', 'K', 'P', 'S', 'T']}


    # see if we need to get rid of anything. This makes a dict that counts 
    # the number of everything then removes from the dict so we don't have to
    # check the list a ton of times. 
    if exclude_aas is not None:
        if isinstance(exclude_aas, list)==False:
            raise Exception('exclude_aas must be a list')
        tdict={}
        for a in set(dis_list_to_use):
            tdict[a]=dis_list_to_use.count(a)
        # get rid of stuff in exclude_aas
        for a in exclude_aas:
            if a in tdict:
                del tdict[a]
            if a in bias_dim_dict['collapse']:
                bias_dim_dict['collapse'].remove(a)
            if a in bias_dim_dict['expand']:
                bias_dim_dict['expand'].remove(a)
            if bias_dim_dict['collapse']==[]:
                raise Exception('Cannot specify W, Y, G, F, Q, and N in exclude_aas')
            if bias_dim_dict['expand']==[]:
                raise Exception('Cannot specify D, E, K, P, S, and T in exclude_aas')
        # make new list
        nl=[]
        for a in tdict:
            for i in range(0,tdict[a]):
                nl.append(a)
        # set dis_list_to_use equal to nl
        dis_list_to_use=nl

    # best seq
    start_seq=input_sequence
    best_seq=start_seq
    # get the starting rg or re.
    if rg_or_re=='rg':
        start_dim = predict_rg(start_seq)
        if allowed_error == 'default_error':
            allowed_error=rg_error
    elif rg_or_re=='re':
        start_dim = predict_re(start_seq)
        if allowed_error == 'default_error':
            allowed_error=re_error
    else:
        raise GooseInputError('rg_or_re must be "rg" or "re"')
    if abs(start_dim - objective_dim)<=allowed_error:
        return start_seq
    else:
        starter=input_sequence
        best_error = abs(start_dim - objective_dim)
        for attempt in range(0,num_attempts):
            if start_dim > objective_dim:
                use_aas = bias_dim_dict['collapse']
            else:
                use_aas = bias_dim_dict['expand']
            # make 32 seq variants. 
            test_seqs=[]
            # change one amino acid at a time. 
            for frac in range(0, 32):
                num_from_dim_change = int(len(input_sequence)/80)
                if num_from_dim_change <1:
                    num_from_dim_change=1
                temp=[]
                for aa in range(0, num_from_dim_change):
                    temp.append(random.choice(use_aas))
                starter=list(starter)
                random.shuffle(starter)
                temp.extend(starter[num_from_dim_change:])
                random.shuffle(temp)
                test_seqs.append(''.join(temp))
            if rg_or_re=='rg':
                seq_dims = batch_predict_rg(test_seqs)
            elif rg_or_re=='re':
                seq_dims = batch_predict_re(test_seqs)
            else:
                raise GooseError('Somehow first exception was not raised, rg_or_re must be "rg" or "re"')
            # get closest rg
            for val in seq_dims:
                cur_seq, cur_dims = seq_dims[val]
                cur_er = abs(cur_dims-objective_dim)
                if cur_er < allowed_error:
                        return cur_seq
                if cur_er < best_error:
                    best_error = cur_er
                    starter = list(cur_seq)
                    start_dim = cur_dims                    
    raise GooseFail(f'Could not generate sequence with specified {rg_or_re}. \nTry increasing number of attempts or with a different rg/re value or length value')



def build_seq_by_dimensions(seq_length, rg_or_re, objective_dim, allowed_error='default_error',
    num_attempts=rg_re_attempt_num, num_sub_attempts=5, reduce_pos_charged=True, exclude_aas=None):
    '''
    Builds a sequence of a given length with a given radius of gyration or Re
    by starting with a sequence that is weighted to be disordered and
    approximately the right dimensions and then subsequently optimizing 
    that sequence.
    
    Parameters
    -----------
    seq_length: int
        Length of sequence to generate
    rg_or_re: str
        'rg' or 're' depending on whether you want to specify radius of gyration or end to end distance
    objective_dim : float
        objective rg or re value
    allowed_error: float
        Allowed error between the specified radius of gyration and the actual radius of gyration
        default is from the backend.parameters module, re_error or rg_error
    num_attempts: int
        Number of attempts to try to generate a sequence with the specified radius of gyration
        default is from the backend.parameters module, rg_re_attempt_num
    num_sub_attempts: int
        Number of sub attempts
        Default is 5
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
    --------
    str
        Sequence with specified radius of gyration or end to end distance
    '''
    # figure out disordered list to use. 
    if rg_or_re.lower()=='re':
        min_possible_value=parameters.get_min_re(seq_length)
        max_possible_value=parameters.get_max_re(seq_length)
    else:
        min_possible_value=parameters.get_min_rg(seq_length)
        max_possible_value=parameters.get_max_rg(seq_length)        
    
    #dis_list_to_use
    dimensions_range = max_possible_value-min_possible_value
    objective_scaled = (objective_dim-min_possible_value)/dimensions_range
    if objective_scaled<=0.1:
        dis_list_to_use=lists.dims_under_10
    elif objective_scaled >0.1 and objective_scaled <=0.22:
        dis_list_to_use=lists.dims_10_22
    elif objective_scaled >0.22 and objective_scaled <=0.32:
        dis_list_to_use=lists.dims_22_32
    elif objective_scaled >0.32 and objective_scaled <=0.42:
        dis_list_to_use=lists.dims_32_42
    elif objective_scaled >0.42 and objective_scaled <=0.53:
        dis_list_to_use=lists.dims_42_53
    elif objective_scaled >0.53 and objective_scaled <=0.70:
        dis_list_to_use=lists.dims_53_70
    elif objective_scaled >0.7 and objective_scaled <=0.80:
        dis_list_to_use=lists.dims_70_80
    else:
        dis_list_to_use=lists.dims_over_80


    # see if we need to get rid of anything. This makes a dict that counts 
    # the number of everything then removes from the dict so we don't have to
    # check the list a ton of times. 
    if exclude_aas is not None:
        if isinstance(exclude_aas, list)==False:
            raise Exception('exclude_aas must be a list')
        tdict={}
        for a in set(dis_list_to_use):
            tdict[a]=dis_list_to_use.count(a)
        # get rid of stuff in exclude_aas
        for a in exclude_aas:
            if a in tdict:
                del tdict[a]
            if a in bias_dim_dict['collapse']:
                bias_dim_dict['collapse'].remove(a)
            if a in bias_dim_dict['expand']:
                bias_dim_dict['expand'].remove(a)
            if bias_dim_dict['collapse']==[]:
                raise Exception('Cannot specify W, Y, G, F, Q, and N in exclude_aas')
            if bias_dim_dict['expand']==[]:
                raise Exception('Cannot specify D, E, K, P, S, and T in exclude_aas')
        # make new list
        nl=[]
        for a in tdict:
            for i in range(0,tdict[a]):
                nl.append(a)
        # set dis_list_to_use equal to nl
        dis_list_to_use=nl

    # try to make the sequence. 
    for i in range(0, num_sub_attempts):
        # base sequence
        start_seq=gen_sequence(seq_length, usedlist=dis_list_to_use)
        try:
            finseq=optimize_seq_dims(start_seq, rg_or_re, objective_dim, allowed_error=allowed_error,
                num_attempts=num_attempts, reduce_pos_charged=reduce_pos_charged, exclude_aas=exclude_aas)
            return finseq
        except:
            continue
    raise GooseFail(f'Could not generate sequence with specified {rg_or_re}. \nTry increasing number of attempts or with a different rg/re value or length value')


def make_rg_re_variant(sequence, increase_or_decrease, rg_or_re, 
    numseqs=32):
    """
    Function to make an Rg or Re variant of a sequence. Doesn't let you
    specify the Rg or Re, but you can get a dict of all unique Rg or Re 
    values if you feel like it. Otherwise just returns a decreased or 
    increased Rg / Re for the sequence. 

    Parameters
    sequence : str
        string of your amino acid sequence
    increase_or_decrease : str
        whether you want to increase or decrease the Rg / Re
    rg_or_re : str
        whether you want to increase or decrease the Rg / Re
    numseqs : int
        number of sequences to generate and predict Rg / Re for. 
        default is 32

    Returns
    -------
    dict
        dict of all the sequences by rg in order from
        smallest rg to largest rg. Easier to do it this way. 
    """
    # ceck some things
    if rg_or_re =='rg':
        start_dims=predict_rg(sequence)
    elif rg_or_re=='re':
        start_dims=predict_re(sequence)
    else:
        raise GooseInputError('rg_or_re must be either rg or re')
    # make sure we are increasing or decreasing. 
    if increase_or_decrease not in ['increase', 'decrease']:
        raise GooseInputError('increase_or_decrease must be either increase or decrease')
    # list to hold seqs
    seqs=[]
    # keep original in list in case we don't make anything better
    seqs.append(sequence)
    # make sequence a list so we can shuffle it
    # aas for start seq
    collapse=['Y', 'G', 'F', 'Q', 'N']
    expand=['D', 'E', 'K', 'P', 'S', 'T']
    if increase_or_decrease=='increase':
        startsequence=create_asymmetry_variant(sequence, 'decrease', aa_class=expand)
        startsequence=list(create_asymmetry_variant(sequence, 'increase', aa_class=collapse))
    else:
        startsequence=create_asymmetry_variant(sequence, 'decrease', aa_class=expand)
        startsequence=list(create_asymmetry_variant(sequence, 'decrease', aa_class=collapse))
    seqs.append(''.join(startsequence))
    # keep batch at 32
    numseqs = int(numseqs/32)*32
    if numseqs<32:
        numseqs=32
    # shuffle em'
    for i in range(0,numseqs-2):
        random.shuffle(startsequence)
        seqs.append(''.join(startsequence))
    # predict some things. 
    if rg_or_re=='rg':
        predictions=batch_predict_rg(seqs, return_seq2prediction=True)
    else:
        predictions=batch_predict_re(seqs, return_seq2prediction=True)
    # make seqs a list 
    seqs=list(predictions.keys())
    # make predicted dim a list
    dims=list(predictions.values())

    # make dict of seqs by rg in order. 
    seq_to_dim={}
    for seq, dim in zip(seqs, dims):
        if increase_or_decrease=='increase':
            if dim > start_dims:
                seq_to_dim[seq]=dim
        else:
            if dim < start_dims:
                seq_to_dim[seq]=dim
    if len(seq_to_dim)==0:
        raise GooseFail('Could not generate any sequences with different Rg / Re')
    else:
        return seq_to_dim

