'''
Testing whether it's possible to generate sequences with 
specific radius of gyration or end to end distance...

Very much in beta. 

'''

from matplotlib import pyplot as plt
from sparrow import Protein as pr
from sparrow.predictors import batch_predict
from goose.backend.sequence_generation_backend import gen_sequence
from goose.backend.lists import disordered_list
from goose.goose_exceptions import GooseError
import random
import numpy as np


def predict_rg(sequence): return pr(sequence).predictor.radius_of_gyration()
def predict_scaled_rg(sequence): return pr(sequence).predictor.radius_of_gyration(use_scaled=False)
def predict_e2e(sequence): return pr(sequence).predictor.end_to_end_distance()
def batch_predict_rg(sequences,show_progress_bar=False): 
    seqobs=[pr(seq) for seq in sequences]
    rgs=batch_predict.batch_predict(seqobs, 'rg',show_progress_bar=show_progress_bar)
    return rgs

def build_rg_seq(seq_length, rg, allowed_error=0.2, num_attempts=100):
    '''
    builds a sequence of a given length with a given radius of gyration
    by randomly selecting amino acids until the radius of gyration is
    this is a very inefficient algorithm, but it works without having
    to bring in new dependencies ...
    '''
    rg_change_dict={'collapse':['Y', 'G', 'F', 'Q', 'N'], 'expand':['D', 'E', 'K', 'P', 'S', 'T']}
    start_seq=gen_sequence(seq_length, usedlist=disordered_list)
    start_rg = predict_rg(start_seq)
    if abs(start_rg - rg)<=allowed_error:
        return start_seq
    else:
        best_error = abs(start_rg - rg)
        starter = list(start_seq).copy()
        for attempt in range(0,num_attempts):
            if start_rg > rg:
                use_aas = rg_change_dict['collapse']
            else:
                use_aas = rg_change_dict['expand']
            # make 32 seq variants. 
            test_seqs=[]
            target_fractions = np.linspace(0.01,1,32)
            for frac in target_fractions:
                num_from_rg_change = int(frac*seq_length)
                if num_from_rg_change==0:
                    num_from_rg_change=1
                if num_from_rg_change > seq_length:
                    num_from_rg_change = seq_length
                temp=[]
                for aa in range(0, num_from_rg_change):
                    temp.append(random.choice(use_aas))
                random.shuffle(starter)
                temp.extend(starter[num_from_rg_change:])
                random.shuffle(temp)
                test_seqs.append(''.join(temp))
            seq_rgs = batch_predict_rg(test_seqs)
            # get closest rg
            for val in seq_rgs:
                cur_er = abs(seq_rgs[val][1]-rg)
                if cur_er < best_error:
                    starter = list(seq_rgs[val][0])
                    best_error=cur_er
                    best_seq = "".join(starter)
            if best_error < allowed_error:
                return best_seq
    raise GooseError('Could not generate sequence with specified radius of gyration. Try increasing number of attempts')


