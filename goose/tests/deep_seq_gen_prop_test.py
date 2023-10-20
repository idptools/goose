'''
Thoroughly tests the sequence generation by properties functionality. 
Probably not necessary to do this every time, mainly for debugging
and making sure we have the correct guard rails on everything. 
'''

import random
from goose import create
from goose import analyze
from goose.backend.sequence_generation import check_disorder
import metapredict as meta
from goose.backend import parameters





def vals_in_range(num_vals, val_range, round_val=5):
    '''
    function to get vals in a range between specific values. 
    '''
    interval_val = float((val_range[1]-val_range[0])/(num_vals-1))
    vals = []
    for i in range(0, num_vals):
        vals.append(round(val_range[0]+(interval_val*i), round_val))
    return vals

def random_seq_lengths(number_seqs, min_length, max_length):
    '''
    function to get random sequence lengths
    '''
    seq_lengths = []
    for i in range(0, number_seqs):
        seq_lengths.append(random.randint(min_length, max_length))
    return seq_lengths

def calc_error(objective_val, actual_val, allowed_error):
    if allowed_error==None:
        allowed_error = 0.0001
    if abs(objective_val-actual_val) > allowed_error:
        return False    
    return True

# test creating a sequence of a length, any properties. 
def test_seq_lengths(number_seqs=10, verbose=False, 
    min_length=parameters.MINIMUM_LENGTH, 
    max_length=parameters.MAXIMUM_LENGTH):
    '''
    Test that we can create sequences of a length. 
    '''
    seq_lengths = random_seq_lengths(number_seqs, min_length, max_length)
    for num, length in enumerate(seq_lengths):
        if verbose:
            print(f'On seq {num+1}')
        seq = create.sequence(length)
        assert len(seq) == length

def test_fcr(seq_per_FCR=10, 
    FCR_values=100, verbose=False,
    min_length=parameters.MINIMUM_LENGTH, 
    max_length=1000, min_FCR = parameters.MINIMUM_FCR,
    max_FCR = parameters.MAXIMUM_FCR):
    '''
    Test that we can create sequences with a specific FCR. 
    '''
    # get random FCR values
    FCR_values = vals_in_range(FCR_values, [min_FCR, max_FCR])
    # create sequences
    for FCR_val in FCR_values:
        if verbose:
            print(f'on FCR {FCR_val}')
            # get random sequence lengths
            seq_lengths = random_seq_lengths(seq_per_FCR, min_length, max_length)
            for num, length in enumerate(seq_lengths):
                FCR_by_length = (round(length*FCR_val))/length
                if verbose:
                    print(f'On seq {num+1}, actual objective FCR = {FCR_by_length}, seq length = {length}')
                seq = create.sequence(length, FCR=FCR_by_length)
                assert calc_error(analyze.properties(seq, fractions=False)['FCR'], FCR_by_length, allowed_error=None)==True





