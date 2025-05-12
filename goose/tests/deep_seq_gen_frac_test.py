'''
Testing for generating sequences by specifying fractions of amino acids. 
3 different tests. 
1. only length 100 but goes across broad fraction space
2. 'actual use' - 10 to 500 amino acids and goes over most of the
fraction space. This is where I'd expect most users to be when
using this funciton.
3. 'long seq' - 10 to 10000 amino acids and goes over much
of the fraction space but is more limited for higher values. 

By doing these 3 tests this *fairly quickly* goes over a lot 
of different possible sequences to try to identify any bugs.
'''
import random

from goose import create
from goose.backend.protein import Protein as pr
from goose.backend import parameters
from goose.backend.parameters import MAX_FRACTION_DICT as max_fracs


def vals_in_range(num_vals, val_range, round_val=5, max_val=None, min_val=None):
    '''
    function to get vals in a range between specific values. 
    '''
    interval_val = float((val_range[1]-val_range[0])/(num_vals-1))
    vals = []
    for i in range(0, num_vals):
        if max_val==None:
            vals.append(round(val_range[0]+(interval_val*i), round_val))
        else:
            finval=round(val_range[0]+(interval_val*i), round_val) <= max_val
            if finval <= max_val and finval >=min_val:
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

def calc_error(objective_val, actual_val, allowed_error=None):
    if allowed_error==None:
        allowed_error = 0.0001
    cur_err = abs(objective_val-actual_val)
    if cur_err > allowed_error:
        return False    
    return True


def test_fractions_basic(min_length=100,
                   max_length=100,
                   seqs_per_length=5,
                   fractions_per_aa=19,
                   verbose=True):
    '''
    function to test fraction functionality
    is the basic test because it's what I used to get the max values for the amino acid
    max fractions. 500,000 attempts at each fraction. However, because this is repeated, 
    not going to go to max value. 
    '''

    # iterate through amino acids
    for aa in max_fracs.keys():
        # get current max fraction for the amino acid
        cur_max = max_fracs[aa]
        # get random sequence lengths
        seq_lengths = random_seq_lengths(seqs_per_length, min_length, max_length)

        # get fractions
        fractions = vals_in_range(fractions_per_aa, [0.01, cur_max-0.1], round_val=2)

        # iterate through fractions
        for frac in fractions:
            # iterate through sequence lengths
            for seq_length in seq_lengths:
                fraction_by_length = (round(seq_length*frac))/seq_length
                ind_val = 1/seq_length
                # make sure when rounding we don't go over the max fraction possible
                if fraction_by_length > cur_max:
                    while fraction_by_length > cur_max:
                        fraction_by_length-=ind_val
                if fraction_by_length < 0:
                    fraction_by_length = 0
                #print what we are up to if verbose
                if verbose:
                    print(f'Attempting to make length {seq_length} fraction {fraction_by_length} for {aa}')
                # make the sequence
                seq = create.seq_fractions(seq_length, **{aa:fraction_by_length}, attempts=500000)
                # verify fraction is as expected
                assert calc_error(seq.count(aa)/len(seq),fraction_by_length)==True


def test_fractions_actual_use(min_length=parameters.MINIMUM_LENGTH, 
                   max_length=500,
                   seqs_per_length=5,
                   fractions_per_aa=19,
                   verbose=False):
    '''
    function to test fraction functionality
    is actual use in that I don't expect people to deviate too high
    into extreme values, so going frm 0.01 to the max fraction -0.24 and 
    doing lengths from 10 to 500
    '''

    # iterate through amino acids
    for aa in max_fracs.keys():
        # get current max fraction for the amino acid
        cur_max = max_fracs[aa]
        # get random sequence lengths
        seq_lengths = random_seq_lengths(seqs_per_length, min_length, max_length)

        # get fractions - doing -0.22 from max
        fractions = vals_in_range(fractions_per_aa, [0.01, cur_max-0.24], round_val=2)

        # iterate through fractions
        for frac in fractions:
            # iterate through sequence lengths
            for seq_length in seq_lengths:
                fraction_by_length = (round(seq_length*frac))/seq_length
                ind_val = 1/seq_length
                # make sure when rounding we don't go over the max fraction possible
                if fraction_by_length > cur_max:
                    while fraction_by_length > cur_max:
                        fraction_by_length-=ind_val
                if fraction_by_length < 0:
                    fraction_by_length = 0
                #print what we are up to if verbose
                if verbose:
                    print(f'Attempting to make length {seq_length} fraction {fraction_by_length} for {aa}')
                # make the sequence
                seq = create.seq_fractions(seq_length, **{aa:fraction_by_length}, attempts=500000)
                # verify fraction is as expected
                assert calc_error(seq.count(aa)/len(seq),fraction_by_length)==True

def test_fractions_long_seq(min_length=parameters.MINIMUM_LENGTH, 
                   max_length=5000,
                   seqs_per_length=5,
                   fractions_per_aa=19,
                   verbose=False):
    '''
    function to test fraction functionality
    is for longer sequence length ranges but doesn't venture
    into extreme values, so going frm 0.01 to the max fraction -0.36 and 
    doing lengths from 10 to 5000
    Exception to this is C. Apparently I can make 100aa seqs with way higher C than
    I can for long sequences.
    '''

    # iterate through amino acids
    for aa in max_fracs.keys():
        # get current max fraction for the amino acid
        cur_max = max_fracs[aa]
        # get random sequence lengths
        seq_lengths = random_seq_lengths(seqs_per_length, min_length, max_length)

        # get fractions - doing -0.36 from max
        if aa=='C':
            cur_max=0.5
        fractions = vals_in_range(fractions_per_aa, [0.01, cur_max-0.36], round_val=2)

        # iterate through fractions
        for frac in fractions:
            # iterate through sequence lengths
            for seq_length in seq_lengths:
                fraction_by_length = (round(seq_length*frac))/seq_length
                ind_val = 1/seq_length
                # make sure when rounding we don't go over the max fraction possible
                if fraction_by_length > cur_max:
                    while fraction_by_length > cur_max:
                        fraction_by_length-=ind_val
                if fraction_by_length < 0:
                    fraction_by_length = 0
                #print what we are up to if verbose
                if verbose:
                    print(f'Attempting to make length {seq_length} fraction {fraction_by_length} for {aa}')
                # make the sequence
                seq = create.seq_fractions(seq_length, **{aa:fraction_by_length}, attempts=500000)
                # verify fraction is as expected
                assert calc_error(seq.count(aa)/len(seq),fraction_by_length)==True





