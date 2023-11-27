'''
Thoroughly tests the sequence generation by properties functionality. 
Probably not necessary to do this every time, mainly for debugging
and making sure we have the correct guard rails on everything. 
'''

import random
from goose import create
from goose.backend.protein import Protein as PR
from goose.backend.sequence_generation import check_disorder
import metapredict as meta
from goose.backend import parameters



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

# test creating a sequence of a length, any properties. 
def test_seq_lengths(number_seqs=10, verbose=False, 
    min_length=parameters.MINIMUM_LENGTH, 
    max_length=parameters.MAXIMUM_LENGTH):
    '''
    Test that we can create sequences of a length. 
    '''
    seqs_made=0
    seq_lengths = random_seq_lengths(number_seqs, min_length, max_length)
    for num, length in enumerate(seq_lengths):
        if verbose:
            print(f'On seq {num+1}')
        seq = create.sequence(length)
        assert len(seq) == length
        seqs_made+=1
    print(f'Successfully made {seqs_made} sequences.')

def test_fcr(seq_per_FCR=10, 
    FCR_values=50, verbose=False,
    min_length=parameters.MINIMUM_LENGTH, 
    max_length=500, min_FCR = parameters.MINIMUM_FCR,
    max_FCR = parameters.MAXIMUM_FCR):
    '''
    Test that we can create sequences with a specific FCR. 
    '''
    # get random FCR values
    seqs_made=0
    FCR_values = vals_in_range(FCR_values, [min_FCR, max_FCR])
    # create sequences
    for FCR_val in FCR_values:
        if verbose:
            print(f'on FCR {FCR_val}')
        # get random sequence lengths
        seq_lengths = random_seq_lengths(seq_per_FCR, min_length, max_length)
        for num, length in enumerate(seq_lengths):
            FCR_by_length = (round(length*FCR_val))/length
            seq = create.sequence(length, FCR=FCR_by_length)
            fin_FCR = PR(seq).FCR
            if verbose:
                print(f'On seq {num+1}, objective FCR = {FCR_by_length}, actual_FCR={fin_FCR}, seq length = {length}')
            assert calc_error(fin_FCR, FCR_by_length, allowed_error=None)==True
            seqs_made+=1
    print(f'Successfully made {seqs_made} sequences.')


def test_ncpr(seq_per_NCPR=10, 
    NCPR_values=50, verbose=False,
    min_length=parameters.MINIMUM_LENGTH, 
    max_length=500, min_NCPR = parameters.MINIMUM_NCPR,
    max_NCPR = parameters.MAXIMUM_NCPR):
    '''
    Test that we can create sequences with a specific NCPR. 
    '''
    # get random FCR values
    NCPR_values = vals_in_range(NCPR_values, [min_NCPR, max_NCPR])
    seqs_made=0
    # create sequences
    for NCPR_val in NCPR_values:
        if verbose:
            print(f'on NCPR {NCPR_val}')
        # get random sequence lengths
        seq_lengths = random_seq_lengths(seq_per_NCPR, min_length, max_length)
        for num, length in enumerate(seq_lengths):
            NCPR_by_length = (round(length*NCPR_val))/length
            seq = create.sequence(length, NCPR=NCPR_by_length)
            fin_NCPR_val = PR(seq).NCPR
            if verbose:
                print(f'On seq {num+1}, objective NCPR = {NCPR_by_length}, final NCPR = {fin_NCPR_val}, seq length = {length}')                
            assert calc_error(NCPR_by_length, fin_NCPR_val, allowed_error=None)==True
            seqs_made+=1
        if verbose:
            print()
    print(f'Successfully made {seqs_made} sequences.')

def test_hydropathy(numseqs=10, 
    prop_val=50, verbose=False,
    min_length=parameters.MINIMUM_LENGTH, 
    max_length=500, min_prop = parameters.MINIMUM_HYDRO,
    max_prop = parameters.MAXIMUM_HYDRO,
    prop_name='hydropathy'):
    '''
    Test that we can create sequences with a specific NCPR. 
    '''
    # get random FCR values
    prop_vals = vals_in_range(prop_val, [min_prop, max_prop], max_val=parameters.MAXIMUM_HYDRO, min_val = parameters.MINIMUM_HYDRO)
    seqs_made=0
    # create sequences
    for prop_val in prop_vals:
        if verbose:
            print(f'on hydropathy {prop_val}')
        # get random sequence lengths
        seq_lengths = random_seq_lengths(numseqs, min_length, max_length)
        for num, length in enumerate(seq_lengths):
            seq = create.sequence(length, hydropathy=prop_val)
            fin_prop_val = PR(seq).hydropathy
            if verbose:
                print(f'On seq {num+1}, objective {prop_name} = {round(prop_val,5)}, final {prop_name} = {round(fin_prop_val,5)}, seq length = {length}')
            assert calc_error(fin_prop_val, prop_val, allowed_error=parameters.HYDRO_ERROR)==True                
            seqs_made+=1
        if verbose:
            print()
    print(f'Successfully made {seqs_made} sequences.')

def test_fcr_ncpr(seq_per_val=5, NCPR_values=25,
    FCR_values=25, verbose=False,
    min_length=parameters.MINIMUM_LENGTH, 
    max_length=500, min_NCPR = parameters.MINIMUM_NCPR,
    max_NCPR = parameters.MAXIMUM_NCPR,
    min_FCR = parameters.MINIMUM_FCR,
    max_FCR = parameters.MAXIMUM_FCR):
    '''
    Test that we can create sequences with a specific NCPR and FCR. 
    '''
    # get random FCR/NCPR values
    NCPR_values = vals_in_range(NCPR_values, [min_NCPR, max_NCPR])
    FCR_values = vals_in_range(FCR_values, [min_FCR, max_FCR])
    seqs_made=0
    # create sequences
    for FCR_val in FCR_values:
        for NCPR_val in NCPR_values:
            if abs(NCPR_val)<=FCR_val:
                if verbose:
                    print(f'on NCPR {NCPR_val}, FCR {FCR_val}')
                # get random sequence lengths
                seq_lengths = random_seq_lengths(seq_per_val, min_length, max_length)
                for num, length in enumerate(seq_lengths):
                    NCPR_by_length = (round(length*NCPR_val))/length
                    FCR_by_length = (round(length*FCR_val))/length
                    if abs(NCPR_by_length)<=FCR_by_length:
                        seq = create.sequence(length, NCPR=NCPR_by_length, FCR=FCR_by_length)
                        fin_NCPR_val = PR(seq).NCPR
                        if verbose:
                            print(f'On seq {num+1}, objective NCPR = {NCPR_by_length}, final NCPR = {fin_NCPR_val}, seq length = {length}')                
                        assert calc_error(NCPR_by_length, fin_NCPR_val, allowed_error=None)==True
                        seqs_made+=1
                    if verbose:
                        print()
    print(f'Successfully made {seqs_made} sequences.')


def test_fcr_ncpr_hydropathy(seq_per_val=5, NCPR_values=10,
    FCR_values=10, hydropathy_values=10, verbose=False,
    min_length=parameters.MINIMUM_LENGTH, 
    max_length=500, min_NCPR = parameters.MINIMUM_NCPR,
    max_NCPR = parameters.MAXIMUM_NCPR,
    min_FCR = parameters.MINIMUM_FCR,
    max_FCR = parameters.MAXIMUM_FCR,
    min_hydro = parameters.MINIMUM_HYDRO,
    max_hydro = parameters.MAXIMUM_HYDRO_CHARGED):
    '''
    Test that we can create sequences with a specific NCPR and FCR. 
    '''
    # get random FCR/NCPR values
    NCPR_values = vals_in_range(NCPR_values, [min_NCPR, max_NCPR])
    FCR_values = vals_in_range(FCR_values, [min_FCR, max_FCR])
    hydro_values = vals_in_range(hydropathy_values, [min_hydro, max_hydro])
    seqs_made=0
    # create sequences
    for FCR_val in FCR_values:
        for NCPR_val in NCPR_values:
            if abs(NCPR_val)<=FCR_val:
                if verbose:
                    print(f'on NCPR {NCPR_val}, FCR {FCR_val}')
                # get random sequence lengths
                seq_lengths = random_seq_lengths(seq_per_val, min_length, max_length)
                for num, length in enumerate(seq_lengths):
                    NCPR_by_length = (round(length*NCPR_val))/length
                    FCR_by_length = (round(length*FCR_val))/length
                    if abs(NCPR_by_length)<=FCR_by_length:
                        seq = create.sequence(length, NCPR=NCPR_by_length, FCR=FCR_by_length)
                        fin_NCPR_val = PR(seq).NCPR
                        if verbose:
                            print(f'On seq {num+1}, objective NCPR = {NCPR_by_length}, final NCPR = {fin_NCPR_val}, seq length = {length}')                
                        assert calc_error(NCPR_by_length, fin_NCPR_val, allowed_error=None)==True
                        seqs_made+=1
                    if verbose:
                        print()
    print(f'Successfully made {seqs_made} sequences.')