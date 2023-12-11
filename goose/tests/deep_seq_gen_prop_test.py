'''
Thoroughly tests the sequence generation by properties functionality. 
By thoroughly, I mean THOROUGHLY. Like figuring out the number of grains 
of sand on a beach by going to the beach and counting the individual grains 
kind of thoroughly. 


Probably not necessary to do this every time, mainly for debugging
and making sure we have the correct guard rails on everything. 
'''

import random
from goose import create
from goose.backend.protein import Protein as PR
from goose.backend.sequence_generation import check_disorder
from goose.backend.sequence_generation_backend import fraction_net_charge
import metapredict as meta
from goose.backend import parameters

#=====================#
# Functions for tests #
#=====================#

def check_hydro_is_possible_testing_only(length, hydropathy, FCR=None, NCPR=None, max_theoretical_val=5.6, print_range=False):
    '''
    Function purely for testing hydropathy, FCR, and NCPR ranges. 
    This function is not used for sequence generation. It's just a way
    to choose which ranges of hydropathy / FCR / NCPR to explore when doing testing to
    check for bugs. 

    Parameters
    ----------
    length : Int
        The length of the sequence as an integer value
    FCR : Float
        The fraction of charged residues wanted for the sequence as a
        decimal value.
    NCPR : Float
        The wanted net charge of the sequence given as a decimal value
    hydropathy : Float
        The wanted mean hydropathy value of the sequence. Uses ajusted
        Kyte-doolittle hydropathy scale that goes from 0 to 9

    max_theoretical_val : float
        A value that is the approximate hydropathy of an amino acid in a sequence. 
        Basically overrides things to not be set to be all isoleucine, which generally doesn't work. 

    Returns
    --------
    Boolean
        Returns True if sequence is possible and False if it is not.
    '''
    # if not specifying FCR or NCPR, just return True. 
    if FCR==None and NCPR==None:
        return True
    # hold charged residues. 
    charges=''
    # account for NCPR
    if NCPR != None:
        added_NCPR_res = round(length*abs(NCPR))
        for i in range(0, added_NCPR_res):
            if NCPR < 0:
                charges+='-'
            else:
                charges+='+'

    # account for FCR
    if FCR != None:
        if NCPR != None:
            # taken from how this is calculated i the backend, which prioritizes NCPR of FCR
            #figure out the remaining FCR fraction
            remaining_FCR_fraction = FCR - abs(NCPR)
            # based on the remamining_FCR_fraction, figure out how many residues to add
            added_FCR_residues = round((length * remaining_FCR_fraction)/2)
            remaining_residues_over_two = round(length - added_NCPR_res/2)

            # figure ot the number of possible residues to add to the sequence
            remaining_residues = (length - added_NCPR_res)

            # if the number of added_FCR_residues are less than or equal to 
            # remaining residues over 2, added_residues = added_FCR_residues        
            if added_FCR_residues <= remaining_residues_over_two:
                added_residues = added_FCR_residues 

            # otherwise, added_residues are qualt to remaining_residues_over_two       
            else:
                added_residues = remaining_residues_over_two

            # if the number of residues to add is greater than the remaining number of residues...
            if (added_residues * 2) > remaining_residues:
                added_residues = added_residues - 1
            # now add the residues needed to get FCR to be as close as possible.
            for res in range(0, added_residues):
                charges+='+-'
        else:
            # if NCPR = None, just get the FCR correct. 
            FCR_residues = round(length*FCR)
            for i in range(0, FCR_residues):
                # 0 can be either. 
                charges+='0'

    # now get the min and max values for the hydropathy based on charges alone. 
    min_hydro_by_charge=0
    max_hydro_by_charge=0
    for charge in charges:
        if charge == '-':
            min_hydro_by_charge+=1
            max_hydro_by_charge+=1
        elif charge == '+':
            min_hydro_by_charge+=0
            max_hydro_by_charge+=0.6
        else:
            min_hydro_by_charge+=0
            max_hydro_by_charge+=1
    # add in 1 for remaning for min (Q/N) and theoretical val for other
    for i in range(0, length-len(charges)):
        min_hydro_by_charge+=1
        max_hydro_by_charge+=max_theoretical_val

    # divide by length
    minval = min_hydro_by_charge/length
    maxval = max_hydro_by_charge/length
    if print_range==True:
        print(f'Min: {minval}, Max: {maxval}')
    if hydropathy >= minval and hydropathy <= maxval:
        return True
    return False




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

#=================#
#  No Properties  #
#=================#

# test creating a sequence of a length, any properties. 
def test_seq_lengths(number_seqs=100, verbose=False, 
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


#=================#
#    1 Property   #
#=================#

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
            assert len(seq) == length
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
            assert len(seq) == length
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
            assert len(seq) == length
            seqs_made+=1
        if verbose:
            print()
    print(f'Successfully made {seqs_made} sequences.')


# setting min length to 15 for kappa. It can do down to 10 more often than not
# but due to random selection of the FCR / NCPR, super short sequences have 
# an occassional fail. 

# making kappa minimum 0.03 for testing. GOOSE can go down to 0, but it takes a very long time
# and can be pretty tricky. 

def test_kappa(seq_per_kappa=10, 
    kappa_values=40, verbose=False,
    min_length=15,
    max_length=500, min_kappa = 0.03,
    max_kappa = parameters.MAXIMUM_KAPPA):
    '''
    Test that we can create sequences with a specific NCPR. 
    '''
    # get random kappa values
    kappa_values = vals_in_range(kappa_values, [min_kappa, max_kappa])
    seqs_made=0
    # create sequences
    for kappa_val in kappa_values:
        if verbose:
            print(f'on kappa {kappa_val}')
        # get random sequence lengths
        seq_lengths = random_seq_lengths(seq_per_kappa, min_length, max_length)
        for num, length in enumerate(seq_lengths):
            if verbose == True:
                print(f'Trying to make seq length={length}, kappa={kappa_val}')
            seq = create.sequence(length, kappa=kappa_val)
            fin_kappa_val = PR(seq).kappa
            if verbose:
                print(f'On seq {num+1}, objective kappa = {kappa_val}, final kappa = {fin_kappa_val}, seq length = {length}')                
            assert calc_error(kappa_val, fin_kappa_val, allowed_error=parameters.MAXIMUM_KAPPA_ERROR)==True
            assert len(seq) == length
            seqs_made+=1
        if verbose:
            print()
    print(f'Successfully made {seqs_made} sequences.')


#=================#
#  2 Properties   #
#                 #
# Notes: For 2    #
# properties, the #
# max length is   #
# set to 300.     #
#=================#

def test_fcr_ncpr(seq_per_val=5, NCPR_values=25,
    FCR_values=25, verbose=False,
    min_length=parameters.MINIMUM_LENGTH, 
    max_length=300, min_NCPR = parameters.MINIMUM_NCPR,
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
                        fin_FCR_val = PR(seq).FCR
                        if verbose:
                            print(f'Seq_{num+1}: obj.NCPR={NCPR_by_length}, fin.NCPR={fin_NCPR_val}: obj.FCR={FCR_by_length}, fin.FCR={fin_FCR_val}: obj.len={length}, fin.len={len(seq)}')
                        assert calc_error(NCPR_by_length, fin_NCPR_val, allowed_error=None)==True
                        # priortizes geting NCPR correct instead of FCR if it is not possible to get both. 
                        assert calc_error(FCR_by_length, fin_FCR_val, allowed_error=(1/length)+0.000001)
                        assert len(seq) == length
                        seqs_made+=1
                    if verbose:
                        print()
    print(f'Successfully made {seqs_made} sequences.')


def test_fcr_hydropathy(seq_per_val=5, hydro_values=25,
    # note: This doesn't go to the in theory max hydropatyh value
    # possible for every charge. It get's pretty close though. 
    FCR_values=25, verbose=False,
    min_length=parameters.MINIMUM_LENGTH, 
    max_length=300, min_FCR = parameters.MINIMUM_FCR,
    max_FCR = parameters.MAXIMUM_FCR,
    min_hydro = parameters.MINIMUM_HYDRO_CHARGED,
    max_hydro = parameters.MAXIMUM_HYDRO_CHARGED):
    '''
    Test that we can create sequences with a specific hydropathy and FCR. 
    '''
    # get random FCR/NCPR values
    hydro_values = vals_in_range(hydro_values, [min_hydro, max_hydro])
    FCR_values = vals_in_range(FCR_values, [min_FCR, max_FCR])
    seqs_made=0
    # create sequences
    for FCR_val in FCR_values:
        for hydro_val in hydro_values:
            if verbose:
                print(f'on hydro {hydro_val}, FCR {FCR_val}')
            # get random sequence lengths
            seq_lengths = random_seq_lengths(seq_per_val, min_length, max_length)
            for num, length in enumerate(seq_lengths):
                FCR_by_length = (round(length*FCR_val))/length
                if check_hydro_is_possible_testing_only(length, hydro_val, FCR=FCR_by_length):
                    seq = create.sequence(length, FCR=FCR_by_length, hydropathy=hydro_val, attempts=1000)
                    fin_FCR_val = PR(seq).FCR
                    fin_hydro_val = PR(seq).hydropathy
                    if verbose:
                        print(f'On seq {num+1}, objective_hydro = {hydro_val}, fin_hydro_val = {fin_hydro_val}, objective FCR = {FCR_by_length}, final FCR = {fin_FCR_val}, seq length = {length}')                
                    assert calc_error(FCR_by_length, fin_FCR_val, allowed_error=None)==True
                    assert calc_error(fin_hydro_val, hydro_val, allowed_error=parameters.HYDRO_ERROR+0.00001)==True
                    assert len(seq) == length
                    seqs_made+=1
            if verbose:
                print()
                
    print(f'Successfully made {seqs_made} sequences.')




def test_ncpr_hydropathy(seq_per_val=10, hydro_values=25,
    NCPR_values=25, verbose=False,
    min_length=parameters.MINIMUM_LENGTH, 
    max_length=300, min_NCPR = parameters.MINIMUM_NCPR,
    max_NCPR = parameters.MAXIMUM_NCPR,
    min_hydro = parameters.MINIMUM_HYDRO_CHARGED,
    max_hydro = parameters.MAXIMUM_HYDRO_CHARGED):
    # note: This doesn't go to the in theory max hydropathy value
    # possible for every charge. It get's pretty close though. 
    # this is mainly because of edge cases. 
    '''
    Test that we can create sequences with a specific hydropathy and NCPR. 
    '''
    # get random NCPR/NCPR values
    hydro_values = vals_in_range(hydro_values, [min_hydro, max_hydro])
    NCPR_values = vals_in_range(NCPR_values, [min_NCPR, max_NCPR])
    seqs_made=0
    # create sequences
    for NCPR_val in NCPR_values:
        for hydro_val in hydro_values:
            # get random sequence lengths
            seq_lengths = random_seq_lengths(seq_per_val, min_length, max_length)
            for num, length in enumerate(seq_lengths):                
                NCPR_by_length = (round(length*NCPR_val))/length
                if check_hydro_is_possible_testing_only(length, hydro_val, NCPR=NCPR_by_length)==True:
                    if verbose:
                        print(f'Trying to make length={length}, NCPR={NCPR_by_length}, hydropathy={hydro_val}')
                    seq = create.sequence(length, NCPR=NCPR_by_length, hydropathy=hydro_val, attempts=1000)
                    fin_NCPR_val = PR(seq).NCPR
                    fin_hydro_val = PR(seq).hydropathy
                    if verbose:
                        print(f'On seq {num+1}, objective_hydro = {hydro_val}, fin_hydro_val = {fin_hydro_val}, objective NCPR = {NCPR_by_length}, final NCPR = {fin_NCPR_val}, seq length = {length}')                
                    assert calc_error(NCPR_by_length, fin_NCPR_val, allowed_error=None)==True
                    assert calc_error(fin_hydro_val, hydro_val, allowed_error=parameters.HYDRO_ERROR+0.00001)==True
                    assert len(seq) == length
                    seqs_made+=1
                
    print(f'Successfully made {seqs_made} sequences.')

# making kappa minimum 0.1 for testing this. GOOSE can go down to 0, but it takes a very long time
# and can be pretty tricky. 

def test_ncpr_kappa(seq_per_val=5, kappa_values=19,
    NCPR_values=19, verbose=False,
    min_length=parameters.MINIMUM_LENGTH, 
    max_length=300, min_NCPR = -0.5,
    max_NCPR = 0.5,
    min_kappa = 0.1,
    max_kappa = 0.95):
    # note: This doesn't go to the in theory min kappa value
    # or the full range of NCPR. Reason for this is you can't really do
    # the full range of kappa for super positive or negatively charged seqs. 
    # or if you can... GOOSE can't. 
    # It get's pretty close though. 
    '''
    Test that we can create sequences with a specific kappa and NCPR. 
    '''
    # get random NCPR/NCPR values
    kappa_values = vals_in_range(kappa_values, [min_kappa, max_kappa])
    NCPR_values = vals_in_range(NCPR_values, [min_NCPR, max_NCPR])
    seqs_made=0
    # create sequences
    for NCPR_val in NCPR_values:
        for kappa_val in kappa_values:
            if verbose:
                print(f'on kappa {kappa_val}, NCPR {NCPR_val}')
            # get random sequence lengths
            seq_lengths = random_seq_lengths(seq_per_val, min_length, max_length)
            for num, length in enumerate(seq_lengths):
                NCPR_by_length = (round(length*NCPR_val))/length
                seq = create.sequence(length, NCPR=NCPR_by_length, kappa=kappa_val, attempts=1000)
                fin_NCPR_val = PR(seq).NCPR
                fin_kappa_val = PR(seq).kappa
                fin_len=len(seq)
                if verbose:
                    print(f'On seq {num+1}, obj.kappa = {kappa_val}, fin.kappa = {fin_kappa_val}: obj.NCPR = {NCPR_by_length}, fin.NCPR = {fin_NCPR_val}: obj.len = {length}, fin.len={len(seq)}')                
                assert calc_error(NCPR_by_length, fin_NCPR_val, allowed_error=None)==True
                assert calc_error(fin_kappa_val, kappa_val, allowed_error=parameters.MAXIMUM_KAPPA_ERROR+0.00001)==True
                assert len(seq) == length
                seqs_made+=1
            if verbose:
                print()
                
    print(f'Successfully made {seqs_made} sequences.')

# making kappa minimum 0.1 for testing this. GOOSE can go down to 0, but it takes a very long time
# and can be pretty tricky. 

def test_FCR_kappa(seq_per_val=5, kappa_values=19,
    FCR_values=19, verbose=False,
    min_length=30, 
    max_length=300, min_FCR = 0.2,
    max_FCR = 1,
    min_kappa = 0.1,
    max_kappa = parameters.MAXIMUM_KAPPA):
    # note: This doesn't go to the in theory min kappa value
    # or the full range of FCR. Reason for this is you can't really do
    # the full range of kappa if you have a short sequence with a low FCR
    # because there' won't be one of each charged residues. 
    # or if you can... GOOSE can't. 
    # It get's pretty close though. 
    '''
    Test that we can create sequences with a specific kappa and FCR. 
    '''
    # get random FCR/FCR values
    kappa_values = vals_in_range(kappa_values, [min_kappa, max_kappa])
    FCR_values = vals_in_range(FCR_values, [min_FCR, max_FCR])
    seqs_made=0
    # create sequences
    for FCR_val in FCR_values:
        for kappa_val in kappa_values:
            if verbose:
                print(f'on kappa {kappa_val}, FCR {FCR_val}')
            # get random sequence lengths
            seq_lengths = random_seq_lengths(seq_per_val, min_length, max_length)
            for num, length in enumerate(seq_lengths):
                FCR_by_length = (round(length*FCR_val))/length
                if verbose:
                    print(f'Trying to make length {length}, FCR {FCR_by_length}, kappa {kappa_val}')
                seq = create.sequence(length, FCR=FCR_by_length, kappa=kappa_val, attempts=1000)
                fin_FCR_val = PR(seq).FCR
                fin_kappa_val = PR(seq).kappa
                fin_len=len(seq)
                if verbose:
                    print(f'On seq {num+1}, obj.kappa = {round(kappa_val,4)}, fin.kappa = {round(fin_kappa_val,4)}: obj.FCR = {round(FCR_by_length,4)}, fin.FCR = {round(fin_FCR_val,4)}: obj.len = {length}, fin.len={len(seq)}, NCPR = {round(PR(seq).NCPR,4)}\n')                
                assert calc_error(FCR_by_length, fin_FCR_val, allowed_error=None)==True
                assert calc_error(fin_kappa_val, kappa_val, allowed_error=parameters.MAXIMUM_KAPPA_ERROR+0.00001)==True
                assert len(seq) == length
                seqs_made+=1
            if verbose:
                print()
                
    print(f'Successfully made {seqs_made} sequences.')

# making kappa minimum 0.1 for testing this. GOOSE can go down to 0, but it takes a very long time
# and can be pretty tricky. 

# This sequence generation is pretty slow and I can't get it much 
# faster so limiting the scope of what is tested a little here.
def test_kappa_hydropathy(seq_per_val=5, hydro_values=10,
    kappa_values=10, verbose=False,
    min_length=30,
    max_length=300, min_kappa = 0.1,
    max_kappa = 0.9,
    min_hydro = 1,
    max_hydro = 5):
    # note: This doesn't go to the in theory max hydropatyh value
    # possible for every charge. It get's pretty close though. 
    # min length set to 25 because kappa is hard to do with 
    # a sequence that's too short. 
    '''
    Test that we can create sequences with a specific hydropathy and kappa. 
    '''
    # get random kappa/NCPR values
    hydro_values = vals_in_range(hydro_values, [min_hydro, max_hydro])
    kappa_values = vals_in_range(kappa_values, [min_kappa, max_kappa])
    seqs_made=0
    # create sequences
    for kappa_val in kappa_values:
        for hydro_val in hydro_values:
            if verbose:
                print(f'on hydro {hydro_val}, kappa {kappa_val}')
            # get random sequence lengths
            seq_lengths = random_seq_lengths(seq_per_val, min_length, max_length)
            for num, length in enumerate(seq_lengths):
                if verbose == True:
                    print(f'Trying to make seq length={length}, kappa={kappa_val}, hydro={hydro_val}')
                seq = create.sequence(length, kappa=kappa_val, hydropathy=hydro_val, attempts=1000)
                fin_kappa_val = PR(seq).kappa
                fin_hydro_val = PR(seq).hydropathy
                fin_FCR_val = PR(seq).FCR
                if verbose:
                    print(f'On seq {num+1}, obj.kappa = {round(kappa_val,4)}, fin.kappa = {round(fin_kappa_val,4)}: obj.hydro = {round(hydro_val,4)}, fin.hydro = {round(fin_hydro_val,4)}: obj.len = {length}, fin.len={len(seq)}, NCPR = {round(PR(seq).NCPR,4)}, FCR = {round(PR(seq).FCR,4)}\n')                
                assert calc_error(kappa_val, fin_kappa_val, allowed_error=parameters.MAXIMUM_KAPPA_ERROR+0.00001)==True
                assert calc_error(fin_hydro_val, hydro_val, allowed_error=parameters.HYDRO_ERROR+0.00001)==True
                assert len(seq) == length
                seqs_made+=1
            if verbose:
                print()
                
    print(f'Successfully made {seqs_made} sequences.')


#=================#
#  3 Properties   #
#=================#

def test_fcr_ncpr_hydropathy(seq_per_val=10, 
    NCPR_values=20,
    FCR_values=20, 
    hydropathy_values=20, 
    verbose=False,
    min_length=15,
    max_length=300, 
    min_NCPR = parameters.MINIMUM_NCPR,
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
                for hydro_val in hydro_values:
                    # get random sequence lengths
                    seq_lengths = random_seq_lengths(seq_per_val, min_length, max_length)
                    for num, length in enumerate(seq_lengths):
                        NCPR_by_length = (round(length*NCPR_val))/length
                        FCR_by_length = (round(length*FCR_val))/length
                        if abs(NCPR_by_length)<=FCR_by_length:
                            if check_hydro_is_possible_testing_only(length, hydro_val, 
                                                                    NCPR=NCPR_by_length, 
                                                                    FCR=FCR_by_length, 
                                                                    max_theoretical_val=5.6):
                                if verbose:
                                    print(f'Making seq - Length:{length}, NCPR:{NCPR_by_length}, FCR:{FCR_by_length}, hydro:{hydro_val}')
                                seq = create.sequence(length, NCPR=NCPR_by_length, FCR=FCR_by_length, hydropathy=hydro_val)
                                fin_NCPR_val = PR(seq).NCPR
                                fin_FCR_val = PR(seq).FCR
                                fin_hydro_val = PR(seq).hydropathy
                                if verbose:
                                    print(f'obj.NCPR = {round(NCPR_by_length,4)}, fin.NCPR = {round(fin_NCPR_val,4)}: obj.hydro = {round(hydro_val,4)}, fin.hydro = {round(fin_hydro_val,4)}: obj.NCPR = {NCPR_by_length}, fin.NCPR = {fin_NCPR_val}: obj.len = {length}, fin.len={len(seq)}\n')                
                                assert calc_error(NCPR_by_length, fin_NCPR_val, allowed_error=None)==True
                                assert calc_error(FCR_by_length, fin_FCR_val, allowed_error=(1/length)+0.00001)==True
                                assert calc_error(hydro_val, fin_hydro_val, allowed_error = parameters.HYDRO_ERROR+0.00001)
                                assert len(seq) == length
                                seqs_made+=1
    print(f'Successfully made {seqs_made} sequences.')


''' 
TO DO
FCR, NCPR, Hydropathy
FCR, NCPR, kappa
FCR, kappa, hydropathy
NCPR, kappa, hydropathy
'''


"""
#=================#
#  4 Properties   #
#=================#

# to do = all 4 properties. 
"""