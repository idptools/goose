'''
Code to test the perforance of GOOSE sequence generation 
'''
# main imports
import random
import time
import numpy as np

from goose import create, analyze
from goose.backend.variant_generation_backend import hydro_range_constant_classes
from goose.backend import parameters

def test_goose_just_length(seq_length=100, num_seqs=100):                            
    '''
    Function that tests how long it takes GOOSE to make sequences of
    a specified length.

    Parameters
    ----------
    seq_length: int
        Length of the sequences to be generated
    num_seqs: int
        Number of sequences to generate 

    Returns
    -------
    durations
        Time in seconds it to generate each sequence
    '''
    # track durations per seq
    durations=[]
    for i in range(num_seqs):
        # start timer
        start = time.time()
        create.sequence(seq_length)
        durations.append(round(time.time()-start, 5))
    
    # return how long it took
    return durations


def test_goose_single_property_generation(seq_length=100, seq_per_prop=25,
                                    hydropathy_range=[1.5, 3.5],
                                    FCR_range = [0, 0.5],
                                    NCPR_range=[-0.5, 0.5],
                                    kappa_range=[0.1, 0.8]):
    '''
    Function that tests how long it takes GOOSE to make sequences of
    a specified length that have randomly selected properties within
    the above specified ranges. 

    Parameters
    ----------
    seq_length: int
        Length of the sequences to be generated
    seq_per_prop: int
        Number of sequences to generate per property
    hydropathy_range: list
        List of two floats representing the range of hydropathy values
        to be generated
    FCR_range: list
        List of two floats representing the range of FCR values to be
        generated
    NCPR_range: list
        List of two floats representing the range of NCPR values to be
        generated
    kappa_range: list
        List of two floats representing the range of kappa values to be
        generated
        
    Returns
    -------
    durations
        time it took to make each sequence
    '''
    hydropathy_vals=np.arange(hydropathy_range[0], hydropathy_range[1], 0.1)
    FCR_vals=np.arange(FCR_range[0], FCR_range[1], 0.01)
    NCPR_vals=np.arange(NCPR_range[0], NCPR_range[1], 0.01)
    kappa_vals=np.arange(kappa_range[0], kappa_range[1], 0.01)
    
    # lists to append randomly selected values to. 
    list_hydro_vals=[]
    list_fcr_vals=[]
    list_ncpr_vals=[]
    list_kappa_vals=[]
    
    for i in range(seq_per_prop):
        list_hydro_vals.append(random.choice(hydropathy_vals))
        list_fcr_vals.append(random.choice(FCR_vals))
        list_ncpr_vals.append(random.choice(NCPR_vals))
        list_kappa_vals.append(random.choice(kappa_vals))

    # track durations
    durations=[]

    for i in range(seq_per_prop):
        start = time.time()
        create.sequence(seq_length, hydropathy=list_hydro_vals[i])
        durations.append(round(time.time()-start, 5))
        start = time.time()
        create.sequence(seq_length, FCR=list_fcr_vals[i])
        durations.append(round(time.time()-start, 5))
        start = time.time()
        create.sequence(seq_length, NCPR=list_ncpr_vals[i])
        durations.append(round(time.time()-start, 5))
        start = time.time()
        create.sequence(seq_length, kappa=list_kappa_vals[i])
        durations.append(round(time.time()-start, 5))
    
    # end time, return how long it took
    return durations


# for seq fractions
def test_generation_by_fractions(seq_length=100, 
                                seq_per_frac=5,
                                max_A=0.55, max_R=0.6,
                                max_N=0.6, max_D=0.6,
                                max_C=0.6, max_Q=0.6,
                                max_E=0.6, max_G=0.6,
                                max_H=0.6, max_I=0.13,
                                max_L=0.02, max_K=0.6,
                                max_M=0.22, max_F=0.6,
                                max_P=0.6, max_S=0.6,
                                max_T=0.6, max_W=0.15,
                                max_Y=0.59, max_V=0.31):
    '''
    Function that tests how long it takes GOOSE to make sequences of
    a specified length that have randomly selected fractions of individual
    amino acids between 0 and the above specified max value.

    Parameters
    ----------
    seq_length: int
        Length of the sequences to be generated
    seq_per_frac: int
        Number of sequences to generate per fraction
    max_A: float
        Maximum fraction of A to be generated
    max_R: float
        Maximum fraction of R to be generated
    max_N: float
        Maximum fraction of N to be generated
    max_D: float
        Maximum fraction of D to be generated
    max_C: float
        Maximum fraction of C to be generated
    max_Q: float
        Maximum fraction of Q to be generated
    max_E: float
        Maximum fraction of E to be generated
    max_G: float
        Maximum fraction of G to be generated
    max_H: float
        Maximum fraction of H to be generated
    max_I: float
        Maximum fraction of I to be generated
    max_L: float
        Maximum fraction of L to be generated
    max_K: float
        Maximum fraction of K to be generated
    max_M: float
        Maximum fraction of M to be generated
    max_F: float
        Maximum fraction of F to be generated
    max_P: float
        Maximum fraction of P to be generated
    max_S: float
        Maximum fraction of S to be generated
    max_T: float
        Maximum fraction of T to be generated
    max_W: float
        Maximum fraction of W to be generated
    max_Y: float
        Maximum fraction of Y to be generated
    max_V: float
        Maximum fraction of V to be generated

    Returns
    -------
    durations
        time it took to make each sequence
    '''
    # make amino acid values. 
    A_vals=np.arange(0, max_A, 0.01)
    R_vals=np.arange(0, max_R, 0.01)
    N_vals=np.arange(0, max_N, 0.01)
    D_vals=np.arange(0, max_D, 0.01)
    C_vals=np.arange(0, max_C, 0.01)
    Q_vals=np.arange(0, max_Q, 0.01)
    E_vals=np.arange(0, max_E, 0.01)
    G_vals=np.arange(0, max_G, 0.01)
    H_vals=np.arange(0, max_H, 0.01)
    I_vals=np.arange(0, max_I, 0.01)
    L_vals=np.arange(0, max_L, 0.01)
    K_vals=np.arange(0, max_K, 0.01)
    M_vals=np.arange(0, max_M, 0.01)
    F_vals=np.arange(0, max_F, 0.01)
    P_vals=np.arange(0, max_P, 0.01)
    S_vals=np.arange(0, max_S, 0.01)
    T_vals=np.arange(0, max_T, 0.01)
    W_vals=np.arange(0, max_W, 0.01)
    Y_vals=np.arange(0, max_Y, 0.01)
    V_vals=np.arange(0, max_V, 0.01)
    
    # lists to append randomly selected values to. 
    # we do this to take out choosing random values from the 
    # timing of the function. 
    list_A_vals=[]
    list_R_vals=[]
    list_N_vals=[]
    list_D_vals=[]
    list_C_vals=[]
    list_Q_vals=[]
    list_E_vals=[]
    list_G_vals=[]
    list_H_vals=[]
    list_I_vals=[]
    list_L_vals=[]
    list_K_vals=[]
    list_M_vals=[]
    list_F_vals=[]
    list_P_vals=[]
    list_S_vals=[]
    list_T_vals=[]
    list_W_vals=[]
    list_Y_vals=[]
    list_V_vals=[]

    for i in range(seq_per_frac):
        list_A_vals.append(random.choice(A_vals))
        list_R_vals.append(random.choice(R_vals))
        list_N_vals.append(random.choice(N_vals))
        list_D_vals.append(random.choice(D_vals))
        list_C_vals.append(random.choice(C_vals))
        list_Q_vals.append(random.choice(Q_vals))
        list_E_vals.append(random.choice(E_vals))
        list_G_vals.append(random.choice(G_vals))
        list_H_vals.append(random.choice(H_vals))
        list_I_vals.append(random.choice(I_vals))
        list_L_vals.append(random.choice(L_vals))
        list_K_vals.append(random.choice(K_vals))
        list_M_vals.append(random.choice(M_vals))
        list_F_vals.append(random.choice(F_vals))
        list_P_vals.append(random.choice(P_vals))
        list_S_vals.append(random.choice(S_vals))
        list_T_vals.append(random.choice(T_vals))
        list_W_vals.append(random.choice(W_vals))
        list_Y_vals.append(random.choice(Y_vals))
        list_V_vals.append(random.choice(V_vals))

    # track durations per seq
    durations=[]

    for i in range(seq_per_frac):
        start=time.time()
        create.seq_fractions(seq_length, A=list_A_vals[i])
        durations.append(round(time.time()-start, 5))
        start=time.time()
        create.seq_fractions(seq_length, R=list_R_vals[i])
        durations.append(round(time.time()-start, 5))
        start=time.time()
        create.seq_fractions(seq_length, N=list_N_vals[i])
        durations.append(round(time.time()-start, 5))
        start=time.time()
        create.seq_fractions(seq_length, D=list_D_vals[i])
        durations.append(round(time.time()-start, 5))
        start=time.time()
        create.seq_fractions(seq_length, C=list_C_vals[i])
        durations.append(round(time.time()-start, 5))
        start=time.time()
        create.seq_fractions(seq_length, Q=list_Q_vals[i])
        durations.append(round(time.time()-start, 5))
        start=time.time()
        create.seq_fractions(seq_length, E=list_E_vals[i])
        durations.append(round(time.time()-start, 5))
        start=time.time()
        create.seq_fractions(seq_length, G=list_G_vals[i])
        durations.append(round(time.time()-start, 5))
        start=time.time()
        create.seq_fractions(seq_length, H=list_H_vals[i])
        durations.append(round(time.time()-start, 5))
        start=time.time()
        create.seq_fractions(seq_length, I=list_I_vals[i])
        durations.append(round(time.time()-start, 5))
        start=time.time()
        create.seq_fractions(seq_length, L=list_L_vals[i])
        durations.append(round(time.time()-start, 5))
        start=time.time()
        create.seq_fractions(seq_length, K=list_K_vals[i])
        durations.append(round(time.time()-start, 5))
        start=time.time()
        create.seq_fractions(seq_length, M=list_M_vals[i])
        durations.append(round(time.time()-start, 5))
        start=time.time()
        create.seq_fractions(seq_length, F=list_F_vals[i])
        durations.append(round(time.time()-start, 5))
        start=time.time()
        create.seq_fractions(seq_length, P=list_P_vals[i])
        durations.append(round(time.time()-start, 5))
        start=time.time()
        create.seq_fractions(seq_length, S=list_S_vals[i])
        durations.append(round(time.time()-start, 5))
        start=time.time()
        create.seq_fractions(seq_length, T=list_T_vals[i])
        durations.append(round(time.time()-start, 5))
        start=time.time()
        create.seq_fractions(seq_length, W=list_W_vals[i])
        durations.append(round(time.time()-start, 5))
        start=time.time()
        create.seq_fractions(seq_length, Y=list_Y_vals[i])
        durations.append(round(time.time()-start, 5))
        start=time.time()
        create.seq_fractions(seq_length, V=list_V_vals[i])
    
    return durations


def test_generation_by_re(seq_length=100, 
                          dimensions=32,
                          number_seq=100):
    '''
    Function to test sequence generation when specifying 
    dimensions and length. 

    Parameters
    ----------
    seq_length: int
        Length of the sequences to be generated
        default 100
    dimensions: float
        dimensions in angstroms for the sequence
        default 32
    number_seq: int
        Number of sequences to generate    
        default 1000

    Returns
    -------
    durations
        duration it took to make each sequence in seconds.
    '''

    # track durations per generation. This technically adds time but it's pretty negligible. 
    durations=[]

    for i in range(number_seq):
        # start timer
        start = time.time()
        create.seq_re(seq_length, dimensions, allowed_error=1)
        durations.append(round(time.time()-start, 5))
    return durations

def test_constant_class_var_generation(seq_length=100,
                                    num_vars=100):
    '''
    Function to test the generation of sequences with constant
    amino acids by class. 

    Parameters
    ----------
    seq_length: int
        Length of the sequences to be generated
        default 100

    num_vars: int
        Number of sequences to generate
        default 100

    Returns
    -------
    durations
        duration it took to make each sequence in seconds.
    '''
    # make a sequence with various classes of amino acids
    starter_seq=create.seq_fractions(seq_length, W=0.05,
                                     D=0.1, K=0.1, Q=0.1,
                                     A=0.1, P=0.05, H=0.05)

    # track durations per generation. This technically adds time but it's pretty negligible. 
    durations=[]

    for i in range(num_vars):
        # start timer
        start = time.time()
        create.constant_class_var(starter_seq)
        durations.append(round(time.time()-start, 5))
    return durations

def test_new_seq_constant_class_var_generation(seq_length=100,
                                                num_vars=100):
    '''
    Function to test the generation of the new_seq_constant_class_var(). 

    Parameters
    ----------
    seq_length: int
        Length of the sequences to be generated
        default 100

    num_vars: int
        Number of sequences to generate
        default 100

    Returns
    -------
    durations
        duration it took to make each sequence in seconds.
    '''
    # make a sequence with various classes of amino acids
    starter_seq=create.seq_fractions(seq_length, W=0.05,
                                     D=0.1, K=0.1, Q=0.1,
                                     A=0.1, P=0.05, H=0.05)


    # track durations per generation. This technically adds time but it's pretty negligible. 
    durations=[]

    for i in range(num_vars):
        # start timer
        start = time.time()
        create.new_seq_constant_class_var(starter_seq)
        durations.append(round(time.time()-start, 5))
    return durations


def test_constant_properties_var_generation(seq_length=100,
                                            num_vars=100):
    '''
    Function to test the generation of the constant_properties_var(). 

    Parameters
    ----------
    seq_length: int
        Length of the sequences to be generated
        default 100

    num_vars: int
        Number of sequences to generate
        default 100

    Returns
    -------
    durations
        duration it took to make each sequence in seconds.
    '''
    # make a sequence with various classes of amino acids
    starter_seq=create.sequence(seq_length, hydropathy=2.5, FCR=0.4, NCPR=0, kappa=0.3)

    # track durations per generation. This technically adds time but it's pretty negligible. 
    durations=[]
    
    for i in range(num_vars):
        # start timer
        start = time.time()
        create.constant_properties_var(starter_seq)
        durations.append(round(time.time()-start, 5))
    
    return durations


def test_constant_residue_var_generation(seq_length=100,
                                            num_vars=100):
    '''
    Function to test the generation of the constant_residue_var(). 

    Parameters
    ----------
    seq_length: int
        Length of the sequences to be generated
        default 100

    num_vars: int
        Number of sequences to generate
        default 100

    Returns
    -------
    durations
        duration it took to make each sequence in seconds.
    '''
    # make a sequence with various classes of amino acids
    starter_seq=create.seq_fractions(seq_length, W=0.05,
                                     D=0.1, K=0.1, Q=0.1,
                                     A=0.1, P=0.05, H=0.05)


    # track durations per generation. This technically adds time but it's pretty negligible. 
    durations=[]

    for i in range(num_vars):
        # start timer
        start = time.time()
        create.constant_residue_var(starter_seq, constant=['Q', 'W'])
        durations.append(round(time.time()-start, 5))

    return durations


def test_region_shuffle_var_generation(seq_length=100,
                                        num_vars=100):
    '''
    Function to test the generation of the region_shuffle_var(). 

    Parameters
    ----------
    seq_length: int
        Length of the sequences to be generated
        default 100

    num_vars: int
        Number of sequences to generate
        default 100

    Returns
    -------
    time
        Time in seconds it takes to generate all the sequences
    '''
    # make a sequence with various classes of amino acids
    starter_seq=create.sequence(seq_length)
    # list of residues in seq
    residues_start = [a for a in range(1, len(starter_seq)-20)]
    # grab random regions in the sequence at least 10 residues apart
    regions=[]
    for i in range(0, num_vars):
        start=random.choice(residues_start)
        end=start+random.randint(10, 20)
        regions.append([start, end])

    # track durations per generation. This technically adds time but it's pretty negligible. 
    durations=[]

    for i in range(num_vars):
        create.region_shuffle_var(starter_seq, shuffle=regions[i])
        # start timer
        start = time.time()    
        durations.append(round(time.time()-start, 5))
    return durations

def test_excluded_shuffle_var_generation(seq_length=100,
                                            num_vars=100):
    '''
    Function to test the generation of the excluded_shuffle_var(). 

    Parameters
    ----------
    seq_length: int
        Length of the sequences to be generated
        default 100

    num_vars: int
        Number of sequences to generate
        default 100

    Returns
    -------
    time
        Time in seconds it takes to generate all the sequences
    '''
    # make a sequence that has at least 5 different amino acids
    starter_seq=create.seq_fractions(seq_length, W=0.05, Q=0.1, S=0.1,
                                        E=0.1)
    # get amino acids in the sequence
    aas_in_sequence=list(set(starter_seq))
    # get amino acids to not_shuffle
    exclude_aas_list=[]
    for i in range(0, num_vars):
        t=[]
        for k in range(0, random.randint(1, len(aas_in_sequence)-1)):
            t.append(random.choice(aas_in_sequence))
        exclude_aas_list.append(list(set(t)))

    # track durations per generation. This technically adds time but it's pretty negligible. 
    durations=[]

    for i in range(num_vars):
        # start timer
        start = time.time()
        create.excluded_shuffle_var(starter_seq, exclude_aas=exclude_aas_list[i])
        durations.append(round(time.time()-start, 5))
    return durations


def test_targeted_shuffle_var_generation(seq_length=100,
                                        num_vars=100):
    '''
    Function to test the generation of the targeted_shuffle_var(). 

    Parameters
    ----------
    seq_length: int
        Length of the sequences to be generated
        default 100

    num_vars: int
        Number of sequences to generate
        default 100

    Returns
    -------
    durations
        duration it took to make each sequence in seconds.
    '''
    # make a sequence that has at least 5 different amino acids
    starter_seq=create.seq_fractions(seq_length, W=0.05, Q=0.1, S=0.1,
                                        E=0.1)
    # get amino acids in the sequence
    aas_in_sequence=list(set(starter_seq))
    # get amino acids to not_shuffle
    target_aas_list=[]
    for i in range(0, num_vars):
        t=[]
        for k in range(0, random.randint(2, len(aas_in_sequence)-1)):
            t.append(random.choice(aas_in_sequence))
        target_aas_list.append(list(set(t)))

    # track durations per generation. This technically adds time but it's pretty negligible. 
    durations=[]

    for i in range(num_vars):
        # start timer
        start = time.time()
        create.targeted_shuffle_var(starter_seq, target_aas=target_aas_list[i])
        durations.append(round(time.time()-start, 5))
    return durations

def test_asymmetry_var_generation(seq_length=100,
                                num_vars=100):
    '''
    Function to test the generation of the targeted_shuffle_var(). 

    Parameters
    ----------
    seq_length: int
        Length of the sequences to be generated
        default 100

    num_vars: int
        Number of sequences to generate
        default 100


    Returns
    -------
    durations
        duration it took to make each sequence in seconds.
    '''
    # make a sequence that has at least all of the different amino acids we can change asymmetry of
    starter_seq=create.seq_fractions(seq_length, W=0.05, Q=0.1, S=0.1,
                                        E=0.1, A=0.1, P=0.1, K=0.1)
    # get amino acids in the sequence
    aas_in_sequence=list(set(starter_seq))
    # possible cluasses to change
    possible_classes = ['negative', 'positive', 'proline', 'aromatic', 'aliphatic', 'polar', 'custom']
    # get classes to change asymmetry of, also make a list of if increasing or decreasing and of the
    # number of changers to make in the sequence
    class_list=[]
    increase_decrease_list=[]
    num_changes_list=[]
    for i in range(0, num_vars):
        num_changes_list.append(random.randint(1, int(seq_length/2)))
        cur_class=random.choice(possible_classes)
        if cur_class=='custom':
            t=[]
            for k in range(0, random.randint(1, len(aas_in_sequence)-1)):
                t.append(random.choice(aas_in_sequence))
            class_list.append(t)
        else:
            class_list.append(cur_class)
        increase_decrease_list.append(random.choice(['increase', 'decrease']))

    # track durations per generation. This technically adds time but it's pretty negligible. 
    durations=[]     

    for i in range(num_vars):
        # start timer
        start = time.time()
        create.asymmetry_var(starter_seq, increase_decrease_list[i], class_list[i], num_changes_list[i])
        durations.append(round(time.time()-start, 5))

    return durations



def test_hydro_class_var_generation(seq_length=100,
                                        num_vars=100):
    '''
    Function to test the generation of the hydro_class_var(). 

    Parameters
    ----------
    seq_length: int
        Length of the sequences to be generated
        default 100

    num_vars: int
        Number of sequences to generate
        default 100

    Returns
    -------
    durations
        duration it took to make each sequence in seconds.
    '''
    # make a sequence we should be able to modulate the hydropathy of
    starter_seq=create.seq_fractions(seq_length, W=0.05, Q=0.1, S=0.1,
                                        E=0.1, A=0.1, K=0.1)

    # get the min and max hydro we can make
    hydro_range=hydro_range_constant_classes(starter_seq)

    # get min, max hydro
    min_hydro=hydro_range[0]
    max_hydro=hydro_range[1]
    # make list of hydro ranges separated by 0.1
    hydro_ranges=np.arange(min_hydro+0.1, max_hydro-0.1, 0.1)
    hydro_change_list=[]

    # choose final hydro value
    for i in range(0, num_vars):
        hydro_change_list.append(random.choice(hydro_ranges))


    # track durations per generation. This technically adds time but it's pretty negligible. 
    durations=[] 

    # make seqs
    for i in range(num_vars):
        start = time.time()
        create.hydro_class_var(starter_seq, hydropathy=hydro_change_list[i])
        durations.append(round(time.time()-start, 5))
    return durations


def test_fcr_class_var_generation(seq_length=100,
                                num_vars=100):
    '''
    Function to test the generation of the fcr_class(). 

    Parameters
    ----------
    seq_length: int
        Length of the sequences to be generated
        default 100

    num_vars: int
        Number of sequences to generate
        default 100

    Returns
    -------
    durations
        duration it took to make each sequence in seconds.
    '''
    # make a sequence we should be able to modulate the FCR of
    starter_seq=create.seq_fractions(seq_length, N=0.1, Q=0.1, S=0.1, T=0.1, W=0.05,
                                        E=0.1, K=0.1, R=0.1, D=0.1)

    # get ranges of fractions to change FCR to
    FCR_frac_interval=round(1/seq_length, 4)

    # make list of FCR vals separted by FCR_frac_interval
    # use between 0.1 and 0.6 because otherwise kappa will be hard to hold constant.
    FCR_vals_possible=np.arange(0.1, 0.6, FCR_frac_interval)
    FCR_changes=[]

    # choose final hydro value
    for i in range(0, num_vars):
        FCR_changes.append(random.choice(FCR_vals_possible))

    # track durations per generation. This technically adds time but it's pretty negligible. 
    durations=[] 

    # start timer
    for i in range(num_vars):
        start = time.time()
        create.fcr_class_var(starter_seq, FCR=FCR_changes[i])
        durations.append(round(time.time()-start, 5))
    
    return durations


def test_ncpr_class_var_generation(seq_length=100,
                                num_vars=100):
    '''
    Function to test the generation of the ncpr_class_var(). 

    Parameters
    ----------
    seq_length: int
        Length of the sequences to be generated
        default 100

    num_vars: int
        Number of sequences to generate
        default 100        

    Returns
    -------
    durations
        duration it took to make each sequence in seconds.
    '''
    # make a sequence we should be able to modulate the FCR of
    starter_seq=create.seq_fractions(seq_length, N=0.1, Q=0.1, S=0.1, T=0.1, W=0.05,
                                        E=0.1, K=0.1, R=0.1, D=0.1)

    # get ranges of fractions to change FCR to
    NCPR_frac_interval=round(1/seq_length, 4)*2

    # make list of NCPR vals separted by NCPR_frac_interval
    # use between -0.3 and 0.3 because otherwise kappa will be hard to hold constant.
    NCPR_vals_possible=np.arange(-0.3, 0.3, NCPR_frac_interval)
    NCPR_changes=[]

    # choose final ncpr value
    for i in range(0, num_vars):
        NCPR_changes.append(random.choice(NCPR_vals_possible))

    # track durations per generation. This technically adds time but it's pretty negligible. 
    durations=[]    

    # start timer
    for i in range(num_vars):
        start = time.time()
        create.ncpr_class_var(starter_seq, NCPR=NCPR_changes[i])
        durations.append(round(time.time()-start, 5))
    return durations


def test_kappa_var_generation(seq_length=100,
                                num_vars=100):
    '''
    Function to test the generation of the kappa_var(). 

    Parameters
    ----------
    seq_length: int
        Length of the sequences to be generated
        default 100

    num_vars: int
        Number of sequences to generate
        default 100

    Returns
    -------
    durations
        duration it took to make each sequence in seconds.
    '''
    # make a sequence we should be able to modulate the FCR of
    starter_seq=create.seq_fractions(seq_length, N=0.1, Q=0.1, S=0.1, T=0.1, W=0.05,
                                        E=0.1, K=0.1, R=0.1, D=0.1)

    # possible kappa values, separated by 0.01
    possible_kappa_vals=np.arange(0.06, 0.95, 0.01)
    kappa_changes=[]

    # choose final kappa value
    for i in range(0, num_vars):
        kappa_changes.append(random.choice(possible_kappa_vals))

    # track durations per generation. This technically adds time but it's pretty negligible. 
    durations=[]

    # start timer
    for i in range(num_vars):
        start = time.time()
        create.kappa_var(starter_seq, kappa=kappa_changes[i])
        durations.append(round(time.time()-start, 5))
    return durations

def test_re_var_generation(seq_length=100,
                            num_vars=100):
    '''
    Function to test the generation of the re_var(). 

    Parameters
    ----------
    seq_length: int
        Length of the sequences to be generated
        default 100

    num_vars: int
        Number of sequences to generate
        default 100

    Returns
    -------
    durations
        duration it took to make each sequence in seconds.
    '''
    # make a sequence we should be able to modulate the FCR of
    starter_seq=create.seq_fractions(seq_length, N=0.1, Q=0.1, S=0.1, T=0.1, W=0.05,
                                        E=0.1, K=0.1, R=0.1, D=0.1)
    
    # choose to increase or decrease
    increase_decrease_list=[]
    for i in range(0, num_vars):
        increase_decrease_list.append(random.choice(['increase', 'decrease']))
    
    # track durations per generation. This technically adds time but it's pretty negligible. 
    durations=[]

    # make the sequences
    for i in range(num_vars):
        start = time.time()
        create.re_var(starter_seq, increase_decrease_list[i])
        durations.append(round(time.time()-start, 5))
    return durations



def test_rg_var_generation(seq_length=100,
                            num_vars=100):
    '''
    Function to test the generation of the rg_var(). 

    Parameters
    ----------
    seq_length: int
        Length of the sequences to be generated
        default 100

    num_vars: int
        Number of sequences to generate
        default 100

    Returns
    -------
    durations
        time it took to make each sequence
    '''
    # make a sequence we should be able to modulate the FCR of
    starter_seq=create.seq_fractions(seq_length, Y=0.04, G=0.1, F=0.03, Q=0.1, 
                                    N=0.1, D=0.1, E=0.1, K=0.1, S=0.1, T=0.1)
    
    # choose to increase or decrease
    increase_decrease_list=[]
    for i in range(0, num_vars):
        increase_decrease_list.append(random.choice(['increase', 'decrease']))

    # track durations per generation. This technically adds time but it's pretty negligible. 
    durations=[]

    # make the sequences
    for i in range(num_vars):
        start = time.time()
        create.rg_var(starter_seq, increase_decrease_list[i])
        durations.append(round(time.time()-start, 5))
    return durations

def test_re_var_generation(seq_length=100,
                            num_vars=100):
    '''
    Function to test the generation of the re_var(). 

    Parameters
    ----------
    seq_length: int
        Length of the sequences to be generated
        default 100

    num_vars: int
        Number of sequences to generate
        default 100

    Returns
    -------
    durations
        duration it took to make each sequence in seconds. 
    '''
    # make a sequence we should be able to modulate the FCR of
    starter_seq=create.seq_fractions(seq_length, Y=0.04, G=0.1, F=0.03, Q=0.1, 
                                    N=0.1, D=0.1, E=0.1, K=0.1, S=0.1, T=0.1)
    
    # choose to increase or decrease
    increase_decrease_list=[]
    for i in range(0, num_vars):
        increase_decrease_list.append(random.choice(['increase', 'decrease']))

    # track durations per generation. This technically adds time but it's pretty negligible. 
    durations=[]

    # make the sequences
    for i in range(num_vars):
        start = time.time()
        create.re_var(starter_seq, increase_decrease_list[i])
        durations.append(round(time.time()-start, 5))
    return durations


def get_goose_performance_numbers(seqs_per_function=500,
    print_function_progress=True):
    '''
    Function that runs and tests everything.
    Returns a dictionary for times for everything that was tested in seconds. 
    
    Parameters
    ----------
    seqs_per_function: int
        Number of sequences to generate per function. 
        default 500

    print_function_progress : bool
        whether to print which function i sbeing tested, default=True

    Returns
    -------
    times
        dictionary of times for each function tested. 
    '''
    function_performance_dict={}

    if print_function_progress:
        print('\nTesting sequence generation by specifying only length')
    function_performance_dict['length']=test_goose_just_length(num_seqs=seqs_per_function)

    if print_function_progress:
        print('\nTesting sequence generation by specifying a single property')    
    function_performance_dict['single_property']=test_goose_single_property_generation(seq_per_prop=int(seqs_per_function/4))

    if print_function_progress:
        print('\nTesting sequence generation by specifying a single fraction')    
    function_performance_dict['fractions']=test_generation_by_fractions(seq_per_frac=int(seqs_per_function/20))

    if print_function_progress:
        print('\nTesting sequence generation by specifying Re')    
    function_performance_dict['Re']=test_generation_by_re(number_seq=seqs_per_function)

    if print_function_progress:
        print('\nTesting sequence generation by making constant class vars')
    function_performance_dict['constant_class_vars']=test_constant_class_var_generation(num_vars=seqs_per_function)

    if print_function_progress:
        print('\nTesting sequence generation by making new constant class vars')
    function_performance_dict['new_seq_constant_class_vars']=test_new_seq_constant_class_var_generation(num_vars=seqs_per_function)

    if print_function_progress:
        print('\nTesting sequence generation by making constant properties vars')
    function_performance_dict['constant_properties_vars']=test_constant_properties_var_generation(num_vars=seqs_per_function)

    if print_function_progress:
        print('\nTesting sequence generation by making constant residue vars')
    function_performance_dict['constant_residue_vars']=test_constant_residue_var_generation(num_vars=seqs_per_function)

    if print_function_progress:
        print('\nTesting sequence generation by making region shuffle vars')
    function_performance_dict['region_shuffle_vars']=test_region_shuffle_var_generation(num_vars=seqs_per_function)

    if print_function_progress:
        print('\nTesting sequence generation by making excluded shuffle vars')
    function_performance_dict['excluded_shuffle_vars']=test_excluded_shuffle_var_generation(num_vars=seqs_per_function)

    if print_function_progress:
        print('\nTesting sequence generation by making targeted shuffle vars')
    function_performance_dict['targeted_shuffle_vars']=test_targeted_shuffle_var_generation(num_vars=seqs_per_function)

    if print_function_progress:
        print('\nTesting sequence generation by making asymmetry vars')
    function_performance_dict['asymmetry_vars']=test_asymmetry_var_generation(num_vars=seqs_per_function)

    if print_function_progress:
        print('\nTesting sequence generation by making hydro class vars')
    function_performance_dict['hydro_class_vars']=test_hydro_class_var_generation(num_vars=seqs_per_function)

    if print_function_progress:
        print('\nTesting sequence generation by making fcr class vars')
    function_performance_dict['fcr_class_vars']=test_fcr_class_var_generation(num_vars=seqs_per_function)

    if print_function_progress:
        print('\nTesting sequence generation by making ncpr class vars')
    function_performance_dict['ncpr_class_vars']=test_ncpr_class_var_generation(num_vars=seqs_per_function)

    if print_function_progress:
        print('\nTesting sequence generation by making kappa vars')
    function_performance_dict['kappa_vars']=test_kappa_var_generation(num_vars=seqs_per_function)

    if print_function_progress:
        print('\nTesting sequence generation by making re vars')
    function_performance_dict['re_vars']=test_re_var_generation(num_vars=seqs_per_function)

    if print_function_progress:
        print('\nTesting sequence generation by making rg vars')
    function_performance_dict['rg_vars']=test_rg_var_generation(num_vars=seqs_per_function)

    return function_performance_dict

def gaussian_random_between(min_val, max_val, mean, std_dev):
    """
    Generates a random number based on a Gaussian distribution between min_val and max_val.
    
    Parameters:
    - min_val: The minimum value of the range.
    - max_val: The maximum value of the range.
    - mean: The mean (mu) of the Gaussian distribution.
    - std_dev: The standard deviation (sigma) of the Gaussian distribution.
    
    Returns:
    - A random number between min_val and max_val based on the Gaussian distribution.
    """
    # Ensure that the mean is within the specified range
    if not (min_val <= mean <= max_val):
        raise ValueError("Mean must be within the range [min_val, max_val]")
    
    # Generate a random number based on the Gaussian distribution
    while True:
        rand_num = np.random.normal(mean, std_dev)
        if min_val <= rand_num <= max_val:
            return rand_num

