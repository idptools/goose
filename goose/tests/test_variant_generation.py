'''
Code to test the variant generation functionality. 
'''

import random

from goose import create
from goose.backend.protein import Protein as pr
from goose.backend import parameters

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

def test_basic(verbose=False):
    '''
    very basic test for each function just to make sure it works at all.
    '''
    # using a fairly flexible sequence. 
    start_seq='KRQRRQDRPRRSNRNNSDQPNRKGSGQDPSDETPGKDKNDDNEPNEQKRRGQESQSRRDGRRDDDNTGEEKGDEDKRETPPSQTDNNQKQQKEDDSENRQKPGGNDKKKGDNRNNPPNNTTQEPQGRKREGSGGPTEKGKGDQEGDTKTNQKSKPQREPPEPGPREKEKEPNQKDEDEQENETDDRKPRKEGPEQKRGDR'
    start_seq_props={'length': 200, 'FCR': 0.5, 'NCPR': 0.0, 'hydropathy': 1.6, 'kappa': 0.155, 'fractions': {'A': 0.0, 'C': 0.0, 'D': 0.125, 'E': 0.125, 'F': 0.0, 'G': 0.1, 'H': 0.0, 'I': 0.0, 'K': 0.125, 'L': 0.0, 'M': 0.0, 'N': 0.1, 'P': 0.1, 'Q': 0.1, 'R': 0.125, 'S': 0.05, 'T': 0.05, 'V': 0.0, 'W': 0.0, 'Y': 0.0}, 'helical regions': [], 'predicted phosphosites': {'S': [24, 29, 53, 55, 81, 95, 131], 'T': [32, 66, 78, 83, 119, 120, 135, 146, 148, 182], 'Y': []}, 'predicted cellular localization': {'mitochondria': 'No mitochondrial targeting sequences predicted.', 'NES': 'No NES sequences predicted.', 'NLS': {'KRQRRQ': [1, 7]}}, 'predicted transcriptional activation': {'Predicted TADs': 'No transcriptional activation domains predicted.'}, 'predicted polymer properties': {'Rg': 43.6639, 'Re': 100.4427}, 'fraction aromatic': 0.0, 'fraction polar': 0.3, 'fraction aliphatic': 0.0}
    if verbose:
        print('on min_var')
    min_var=create.minimal_var(start_seq, hydropathy = 1.5, FCR = 0.5, NCPR = 0.2)
    if verbose:
        print('on new_seq_constant_class_var')
    new_seq_constant_class_var=create.new_seq_constant_class_var(start_seq, attempts=5)
    if verbose:
        print('on constant_class_var')     
    constant_class_var=create.constant_class_var(start_seq, attempts=5)
    if verbose:
        print('on constant_properties_var')    
    constant_prop_var=create.constant_properties_var(start_seq, attempts=5)
    if verbose:
        print('on hydro_class_var')    
    hydro_class_var=create.hydro_class_var(start_seq, hydropathy=2)
    if verbose:
        print('on constant_residue_var')    
    constant_res_var=create.constant_residue_var(start_seq, constant=['R'])
    if verbose:
        print('on region_shuffle_var')    
    region_shuffle_var=create.region_shuffle_var(start_seq, shuffle=[50,100], attempts=5)
    if verbose:
        print('on kappa_var')    
    kappa_var=create.kappa_var(start_seq, 0.3)
    if verbose:
        print('on asymmetry_var')    
    asymmetry_var=create.asymmetry_var(start_seq, 'increase', ['K', 'R', 'D', 'E'], number_changes=None)
    if verbose:
        print('on fcr_class_var')    
    fcr_class_var=create.fcr_class_var(start_seq, 0.4, attempts=10)
    if verbose:
        print('on ncpr_class_var')    
    ncpr_class_var=create.ncpr_class_var(start_seq, 0.1, attempts=10)
    if verbose:
        print('on all_props_class_var')    
    all_props_class_var=create.all_props_class_var(start_seq, hydropathy=2, FCR=0.4, NCPR=-0.1, kappa=0.3)
    if verbose:
        print('on targeted_shuffle_var')    
    target_shuff_var=create.targeted_shuffle_var(start_seq, ['P', 'D'], attempts=10)
    if verbose:
        print('on excluded_shuffle_var')    
    exclude_shuff_var=create.excluded_shuffle_var(start_seq, ['P', 'D'], attempts=10)
    if verbose:
        print('on re_var')    
    re_var=create.re_var(start_seq, 'increase', return_all=False, return_all_interval=0.2)
    if verbose:
        print('on rg_var')    
    rg_var=create.rg_var(start_seq, 'decrease', return_all=False, return_all_interval=0.2)
    if verbose:
        print('on Done!')    

