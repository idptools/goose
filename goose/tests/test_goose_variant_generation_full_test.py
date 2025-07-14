"""
A much more thorough test of the variant generation features. 
"""
import pytest
import random
import numpy as np
from sparrow.protein import Protein
from goose import create
from goose import goose_exceptions
from goose.backend import parameters
from goose.backend_variant_generation.helper_functions import find_hydro_range_constant_class
from sparrow.protein import Protein
from goose.data.defined_aa_classes import aa_classes
from goose.backend_variant_generation.verify_variant_functions import (
    verify_same_number_by_class,
    verify_same_order_by_class,
    verify_constant_properties,
    verify_same_residue_count,
    verify_constant_hydropathy,
    verify_constant_kappa,
    verify_constant_FCR,
    verify_constant_NCPR,
    verify_target_hydropathy,
    verify_target_kappa,
    verify_target_FCR,
    verify_target_NCPR,
    verify_constant_region,
    verify_region_changed,
    verify_constant_residue_positions,
    verify_disorder, 
    verify_changed_iwd,
    verify_dimensions,
    verify_same_length)

'''
For each of the functions, we need to verify specific things, going to list that here:
- ``'shuffle_specific_regions'`` - verify_constant_region, verify_region_changed, verify_residue_count
- ``'shuffle_except_specific_regions'`` - verify_constant_region, verify_region_changed, verify_residue_count
- ``'shuffle_specific_residues'`` - verify_constant_residue_positions, verify_region_changed, verify_residue_count
- ``'shuffle_except_specific_residues'`` - verify_constant_residue_positions, verify_region_changed, verify_residue_count
- ``'weighted_shuffle_specific_residues'`` - verify_constant_residue_positions, verify_region_changed, verify_residue_count
- ``'targeted_reposition_specific_residues'`` - verify_constant_residue_positions, verify_region_changed, verify_residue_count
- ``'change_residue_asymmetry'`` - verify_changed_iwd, verify_residue_count
- ``'constant_properties'`` - verify_constant_properties
- ``'constant_residues_and_properties'`` - verify_constant_properties, verify_constant_residue_positions
- ``'constant_properties_and_class'`` - verify_constant_properties, verify_same_number_by_class
- ``'constant_properties_and_class_by_order'`` -  verify_constant_properties, verify_same_number_by_class, verify_same_order_by_class
- ``'change_hydropathy_constant_class'`` - verify_target_hydropathy, verify_same_number_by_class
- ``'change_fcr_minimize_class_changes'`` - verify_target_FCR
- ``'change_ncpr_constant_class'`` - verify_target_NCPR
- ``'change_kappa'`` - verify_target_kappa, verify_residue_count
- ``'change_properties_minimize_differences'`` - verify_target_hydropathy, verify_target_kappa, verify_target_FCR, verify_target_NCPR
    Note: only verify things that are changed, otherwise make sure they are within tolerance
- ``'change_any_properties'`` -verify_target_hydropathy, verify_target_kappa, verify_target_FCR, verify_target_NCPR
    Note: only verify things that are changed, otherwise make sure they are within tolerance
- ``'change_dimensions'`` - verify_dimensions

Note: For all of them, we need to use verify_disorder and verify_same_length
'''

# =============================================================================
# Test sequences - create a variety of test sequences with different properties
@pytest.fixture
def test_sequences():
    """Create test sequences with different properties for variant testing."""
    return {
        'shuffle_test':'GSGSGSGSGS'+create.sequence(90, FCR=0.4, hydropathy=2.3),
        'no_charge': create.sequence(100, FCR=0.0),
        'only_negative': create.sequence(100, FCR=0.2, NCPR=-0.2),
        'only_positive': create.sequence(100, FCR=0.2, NCPR=0.2),
        'balanced_charge': create.sequence(100, FCR=0.2, NCPR=0.0),
        'only_charged': create.sequence(100, FCR=1.0, NCPR=0),
        'high_kappa': create.sequence(100, FCR=0.2, NCPR=0.0, hydropathy=2.0, kappa=0.8),
        'low_kappa': create.sequence(100, FCR=0.2, NCPR=0.0, hydropathy=2.0, kappa=0.2),
        'hydro_3': create.sequence(100, FCR=0.2, NCPR=0.0, hydropathy=3.0)
    }

# =============================================================================
# SHUFFLING VARIANT TESTS
# =============================================================================

def test_shuffle_specific_regions(test_sequences):
    """Test shuffle_specific_regions variant"""
    seq = test_sequences['shuffle_test']
    target_regions= (0, 40)
    constant_regions = [41, 100]
    variant_seq = create.variant(seq, 'shuffle_specific_regions',
                                 shuffle_regions=[target_regions])
    assert verify_constant_region(seq, variant_seq, constant_regions)
    assert verify_region_changed(seq, variant_seq, target_regions)
    assert verify_same_residue_count(seq, variant_seq)
    assert verify_disorder(seq, variant_seq)
    assert verify_same_length(seq, variant_seq)
def test_shuffle_except_specific_regions(test_sequences):
    """Test shuffle_except_specific_regions variant."""
    seq = test_sequences['shuffle_test']
    excluded_regions= (0, 40)
    shuffled_regions = [41, 100]
    variant_seq = create.variant(seq, 'shuffle_except_specific_regions',
                                 excluded_regions=[excluded_regions])
    assert verify_constant_region(seq, variant_seq, list(excluded_regions))
    assert verify_region_changed(seq, variant_seq, shuffled_regions)
    assert verify_same_residue_count(seq, variant_seq)
    assert verify_disorder(seq, variant_seq)
    assert verify_same_length(seq, variant_seq)

def test_shuffle_specific_residues(test_sequences):
    """Test shuffle_specific_residues variant."""
    seq = test_sequences['shuffle_test']
    target_residues = ['G', 'S']
    constant_residues = [a for a in seq if a not in target_residues]
    variant_seq = create.variant(seq, 'shuffle_specific_residues',
                                 target_residues=target_residues)
    assert verify_constant_residue_positions(seq, variant_seq, constant_residues)
    assert verify_region_changed(seq, variant_seq, [0, 10])
    assert verify_same_residue_count(seq, variant_seq)
    assert verify_disorder(seq, variant_seq)
    assert verify_same_length(seq, variant_seq)

def test_shuffle_except_specific_residues(test_sequences):
    """Test shuffle_except_specific_residues variant."""
    seq = test_sequences['shuffle_test']
    excluded_residues = ['G', 'S']
    variant_seq = create.variant(seq, 'shuffle_except_specific_residues',
                                 excluded_residues=excluded_residues)
    assert verify_constant_residue_positions(seq, variant_seq, excluded_residues)
    assert verify_region_changed(seq, variant_seq, [10, len(seq)])
    assert verify_same_residue_count(seq, variant_seq)
    assert verify_disorder(seq, variant_seq)
    assert verify_same_length(seq, variant_seq)

def test_weighted_shuffle_specific_residues(test_sequences):
    """Test shuffle_specific_residues variant."""
    seq = test_sequences['shuffle_test']
    target_residues = ['G', 'S']
    constant_residues = [a for a in seq if a not in target_residues]
    variant_seq = create.variant(seq, 'weighted_shuffle_specific_residues',
                                 target_residues=target_residues,
                                 shuffle_weight=1)
    assert verify_constant_residue_positions(seq, variant_seq, constant_residues)
    assert verify_region_changed(seq, variant_seq, [0, 10])
    assert verify_same_residue_count(seq, variant_seq)
    assert verify_disorder(seq, variant_seq)
    assert verify_same_length(seq, variant_seq)

def test_targeted_reposition_specific_residues(test_sequences):
    """Test targeted_reposition_specific_residues variant."""
    seq = test_sequences['shuffle_test']
    target_residues = ['G', 'S']
    variant_seq = create.variant(seq, 'targeted_reposition_specific_residues',
                                 target_residues=target_residues)
    assert verify_region_changed(seq, variant_seq, [0, 10])
    assert verify_same_residue_count(seq, variant_seq)
    assert verify_disorder(seq, variant_seq)
    assert verify_same_length(seq, variant_seq)

def test_change_residue_asymmetry(test_sequences):
    """Test change_residue_asymmetry variant."""
    seq = test_sequences['shuffle_test']
    seq = list(seq)
    random.shuffle(seq)
    seq = ''.join(seq)
    targets = ['G', 'S']
    variant_seq = create.variant(seq, 'change_residue_asymmetry', target_residues=targets,
                                 increase_or_decrease='decrease',
                                 num_changes=2)
    assert verify_same_residue_count(seq, variant_seq)
    assert verify_changed_iwd(seq, variant_seq,
                              'decrease', residues=targets)
    assert verify_disorder(seq, variant_seq)
    assert verify_same_length(seq, variant_seq)

    # now take the sequence with decreased asymmetry and increase it
    seq=variant_seq
    variant_seq = create.variant(seq, 'change_residue_asymmetry', target_residues=targets,
                                 increase_or_decrease='increase',
                                 num_changes=2)
    assert verify_same_residue_count(seq, variant_seq)
    assert verify_changed_iwd(seq, variant_seq,
                              'increase', residues=targets)
    assert verify_disorder(seq, variant_seq)
    assert verify_same_length(seq, variant_seq)


# =============================================================================
# CONSTANT PROPERTIES VARIANT TESTS
# =============================================================================
def test_constant_properties(test_sequences):
    seqs = [test_sequences['only_negative'],
            test_sequences['only_positive'],
            test_sequences['balanced_charge'],
            test_sequences['no_charge'],
            test_sequences['high_kappa'],
            test_sequences['low_kappa'],
            test_sequences['hydro_3']]
    for seq in seqs:
        variant_seq = create.variant(seq, 'constant_properties')
        assert verify_constant_properties(seq, variant_seq)
        assert verify_disorder(seq, variant_seq)
        assert verify_same_length(seq, variant_seq)

def test_constant_residues_and_properties(test_sequences):
    seqs = [test_sequences['only_negative'],
            test_sequences['only_positive'],
            test_sequences['balanced_charge'],
            test_sequences['no_charge'],
            test_sequences['high_kappa'],
            test_sequences['low_kappa'],
            test_sequences['hydro_3']]
    for seq in seqs:
        # choose a residue to keep constant
        residue_options = [a for a in seq]
        constant_residue = random.choice(residue_options)
        variant_seq = create.variant(seq, 'constant_residues_and_properties',
                                      constant_residues=[constant_residue])
        assert verify_constant_properties(seq, variant_seq)
        assert verify_constant_residue_positions(seq, variant_seq, residues=[constant_residue])
        assert verify_disorder(seq, variant_seq)
        assert verify_same_length(seq, variant_seq)

def test_constant_residues_and_class(test_sequences):
    """Test constant_properties_and_class variant."""
    seqs = [test_sequences['only_negative'],
            test_sequences['only_positive'],
            test_sequences['balanced_charge'],
            test_sequences['no_charge'],
            test_sequences['high_kappa'],
            test_sequences['low_kappa'],
            test_sequences['hydro_3']]
    for seq in seqs:
        # make variant
        variant_seq = create.variant(seq, 'constant_properties_and_class')
        assert verify_constant_properties(seq, variant_seq)
        assert verify_same_number_by_class(seq, variant_seq)
        assert verify_disorder(seq, variant_seq)
        assert verify_same_length(seq, variant_seq)

def test_constant_properties_and_class_by_order(test_sequences):
    """Test constant_properties_and_class variant."""
    seqs = [test_sequences['only_negative'],
            test_sequences['only_positive'],
            test_sequences['balanced_charge'],
            test_sequences['no_charge'],
            test_sequences['high_kappa'],
            test_sequences['low_kappa'],
            test_sequences['hydro_3']]
    for seq in seqs:
        # make variant
        variant_seq = create.variant(seq, 'constant_properties_and_class_by_order')
        assert verify_constant_properties(seq, variant_seq)
        assert verify_same_number_by_class(seq, variant_seq)
        assert verify_same_order_by_class(seq, variant_seq)
        assert verify_disorder(seq, variant_seq)
        assert verify_same_length(seq, variant_seq)

def test_change_hydropathy_constant_class(test_sequences):
    """Test change_hydropathy_constant_class variant."""
    seqs = [test_sequences['only_negative'],
            test_sequences['only_positive'],
            test_sequences['balanced_charge'],
            test_sequences['no_charge'],
            test_sequences['hydro_3']]
    for seq in seqs:
        # figure out hydro
        min_hydro, max_hydro = find_hydro_range_constant_class(seq)
        # get current hydro
        cur_hydro = Protein(seq).hydrophobicity
        # get two target hydro
        target_low = (cur_hydro+min_hydro)/2
        target_high = (cur_hydro+max_hydro)/2
        for target_hydro in [target_low, target_high]:
            # make variant
            variant_seq = create.variant(seq, 'change_hydropathy_constant_class',
                                         target_hydropathy=target_hydro)
            assert verify_target_hydropathy(variant_seq, target_hydro)
            assert verify_same_number_by_class(seq, variant_seq)
            assert verify_disorder(seq, variant_seq)
            assert verify_same_length(seq, variant_seq)

def test_change_fcr_minimize_class_changes(test_sequences):
    """Test change_fcr_minimize_class_changes variant."""
    seqs = [test_sequences['only_negative'],
            test_sequences['only_positive'],
            test_sequences['balanced_charge'],
            test_sequences['no_charge'],
            test_sequences['hydro_3']]
    for seq in seqs:
        # set target
        cur_FCR = Protein(seq).FCR
        target_1 = cur_FCR + random.uniform(0.05, 0.15)
        target_2 = cur_FCR - random.uniform(0.05, 0.15)
        target_1 = round(target_1 * len(seq)) / len(seq)  # round to nearest residue fraction
        target_2 = round(target_2 * len(seq)) / len(seq)
        if target_1 < 1:
            # make variant
            variant_seq = create.variant(seq, 'change_fcr_minimize_class_changes',
                                        target_FCR=target_1)
            assert verify_target_FCR(variant_seq, target_1)
            assert verify_disorder(seq, variant_seq)
            assert verify_same_length(seq, variant_seq)
        
        if target_2 > 0:
            # make variant
            variant_seq = create.variant(seq, 'change_fcr_minimize_class_changes',
                                        target_FCR=target_2)
            assert verify_target_FCR(variant_seq, target_2)
            assert verify_disorder(seq, variant_seq)
            assert verify_same_length(seq, variant_seq)


def test_change_ncpr_constant_class(test_sequences):
    """Test change_ncpr_constant_class variant."""
    seqs = [test_sequences['only_negative'],
            test_sequences['only_positive'],
            test_sequences['balanced_charge'],
            test_sequences['high_kappa'],
            test_sequences['low_kappa'],
            test_sequences['hydro_3']]
    for seq in seqs:
        # get current FCR
        cur_FCR = Protein(seq).FCR
        # set target
        target_NCPR = random.uniform(-cur_FCR, cur_FCR)
        target_NCPR = round(target_NCPR * len(seq)) / len(seq)  # round to nearest residue fraction
        # make sure we don't set NCPR greater than FCR
        if target_NCPR > cur_FCR:
            target_NCPR = cur_FCR * np.sign(target_NCPR)
        # make variant
        variant_seq = create.variant(seq, 'change_ncpr_constant_class',
                                     target_NCPR=target_NCPR)
        assert verify_target_NCPR(variant_seq, target_NCPR)
        assert (Protein(variant_seq).FCR - cur_FCR) < (1/len(seq)) + 1e-5
        assert verify_disorder(seq, variant_seq)
        assert verify_same_length(seq, variant_seq)


def test_change_kappa(test_sequences):
    """Test change_kappa variant."""
    seqs = [test_sequences['balanced_charge'],
            test_sequences['high_kappa'],
            test_sequences['low_kappa']]
    for seq in seqs:
        # set target
        target_kappa = random.uniform(0.1, 0.9)
        # make variant
        variant_seq = create.variant(seq, 'change_kappa',
                                     target_kappa=target_kappa)
        assert verify_target_kappa(variant_seq, target_kappa)
        assert verify_disorder(seq, variant_seq)
        assert verify_same_length(seq, variant_seq)

def test_change_properties_minimize_differences(test_sequences):
    """Test change_properties_minimize_differences variant."""
    seqs = [test_sequences['only_negative'],
            test_sequences['only_positive'],
            test_sequences['balanced_charge'],
            test_sequences['hydro_3']]
    for seq in seqs:
        # set targets
        target_hydro = random.uniform(2, 4.0)
        target_FCR = random.uniform(0.1, 0.3)
        target_FCR = round(target_FCR * len(seq)) / len(seq)  # round to nearest residue fraction
        target_NCPR = random.uniform(-target_FCR, target_FCR)
        target_NCPR = round(target_NCPR * len(seq)) / len(seq)  # round to nearest residue fraction
        # only set kappa if we can actually change it.
        if target_FCR > 0 and target_NCPR != target_FCR:
            target_kappa = random.uniform(0.1, 0.5)
            check_kappa=True
        else:
            target_kappa=None
            check_kappa=False
        # make variant
        variant_seq = create.variant(seq, 'change_properties_minimize_differences',
                                     target_hydropathy=target_hydro,
                                     target_FCR=target_FCR,
                                     target_NCPR=target_NCPR,
                                     target_kappa=target_kappa)
        assert verify_target_hydropathy(variant_seq, target_hydro)
        assert verify_target_FCR(variant_seq, target_FCR)
        assert verify_target_NCPR(variant_seq, target_NCPR)
        if check_kappa:
            assert verify_target_kappa(variant_seq, target_kappa)
        assert verify_disorder(seq, variant_seq)
        assert verify_same_length(seq, variant_seq)

def test_change_any_properties(test_sequences):
    """Test change_any_properties variant."""
    seqs = [test_sequences['only_negative'],
            test_sequences['only_positive'],
            test_sequences['balanced_charge'],
            test_sequences['hydro_3']]
    for seq in seqs:
        # set targets
        starting_hydro = Protein(seq).hydrophobicity
        target_hydro = starting_hydro + random.uniform(-1.0, 1.0)
        if target_hydro < 1:
            target_hydro = 1.5
        target_FCR = random.uniform(0.1, 0.3)
        target_FCR = round(target_FCR * len(seq)) / len(seq)  # round to nearest residue fraction
        target_NCPR = random.uniform(-target_FCR, target_FCR)
        target_NCPR = round(target_NCPR * len(seq)) / len(seq)  # round to nearest residue fracti
        # only set kappa if we can actually change it.
        if target_FCR > 0 and target_NCPR != target_FCR:
            target_kappa = random.uniform(0.1, 0.5)
            check_kappa=True
        else:
            check_kappa=False
            target_kappa=None
        # make variant
        variant_seq = create.variant(seq, 'change_any_properties',
                                     target_hydropathy=target_hydro,
                                     target_FCR=target_FCR,
                                     target_NCPR=target_NCPR,
                                     target_kappa=target_kappa)
        assert verify_target_hydropathy(variant_seq, target_hydro)
        assert verify_target_FCR(variant_seq, target_FCR)
        assert verify_target_NCPR(variant_seq, target_NCPR)
        if check_kappa:
            assert verify_target_kappa(variant_seq, target_kappa)
        assert verify_disorder(seq, variant_seq)
        assert verify_same_length(seq, variant_seq)

def test_change_dimensions(test_sequences):
    """Test change_dimensions variant."""
    # doing fewer sequences because this is rather slow.
    seqs = [test_sequences['balanced_charge']]
    for dim_type in ['re', 'rg']:
        for change_type in ['increase', 'decrease']:
            for seq in seqs:
                variant_seq = create.variant(seq, 'change_dimensions',
                                             increase_or_decrease=change_type,
                                             rg_or_re=dim_type)
                assert verify_dimensions(seq, variant_seq, increase_or_decrease=change_type, rg_or_re=dim_type)
                assert verify_disorder(seq, variant_seq)
                assert verify_same_length(seq, variant_seq)