"""
Basic pytest tests for create.sequence variant functions.

This module contains Basic tests for the main sequence variant functions.
Tests should ensure that the functions behave as expected. 
"""
import random
import pytest
from goose import create
from goose import goose_exceptions
from goose.backend import parameters
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

# =============================================================================
# Test sequences - create a variety of test sequences with different properties
@pytest.fixture
def test_sequences():
    """Create test sequences with different properties for variant testing."""
    return {
        'basic': create.sequence(50, FCR=0.3, NCPR=0.1, hydropathy=2.5),
        'charged': create.sequence(60, FCR=0.4, NCPR=0.2, kappa=0.3),
        'hydrophobic': create.sequence(40, hydropathy=3.5, FCR=0.2),
        'long': create.sequence(100, FCR=0.25, NCPR=0.05, hydropathy=2.0),
        'mixed': create.sequence(80, FCR=0.35, NCPR=-0.1, hydropathy=2.8, kappa=0.4)
    }

# =============================================================================
# SHUFFLING VARIANT TESTS
# =============================================================================

def test_shuffle_specific_regions(test_sequences):
    """Test shuffle_specific_regions variant type."""
    seq = test_sequences['basic']
    shuffle_regions = [(10, 20), (30, 40)]
    
    variant_seq = create.variant(seq, 'shuffle_specific_regions', 
                                shuffle_regions=shuffle_regions)
    
    # Verify basic properties
    assert verify_same_length(seq, variant_seq)
    assert verify_disorder(seq, variant_seq)
    assert verify_same_residue_count(seq, variant_seq)
    
    # Verify constant regions (regions not shuffled)
    assert verify_constant_region(seq, variant_seq, (0, 10))
    assert verify_constant_region(seq, variant_seq, (20, 30))
    assert verify_constant_region(seq, variant_seq, (40, len(seq)))
    
    # Verify shuffled regions changed (if they have diverse amino acids)
    # Note: verify_region_changed handles cases where regions might not change
    assert verify_region_changed(seq, variant_seq, (10, 20))
    assert verify_region_changed(seq, variant_seq, (30, 40))

def test_shuffle_except_specific_regions(test_sequences):
    """Test shuffle_except_specific_regions variant type."""
    seq = test_sequences['long']
    excluded_regions = [(20, 30), (60, 70)]
    
    variant_seq = create.variant(seq, 'shuffle_except_specific_regions',
                                excluded_regions=excluded_regions)
    
    # Verify basic properties
    assert verify_same_length(seq, variant_seq)
    assert verify_disorder(seq, variant_seq)
    assert verify_same_residue_count(seq, variant_seq)
    
    # Verify excluded regions remain constant
    assert verify_constant_region(seq, variant_seq, (20, 30))
    assert verify_constant_region(seq, variant_seq, (60, 70))
    
    # Verify other regions were shuffled
    assert verify_region_changed(seq, variant_seq, (0, 20))
    assert verify_region_changed(seq, variant_seq, (70, len(seq)))

def test_shuffle_specific_residues(test_sequences):
    """Test shuffle_specific_residues variant type."""
    seq = test_sequences['charged']
    target_residues = ['K', 'R', 'E', 'D']
    
    variant_seq = create.variant(seq, 'shuffle_specific_residues',
                                target_residues=target_residues)
    
    # Verify basic properties
    assert verify_same_length(seq, variant_seq)
    assert verify_disorder(seq, variant_seq)
    assert verify_same_residue_count(seq, variant_seq)
    
    # For residues not in target_residues, positions should remain the same
    non_target_residues = [aa for aa in set(seq) if aa not in target_residues]
    if non_target_residues:
        assert verify_constant_residue_positions(seq, variant_seq, non_target_residues)

def test_shuffle_except_specific_residues(test_sequences):
    """Test shuffle_except_specific_residues variant type."""
    seq = test_sequences['mixed']
    excluded_residues = ['G', 'P', 'C']
    
    variant_seq = create.variant(seq, 'shuffle_except_specific_residues',
                                excluded_residues=excluded_residues)
    
    # Verify basic properties
    assert verify_same_length(seq, variant_seq)
    assert verify_disorder(seq, variant_seq)
    assert verify_same_residue_count(seq, variant_seq)
    
    # Verify excluded residues maintain their positions
    present_excluded = [aa for aa in excluded_residues if aa in seq]
    if present_excluded:
        assert verify_constant_residue_positions(seq, variant_seq, present_excluded)

def test_weighted_shuffle_specific_residues(test_sequences):
    """Test weighted_shuffle_specific_residues variant type."""
    seq = test_sequences['hydrophobic']
    target_residues = ['A', 'L', 'I', 'V']
    
    variant_seq = create.variant(seq, 'weighted_shuffle_specific_residues',
                                target_residues=target_residues,
                                shuffle_weight=0.7)
    
    # Verify basic properties
    assert verify_same_length(seq, variant_seq)
    assert verify_disorder(seq, variant_seq)
    assert verify_same_residue_count(seq, variant_seq)

def test_targeted_reposition_specific_residues(test_sequences):
    """Test targeted_reposition_specific_residues variant type."""
    seq = test_sequences['basic']
    possible_residues = list(set([a for a in seq]))
    target_residue = random.choice(possible_residues)
    
    variant_seq = create.variant(seq, 'targeted_reposition_specific_residues',
                                target_residues=target_residue)
    
    # Verify basic properties
    assert verify_same_length(seq, variant_seq)
    assert verify_disorder(seq, variant_seq)
    assert verify_same_residue_count(seq, variant_seq)

# =============================================================================
# ASYMMETRY VARIANT TESTS
# =============================================================================

def test_change_residue_asymmetry(test_sequences):
    """Test change_residue_asymmetry variant type."""
    seq = test_sequences['charged']
    target_residues = ['polar']
    
    variant_seq = create.variant(seq, 'change_residue_asymmetry',
                                target_residues=target_residues,
                                increase_or_decrease='increase')
    
    # Verify basic properties
    assert verify_same_length(seq, variant_seq)
    assert verify_disorder(seq, variant_seq)
    assert verify_same_residue_count(seq, variant_seq)
    
    # Verify asymmetry changed
    assert verify_changed_iwd(seq, variant_seq, 'increase', target_residues)

# =============================================================================
# PROPERTY-BASED VARIANT TESTS
# =============================================================================

def test_constant_properties(test_sequences):
    """Test constant_properties variant type."""
    seq = test_sequences['mixed']
    
    variant_seq = create.variant(seq, 'constant_properties')
    
    # Verify basic properties
    assert verify_same_length(seq, variant_seq)
    assert verify_disorder(seq, variant_seq)
    
    # Verify properties remain constant
    assert verify_constant_properties(seq, variant_seq)

def test_constant_residues_and_properties(test_sequences):
    """Test constant_residues_and_properties variant type."""
    seq = test_sequences['long']
    constant_residues = ['G', 'P']
    
    variant_seq = create.variant(seq, 'constant_residues_and_properties',
                                constant_residues=constant_residues)
    
    # Verify basic properties
    assert verify_same_length(seq, variant_seq)
    assert verify_disorder(seq, variant_seq)
    
    # Verify properties remain constant
    assert verify_constant_properties(seq, variant_seq)
    
    # Verify constant residues maintain positions
    present_constant = [aa for aa in constant_residues if aa in seq]
    if present_constant:
        assert verify_constant_residue_positions(seq, variant_seq, present_constant)

def test_constant_properties_and_class(test_sequences):
    """Test constant_properties_and_class variant type."""
    seq = test_sequences['basic']
    
    variant_seq = create.variant(seq, 'constant_properties_and_class')
    
    # Verify basic properties
    assert verify_same_length(seq, variant_seq)
    assert verify_disorder(seq, variant_seq)
    
    # Verify properties and class distribution remain constant
    assert verify_constant_properties(seq, variant_seq)
    assert verify_same_number_by_class(seq, variant_seq)

def test_constant_properties_and_class_by_order(test_sequences):
    """Test constant_properties_and_class_by_order variant type."""
    seq = test_sequences['hydrophobic']
    
    variant_seq = create.variant(seq, 'constant_properties_and_class_by_order')
    
    # Verify basic properties
    assert verify_same_length(seq, variant_seq)
    assert verify_disorder(seq, variant_seq)
    
    # Verify properties, class distribution, and order remain constant
    assert verify_constant_properties(seq, variant_seq)
    assert verify_same_number_by_class(seq, variant_seq)
    assert verify_same_order_by_class(seq, variant_seq)

# =============================================================================
# PROPERTY MODIFICATION VARIANT TESTS
# =============================================================================

def test_change_hydropathy_constant_class(test_sequences):
    """Test change_hydropathy_constant_class variant type."""
    seq = test_sequences['basic']
    original_protein = Protein(seq)
    target_hydropathy = original_protein.hydrophobicity + 0.5
    
    variant_seq = create.variant(seq, 'change_hydropathy_constant_class',
                                target_hydropathy=target_hydropathy)
    
    # Verify basic properties
    assert verify_same_length(seq, variant_seq)
    assert verify_disorder(seq, variant_seq)
    
    # Verify target hydropathy achieved and class distribution maintained
    assert verify_target_hydropathy(variant_seq, target_hydropathy)
    assert verify_same_number_by_class(seq, variant_seq)

def test_change_fcr_minimize_class_changes(test_sequences):
    """Test change_fcr_minimize_class_changes variant type."""
    seq = test_sequences['charged']
    original_protein = Protein(seq)
    target_FCR = min(original_protein.FCR + 0.1, 0.8)  # Ensure reasonable target
    
    variant_seq = create.variant(seq, 'change_fcr_minimize_class_changes',
                                target_FCR=target_FCR)
    
    # Verify basic properties
    assert verify_same_length(seq, variant_seq)
    assert verify_disorder(seq, variant_seq)
    
    # Verify target FCR achieved
    assert verify_target_FCR(variant_seq, target_FCR)

def test_change_ncpr_constant_class(test_sequences):
    """Test change_ncpr_constant_class variant type."""
    seq = test_sequences['mixed']
    original_protein = Protein(seq)
    target_NCPR = original_protein.NCPR + 0.05
    
    variant_seq = create.variant(seq, 'change_ncpr_constant_class',
                                target_NCPR=target_NCPR)
    
    # Verify basic properties
    assert verify_same_length(seq, variant_seq)
    assert verify_disorder(seq, variant_seq)
    
    # Verify target NCPR achieved
    assert verify_target_NCPR(variant_seq, target_NCPR)

def test_change_kappa(test_sequences):
    """Test change_kappa variant type."""
    seq = test_sequences['charged']
    target_kappa = 0.5
    
    variant_seq = create.variant(seq, 'change_kappa',
                                target_kappa=target_kappa)
    
    # Verify basic properties
    assert verify_same_length(seq, variant_seq)
    assert verify_disorder(seq, variant_seq)
    assert verify_same_residue_count(seq, variant_seq)
    
    # Verify target kappa achieved
    assert verify_target_kappa(variant_seq, target_kappa)

def test_change_properties_minimize_differences(test_sequences):
    """Test change_properties_minimize_differences variant type."""
    seq = test_sequences['long']
    original_protein = Protein(seq)
    
    # Only change a subset of properties
    target_hydropathy = original_protein.hydrophobicity + 0.3
    target_kappa = 0.6
    
    variant_seq = create.variant(seq, 'change_properties_minimize_differences',
                                target_hydropathy=target_hydropathy,
                                target_kappa=target_kappa)
    
    # Verify basic properties
    assert verify_same_length(seq, variant_seq)
    assert verify_disorder(seq, variant_seq)
    
    # Verify only changed properties
    assert verify_target_hydropathy(variant_seq, target_hydropathy)
    assert verify_target_kappa(variant_seq, target_kappa)
    
    # Verify unchanged properties remain within tolerance
    assert verify_constant_FCR(seq, variant_seq)
    assert verify_constant_NCPR(seq, variant_seq)

def test_change_any_properties(test_sequences):
    """Test change_any_properties variant type."""
    seq = test_sequences['basic']
    original_protein = Protein(seq)
    
    # Change multiple properties
    target_hydropathy = original_protein.hydrophobicity + 0.1
    target_FCR = min(original_protein.FCR + 0.1, 0.2)
    target_kappa = 0.3
    
    variant_seq = create.variant(seq, 'change_any_properties',
                                target_hydropathy=target_hydropathy,
                                target_FCR=target_FCR,
                                target_kappa=target_kappa)
    
    # Verify basic properties
    assert verify_same_length(seq, variant_seq)
    assert verify_disorder(seq, variant_seq)
    
    # Verify changed properties
    assert verify_target_hydropathy(variant_seq, target_hydropathy)
    assert verify_target_FCR(variant_seq, target_FCR)
    assert verify_target_kappa(variant_seq, target_kappa)

# =============================================================================
# DIMENSIONAL VARIANT TESTS
# =============================================================================

def test_change_dimensions_increase_re(test_sequences):
    """Test change_dimensions variant type - increase Re."""
    seq = test_sequences['charged']
    
    variant_seq = create.variant(seq, 'change_dimensions',
                                increase_or_decrease='increase',
                                rg_or_re='re')
    
    # Verify basic properties
    assert verify_same_length(seq, variant_seq)
    assert verify_disorder(seq, variant_seq)
    
    # Verify dimensions changed as expected
    assert verify_dimensions(seq, variant_seq, 're', 'increase')

def test_change_dimensions_decrease_rg(test_sequences):
    """Test change_dimensions variant type - decrease Rg."""
    seq = test_sequences['mixed']
    
    variant_seq = create.variant(seq, 'change_dimensions',
                                increase_or_decrease='decrease',
                                rg_or_re='rg')
    
    # Verify basic properties
    assert verify_same_length(seq, variant_seq)
    assert verify_disorder(seq, variant_seq)
    
    # Verify dimensions changed as expected
    assert verify_dimensions(seq, variant_seq, 'rg', 'decrease')

# =============================================================================
# ERROR HANDLING TESTS
# =============================================================================

def test_invalid_variant_type(test_sequences):
    """Test that invalid variant types raise appropriate errors."""
    seq = test_sequences['basic']
    
    with pytest.raises(goose_exceptions.GooseInputError):
        create.variant(seq, 'invalid_variant_type')

def test_missing_required_parameters(test_sequences):
    """Test that missing required parameters raise appropriate errors."""
    seq = test_sequences['basic']
    
    # Test missing shuffle_regions parameter
    with pytest.raises(goose_exceptions.GooseInputError):
        create.variant(seq, 'shuffle_specific_regions')
    
    # Test missing target_hydropathy parameter
    with pytest.raises(goose_exceptions.GooseInputError):
        create.variant(seq, 'change_hydropathy_constant_class')

def test_invalid_sequence_input():
    """Test that invalid sequence inputs raise appropriate errors."""
    # Test with invalid sequence
    with pytest.raises(goose_exceptions.GooseInputError):
        create.variant("INVALID123", 'constant_properties')
    
    # Test with empty sequence
    with pytest.raises(goose_exceptions.GooseInputError):
        create.variant("", 'constant_properties')

# =============================================================================
# EDGE CASES
# =============================================================================

def test_variant_with_single_residue_type():
    """Test variant generation with sequences containing limited residue diversity."""
    # Create a sequence with limited diversity
    seq = "A" * 30 + "G" * 20  # Simple sequence with only A and G
    
    # Test shuffle - should still work even with limited diversity
    variant_seq = create.variant(seq, 'shuffle_specific_regions',
                                shuffle_regions=[(10, 20)])
    
    assert verify_same_length(seq, variant_seq)
    assert verify_same_residue_count(seq, variant_seq)

def test_variant_with_small_regions():
    """Test variant generation with very small regions."""
    seq = create.sequence(30, FCR=0.3)
    
    # Test with small regions
    variant_seq = create.variant(seq, 'shuffle_specific_regions',
                                shuffle_regions=[(5, 7), (15, 17)])
    
    assert verify_same_length(seq, variant_seq)
    assert verify_same_residue_count(seq, variant_seq)

def test_variant_properties_tolerance():
    """Test that property changes are within expected tolerances."""
    seq = create.sequence(80, FCR=0.4, NCPR=0.1, hydropathy=3.0, kappa=0.3)
    
    # Test constant properties variant
    variant_seq = create.variant(seq, 'constant_properties')
    
    # Properties should be maintained within tolerances
    original_protein = Protein(seq)
    variant_protein = Protein(variant_seq)
    
    assert abs(original_protein.FCR - variant_protein.FCR) <= 0.001
    assert abs(original_protein.NCPR - variant_protein.NCPR) <= (1/len(seq)) + 1e-6
    assert abs(original_protein.hydrophobicity - variant_protein.hydrophobicity) <= parameters.MAXIMUM_HYDRO_ERROR
    assert abs(original_protein.kappa - variant_protein.kappa) <= parameters.MAXIMUM_KAPPA_ERROR

