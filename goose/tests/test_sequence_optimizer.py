"""
Comprehensive pytest tests for SequenceOptimizer.

This module contains thorough tests for the SequenceOptimizer class in goose.optimize.
Tests cover initialization, property addition, optimization functionality, 
and various optimization scenarios.
"""

import pytest
import numpy as np
from goose.optimize import SequenceOptimizer
from goose.backend.optimizer_properties import (
    FCR, NCPR, Hydrophobicity, Kappa, RadiusOfGyration, 
    EndToEndDistance, AminoAcidFractions, SCD, SHD, 
    Complexity, FractionDisorder, CustomProperty, ComputeIWD
)
from sparrow.protein import Protein


class TestSequenceOptimizerInitialization:
    """Test suite for SequenceOptimizer initialization."""

    def test_basic_initialization(self):
        """Test basic optimizer initialization with minimal parameters."""
        optimizer = SequenceOptimizer(target_length=50, verbose=False)
        assert optimizer.target_length == 50
        assert optimizer.max_iterations == 1000
        assert optimizer.verbose == False
        assert len(optimizer.properties) == 0

    def test_initialization_with_custom_parameters(self):
        """Test optimizer initialization with custom parameters."""
        optimizer = SequenceOptimizer(
            target_length=100,
            max_iterations=5000,
            verbose=True,
            num_candidates=10,
            enable_shuffling=False
        )
        assert optimizer.target_length == 100
        assert optimizer.max_iterations == 5000
        assert optimizer.verbose == True
        assert optimizer.num_candidates == 10
        assert optimizer.enable_shuffling == False

    def test_initialization_with_adaptive_scaling_disabled(self):
        """Test optimizer initialization with adaptive scaling disabled."""
        optimizer = SequenceOptimizer(
            target_length=50,
            enable_adaptive_scaling=False,
            verbose=False
        )
        assert optimizer.enable_adaptive_scaling == False

    def test_initialization_with_error_tolerance(self):
        """Test optimizer initialization with error tolerance."""
        optimizer = SequenceOptimizer(
            target_length=50,
            error_tolerance=1e-4,
            enable_error_tolerance=True,
            verbose=False
        )
        assert optimizer.error_tolerance == 1e-4
        assert optimizer.enable_error_tolerance == True

    def test_initialization_with_scalar_aa_fraction_targets(self):
        """Scalar aa_fraction_ranges values should be accepted as exact targets."""
        optimizer = SequenceOptimizer(
            target_length=100,
            aa_fraction_ranges={
                'A': 0.10,
                ('W', 'F', 'Y'): 0.05,
            },
            verbose=False,
        )

        assert optimizer._aa_count_bounds['A'] == (10, 10)
        assert (('W', 'F', 'Y'), 5, 5) in optimizer._aa_group_bounds

    def test_initialization_with_bound_aa_fraction_targets(self):
        """Range-based aa_fraction_ranges values should keep existing behavior."""
        optimizer = SequenceOptimizer(
            target_length=100,
            aa_fraction_ranges={
                'A': (0.05, 0.15),
                'WFY': (0.03, 0.12),
            },
            verbose=False,
        )

        assert optimizer._aa_count_bounds['A'] == (5, 15)
        assert (('W', 'F', 'Y'), 3, 12) in optimizer._aa_group_bounds


class TestPropertyAddition:
    """Test suite for adding properties to SequenceOptimizer."""

    def test_add_single_property(self):
        """Test adding a single property."""
        optimizer = SequenceOptimizer(target_length=50, verbose=False)
        optimizer.add_property(FCR, target_value=0.3, weight=1.0)
        assert len(optimizer.properties) == 1
        assert isinstance(optimizer.properties[0], FCR)
        assert optimizer.properties[0].target_value == 0.3
        assert optimizer.properties[0].weight == 1.0

    def test_add_multiple_properties(self):
        """Test adding multiple properties."""
        optimizer = SequenceOptimizer(target_length=50, verbose=False)
        optimizer.add_property(FCR, target_value=0.3, weight=1.0)
        optimizer.add_property(NCPR, target_value=0.2, weight=1.0)
        optimizer.add_property(Hydrophobicity, target_value=3.0, weight=1.0)
        assert len(optimizer.properties) == 3

    def test_add_property_with_tolerance(self):
        """Test adding a property with tolerance."""
        optimizer = SequenceOptimizer(target_length=50, verbose=False)
        optimizer.add_property(FCR, target_value=0.3, tolerance=0.01, weight=1.0)
        assert len(optimizer.properties) == 1
        assert optimizer.properties[0]._tolerance == 0.01

    def test_add_property_with_custom_weight(self):
        """Test adding a property with custom weight."""
        optimizer = SequenceOptimizer(target_length=50, verbose=False)
        optimizer.add_property(FCR, target_value=0.3, weight=2.0)
        assert optimizer.properties[0].weight == 2.0

    def test_add_kappa_property(self):
        """Test adding kappa property."""
        optimizer = SequenceOptimizer(target_length=50, verbose=False)
        optimizer.add_property(Kappa, target_value=0.5, weight=1.0)
        assert len(optimizer.properties) == 1
        assert isinstance(optimizer.properties[0], Kappa)

    def test_add_radius_of_gyration_property(self):
        """Test adding radius of gyration property."""
        optimizer = SequenceOptimizer(target_length=50, verbose=False)
        optimizer.add_property(RadiusOfGyration, target_value=25.0, weight=1.0)
        assert len(optimizer.properties) == 1
        assert isinstance(optimizer.properties[0], RadiusOfGyration)

    def test_add_amino_acid_fractions_property(self):
        """Test adding amino acid fractions property."""
        optimizer = SequenceOptimizer(target_length=50, verbose=False)
        fractions = {'A': 0.2, 'G': 0.3, 'S': 0.5}
        optimizer.add_property(AminoAcidFractions, target_fractions=fractions, weight=1.0)
        assert len(optimizer.properties) == 1
        assert isinstance(optimizer.properties[0], AminoAcidFractions)


class TestBasicOptimization:
    """Test suite for basic optimization functionality."""

    def test_optimize_single_property_fcr(self):
        """Test optimization with single FCR property."""
        optimizer = SequenceOptimizer(
            target_length=50, 
            max_iterations=5000,
            verbose=False
        )
        optimizer.add_property(FCR, target_value=0.3, weight=1.0)
        result = optimizer.run()
        
        assert isinstance(result, str)
        assert len(result) == 50
        
        # Check that FCR is close to target
        protein = Protein(result)
        assert abs(protein.FCR - 0.3) < 0.05

    def test_optimize_single_property_ncpr(self):
        """Test optimization with single NCPR property."""
        optimizer = SequenceOptimizer(
            target_length=50, 
            max_iterations=5000,
            verbose=False
        )
        optimizer.add_property(NCPR, target_value=0.2, weight=1.0)
        result = optimizer.run()
        
        assert isinstance(result, str)
        assert len(result) == 50
        
        # Check that NCPR is close to target
        protein = Protein(result)
        assert abs(protein.NCPR - 0.2) < 0.05

    def test_optimize_single_property_hydrophobicity(self):
        """Test optimization with single hydrophobicity property."""
        optimizer = SequenceOptimizer(
            target_length=50, 
            max_iterations=5000,
            verbose=False
        )
        optimizer.add_property(Hydrophobicity, target_value=3.5, weight=1.0)
        result = optimizer.run()
        
        assert isinstance(result, str)
        assert len(result) == 50
        
        # Check that hydrophobicity is close to target
        protein = Protein(result)
        assert abs(protein.hydrophobicity - 3.5) < 0.3

    def test_optimize_single_property_kappa(self):
        """Test optimization with single kappa property."""
        optimizer = SequenceOptimizer(
            target_length=50, 
            max_iterations=5000,
            verbose=False
        )
        optimizer.add_property(Kappa, target_value=0.3, weight=1.0)
        result = optimizer.run()
        
        assert isinstance(result, str)
        assert len(result) == 50
        
        # Check that kappa is close to target
        protein = Protein(result)
        assert abs(protein.kappa - 0.3) < 0.1


class TestMultiPropertyOptimization:
    """Test suite for multi-property optimization."""

    def test_optimize_fcr_and_ncpr(self):
        """Test optimization with FCR and NCPR properties."""
        optimizer = SequenceOptimizer(
            target_length=50, 
            max_iterations=5000,
            verbose=False
        )
        optimizer.add_property(FCR, target_value=0.3, weight=1.0)
        optimizer.add_property(NCPR, target_value=0.15, weight=1.0)
        result = optimizer.run()
        
        assert isinstance(result, str)
        assert len(result) == 50
        
        protein = Protein(result)
        assert abs(protein.FCR - 0.3) < 0.05
        assert abs(protein.NCPR - 0.15) < 0.05

    def test_optimize_fcr_ncpr_hydrophobicity(self):
        """Test optimization with FCR, NCPR, and hydrophobicity properties."""
        optimizer = SequenceOptimizer(
            target_length=50, 
            max_iterations=5000,
            verbose=False
        )
        optimizer.add_property(FCR, target_value=0.25, weight=1.0)
        optimizer.add_property(NCPR, target_value=0.1, weight=1.0)
        optimizer.add_property(Hydrophobicity, target_value=3.0, weight=1.0)
        result = optimizer.run()
        
        assert isinstance(result, str)
        assert len(result) == 50
        
        protein = Protein(result)
        assert abs(protein.FCR - 0.25) < 0.05
        assert abs(protein.NCPR - 0.1) < 0.05
        assert abs(protein.hydrophobicity - 3.0) < 0.3

    def test_optimize_with_different_weights(self):
        """Test optimization with different property weights."""
        optimizer = SequenceOptimizer(
            target_length=50, 
            max_iterations=5000,
            verbose=False
        )
        optimizer.add_property(FCR, target_value=0.3, weight=2.0)  # Higher weight
        optimizer.add_property(NCPR, target_value=0.2, weight=0.5)  # Lower weight
        result = optimizer.run()
        
        assert isinstance(result, str)
        assert len(result) == 50
        
        # FCR should be closer to target due to higher weight
        protein = Protein(result)
        fcr_error = abs(protein.FCR - 0.3)
        ncpr_error = abs(protein.NCPR - 0.2)
        
        # FCR should generally be closer to target (though not guaranteed)
        assert fcr_error < 0.1

    def test_optimize_four_properties(self):
        """Test optimization with four properties."""
        optimizer = SequenceOptimizer(
            target_length=100, 
            max_iterations=5000,
            verbose=False
        )
        optimizer.add_property(FCR, target_value=0.25, weight=1.0)
        optimizer.add_property(NCPR, target_value=0.1, weight=1.0)
        optimizer.add_property(Hydrophobicity, target_value=3.2, weight=1.0)
        optimizer.add_property(Kappa, target_value=0.4, weight=1.0)
        result = optimizer.run()
        
        assert isinstance(result, str)
        assert len(result) == 100
        
        protein = Protein(result)
        assert abs(protein.FCR - 0.25) < 0.06
        assert abs(protein.NCPR - 0.1) < 0.06
        assert abs(protein.hydrophobicity - 3.2) < 0.4
        assert abs(protein.kappa - 0.4) < 0.15


class TestAminoAcidFractions:
    """Test suite for amino acid fraction optimization."""

    def test_optimize_simple_fractions(self):
        """Test optimization with simple amino acid fractions."""
        optimizer = SequenceOptimizer(
            target_length=50, 
            max_iterations=5000,
            verbose=False
        )
        fractions = {'A': 0.3, 'G': 0.3, 'S': 0.4}
        optimizer.add_property(AminoAcidFractions, target_fractions=fractions, weight=1.0)
        result = optimizer.run()
        
        assert isinstance(result, str)
        assert len(result) == 50
        
        # Check that fractions are reasonably close
        a_count = result.count('A') / len(result)
        g_count = result.count('G') / len(result)
        s_count = result.count('S') / len(result)
        
        assert abs(a_count - 0.3) < 0.1
        assert abs(g_count - 0.3) < 0.1
        assert abs(s_count - 0.4) < 0.1

    def test_optimize_fractions_with_fcr(self):
        """Test optimization with amino acid fractions and FCR."""
        optimizer = SequenceOptimizer(
            target_length=60, 
            max_iterations=5000,
            verbose=False
        )
        fractions = {'D': 0.2, 'K': 0.2, 'G': 0.3, 'S': 0.3}
        optimizer.add_property(AminoAcidFractions, target_fractions=fractions, weight=1.0)
        optimizer.add_property(FCR, target_value=0.4, weight=1.0)
        result = optimizer.run()
        
        assert isinstance(result, str)
        assert len(result) == 60
        
        protein = Protein(result)
        assert abs(protein.FCR - 0.4) < 0.08


class TestEdgeCases:
    """Test suite for edge cases and error conditions."""

    def test_no_properties_raises_error(self):
        """Test that running optimizer without properties raises an error."""
        optimizer = SequenceOptimizer(target_length=50, verbose=False)
        with pytest.raises(ValueError, match="No properties defined"):
            optimizer.run()

    def test_very_short_sequence(self):
        """Test optimization with very short sequence."""
        optimizer = SequenceOptimizer(
            target_length=15, 
            max_iterations=5000,
            verbose=False
        )
        optimizer.add_property(FCR, target_value=0.3, weight=1.0)
        result = optimizer.run()
        
        assert isinstance(result, str)
        assert len(result) == 15

    def test_longer_sequence(self):
        """Test optimization with longer sequence."""
        optimizer = SequenceOptimizer(
            target_length=200, 
            max_iterations=5000,
            verbose=False
        )
        optimizer.add_property(FCR, target_value=0.25, weight=1.0)
        optimizer.add_property(NCPR, target_value=0.1, weight=1.0)
        result = optimizer.run()
        
        assert isinstance(result, str)
        assert len(result) == 200

    def test_extreme_fcr_value(self):
        """Test optimization with extreme FCR value."""
        optimizer = SequenceOptimizer(
            target_length=50, 
            max_iterations=5000,
            verbose=False
        )
        optimizer.add_property(FCR, target_value=0.8, weight=1.0)
        result = optimizer.run()
        
        assert isinstance(result, str)
        assert len(result) == 50
        
        protein = Protein(result)
        # Should get reasonably close even with extreme value
        assert protein.FCR > 0.6

    def test_zero_fcr_value(self):
        """Test optimization with zero FCR."""
        optimizer = SequenceOptimizer(
            target_length=50, 
            max_iterations=5000,
            verbose=False
        )
        optimizer.add_property(FCR, target_value=0.0, weight=1.0)
        result = optimizer.run()
        
        assert isinstance(result, str)
        assert len(result) == 50
        
        protein = Protein(result)
        # Should get very close to zero
        assert protein.FCR < 0.05


class TestToleranceAndConvergence:
    """Test suite for tolerance and convergence features."""

    def test_tolerance_stops_early(self):
        """Test that tolerance allows early stopping."""
        optimizer = SequenceOptimizer(
            target_length=50, 
            max_iterations=5000,
            error_tolerance=0.01,
            enable_error_tolerance=True,
            verbose=False
        )
        optimizer.add_property(FCR, target_value=0.3, tolerance=0.02, weight=1.0)
        result = optimizer.run()
        
        assert isinstance(result, str)
        assert len(result) == 50
        # Should converge before max iterations
        assert optimizer.iteration < 500

    def test_property_specific_tolerance(self):
        """Test that property-specific tolerance works."""
        optimizer = SequenceOptimizer(
            target_length=50, 
            max_iterations=5000,
            verbose=False
        )
        # Tight tolerance for FCR, loose for NCPR
        optimizer.add_property(FCR, target_value=0.3, tolerance=0.005, weight=1.0)
        optimizer.add_property(NCPR, target_value=0.2, tolerance=0.05, weight=1.0)
        result = optimizer.run()
        
        assert isinstance(result, str)
        assert len(result) == 50


class TestAdvancedFeatures:
    """Test suite for advanced optimizer features."""

    def test_shuffling_disabled(self):
        """Test optimization with shuffling disabled."""
        optimizer = SequenceOptimizer(
            target_length=50, 
            max_iterations=5000,
            enable_shuffling=False,
            verbose=False
        )
        optimizer.add_property(FCR, target_value=0.3, weight=1.0)
        result = optimizer.run()
        
        assert isinstance(result, str)
        assert len(result) == 50

    def test_adaptive_scaling_disabled(self):
        """Test optimization with adaptive scaling disabled."""
        optimizer = SequenceOptimizer(
            target_length=50, 
            max_iterations=5000,
            enable_adaptive_scaling=False,
            verbose=False
        )
        optimizer.add_property(FCR, target_value=0.3, weight=1.0)
        result = optimizer.run()
        
        assert isinstance(result, str)
        assert len(result) == 50

    def test_custom_starting_sequence(self):
        """Test optimization with custom starting sequence."""
        starting_seq = "G" * 50
        optimizer = SequenceOptimizer(
            target_length=50, 
            max_iterations=5000,
            verbose=False
        )
        optimizer.starting_sequence = starting_seq
        optimizer.add_property(FCR, target_value=0.3, weight=1.0)
        result = optimizer.run()
        
        assert isinstance(result, str)
        assert len(result) == 50
        # Result should be different from starting sequence
        assert result != starting_seq

    def test_fewer_candidates(self):
        """Test optimization with fewer candidates per iteration."""
        optimizer = SequenceOptimizer(
            target_length=50, 
            max_iterations=5000,
            num_candidates=2,
            verbose=False
        )
        optimizer.add_property(FCR, target_value=0.3, weight=1.0)
        result = optimizer.run()
        
        assert isinstance(result, str)
        assert len(result) == 50

    def test_more_candidates(self):
        """Test optimization with more candidates per iteration."""
        optimizer = SequenceOptimizer(
            target_length=50, 
            max_iterations=5000,
            num_candidates=10,
            verbose=False
        )
        optimizer.add_property(FCR, target_value=0.3, weight=1.0)
        result = optimizer.run()
        
        assert isinstance(result, str)
        assert len(result) == 50


class TestComplexProperties:
    """Test suite for complex property optimization."""

    def test_optimize_scd(self):
        """Test optimization with SCD property."""
        optimizer = SequenceOptimizer(
            target_length=50, 
            max_iterations=5000,
            verbose=False
        )
        optimizer.add_property(SCD, target_value=20.0, weight=1.0)
        result = optimizer.run()
        
        assert isinstance(result, str)
        assert len(result) == 50

    def test_optimize_complexity(self):
        """Test optimization with complexity property."""
        optimizer = SequenceOptimizer(
            target_length=50, 
            max_iterations=5000,
            verbose=False
        )
        optimizer.add_property(Complexity, target_value=0.5, weight=1.0)
        result = optimizer.run()
        
        assert isinstance(result, str)
        assert len(result) == 50

    def test_optimize_fraction_disorder(self):
        """Test optimization with fraction disorder property."""
        optimizer = SequenceOptimizer(
            target_length=50, 
            max_iterations=5000,
            verbose=False
        )
        optimizer.add_property(FractionDisorder, target_value=0.8, weight=1.0)
        result = optimizer.run()
        
        assert isinstance(result, str)
        assert len(result) == 50

# custom properties to use for testing
class AlanineFraction(CustomProperty):
    def __init__(self, target_value, weight=1.0):
        super().__init__(target_value=target_value, weight=weight)
    def calculate_raw_value(self, protein_obj):        
        sequence = protein_obj.sequence
        alanine_fraction = sequence.count('A') / len(sequence)
        return alanine_fraction
    
class GlycineContent(CustomProperty):
    def __init__(self, target_value, weight=1.0):
        super().__init__(target_value=target_value, weight=weight)
    def calculate_raw_value(self, protein_obj):
        sequence = protein_obj.sequence
        return sequence.count('G') / len(sequence)

class TestCustomProperty:
    """Test suite for custom property optimization."""


    def test_optimize_custom_property(self):
        """Test optimization with custom property."""
        
        optimizer = SequenceOptimizer(
            target_length=50, 
            max_iterations=5000,
            verbose=False
        )
        optimizer.add_property(
            AlanineFraction, 
            target_value=0.3,
            weight=1.0
        )
        result = optimizer.run()
        
        assert isinstance(result, str)
        assert len(result) == 50
        
        # Check that alanine fraction is close to target
        alanine_fraction = result.count('A') / len(result)
        assert abs(alanine_fraction - 0.3) < 0.1

    def test_optimize_custom_with_standard_properties(self):
        """Test optimization with custom property combined with standard properties."""
        optimizer = SequenceOptimizer(
            target_length=60, 
            max_iterations=5000,
            verbose=False
        )
        optimizer.add_property(
            GlycineContent, 
            target_value=0.25,
            weight=1.0
        )
        optimizer.add_property(FCR, target_value=0.3, weight=1.0)
        result = optimizer.run()
        
        assert isinstance(result, str)
        assert len(result) == 60


class TestReproducibility:
    """Test suite for reproducibility and consistency."""

    def test_multiple_runs_different_results(self):
        """Test that multiple runs produce different results (stochastic)."""
        results = []
        for _ in range(3):
            optimizer = SequenceOptimizer(
                target_length=50, 
                max_iterations=5000,
                verbose=False
            )
            optimizer.add_property(FCR, target_value=0.3, weight=1.0)
            result = optimizer.run()
            results.append(result)
        
        # Results should be different (stochastic optimization)
        assert len(set(results)) > 1

    def test_all_results_meet_length_requirement(self):
        """Test that all optimization results meet length requirement."""
        for _ in range(3):
            optimizer = SequenceOptimizer(
                target_length=50, 
                max_iterations=5000,
                verbose=False
            )
            optimizer.add_property(FCR, target_value=0.3, weight=1.0)
            result = optimizer.run()
            assert len(result) == 50


class TestBestSequenceTracking:
    """Test suite for best sequence tracking during optimization."""

    def test_best_sequence_is_stored(self):
        """Test that best sequence is stored after optimization."""
        optimizer = SequenceOptimizer(
            target_length=50, 
            max_iterations=5000,
            verbose=False
        )
        optimizer.add_property(FCR, target_value=0.3, weight=1.0)
        result = optimizer.run()
        
        assert optimizer.best_sequence == result
        assert len(optimizer.best_sequence) == 50

    def test_best_error_is_tracked(self):
        """Test that best error is tracked during optimization."""
        optimizer = SequenceOptimizer(
            target_length=50, 
            max_iterations=5000,
            verbose=False
        )
        optimizer.add_property(FCR, target_value=0.3, weight=1.0)
        optimizer.run()
        
        assert hasattr(optimizer, 'best_error')
        assert isinstance(optimizer.best_error, float)
        assert optimizer.best_error >= 0

    def test_iteration_count_is_tracked(self):
        """Test that iteration count is tracked."""
        optimizer = SequenceOptimizer(
            target_length=50, 
            max_iterations=5000,
            verbose=False
        )
        optimizer.add_property(FCR, target_value=0.3, weight=1.0)
        optimizer.run()
        
        assert hasattr(optimizer, 'iteration')
        assert isinstance(optimizer.iteration, int)
        assert optimizer.iteration > 0
        assert optimizer.iteration <= 5000


class TestSetInitialSequence:
    """Tests for set_initial_sequence optimization convergence."""

    def test_set_initial_sequence_converges_basic_properties(self):
        """Test that set_initial_sequence converges for NCPR + Hydrophobicity."""
        import random
        random.seed(99)
        start_seq = ''.join(random.choices('ACDEFGHIKLMNPQRSTVWY', k=50))
        
        optimizer = SequenceOptimizer(target_length=50, max_iterations=500, verbose=False)
        optimizer.add_property(NCPR, target_value=0.15, tolerance=0.03)
        optimizer.add_property(Hydrophobicity, target_value=3.8, tolerance=0.15)
        optimizer.set_initial_sequence(start_seq)
        result = optimizer.run()
        
        p = Protein(result)
        assert abs(p.NCPR - 0.15) <= 0.035
        assert abs(p.hydrophobicity - 3.8) <= 0.2

    def test_set_initial_sequence_does_not_precompute_normalization(self):
        """Test that set_initial_sequence defers normalization to run()."""
        import random
        random.seed(7)
        start_seq = ''.join(random.choices('ACDEFGHIKLMNPQRSTVWY', k=50))
        
        optimizer = SequenceOptimizer(target_length=50, max_iterations=100, verbose=False)
        optimizer.add_property(NCPR, target_value=0.1)
        optimizer.set_initial_sequence(start_seq)
        
        # Normalization should NOT be initialized yet
        assert not optimizer._normalization_initialized
        # best_error should still be inf (not evaluated yet)
        assert optimizer.best_error == float('inf')

    def test_set_initial_sequence_validates_length(self):
        """Test that set_initial_sequence rejects wrong-length sequences."""
        optimizer = SequenceOptimizer(target_length=50, verbose=False)
        optimizer.add_property(NCPR, target_value=0.1)
        
        with pytest.raises(ValueError, match="does not match target length"):
            optimizer.set_initial_sequence("AAAA")

    def test_set_initial_sequence_comparable_to_normal_path(self):
        """Test that set_initial_sequence achieves similar results to normal path."""
        import random
        random.seed(123)
        start_seq = ''.join(random.choices('ACDEFGHIKLMNPQRSTVWY', k=60))
        
        # With set_initial_sequence
        opt1 = SequenceOptimizer(target_length=60, max_iterations=300, verbose=False)
        opt1.add_property(FCR, target_value=0.25, tolerance=0.03)
        opt1.set_initial_sequence(start_seq)
        result1 = opt1.run()
        
        # Without (normal path)
        opt2 = SequenceOptimizer(target_length=60, max_iterations=300, verbose=False)
        opt2.add_property(FCR, target_value=0.25, tolerance=0.03)
        result2 = opt2.run()
        
        p1 = Protein(result1)
        p2 = Protein(result2)
        # Both should achieve the target
        assert abs(p1.FCR - 0.25) <= 0.05
        assert abs(p2.FCR - 0.25) <= 0.05


class TestFixedPositions:
    """Tests for the fixed_positions constraint."""

    def test_mapping_form_resolves_at_init(self):
        """A {index: residue} mapping is usable without a starting sequence."""
        optimizer = SequenceOptimizer(
            target_length=30, fixed_positions={0: 'M', 10: 'W'}, verbose=False
        )
        assert optimizer.get_fixed_positions() == {0: 'M', 10: 'W'}
        assert 0 not in optimizer._free_indices
        assert 10 not in optimizer._free_indices
        assert len(optimizer._free_indices) == 28

    def test_negative_indices_are_normalized(self):
        """Negative indices count from the end, Python-style."""
        optimizer = SequenceOptimizer(
            target_length=30, fixed_positions={-1: 'K', -2: 'R'}, verbose=False
        )
        assert optimizer.get_fixed_positions() == {29: 'K', 28: 'R'}

    def test_index_form_resolves_from_initial_sequence(self):
        """Bare indices take their residues from the starting sequence."""
        start = 'M' + 'A' * 28 + 'K'
        optimizer = SequenceOptimizer(
            target_length=30, fixed_positions=[0, 29], verbose=False
        )
        assert optimizer.get_fixed_positions() is None
        optimizer.set_initial_sequence(start)
        assert optimizer.get_fixed_positions() == {0: 'M', 29: 'K'}

    def test_index_form_reresolves_on_second_initial_sequence(self):
        """A second set_initial_sequence() re-reads the fixed residues."""
        optimizer = SequenceOptimizer(
            target_length=10, fixed_positions=[0, 9], verbose=False
        )
        optimizer.set_initial_sequence('MAAAAAAAAK')
        assert optimizer.get_fixed_positions() == {0: 'M', 9: 'K'}
        optimizer.set_initial_sequence('WAAAAAAAAR')
        assert optimizer.get_fixed_positions() == {0: 'W', 9: 'R'}

    def test_mapping_overrides_conflicting_initial_sequence(self):
        """An explicit mapping wins over the starting sequence's residues."""
        optimizer = SequenceOptimizer(
            target_length=20, fixed_positions={0: 'M', 5: 'W'}, verbose=False
        )
        optimizer.set_initial_sequence('A' * 20)
        assert optimizer.starting_sequence[0] == 'M'
        assert optimizer.starting_sequence[5] == 'W'

    def test_no_fixed_positions_leaves_state_empty(self):
        """Default behavior is unchanged when fixed_positions is not used."""
        optimizer = SequenceOptimizer(target_length=30, verbose=False)
        assert optimizer.get_fixed_positions() is None
        assert optimizer._fixed_map is None
        assert optimizer._free_indices == []

    def test_assignment_after_construction_takes_effect(self):
        """Setting .fixed_positions post-construction parses and validates."""
        optimizer = SequenceOptimizer(target_length=40, verbose=False)
        optimizer.fixed_positions = {0: 'M', 39: 'K'}
        assert optimizer.get_fixed_positions() == {0: 'M', 39: 'K'}
        assert 0 not in optimizer._free_indices

        with pytest.raises(ValueError, match="out of range"):
            optimizer.fixed_positions = {99: 'M'}

    def test_assignment_resolves_against_existing_initial_sequence(self):
        """Index form assigned after set_initial_sequence resolves immediately."""
        optimizer = SequenceOptimizer(target_length=20, verbose=False)
        optimizer.set_initial_sequence('M' + 'A' * 18 + 'K')
        optimizer.fixed_positions = [0, 19]
        assert optimizer.get_fixed_positions() == {0: 'M', 19: 'K'}

    def test_assignment_revalidates_starting_sequence_composition(self):
        """A mapping assigned later cannot silently break aa_fraction_ranges."""
        optimizer = SequenceOptimizer(
            target_length=20, verbose=False, aa_fraction_ranges={'W': (0.0, 0.05)}
        )
        optimizer.set_initial_sequence('W' + 'A' * 19)  # already at the W cap
        with pytest.raises(ValueError, match="violates aa_fraction_ranges"):
            optimizer.fixed_positions = {5: 'W'}

    def test_assignment_stamps_existing_starting_sequence(self):
        """A mapping assigned later is stamped onto the starting sequence."""
        optimizer = SequenceOptimizer(target_length=10, verbose=False)
        optimizer.set_initial_sequence('A' * 10)
        optimizer.fixed_positions = {0: 'M', 9: 'K'}
        assert optimizer.starting_sequence == 'M' + 'A' * 8 + 'K'

    def test_rejects_non_standard_residue_from_initial_sequence(self):
        """Index form must reject a non-standard residue at a fixed index."""
        for ranges in (None, {'G': (0.0, 0.5)}):
            optimizer = SequenceOptimizer(
                target_length=10, fixed_positions=[0], verbose=False,
                aa_fraction_ranges=ranges,
            )
            with pytest.raises(ValueError, match="invalid amino acid"):
                optimizer.set_initial_sequence('X' + 'A' * 9)

    def test_rejects_lowercase_residue_from_initial_sequence(self):
        optimizer = SequenceOptimizer(
            target_length=10, fixed_positions=[0], verbose=False
        )
        with pytest.raises(ValueError, match="invalid amino acid"):
            optimizer.set_initial_sequence('m' + 'A' * 9)

    def test_failed_assignment_leaves_state_untouched(self):
        """A rejected assignment must not drop the existing constraint."""
        optimizer = SequenceOptimizer(
            target_length=20, fixed_positions={0: 'M'}, verbose=False
        )
        optimizer.set_initial_sequence('M' + 'A' * 19)
        before = (
            optimizer.fixed_positions,
            optimizer.get_fixed_positions(),
            optimizer.starting_sequence,
        )
        for bad in ({99: 'W'}, {0: 'Z'}, ['M']):
            with pytest.raises(ValueError):
                optimizer.fixed_positions = bad
        assert (
            optimizer.fixed_positions,
            optimizer.get_fixed_positions(),
            optimizer.starting_sequence,
        ) == before

    def test_failed_assignment_rolls_back_composition_revalidation(self):
        optimizer = SequenceOptimizer(
            target_length=20, verbose=False, aa_fraction_ranges={'W': (0.0, 0.05)}
        )
        optimizer.set_initial_sequence('W' + 'A' * 19)
        before = (optimizer.get_fixed_positions(), optimizer.starting_sequence)
        with pytest.raises(ValueError):
            optimizer.fixed_positions = {5: 'W'}
        assert (optimizer.get_fixed_positions(), optimizer.starting_sequence) == before

    def test_rejected_initial_sequence_leaves_fixed_map_consistent(self):
        """A rejected sequence must not resolve the map against its residues."""
        optimizer = SequenceOptimizer(
            target_length=10, fixed_positions=[0], verbose=False,
            aa_fraction_ranges={'G': (0.5, 0.6)},
        )
        optimizer.set_initial_sequence('GGGGGAAAAA')
        before = (optimizer.get_fixed_positions(), optimizer.starting_sequence)
        with pytest.raises(ValueError, match="violates aa_fraction_ranges"):
            optimizer.set_initial_sequence('A' * 10)
        assert (optimizer.get_fixed_positions(), optimizer.starting_sequence) == before
        # The map still agrees with the sequence actually in use.
        assert optimizer.get_fixed_positions()[0] == optimizer.starting_sequence[0]

    def test_unresolved_index_form_cannot_generate_candidates(self):
        """Mutating before resolution would silently violate the constraint."""
        optimizer = SequenceOptimizer(
            target_length=20, fixed_positions=[0, 19], verbose=False
        )
        with pytest.raises(ValueError, match="no starting sequence was set"):
            optimizer._generate_candidates('M' + 'A' * 18 + 'K', 50, 0)

    def test_assignment_of_none_clears_the_constraint(self):
        optimizer = SequenceOptimizer(
            target_length=20, fixed_positions={0: 'M'}, verbose=False
        )
        optimizer.fixed_positions = None
        assert optimizer.get_fixed_positions() is None
        assert optimizer._free_indices == []

    def test_mutation_primitives_never_touch_fixed_positions(self):
        """All fixed-aware mutators leave the pinned residues in place."""
        optimizer = SequenceOptimizer(
            target_length=40, fixed_positions={0: 'M', 7: 'W', 39: 'K'}, verbose=False
        )
        seq = optimizer._apply_fixed('A' * 40)
        for _ in range(200):
            for mutant in (
                optimizer._fixed_point_mutation(seq),
                optimizer._fixed_multi_mutations(seq, 10),
                optimizer._fixed_aware_shuffle(seq),
                optimizer._swap_positions_mutation(seq, 5),
            ):
                assert len(mutant) == 40
                assert mutant[0] == 'M'
                assert mutant[7] == 'W'
                assert mutant[39] == 'K'

    def test_fixed_aware_shuffle_preserves_composition(self):
        """Shuffling only permutes the free positions."""
        import random
        random.seed(11)
        optimizer = SequenceOptimizer(
            target_length=40, fixed_positions={0: 'M', 39: 'K'}, verbose=False
        )
        seq = optimizer._apply_fixed(
            ''.join(random.choices('ACDEFGHIKLMNPQRSTVWY', k=40))
        )
        for _ in range(50):
            mutant = optimizer._fixed_aware_shuffle(seq)
            assert sorted(mutant) == sorted(seq)
            assert mutant[0] == 'M' and mutant[39] == 'K'

    def test_generated_candidates_respect_fixed_positions(self):
        """_generate_candidates honors the constraint, shuffle iterations included."""
        optimizer = SequenceOptimizer(
            target_length=40, fixed_positions={0: 'M', 7: 'W', 39: 'K'}, verbose=False
        )
        seq = optimizer._apply_fixed('A' * 40)
        # iteration 50 triggers the shuffle branch (shuffle_frequency default 50)
        for iteration in (1, 50, 100):
            for _ in range(50):
                for candidate in optimizer._generate_candidates(seq, iteration, 0):
                    assert len(candidate) == 40
                    assert candidate[0] == 'M'
                    assert candidate[7] == 'W'
                    assert candidate[39] == 'K'

    def test_run_preserves_fixed_positions_without_initial_sequence(self):
        """A full run from scratch returns a sequence carrying the fixed residues."""
        optimizer = SequenceOptimizer(
            target_length=50, max_iterations=150, verbose=False,
            fixed_positions={0: 'M', 10: 'W', 49: 'K'}, seed=7,
        )
        optimizer.add_property(FCR, target_value=0.25)
        result = optimizer.run()

        assert len(result) == 50
        assert result[0] == 'M'
        assert result[10] == 'W'
        assert result[49] == 'K'

    def test_run_preserves_fixed_positions_and_optimizes_the_rest(self):
        """The pinned residues survive while free positions move toward target."""
        import random
        random.seed(3)
        # Uncharged start (FCR 0) so reaching the target demands many edits.
        start = ''.join(random.choices('GSQNAPT', k=60))
        keep = [0, 5, 12, 19, 33, 59]

        optimizer = SequenceOptimizer(
            target_length=60, max_iterations=300, verbose=False,
            fixed_positions=keep, seed=3,
        )
        optimizer.set_initial_sequence(start)
        optimizer.add_property(FCR, target_value=0.5, tolerance=0.02)
        result = optimizer.run()

        # Fixed positions unchanged...
        for i in keep:
            assert result[i] == start[i]
        # ...while the rest of the sequence was substantially rewritten to
        # reach a target the starting sequence was far from.
        assert sum(1 for i in range(60) if result[i] != start[i]) > 10
        assert abs(Protein(result).FCR - 0.5) <= 0.05

    def test_run_with_range_object(self):
        """A contiguous motif can be pinned with range()."""
        import random
        random.seed(6)
        start = ''.join(random.choices('GSQNAPTED', k=40))

        optimizer = SequenceOptimizer(
            target_length=40, max_iterations=150, verbose=False,
            fixed_positions=range(10, 21), seed=6,
        )
        optimizer.set_initial_sequence(start)
        optimizer.add_property(FCR, target_value=0.4)
        result = optimizer.run()

        assert result[10:21] == start[10:21]

    def test_combines_with_aa_fraction_ranges(self):
        """Fixed residues count toward the composition bounds and both hold."""
        optimizer = SequenceOptimizer(
            target_length=50, max_iterations=150, verbose=False,
            fixed_positions={0: 'M', 25: 'W'},
            aa_fraction_ranges={'W': (0.02, 0.10), 'E': (0.1, 0.3)},
            seed=11,
        )
        optimizer.add_property(FCR, target_value=0.3)
        result = optimizer.run()

        assert result[0] == 'M' and result[25] == 'W'
        assert 0.02 <= result.count('W') / 50 <= 0.10
        assert 0.1 <= result.count('E') / 50 <= 0.3

    def test_constrained_builder_honors_fixed_positions(self):
        """The composition-constrained seed builder places fixed residues."""
        optimizer = SequenceOptimizer(
            target_length=40, verbose=False,
            fixed_positions={0: 'M', 3: 'C', 39: 'W'},
            aa_fraction_ranges={'C': (0.0, 0.1), 'G': (0.1, 0.2)},
        )
        for _ in range(100):
            seq = optimizer._build_constrained_random_sequence(40)
            assert len(seq) == 40
            assert seq[0] == 'M' and seq[3] == 'C' and seq[39] == 'W'
            assert seq.count('C') <= 4
            assert 4 <= seq.count('G') <= 8

    def test_fully_constrained_composition_with_fixed_positions(self):
        """The swap-only mutation path also respects fixed positions."""
        fractions = {aa: 2 / 40 for aa in 'ACDEFGHIKLMNPQRSTVWY'}
        optimizer = SequenceOptimizer(
            target_length=40, max_iterations=100, verbose=False,
            fixed_positions={0: 'W', 20: 'C', 39: 'K'},
            aa_fraction_ranges=fractions, seed=8,
        )
        assert optimizer._composition_fully_constrained
        optimizer.add_property(Kappa, target_value=0.25)
        result = optimizer.run()

        assert result[0] == 'W' and result[20] == 'C' and result[39] == 'K'
        assert all(result.count(aa) == 2 for aa in 'ACDEFGHIKLMNPQRSTVWY')

    def test_emergency_injection_preserves_fixed_positions(self):
        """The stagnation escape hatch cannot discard the fixed residues."""
        optimizer = SequenceOptimizer(
            target_length=40, verbose=False,
            fixed_positions={0: 'M', 39: 'K'}, seed=1,
        )
        optimizer.add_property(FCR, target_value=0.3)
        optimizer.best_sequence = optimizer._apply_fixed('A' * 40)
        optimizer.best_error = float('inf')
        for _ in range(10):
            optimizer._emergency_diversity_injection(float('inf'))
            assert optimizer.best_sequence[0] == 'M'
            assert optimizer.best_sequence[39] == 'K'

    def test_single_free_position(self):
        """Pinning all but one position still runs."""
        optimizer = SequenceOptimizer(
            target_length=5, max_iterations=25, verbose=False,
            fixed_positions={0: 'M', 1: 'W', 2: 'K', 4: 'E'}, seed=4,
        )
        optimizer.add_property(FCR, target_value=0.5)
        result = optimizer.run()

        assert len(result) == 5
        assert result[0] == 'M' and result[1] == 'W'
        assert result[2] == 'K' and result[4] == 'E'

    def test_reproducible_with_seed(self):
        """The pure-Python fixed path honors ``seed``."""
        def build_and_run(seed):
            optimizer = SequenceOptimizer(
                target_length=40, max_iterations=100, verbose=False,
                fixed_positions={0: 'M', 39: 'K'}, seed=seed,
            )
            optimizer.add_property(FCR, target_value=0.3)
            return optimizer.run()

        assert build_and_run(42) == build_and_run(42)

    def test_rejects_out_of_range_index(self):
        with pytest.raises(ValueError, match="out of range"):
            SequenceOptimizer(target_length=10, fixed_positions={20: 'M'}, verbose=False)

    def test_rejects_invalid_residue(self):
        with pytest.raises(ValueError, match="not a valid amino acid"):
            SequenceOptimizer(target_length=10, fixed_positions={0: 'Z'}, verbose=False)

    def test_rejects_pinning_every_position(self):
        with pytest.raises(ValueError, match="at least one position must remain free"):
            SequenceOptimizer(
                target_length=3, fixed_positions={0: 'M', 1: 'W', 2: 'K'}, verbose=False
            )

    def test_rejects_empty_specification(self):
        with pytest.raises(ValueError, match="empty"):
            SequenceOptimizer(target_length=10, fixed_positions={}, verbose=False)

    def test_rejects_residue_letters_in_iterable_form(self):
        with pytest.raises(ValueError, match="integer indices"):
            SequenceOptimizer(
                target_length=10, fixed_positions=['M', 'W'], verbose=False
            )

    def test_rejects_fixed_residues_exceeding_aa_maximum(self):
        with pytest.raises(ValueError, match="maximum allowed count"):
            SequenceOptimizer(
                target_length=100,
                fixed_positions={i: 'W' for i in range(10)},
                aa_fraction_ranges={'W': (0.0, 0.05)},
                verbose=False,
            )

    def test_rejects_unreachable_aa_minimum(self):
        with pytest.raises(ValueError, match="free position"):
            SequenceOptimizer(
                target_length=10,
                fixed_positions={i: 'A' for i in range(9)},
                aa_fraction_ranges={'G': (0.5, 0.6)},
                verbose=False,
            )

    def test_index_form_without_initial_sequence_raises_on_run(self):
        optimizer = SequenceOptimizer(
            target_length=20, fixed_positions=[0, 1], verbose=False
        )
        optimizer.add_property(FCR, target_value=0.3)
        with pytest.raises(ValueError, match="no starting sequence was set"):
            optimizer.run()


if __name__ == '__main__':
    pytest.main([__file__, '-v'])
