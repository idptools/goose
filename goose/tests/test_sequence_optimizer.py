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
    def __init__(self, target_value, weight=1.0, tolerance=0.01):
        super().__init__(target_value, weight, tolerance)
    def calculate_raw_value(self, protein_obj):        
        sequence = protein_obj.sequence
        alanine_fraction = sequence.count('A') / len(sequence)
        return alanine_fraction
    
class GlycineContent(CustomProperty):
    def __init__(self, target_value, weight=1.0, tolerance=0.01):
        super().__init__(target_value, weight, tolerance)
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


if __name__ == '__main__':
    pytest.main([__file__, '-v'])
