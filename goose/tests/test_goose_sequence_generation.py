"""
Comprehensive pytest tests for GOOSE sequence generation functions.

This module contains thorough tests for the main sequence generation functions
in goose.create: sequence, seq_by_fractions, and seq_by_classes.
Tests cover all parameters, edge cases, and error conditions.
"""

import pytest
import goose
from goose import goose_exceptions
from goose.backend import parameters
import re


class TestSequenceGeneration:
    """Test suite for the main sequence() function."""

    def test_basic_sequence_generation(self):
        """Test basic sequence generation with default parameters."""
        seq = goose.sequence(50)
        assert isinstance(seq, str)
        assert len(seq) == 50
        assert all(aa in 'ACDEFGHIKLMNPQRSTVWY' for aa in seq)

    def test_sequence_with_fcr(self):
        """Test sequence generation with FCR parameter."""
        seq = goose.sequence(50, FCR=0.3)
        assert isinstance(seq, str)
        assert len(seq) == 50
        # Count charged residues (D, E, K, R)
        charged_count = sum(1 for aa in seq if aa in 'DEKR')
        fcr_actual = charged_count / len(seq)
        # Allow some tolerance for FCR matching
        assert abs(fcr_actual - 0.3) < 0.1

    def test_sequence_with_ncpr(self):
        """Test sequence generation with NCPR parameter."""
        seq = goose.sequence(50, NCPR=0.2)
        assert isinstance(seq, str)
        assert len(seq) == 50
        # Count positive and negative residues
        pos_count = sum(1 for aa in seq if aa in 'KR')
        neg_count = sum(1 for aa in seq if aa in 'DE')
        ncpr_actual = (pos_count - neg_count) / len(seq)
        # Allow some tolerance for NCPR matching
        assert abs(ncpr_actual - 0.2) < 0.15

    def test_sequence_with_hydropathy(self):
        """Test sequence generation with hydropathy parameter."""
        seq = goose.sequence(50, hydropathy=3.0)
        assert isinstance(seq, str)
        assert len(seq) == 50

    def test_sequence_with_kappa(self):
        """Test sequence generation with kappa parameter."""
        seq = goose.sequence(50, kappa=0.5)
        assert isinstance(seq, str)
        assert len(seq) == 50

    def test_sequence_with_multiple_properties(self):
        """Test sequence generation with multiple properties."""
        seq = goose.sequence(50, FCR=0.3, NCPR=0.1, hydropathy=3.0)
        assert isinstance(seq, str)
        assert len(seq) == 50

    def test_sequence_with_exclude_residues(self):
        """Test sequence generation with excluded residues."""
        seq = goose.sequence(50, exclude=['C', 'M'])
        assert isinstance(seq, str)
        assert len(seq) == 50
        assert 'C' not in seq
        assert 'M' not in seq

    def test_sequence_with_custom_probabilities_dict(self):
        """Test sequence generation with custom probabilities dictionary."""
        custom_probs = {'A': 0.5, 'G': 0.3, 'S': 0.2}
        seq = goose.sequence(50, custom_probabilities=custom_probs)
        assert isinstance(seq, str)
        assert len(seq) == 50

    def test_sequence_with_custom_probabilities_string(self):
        """Test sequence generation with custom probabilities string."""
        seq = goose.sequence(50, custom_probabilities='human')
        assert isinstance(seq, str)
        assert len(seq) == 50

    def test_sequence_with_attempts(self):
        """Test sequence generation with different attempt numbers."""
        seq = goose.sequence(50, attempts=5)
        assert isinstance(seq, str)
        assert len(seq) == 50

    def test_sequence_with_disorder_cutoff(self):
        """Test sequence generation with custom disorder cutoff."""
        seq = goose.sequence(50, disorder_cutoff=0.3)
        assert isinstance(seq, str)
        assert len(seq) == 50

    def test_sequence_with_use_weighted_probabilities(self):
        """Test sequence generation with weighted probabilities."""
        seq = goose.sequence(50, use_weighted_probabilities=True)
        assert isinstance(seq, str)
        assert len(seq) == 50

    def test_sequence_with_strict_disorder(self):
        """Test sequence generation with strict disorder checking."""
        seq = goose.sequence(50, strict_disorder=True)
        assert isinstance(seq, str)
        assert len(seq) == 50

    def test_sequence_with_return_all_sequences(self):
        """Test sequence generation returning all sequences."""
        seqs = goose.sequence(30, return_all_sequences=True, attempts=3)
        assert isinstance(seqs, list)
        assert len(seqs) >= 1
        for seq in seqs:
            assert isinstance(seq, str)
            assert len(seq) == 30

    def test_sequence_with_metapredict_version(self):
        """Test sequence generation with different MetaPredict versions."""
        for version in [1, 2, 3]:
            seq = goose.sequence(50, metapredict_version=version)
            assert isinstance(seq, str)
            assert len(seq) == 50

    def test_sequence_with_max_consecutive_ordered(self):
        """Test sequence generation with max consecutive ordered residues."""
        seq = goose.sequence(50, max_consecutive_ordered=5)
        assert isinstance(seq, str)
        assert len(seq) == 50

    def test_sequence_with_max_total_ordered(self):
        """Test sequence generation with max total ordered fraction."""
        seq = goose.sequence(50, max_total_ordered=0.1)
        assert isinstance(seq, str)
        assert len(seq) == 50

    def test_sequence_with_batch_size(self):
        """Test sequence generation with custom batch size."""
        seq = goose.sequence(50, batch_size=10)
        assert isinstance(seq, str)
        assert len(seq) == 50

    def test_sequence_with_tolerance_parameters(self):
        """Test sequence generation with custom tolerance parameters."""
        seq = goose.sequence(50, hydropathy=3.0, hydropathy_tolerance=0.1, 
                           kappa=0.5, kappa_tolerance=0.05)
        assert isinstance(seq, str)
        assert len(seq) == 50

    def test_sequence_legacy_parameter_names(self):
        """Test sequence generation with legacy parameter names."""
        # Test hydrophobicity -> hydropathy
        seq = goose.sequence(50, hydrophobicity=3.0)
        assert isinstance(seq, str)
        assert len(seq) == 50

        # Test cutoff -> disorder_cutoff
        seq = goose.sequence(50, cutoff=0.4)
        assert isinstance(seq, str)
        assert len(seq) == 50

    def test_sequence_length_validation(self):
        """Test sequence length validation."""
        # Test minimum length
        with pytest.raises(goose_exceptions.GooseInputError):
            goose.sequence(5)  # Below minimum
        
        # Test maximum length (this might take a while, so use a reasonable upper bound)
        seq = goose.sequence(100)  # Should work
        assert len(seq) == 100

    def test_sequence_invalid_parameters(self):
        """Test sequence generation with invalid parameters."""
        # Test invalid FCR
        with pytest.raises(goose_exceptions.GooseInputError):
            goose.sequence(50, FCR=1.5)
        
        # Test invalid NCPR
        with pytest.raises(goose_exceptions.GooseInputError):
            goose.sequence(50, NCPR=2.0)
        
        # Test invalid hydropathy
        with pytest.raises(goose_exceptions.GooseInputError):
            goose.sequence(50, hydropathy=10.0)
        
        # Test invalid kappa
        with pytest.raises(goose_exceptions.GooseInputError):
            goose.sequence(50, kappa=2.0)

    def test_sequence_unknown_parameters(self):
        """Test sequence generation with unknown parameters."""
        with pytest.raises(goose_exceptions.GooseInputError):
            goose.sequence(50, unknown_param=0.5)

    def test_sequence_extreme_parameters(self):
        """Test sequence generation with extreme but valid parameters."""
        # Test extreme FCR
        seq = goose.sequence(50, FCR=1.0)
        assert isinstance(seq, str)
        assert len(seq) == 50
        
        # Test extreme NCPR
        seq = goose.sequence(50, NCPR=-1.0)
        assert isinstance(seq, str)
        assert len(seq) == 50


class TestSeqByFractions:
    """Test suite for the seq_by_fractions() function."""

    def test_basic_seq_by_fractions(self):
        """Test basic sequence generation by fractions."""
        seq = goose.seq_by_fractions(50)
        assert isinstance(seq, str)
        assert len(seq) == 50
        assert all(aa in 'ACDEFGHIKLMNPQRSTVWY' for aa in seq)

    def test_seq_by_fractions_single_amino_acid(self):
        """Test sequence generation with single amino acid fraction."""
        seq = goose.seq_by_fractions(100, A=0.3)
        assert isinstance(seq, str)
        assert len(seq) == 100
        # Count alanine residues
        a_count = seq.count('A')
        a_fraction = a_count / len(seq)
        # Allow some tolerance for fraction matching
        assert abs(a_fraction - 0.3) < 0.1

    def test_seq_by_fractions_multiple_amino_acids(self):
        """Test sequence generation with multiple amino acid fractions."""
        seq = goose.seq_by_fractions(100, A=0.2, G=0.1, S=0.15)
        assert isinstance(seq, str)
        assert len(seq) == 100
        
        # Check individual fractions
        a_fraction = seq.count('A') / len(seq)
        g_fraction = seq.count('G') / len(seq)
        s_fraction = seq.count('S') / len(seq)
        
        assert abs(a_fraction - 0.2) < 0.1
        assert abs(g_fraction - 0.1) < 0.1
        assert abs(s_fraction - 0.15) < 0.1

    def test_seq_by_fractions_all_amino_acids(self):
        """Test sequence generation with all amino acid fractions specified."""
        fractions = {
            'A': 0.05, 'C': 0.05, 'D': 0.05, 'E': 0.05, 'F': 0.05,
            'G': 0.05, 'H': 0.05, 'I': 0.05, 'K': 0.05, 'L': 0.05,
            'M': 0.05, 'N': 0.05, 'P': 0.05, 'Q': 0.05, 'R': 0.05,
            'S': 0.05, 'T': 0.05, 'V': 0.05, 'W': 0.05, 'Y': 0.05
        }
        seq = goose.seq_by_fractions(100, **fractions)
        assert isinstance(seq, str)
        assert len(seq) == 100

    def test_seq_by_fractions_with_max_aa_fractions(self):
        """Test sequence generation with custom max amino acid fractions."""
        max_fractions = {'A': 0.5, 'G': 0.3}
        seq = goose.seq_by_fractions(100, A=0.4, max_aa_fractions=max_fractions)
        assert isinstance(seq, str)
        assert len(seq) == 100

    def test_seq_by_fractions_with_remaining_probabilities_dict(self):
        """Test sequence generation with remaining probabilities dictionary."""
        remaining_probs = {'A': 0.5, 'G': 0.3, 'S': 0.2}
        seq = goose.seq_by_fractions(50, A=0.2, remaining_probabilities=remaining_probs)
        assert isinstance(seq, str)
        assert len(seq) == 50

    def test_seq_by_fractions_with_remaining_probabilities_string(self):
        """Test sequence generation with remaining probabilities string."""
        seq = goose.seq_by_fractions(50, A=0.2, remaining_probabilities='human')
        assert isinstance(seq, str)
        assert len(seq) == 50

    def test_seq_by_fractions_with_attempts(self):
        """Test sequence generation with different attempt numbers."""
        seq = goose.seq_by_fractions(50, A=0.2, attempts=5)
        assert isinstance(seq, str)
        assert len(seq) == 50

    def test_seq_by_fractions_with_disorder_cutoff(self):
        """Test sequence generation with custom disorder cutoff."""
        seq = goose.seq_by_fractions(50, A=0.2, disorder_cutoff=0.3)
        assert isinstance(seq, str)
        assert len(seq) == 50

    def test_seq_by_fractions_with_strict_disorder(self):
        """Test sequence generation with strict disorder checking."""
        seq = goose.seq_by_fractions(50, A=0.2, strict_disorder=True)
        assert isinstance(seq, str)
        assert len(seq) == 50

    def test_seq_by_fractions_with_return_all_sequences(self):
        """Test sequence generation returning all sequences."""
        seqs = goose.seq_by_fractions(30, A=0.2, return_all_sequences=True, attempts=3)
        assert isinstance(seqs, list)
        assert len(seqs) >= 1
        for seq in seqs:
            assert isinstance(seq, str)
            assert len(seq) == 30

    def test_seq_by_fractions_with_metapredict_version(self):
        """Test sequence generation with different MetaPredict versions."""
        for version in [1, 2, 3]:
            seq = goose.seq_by_fractions(50, A=0.2, metapredict_version=version)
            assert isinstance(seq, str)
            assert len(seq) == 50

    def test_seq_by_fractions_with_ordered_constraints(self):
        """Test sequence generation with ordered residue constraints."""
        seq = goose.seq_by_fractions(50, A=0.2, max_consecutive_ordered=5, 
                                   max_total_ordered=0.1)
        assert isinstance(seq, str)
        assert len(seq) == 50

    def test_seq_by_fractions_with_batch_size(self):
        """Test sequence generation with custom batch size."""
        seq = goose.seq_by_fractions(50, A=0.2, batch_size=10)
        assert isinstance(seq, str)
        assert len(seq) == 50

    def test_seq_by_fractions_legacy_parameter_names(self):
        """Test sequence generation with legacy parameter names."""
        # Test cutoff -> disorder_cutoff
        seq = goose.seq_by_fractions(50, A=0.2, cutoff=0.4)
        assert isinstance(seq, str)
        assert len(seq) == 50

    def test_seq_by_fractions_length_validation(self):
        """Test sequence length validation."""
        # Test minimum length
        with pytest.raises(goose_exceptions.GooseInputError):
            goose.seq_by_fractions(5, A=0.2)

    def test_seq_by_fractions_invalid_fractions(self):
        """Test sequence generation with invalid fractions."""
        # Test fraction > 1
        with pytest.raises(goose_exceptions.GooseInputError):
            goose.seq_by_fractions(50, A=1.5)
        
        # Test negative fraction
        with pytest.raises(goose_exceptions.GooseInputError):
            goose.seq_by_fractions(50, A=-0.1)
        
        # Test sum of fractions > 1
        with pytest.raises(goose_exceptions.GooseInputError):
            goose.seq_by_fractions(50, A=0.6, G=0.6)

    def test_seq_by_fractions_unknown_parameters(self):
        """Test sequence generation with unknown parameters."""
        with pytest.raises(goose_exceptions.GooseInputError):
            goose.seq_by_fractions(50, A=0.2, unknown_param=0.5)

    def test_seq_by_fractions_extreme_fractions(self):
        """Test sequence generation with extreme but valid fractions."""
        # Test very high fraction
        seq = goose.seq_by_fractions(50, A=0.8)
        assert isinstance(seq, str)
        assert len(seq) == 50
        
        # Test very low fraction
        seq = goose.seq_by_fractions(100, A=0.01)
        assert isinstance(seq, str)
        assert len(seq) == 100


class TestSeqByClasses:
    """Test suite for the seq_by_classes() function."""

    def test_basic_seq_by_classes(self):
        """Test basic sequence generation by classes."""
        seq = goose.seq_by_classes(50)
        assert isinstance(seq, str)
        assert len(seq) == 50
        assert all(aa in 'ACDEFGHIKLMNPQRSTVWY' for aa in seq)

    def test_seq_by_classes_single_class(self):
        """Test sequence generation with single class fraction."""
        seq = goose.seq_by_classes(100, aromatic=0.2)
        assert isinstance(seq, str)
        assert len(seq) == 100
        # Count aromatic residues (F, W, Y)
        aromatic_count = sum(1 for aa in seq if aa in 'FWY')
        aromatic_fraction = aromatic_count / len(seq)
        # Allow some tolerance for fraction matching
        assert abs(aromatic_fraction - 0.2) < 0.1

    def test_seq_by_classes_multiple_classes(self):
        """Test sequence generation with multiple class fractions."""
        seq = goose.seq_by_classes(100, aromatic=0.15, aliphatic=0.2, polar=0.1)
        assert isinstance(seq, str)
        assert len(seq) == 100
        
        # Check individual class fractions
        aromatic_count = sum(1 for aa in seq if aa in 'FWY')
        aliphatic_count = sum(1 for aa in seq if aa in 'AILV')
        polar_count = sum(1 for aa in seq if aa in 'NQST')
        
        aromatic_fraction = aromatic_count / len(seq)
        aliphatic_fraction = aliphatic_count / len(seq)
        polar_fraction = polar_count / len(seq)
        
        assert abs(aromatic_fraction - 0.15) < 0.1
        assert abs(aliphatic_fraction - 0.2) < 0.1
        assert abs(polar_fraction - 0.1) < 0.1

    def test_seq_by_classes_all_classes(self):
        """Test sequence generation with all class fractions specified."""
        seq = goose.seq_by_classes(100, aromatic=0.1, aliphatic=0.1, polar=0.1,
                                 positive=0.1, negative=0.1, glycine=0.1,
                                 proline=0.1, cysteine=0.05, histidine=0.05)
        assert isinstance(seq, str)
        assert len(seq) == 100

    def test_seq_by_classes_charged_residues(self):
        """Test sequence generation with charged residue classes."""
        seq = goose.seq_by_classes(100, positive=0.2, negative=0.15)
        assert isinstance(seq, str)
        assert len(seq) == 100
        
        # Check charged residue fractions
        positive_count = sum(1 for aa in seq if aa in 'KR')
        negative_count = sum(1 for aa in seq if aa in 'DE')
        
        positive_fraction = positive_count / len(seq)
        negative_fraction = negative_count / len(seq)
        
        assert abs(positive_fraction - 0.2) < 0.1
        assert abs(negative_fraction - 0.15) < 0.1

    def test_seq_by_classes_special_residues(self):
        """Test sequence generation with special residue classes."""
        seq = goose.seq_by_classes(100, glycine=0.15, proline=0.1, 
                                 cysteine=0.05, histidine=0.05)
        assert isinstance(seq, str)
        assert len(seq) == 100
        
        # Check special residue fractions
        g_fraction = seq.count('G') / len(seq)
        p_fraction = seq.count('P') / len(seq)
        c_fraction = seq.count('C') / len(seq)
        h_fraction = seq.count('H') / len(seq)
        
        assert abs(g_fraction - 0.15) < 0.1
        assert abs(p_fraction - 0.1) < 0.1
        assert abs(c_fraction - 0.05) < 0.1
        assert abs(h_fraction - 0.05) < 0.1

    def test_seq_by_classes_with_max_class_fractions(self):
        """Test sequence generation with custom max class fractions."""
        max_fractions = {'aromatic': 0.5, 'aliphatic': 0.4}
        seq = goose.seq_by_classes(100, aromatic=0.3, aliphatic=0.25, 
                                 max_class_fractions=max_fractions)
        assert isinstance(seq, str)
        assert len(seq) == 100

    def test_seq_by_classes_with_remaining_probabilities_dict(self):
        """Test sequence generation with remaining probabilities dictionary."""
        remaining_probs = {'A': 0.5, 'G': 0.3, 'S': 0.2}
        seq = goose.seq_by_classes(50, aromatic=0.2, remaining_probabilities=remaining_probs)
        assert isinstance(seq, str)
        assert len(seq) == 50

    def test_seq_by_classes_with_remaining_probabilities_string(self):
        """Test sequence generation with remaining probabilities string."""
        seq = goose.seq_by_classes(50, aromatic=0.2, remaining_probabilities='human')
        assert isinstance(seq, str)
        assert len(seq) == 50

    def test_seq_by_classes_with_attempts(self):
        """Test sequence generation with different attempt numbers."""
        seq = goose.seq_by_classes(50, aromatic=0.2, num_attempts=5)
        assert isinstance(seq, str)
        assert len(seq) == 50

    def test_seq_by_classes_with_disorder_cutoff(self):
        """Test sequence generation with custom disorder cutoff."""
        seq = goose.seq_by_classes(50, aromatic=0.2, disorder_cutoff=0.3)
        assert isinstance(seq, str)
        assert len(seq) == 50

    def test_seq_by_classes_with_strict_disorder(self):
        """Test sequence generation with strict disorder checking."""
        seq = goose.seq_by_classes(50, aromatic=0.2, strict_disorder=True)
        assert isinstance(seq, str)
        assert len(seq) == 50

    def test_seq_by_classes_with_metapredict_version(self):
        """Test sequence generation with different MetaPredict versions."""
        for version in [1, 2, 3]:
            seq = goose.seq_by_classes(50, aromatic=0.2, metapredict_version=version)
            assert isinstance(seq, str)
            assert len(seq) == 50

    def test_seq_by_classes_with_ordered_constraints(self):
        """Test sequence generation with ordered residue constraints."""
        seq = goose.seq_by_classes(50, aromatic=0.2, max_consecutive_ordered=5, 
                                 max_total_ordered=0.1)
        assert isinstance(seq, str)
        assert len(seq) == 50

    def test_seq_by_classes_legacy_parameter_names(self):
        """Test sequence generation with legacy parameter names."""
        # Test cutoff -> disorder_cutoff
        seq = goose.seq_by_classes(50, aromatic=0.2, cutoff=0.4)
        assert isinstance(seq, str)
        assert len(seq) == 50

    def test_seq_by_classes_length_validation(self):
        """Test sequence length validation."""
        # Test minimum length
        with pytest.raises(goose_exceptions.GooseInputError):
            goose.seq_by_classes(5, aromatic=0.2)

    def test_seq_by_classes_invalid_fractions(self):
        """Test sequence generation with invalid class fractions."""
        # Test fraction > 1
        with pytest.raises(goose_exceptions.GooseInputError):
            goose.seq_by_classes(50, aromatic=1.5)
        
        # Test negative fraction
        with pytest.raises(goose_exceptions.GooseInputError):
            goose.seq_by_classes(50, aromatic=-0.1)
        
        # Test sum of fractions > 1
        with pytest.raises(goose_exceptions.GooseInputError):
            goose.seq_by_classes(50, aromatic=0.6, aliphatic=0.6)

    def test_seq_by_classes_extreme_fractions(self):
        """Test sequence generation with extreme but valid fractions."""
        # Test very high fraction
        seq = goose.seq_by_classes(50, aromatic=0.8)
        assert isinstance(seq, str)
        assert len(seq) == 50
        
        # Test very low fraction
        seq = goose.seq_by_classes(100, aromatic=0.01)
        assert isinstance(seq, str)
        assert len(seq) == 100


class TestSequenceGenerationIntegration:
    """Integration tests for sequence generation functions."""

    def test_sequence_consistency(self):
        """Test that sequence generation is consistent across calls."""
        # Multiple calls should produce valid sequences
        for _ in range(5):
            seq = goose.sequence(50, FCR=0.3)
            assert isinstance(seq, str)
            assert len(seq) == 50

    def test_seq_by_fractions_consistency(self):
        """Test that seq_by_fractions is consistent across calls."""
        # Multiple calls should produce valid sequences
        for _ in range(5):
            seq = goose.seq_by_fractions(50, A=0.2, G=0.1)
            assert isinstance(seq, str)
            assert len(seq) == 50

    def test_seq_by_classes_consistency(self):
        """Test that seq_by_classes is consistent across calls."""
        # Multiple calls should produce valid sequences
        for _ in range(5):
            seq = goose.seq_by_classes(50, aromatic=0.2, aliphatic=0.1)
            assert isinstance(seq, str)
            assert len(seq) == 50

    def test_function_parameter_compatibility(self):
        """Test parameter compatibility across functions."""
        # Test common parameters work across all functions
        common_params = {
            'disorder_cutoff': 0.4,
            'strict_disorder': False,
            'metapredict_version': 3,
            'max_consecutive_ordered': 5,
            'max_total_ordered': 0.1
        }
        
        # Test with sequence()
        seq1 = goose.sequence(50, **common_params)
        assert isinstance(seq1, str)
        assert len(seq1) == 50
        
        # Test with seq_by_fractions()
        seq2 = goose.seq_by_fractions(50, A=0.2, **common_params)
        assert isinstance(seq2, str)
        assert len(seq2) == 50
        
        # Test with seq_by_classes()
        seq3 = goose.seq_by_classes(50, aromatic=0.2, **common_params)
        assert isinstance(seq3, str)
        assert len(seq3) == 50

    def test_edge_case_lengths(self):
        """Test sequence generation with edge case lengths."""
        # Test minimum allowed length
        seq = goose.sequence(10)
        assert len(seq) == 10
        
        # Test larger lengths
        seq = goose.sequence(200)
        assert len(seq) == 200

    def test_amino_acid_composition_validation(self):
        """Test that generated sequences have valid amino acid composition."""
        valid_aas = set('ACDEFGHIKLMNPQRSTVWY')
        
        # Test sequence()
        seq = goose.sequence(100, FCR=0.3, NCPR=0.1)
        assert set(seq).issubset(valid_aas)
        
        # Test seq_by_fractions()
        seq = goose.seq_by_fractions(100, A=0.2, G=0.1, S=0.1)
        assert set(seq).issubset(valid_aas)
        
        # Test seq_by_classes()
        seq = goose.seq_by_classes(100, aromatic=0.2, aliphatic=0.1)
        assert set(seq).issubset(valid_aas)

    def test_error_handling_consistency(self):
        """Test that error handling is consistent across functions."""
        # Test invalid length across all functions
        with pytest.raises(goose_exceptions.GooseInputError):
            goose.sequence(5)
        
        with pytest.raises(goose_exceptions.GooseInputError):
            goose.seq_by_fractions(5, A=0.2)
        
        with pytest.raises(goose_exceptions.GooseInputError):
            goose.seq_by_classes(5, aromatic=0.2)

    def test_performance_with_difficult_constraints(self):
        """Test performance with difficult but achievable constraints."""
        # Test sequence with multiple constraints
        seq = goose.sequence(50, FCR=0.4, NCPR=0.2, hydropathy=3.0, attempts=5)
        assert isinstance(seq, str)
        assert len(seq) == 50
        
        # Test seq_by_fractions with many constraints
        seq = goose.seq_by_fractions(50, A=0.3, G=0.2, S=0.2, attempts=5)
        assert isinstance(seq, str)
        assert len(seq) == 50
        
        # Test seq_by_classes with many constraints
        seq = goose.seq_by_classes(50, aromatic=0.2, aliphatic=0.2, polar=0.2, 
                                 positive=0.1, num_attempts=5)
        assert isinstance(seq, str)
        assert len(seq) == 50


if __name__ == '__main__':
    pytest.main([__file__])
