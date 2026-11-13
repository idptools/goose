"""Compare optimized profiles with independent per-window Sparrow calculations."""

from itertools import product
import random

import numpy as np
import pytest
import sparrow

from goose.backend.optimizer_properties import (
    ConstraintType, CustomProperty, FCR, Hydrophobicity, Kappa, NCPR,
    _cached_profile_kappa,
)


PROPERTIES = (FCR, NCPR, Hydrophobicity, Kappa)


def reference_profile(prop, sequence):
    """Original window-by-window algorithm, including its end convention."""
    w = prop.window_size
    values = [prop.calculate_raw_value(sparrow.Protein(sequence[i:i + w]))
              for i in range(len(sequence) - w + 1)]
    if prop.end_mode != 'clip':
        left = w // 2
        right = len(sequence) - len(values) - left
        first, last = (values[0], values[-1]) if prop.end_mode == 'extend' else (0, 0)
        values = [first] * left + values + [last] * right
    return np.array(values)


@pytest.mark.parametrize('property_class', PROPERTIES)
@pytest.mark.parametrize('window', [1, 2, 5, 6, 7, 20, 21, 81])
@pytest.mark.parametrize('end_mode', ['clip', 'extend', 'zeros'])
def test_profiles_match_sparrow(property_class, window, end_mode):
    rng = random.Random(105)
    prop = property_class(target_value=0, window_size=window, end_mode=end_mode)
    sequences = [
        ''.join(rng.choices('ACDEFGHIKLMNPQRSTVWY', k=81)),
        'E' * 81, 'K' * 81, 'A' * 81,
        ('KRDE' * 21)[:81], 'K' * 40 + 'E' * 41,
    ]
    for sequence in sequences:
        for text in (sequence, sequence.lower()):
            actual = prop.calculate_linear_profile(text)
            expected = reference_profile(prop, text)
            assert actual.dtype == expected.dtype
            np.testing.assert_allclose(actual, expected, rtol=1e-11, atol=1e-11)
            if property_class is Kappa:
                np.testing.assert_array_equal(actual, expected)


@pytest.mark.parametrize('length', [6, 7])
def test_kappa_exhaustive_charge_patterns(length):
    # Include every positive/negative/neutral ordering at these short lengths.
    prop = Kappa(target_value=0, window_size=length, end_mode='clip')
    for residues in product('KEA', repeat=length):
        sequence = ''.join(residues)
        assert prop.calculate_linear_profile(sequence)[0] == sparrow.Protein(sequence).kappa


@pytest.mark.parametrize('property_class', PROPERTIES)
@pytest.mark.parametrize('constraint', list(ConstraintType))
@pytest.mark.parametrize('end_mode', ['clip', 'extend', 'zeros'])
def test_calculate_preserves_profile_error_and_raw_average(property_class, constraint, end_mode):
    target = 'KKERDASQGPEKKDEQASGQDERAK'
    sequence = target[::-1]
    prop = property_class(
        calculate_as_linear_profile=True, target_sequence=target,
        window_size=6, end_mode=end_mode, constraint_type=constraint,
    )
    prop._set_sequence_length_hint(len(sequence))
    expected_target = reference_profile(prop, target)
    expected_profile = reference_profile(prop, sequence)
    difference = expected_profile - expected_target
    if constraint is ConstraintType.EXACT:
        expected_error = np.sum(np.abs(difference)) / len(sequence)
    elif constraint is ConstraintType.MINIMUM:
        expected_error = np.sum(np.maximum(0, -difference)) / len(sequence)
    else:
        expected_error = np.sum(np.maximum(0, difference)) / len(sequence)
    raw, error = prop.calculate(sparrow.Protein(sequence))
    assert raw == pytest.approx(np.mean(expected_profile), abs=1e-11)
    assert error == pytest.approx(expected_error, abs=1e-11)


@pytest.mark.parametrize('property_class', PROPERTIES)
def test_subclass_raw_calculation_is_respected(property_class):
    class OverriddenProperty(property_class):
        def calculate_raw_value(self, protein):
            return protein.sequence.count('A') + 12.0

    prop = OverriddenProperty(target_value=0, window_size=6)
    sequence = 'ADEKRQASQAPA'
    np.testing.assert_array_equal(prop.calculate_linear_profile(sequence),
                                  reference_profile(prop, sequence))


@pytest.mark.parametrize('property_class', PROPERTIES)
def test_instance_raw_calculation_is_respected(property_class):
    prop = property_class(target_value=0, window_size=6)
    prop.calculate_raw_value = lambda protein: protein.sequence.count('A') + 12.0
    sequence = 'ADEKRQASQAPA'
    np.testing.assert_array_equal(prop.calculate_linear_profile(sequence),
                                  reference_profile(prop, sequence))


def test_custom_property_can_use_window_call_state():
    class StatefulProperty(CustomProperty):
        def calculate_raw_value(self, protein):
            self.calls += 1
            return self.calls

    prop = StatefulProperty(target_value=0, window_size=2, end_mode='clip')
    prop.calls = 0
    np.testing.assert_array_equal(prop.calculate_linear_profile('AAAA'), [1, 2, 3])
    np.testing.assert_array_equal(prop.calculate_linear_profile('AAAA'), [4, 5, 6])


@pytest.mark.parametrize('property_class', PROPERTIES)
def test_scalar_calculation_still_uses_sparrow(property_class):
    protein = sparrow.Protein('DEKRQASGPDEKR')
    prop = property_class(target_value=0)
    assert prop.calculate(protein)[0] == prop.calculate_raw_value(protein)


@pytest.mark.parametrize('window', [5, 6, 21])
def test_optimizer_batch_matches_individual_profile_errors(window):
    from goose.optimize import SequenceOptimizer

    target = 'KKERDASQGPEKKDEQASGQDERAK'
    sequences = [target, target[::-1], 'A' * len(target)]
    optimizer = SequenceOptimizer(target_length=len(target), verbose=False)
    for cls in PROPERTIES:
        optimizer.add_property(cls, calculate_as_linear_profile=True,
                               target_sequence=target, window_size=window)
    batch = optimizer._evaluate_sequences_batch(sequences, use_cache=False)
    for sequence, (_, info) in zip(sequences, batch):
        for prop in optimizer.properties:
            profile = reference_profile(prop, sequence)
            target_profile = reference_profile(prop, target)
            expected_error = np.sum(np.abs(profile - target_profile)) / len(sequence)
            actual = info[prop.tracking_property_name]
            assert actual['raw_value'] == pytest.approx(np.mean(profile), abs=1e-11)
            assert actual['raw_error'] == pytest.approx(expected_error, abs=1e-11)


def test_kappa_cache_reuses_equivalent_charge_patterns():
    _cached_profile_kappa.cache_clear()
    try:
        prop = Kappa(target_value=0, window_size=6, end_mode='clip')
        first = prop.calculate_linear_profile('KKAAEE')
        info = _cached_profile_kappa.cache_info()
        second = prop.calculate_linear_profile('RRGGDD')
        np.testing.assert_array_equal(first, second)
        assert _cached_profile_kappa.cache_info().hits == info.hits + 1
        assert _cached_profile_kappa.cache_info().misses == info.misses
    finally:
        _cached_profile_kappa.cache_clear()


@pytest.mark.parametrize('property_class', PROPERTIES)
def test_unsupported_residue_uses_generic_behavior(property_class):
    prop = property_class(target_value=0, window_size=6)
    # Preserve each Sparrow property's handling of nonstandard residues.
    sequence = 'KDEXQAAKDER'
    if property_class is not Kappa:
        error_type = (sparrow.sparrow_exceptions.CalculationException
                      if property_class is Hydrophobicity else KeyError)
        with pytest.raises(error_type):
            reference_profile(prop, sequence)
        with pytest.raises(error_type):
            prop.calculate_linear_profile(sequence)
    else:
        np.testing.assert_array_equal(prop.calculate_linear_profile(sequence),
                                      reference_profile(prop, sequence))
