"""Regressions for outer batch attempts and opt-in inner optimization limits."""

import numpy as np
import pytest

from goose import create, property_library
from goose.backend_property_calculation.calculate_properties_batch import sequences_to_matrices
from goose.backend_property_optimization.optimize_hydropathy import optimize_hydropathy
from goose.backend_sequence_generation import sequence_generation
from goose.goose_exceptions import GooseFail, GooseInputError


SEQUENCE = 'K' * 5 + 'E' * 5 + 'G' * 15


def instrument_property_pipeline(monkeypatch, failing_stage=None):
    events = {'generate': [], 'hydropathy': [], 'kappa': []}

    def generate(self, length, **kwargs):
        events['generate'].append(kwargs)
        return sequences_to_matrices([SEQUENCE] * kwargs['num_sequences'])

    def hydropathy(sequences, target, **kwargs):
        events['hydropathy'].append(kwargs)
        return None if failing_stage == 'hydropathy' else sequences

    def kappa(sequences, target, **kwargs):
        events['kappa'].append(kwargs)
        return None if failing_stage == 'kappa' else [SEQUENCE]

    monkeypatch.setattr(sequence_generation.SequenceGenerator, 'generate_sequences_vectorized', generate)
    monkeypatch.setattr(sequence_generation, 'optimize_hydropathy', hydropathy)
    monkeypatch.setattr(sequence_generation, 'optimize_kappa', kappa)
    monkeypatch.setattr(sequence_generation, 'check_disorder', lambda sequences, **kwargs: sequences)
    return events


def test_property_library_forwards_outer_and_inner_limits_end_to_end(monkeypatch):
    events = instrument_property_pipeline(monkeypatch)
    library = property_library.generate_property_library(
        25, hydropathy_range=3, FCR_range=0.4, NCPR_range=0, kappa_range=0.5,
        attempts=1, max_hydropathy_iterations=7, max_kappa_iterations=11,
        on_generation_failure='raise', verbose=False)
    assert len(library) == 1
    assert len(events['generate']) == len(events['hydropathy']) == len(events['kappa']) == 1
    # One attempt is a batch of ten, rather than one candidate or mutation.
    assert events['generate'][0]['num_sequences'] == 10
    assert events['hydropathy'][0]['max_iterations'] == 7
    assert events['kappa'][0]['num_iterations'] == 11


@pytest.mark.parametrize('limits', [{}, {'max_hydropathy_iterations': None, 'max_kappa_iterations': None}])
def test_omitted_limits_preserve_optimizer_defaults(monkeypatch, limits):
    events = instrument_property_pipeline(monkeypatch)
    assert create.sequence(25, hydropathy=3, FCR=0.4, NCPR=0, kappa=0.5,
                           attempts=1, **limits) == SEQUENCE
    assert 'max_iterations' not in events['hydropathy'][0]
    assert 'num_iterations' not in events['kappa'][0]


@pytest.mark.parametrize('stage', ['hydropathy', 'kappa'])
@pytest.mark.parametrize('attempts', [1, 3])
def test_exhausted_inner_searches_respect_batch_attempts_and_raise_goosefail(monkeypatch, stage, attempts):
    events = instrument_property_pipeline(monkeypatch, failing_stage=stage)
    with pytest.raises(GooseFail):
        create.sequence(25, hydropathy=3, FCR=0.4, NCPR=0, kappa=0.5,
                        attempts=attempts, max_hydropathy_iterations=1, max_kappa_iterations=1)
    assert len(events['generate']) == attempts
    assert len(events[stage]) == attempts


@pytest.mark.parametrize('name', ['max_hydropathy_iterations', 'max_kappa_iterations'])
@pytest.mark.parametrize('value', [0, -1, 1.5, True, float('nan'), float('inf'), '1'])
def test_invalid_inner_limits_raise_before_generation_or_grid_planning(monkeypatch, name, value):
    def unexpected(*args, **kwargs):
        pytest.fail('Invalid iteration limit reached generation')

    monkeypatch.setattr(sequence_generation, 'by_properties', unexpected)
    with pytest.raises(GooseInputError, match=name):
        create.sequence(25, **{name: value})
    with pytest.raises(GooseInputError, match=name):
        property_library.enumerate_property_combinations(25, FCR_range=0.4, **{name: value})
    with pytest.raises(GooseInputError, match=name):
        property_library.generate_property_library(25, FCR_range=0.4, verbose=False, **{name: value})


@pytest.mark.parametrize('convert_back', [True, False])
def test_hydropathy_optimizer_keeps_a_single_successful_candidate(convert_back):
    sequence = 'G' * 25
    result = optimize_hydropathy([sequence], 4.1, max_iterations=1,
                                 need_to_convert_seqs=True,
                                 convert_back_to_sequences=convert_back)
    if convert_back:
        assert result == [sequence]
    else:
        np.testing.assert_array_equal(result, sequences_to_matrices([sequence]))


def test_batch_size_one_can_generate_a_hydropathy_constrained_sequence():
    # At FCR=0 and with all remaining probability assigned to G, every starter
    # has exactly the requested hydropathy. Exercise the real generation and
    # optimizer code with a single candidate; no stochastic success is needed.
    sequence = create.sequence(
        25, hydropathy=4.1, FCR=0, attempts=1, batch_size=1,
        max_hydropathy_iterations=1, max_kappa_iterations=1, disorder_cutoff=0,
        custom_probabilities={'G': 1.0})
    assert sequence == 'G' * 25
