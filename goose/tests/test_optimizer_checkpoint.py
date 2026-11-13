"""Regressions for returning a solution that satisfies the stopping tolerances."""

import importlib

import numpy as np
import pytest
import sparrow

from goose.backend.optimizer_properties import CustomProperty, FCR, Hydrophobicity, Kappa, NCPR
from goose.optimize import SequenceOptimizer


TARGET = 'VTKQQYSENRYLFSGIVWLMLQKRMGEVHSIQTWEERYCQYMIAKRIGRMNLCYICMGSKLQYCDSKLTNPNESYCVVVTVFTANCEVGT'
FEASIBLE = 'VTKDQYSENRYLFSGIVWLMLQKRMGEVHSIQTGEVRYCQYMIAKRIGRMNLCYICMGSKLQYCDSKFTNPNESYCVSVTVFTYNCEVGT'
INFEASIBLE = 'VTKQQYSENRYLFSGIVWLMLQKRMGTVHSIQTWEERYCQYMIAKRIGRMNLCYICMGSKLQMCDSKLTNANESYCVVVTVFTANCEVGT'


def profile_optimizer(**kwargs):
    optimizer = SequenceOptimizer(target_length=len(TARGET), verbose=False,
                                  enable_adaptive_scaling=False, **kwargs)
    for cls, window, tolerance in [(NCPR, 6, .07), (FCR, 6, .07),
                                   (Kappa, 15, .05), (Hydrophobicity, 8, .5)]:
        optimizer.add_property(cls, calculate_as_linear_profile=True,
                               target_sequence=TARGET, window_size=window,
                               tolerance=tolerance)
    return optimizer


def assert_independent_profile_errors_within_tolerance(optimizer, sequence):
    # Use direct Sparrow window evaluations, independent of the optimized paths.
    for prop in optimizer.properties:
        def profile(text):
            w = prop.window_size
            values = [prop.calculate_raw_value(sparrow.Protein(text[i:i + w]))
                      for i in range(len(text) - w + 1)]
            return np.array([values[0]] * (w // 2) + values
                            + [values[-1]] * (w - 1 - w // 2))
        error = np.sum(np.abs(profile(sequence) - profile(TARGET))) / len(sequence)
        assert error <= prop._tolerance


def test_final_restore_preserves_accepted_profile_solution(monkeypatch):
    optimizer = profile_optimizer(max_iterations=1)
    _, feasible_info = optimizer._evaluate_sequence(FEASIBLE)
    _, infeasible_info = optimizer._evaluate_sequence(INFEASIBLE)
    assert all(info['within_tolerance'] for info in feasible_info.values())
    assert not infeasible_info['Kappa']['within_tolerance']
    assert optimizer._weighted_raw_sum(infeasible_info, False) < optimizer._weighted_raw_sum(feasible_info, False)

    optimizer.set_initial_sequence(INFEASIBLE)
    monkeypatch.setattr(optimizer, '_generate_candidates', lambda *args: [FEASIBLE])
    assert optimizer.run() == FEASIBLE
    assert optimizer._check_error_tolerance()
    assert_independent_profile_errors_within_tolerance(optimizer, optimizer.best_sequence)


@pytest.mark.parametrize('early_stop', [True, False])
def test_initial_candidate_checkpoint_preserves_feasible_solution(monkeypatch, early_stop):
    optimizer = profile_optimizer(max_iterations=0, enable_error_tolerance=early_stop)
    module = importlib.import_module('goose.optimize')
    monkeypatch.setattr(module, 'build_diverse_initial_sequences',
                        lambda *args, **kwargs: [INFEASIBLE, FEASIBLE])
    assert optimizer.run() == FEASIBLE
    assert_independent_profile_errors_within_tolerance(optimizer, optimizer.best_sequence)


class LookupProperty(CustomProperty):
    """A small deterministic error landscape for testing candidate acceptance."""

    def __init__(self, index, **kwargs):
        super().__init__(target_value=0, **kwargs)
        self.index = index

    def calculate_raw_value(self, protein):
        return {'AA': (.06, .01), 'CC': (.04, .49), 'DD': (.055, .01)}[protein.sequence][self.index]


@pytest.mark.parametrize('early_stop', [True, False])
def test_relaxed_acceptance_preserves_a_solution_in_same_batch(monkeypatch, early_stop):
    optimizer = SequenceOptimizer(target_length=2, max_iterations=2, verbose=False,
                                  enable_adaptive_scaling=False, stagnation_threshold=1,
                                  enable_error_tolerance=early_stop)
    optimizer.add_property(LookupProperty, index=0, tolerance=.05)
    optimizer.add_property(LookupProperty, index=1, tolerance=.5)
    optimizer.set_initial_sequence('AA')
    monkeypatch.setattr(optimizer, '_normalization_initialized', True)
    monkeypatch.setattr(optimizer, '_generate_candidates', lambda *args: ['CC', 'DD'])
    assert optimizer.run() == 'CC'
    # Either stop immediately or restore the feasible checkpoint after exploration.
    assert optimizer.iteration == (1 if early_stop else 2)
    assert optimizer._all_properties_within_tolerance(optimizer._get_best_sequence_property_info())


def test_checkpoint_uses_global_tolerance_for_properties_without_local_tolerance(monkeypatch):
    class GlobalToleranceError(LookupProperty):
        def calculate_raw_value(self, protein):
            return {'AA': (.002, .01), 'CC': (.0009, .49)}[protein.sequence][self.index]

    optimizer = SequenceOptimizer(target_length=2, max_iterations=1, verbose=False,
                                  error_tolerance=.001, enable_adaptive_scaling=False)
    optimizer.add_property(GlobalToleranceError, index=0)
    optimizer.add_property(GlobalToleranceError, index=1, tolerance=.5)
    optimizer.set_initial_sequence('AA')
    monkeypatch.setattr(optimizer, '_normalization_initialized', True)
    monkeypatch.setattr(optimizer, '_generate_candidates', lambda *args: ['CC'])
    assert optimizer.run() == 'CC'
    assert optimizer._check_error_tolerance()


def test_profile_report_describes_profile_target(monkeypatch, capsys):
    optimizer = SequenceOptimizer(target_length=len(TARGET), max_iterations=0, verbose=True)
    optimizer.add_property(Kappa, calculate_as_linear_profile=True,
                           target_sequence=TARGET, window_size=15, tolerance=.05)
    optimizer.set_initial_sequence(TARGET)
    monkeypatch.setattr(optimizer, '_calculate_initial_normalization', lambda *args: None)
    optimizer.run()
    report = capsys.readouterr().err
    assert 'target: exact profile' in report
    assert 'mean profile value' in report
    assert 'profile_error: 0.0000' in report
    assert 'target: exact value of 0.000' not in report
