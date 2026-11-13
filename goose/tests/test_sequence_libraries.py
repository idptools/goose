"""Library grids, whole-residue composition targets, and generation behavior."""

from math import floor

import pytest

from goose import create, property_library
from goose.backend import parameters
from goose.class_library import enumerate_class_combinations, generate_class_library
from goose.data.defined_aa_classes import aa_classes
from goose.dimensions_library import enumerate_dimension_combinations, generate_dimensions_library
from goose.fractions_library import enumerate_fraction_combinations, generate_fractions_library
from goose.goose_exceptions import GooseFail, GooseInputError


@pytest.mark.parametrize('length', [10, 25, 37, 100, 101])
@pytest.mark.parametrize('enumerator,limits', [
    (enumerate_fraction_combinations, parameters.MAX_FRACTION_DICT),
    (enumerate_class_combinations, parameters.MAX_CLASS_FRACTIONS),
])
def test_default_composition_sweeps_are_individual_and_attainable(length, enumerator, limits):
    possible, skipped = enumerator(length)
    assert not skipped
    assert {next(iter(target['counts'])) for target in possible} == set(limits)
    assert len({target['name'] for target in possible}) == len(possible)
    for key, maximum in limits.items():
        targets = [target for target in possible if key in target['counts']]
        counts = [target['counts'][key] for target in targets]
        assert counts[0] == 0
        assert counts[-1] == floor(maximum * length + 1e-9)
        assert len(counts) == len(set(counts))
        expected_step = max(1, floor(length * 0.05 + 0.5))
        assert all(right - left == expected_step for left, right in zip(counts[:-2], counts[1:-1]))
        for target in targets:
            assert len(target['counts']) == 1
            assert target[key] == target['counts'][key] / length
            assert target[key] <= maximum


@pytest.mark.parametrize('enumerator,ranges_arg,key', [
    (enumerate_fraction_combinations, 'fraction_ranges', 'E'),
    (enumerate_class_combinations, 'class_ranges', 'positive'),
])
def test_fraction_bounds_move_inward_and_step_rounds_up(enumerator, ranges_arg, key):
    possible, _ = enumerator(25, **{ranges_arg: {key: [0.01, 0.30, 0.1]}})
    assert [target['counts'][key] for target in possible] == [1, 4, 7]
    possible, _ = enumerator(25, **{ranges_arg: {key: [0, 0.20, 0.1]}})
    assert [target['counts'][key] for target in possible] == [0, 3, 5]


def test_ranges_select_axes_and_third_value_overrides_global_increment():
    possible, _ = enumerate_fraction_combinations(
        25, fraction_ranges={'E': [0, 0.12, 0.04]}, fraction_increment=0.2)
    assert [target['E'] for target in possible] == [0, 0.04, 0.08, 0.12]
    assert all(set(target['counts']) == {'E'} for target in possible)
    possible, _ = enumerate_class_combinations(25, classes='polar')
    assert all(set(target['counts']) == {'polar'} for target in possible)
    possible, _ = enumerate_fraction_combinations(25, amino_acids='DE')
    assert {next(iter(target['counts'])) for target in possible} == {'D', 'E'}


def test_maximum_overrides_control_defaults_without_mutating_global_limits():
    overrides = {'I': 0.8}
    original = dict(parameters.MAX_FRACTION_DICT)
    possible, _ = enumerate_fraction_combinations(25, amino_acids='I', max_aa_fractions=overrides)
    assert possible[-1]['counts'] == {'I': 20}
    assert parameters.MAX_FRACTION_DICT == original
    assert overrides == {'I': 0.8}
    possible, _ = enumerate_class_combinations(25, classes='aliphatic', max_class_fractions={'aliphatic': 0.4})
    assert possible[-1]['aliphatic'] == 0.4


def test_fraction_bounds_remain_inside_near_integer_boundaries():
    possible, _ = enumerate_fraction_combinations(25, amino_acids='E', max_aa_fractions={'E': 0.2 - 1e-12})
    assert possible[-1]['E'] == 0.16
    possible, _ = enumerate_fraction_combinations(25, fraction_ranges={'E': [0.2 + 1e-12, 0.28]})
    assert possible[0]['E'] == 0.24
    possible, _ = enumerate_fraction_combinations(100, fraction_ranges={'E': 0.1 + 0.2})
    assert possible[0]['counts'] == {'E': 30}


@pytest.mark.parametrize('enumerator,ranges_arg,keys', [
    (enumerate_fraction_combinations, 'fraction_ranges', ('D', 'E')),
    (enumerate_class_combinations, 'class_ranges', ('positive', 'negative')),
])
def test_combinations_are_opt_in_and_invalid_sums_are_skipped(enumerator, ranges_arg, keys):
    ranges = {key: [0, 1, 0.5] for key in keys}
    independent, _ = enumerator(100, **{ranges_arg: ranges})
    assert len(independent) == 6
    possible, skipped = enumerator(100, combine=True, **{ranges_arg: ranges})
    assert len(possible) == 6
    assert len(skipped) == 3
    assert all(sum(target['counts'].values()) <= 100 for target in possible)
    assert all(target['reason_type'] == 'incompatible_fraction_sum' for target in skipped)


def test_fully_specified_classes_need_to_fill_the_sequence():
    possible, skipped = enumerate_class_combinations(25, {key: 0 for key in aa_classes}, combine=True)
    assert not possible
    assert len(skipped) == 1


@pytest.mark.parametrize('kwargs', [
    {'fraction_increment': 0.001}, {'fraction_increment': 0},
    {'fraction_increment': float('nan')}, {'fraction_increment': True},
    {'fraction_ranges': {'E': 0.1}},  # 2.5 residues
    {'fraction_ranges': {'E': [0.01, 0.03]}},  # no attainable target
    {'fraction_ranges': {'I': [0, 1]}},  # above GOOSE maximum
    {'fraction_ranges': {'E': [0.4, 0.2]}},
    {'fraction_ranges': {'E': [0, 1, 0.001]}},
    {'fraction_ranges': {'E': [0, float('inf')]}},
    {'fraction_ranges': {'Z': [0, 1]}}, {'fraction_ranges': {}},
    {'fraction_ranges': []}, {'amino_acids': []}, {'amino_acids': 'EE'},
    {'amino_acids': 'D', 'fraction_ranges': {'E': 0}}, {'combine': 1}, {'combine': True},
    {'max_aa_fractions': {'Z': 0.5}}, {'max_aa_fractions': {'E': float('nan')}},
])
def test_invalid_fraction_plans_raise(kwargs):
    with pytest.raises(GooseInputError):
        enumerate_fraction_combinations(25, **kwargs)


@pytest.mark.parametrize('kwargs', [
    {'class_increment': 0.001}, {'class_ranges': {'polar': 0.1}},
    {'classes': 'charged'}, {'classes': ['positive', 'positive']},
    {'max_class_fractions': {'positive': -0.1}},
])
def test_invalid_class_plans_raise(kwargs):
    with pytest.raises(GooseInputError):
        enumerate_class_combinations(25, **kwargs)


@pytest.mark.parametrize('dimension', ['rg', 're'])
def test_dimension_defaults_use_exact_existing_bounds(dimension):
    possible, skipped = enumerate_dimension_combinations(25, rg_or_re=dimension, dimension_increment=1.3)
    assert not skipped
    minimum = getattr(parameters, f'get_min_{dimension}')(25)
    maximum = getattr(parameters, f'get_max_{dimension}')(25)
    values = [target[dimension] for target in possible]
    assert values[0] == minimum
    assert values[-1] == maximum
    assert values == sorted(set(values))
    assert all(minimum <= value <= maximum for value in values)
    assert all(right - left == pytest.approx(1.3) for left, right in zip(values[:-2], values[1:-1]))


def test_dimension_bounds_may_be_omitted_individually_or_equal():
    targets, _ = enumerate_dimension_combinations(100, min_rg=20, dimension_increment=5)
    assert targets[0]['rg'] == 20
    assert targets[-1]['rg'] == parameters.get_max_rg(100)
    targets, _ = enumerate_dimension_combinations(100, max_re=60, rg_or_re='Re')
    assert targets[0]['re'] == parameters.get_min_re(100)
    assert targets[-1]['re'] == 60
    targets, _ = enumerate_dimension_combinations(100, min_rg=20, max_rg=20)
    assert targets == [{'rg': 20.0, 'name': 'rg_20.0'}]


def test_dimension_endpoint_is_included_and_close_targets_have_unique_names():
    targets, _ = enumerate_dimension_combinations(100, min_rg=20, max_rg=25, dimension_increment=2)
    assert [target['rg'] for target in targets] == [20, 22, 24, 25]
    targets, _ = enumerate_dimension_combinations(100, min_rg=20, max_rg=20.000000000003,
                                                 dimension_increment=1e-12)
    assert len({target['name'] for target in targets}) == len(targets)


@pytest.mark.parametrize('kwargs', [
    {'rg_or_re': 'both'}, {'rg_or_re': None}, {'min_re': 20},
    {'rg_or_re': 're', 'max_rg': 30}, {'min_rg': 0}, {'max_rg': 100},
    {'min_rg': 30, 'max_rg': 20}, {'min_rg': float('nan')},
    {'dimension_increment': 0}, {'dimension_increment': -1},
    {'dimension_increment': float('inf')}, {'dimension_increment': True},
    {'dimension_increment': 1e-30},
])
def test_invalid_dimension_plans_raise(kwargs):
    with pytest.raises(GooseInputError):
        enumerate_dimension_combinations(100, **kwargs)


@pytest.mark.parametrize('generator', [generate_fractions_library, generate_class_library, generate_dimensions_library])
@pytest.mark.parametrize('length', [False, 9, 10001, 25.5, 'E' * 25])
def test_library_length_requires_an_integer_in_supported_range(generator, length):
    with pytest.raises(GooseInputError):
        generator(length, verbose=False)


@pytest.mark.parametrize('kwargs', [
    {'num_copies': 0}, {'num_copies': True}, {'num_copies': 1.5},
    {'require_unique_copies': 'yes'}, {'on_generation_failure': 'ignore'},
    {'return_details': 1}, {'return_all_sequences': True}, {'verbose': 0},
    {'attempts': 0}, {'disorder_cutoff': float('nan')},
    {'remaining_probabilities': {'E': float('nan')}},
])
def test_invalid_generation_controls_fail_before_any_generation(monkeypatch, kwargs):
    def unexpected(*args, **kw):
        pytest.fail('Invalid input reached sequence generation')
    monkeypatch.setattr(create, 'seq_by_fractions', unexpected)
    with pytest.raises(GooseInputError):
        generate_fractions_library(25, **kwargs)


def test_fraction_generation_passes_only_one_target_and_validates_counts(monkeypatch):
    calls = []

    def generate(length, **kwargs):
        calls.append(kwargs)
        count = round(kwargs['E'] * length)
        return 'E' * count + 'G' * (length - count)

    monkeypatch.setattr(create, 'seq_by_fractions', generate)
    details = generate_fractions_library(25, num_copies=2, fraction_ranges={'E': [0, 0.08]},
                                          require_unique_copies=False, return_details=True, verbose=False)
    assert len(details['library']) == 6
    assert len(calls) == 6
    assert all(set(call) == {'E', 'attempts'} for call in calls)
    assert all(call['attempts'] == parameters.DEFAULT_ATTEMPTS for call in calls)
    assert details['duplicates'] == ['E_0.0_rep_2', 'E_0.04_rep_2', 'E_0.08_rep_2']


def test_class_generation_dispatch_and_passthrough(monkeypatch):
    calls = []

    def generate(length, **kwargs):
        calls.append(kwargs)
        return 'K' * round(kwargs['positive'] * length) + 'G' * (length - round(kwargs['positive'] * length))

    monkeypatch.setattr(create, 'seq_by_classes', generate)
    library = generate_class_library(25, class_ranges={'positive': [0, 0.08]},
                                     num_attempts=7, disorder_cutoff=0.3, verbose=False)
    assert len(library) == 3
    assert all(call['num_attempts'] == 7 and call['disorder_cutoff'] == 0.3 for call in calls)


@pytest.mark.parametrize('dimension', ['rg', 're'])
def test_dimension_generation_dispatch_and_copy_names(monkeypatch, dimension):
    calls = []

    def generate(length, objective, **kwargs):
        calls.append((objective, kwargs))
        return 'E' * length

    monkeypatch.setattr(create, f'seq_by_{dimension}', generate)
    other = 're' if dimension == 'rg' else 'rg'
    monkeypatch.setattr(create, f'seq_by_{other}', lambda *args, **kwargs: pytest.fail('Wrong dimension'))
    bounds = {f'min_{dimension}': 20, f'max_{dimension}': 25}
    library = generate_dimensions_library(100, 2, rg_or_re=dimension, dimension_increment=2,
                                          require_unique_copies=False, allowed_error=0.2,
                                          verbose=False, **bounds)
    assert len(library) == 8
    assert list(library)[-1] == f'{dimension}_25.0_rep_2'
    assert [objective for objective, _ in calls] == [20, 20, 22, 22, 24, 24, 25, 25]
    assert all(options == {'allowed_error': 0.2, 'attempts': parameters.DEFAULT_ATTEMPTS}
               for _, options in calls)


def test_failed_generations_and_mismatched_counts_are_recorded_or_raised(monkeypatch):
    def failing(*args, **kwargs):
        raise GooseFail('no candidate passed disorder filtering')

    monkeypatch.setattr(create, 'seq_by_fractions', failing)
    details = generate_fractions_library(25, fraction_ranges={'E': 0.2},
                                          return_details=True, verbose=False)
    assert not details['library']
    assert details['failed'][0]['counts'] == {'E': 5}
    assert 'no candidate passed' in details['failed'][0]['reason']
    with pytest.raises(GooseFail):
        generate_fractions_library(25, fraction_ranges={'E': 0.2},
                                   on_generation_failure='raise', verbose=False)
    monkeypatch.setattr(create, 'seq_by_fractions', lambda *args, **kwargs: 'G' * 25)
    details = generate_fractions_library(25, fraction_ranges={'E': 0.2},
                                          return_details=True, verbose=False)
    assert not details['library']
    assert 'Expected 5 residues for E, got 0' in details['failed'][0]['reason']


def test_input_errors_from_generators_are_not_silently_skipped(monkeypatch):
    def invalid(*args, **kwargs):
        raise GooseInputError('unsupported request')

    monkeypatch.setattr(create, 'seq_by_classes', invalid)
    with pytest.raises(GooseInputError, match='unsupported request'):
        generate_class_library(25, class_ranges={'polar': 0.2}, verbose=False)


def test_unique_copies_retry_and_report_unavoidable_homopolymer_duplicates(monkeypatch):
    calls = []
    monkeypatch.setattr(create, 'seq_by_fractions', lambda *args, **kwargs: calls.append(1) or 'E' * 25)
    details = generate_fractions_library(25, 2, fraction_ranges={'E': 1},
                                          return_details=True, verbose=False)
    assert len(details['library']) == 2
    assert details['duplicates'] == ['E_1.0_rep_2']
    assert len(calls) == 6


def test_public_functions_are_also_available_from_property_library():
    assert property_library.generate_fractions_library is generate_fractions_library
    assert property_library.generate_dimensions_library is generate_dimensions_library
    assert property_library.generate_class_library is generate_class_library
    possible, skipped = property_library.enumerate_property_combinations(100, FCR_range=0.2)
    assert len(possible) == 1 and not skipped


@pytest.mark.parametrize('length', [25, 37])
def test_real_fraction_and_class_backends_preserve_library_counts(length):
    # Exercise the actual allocation and disorder-prediction paths. A zero
    # cutoff isolates composition correctness from stochastic feasibility.
    step = 1 / length
    fraction_library = generate_fractions_library(
        length, fraction_ranges={'E': [0, 2 * step, step]}, attempts=2,
        disorder_cutoff=0, on_generation_failure='raise', verbose=False)
    class_library = generate_class_library(
        length, class_ranges={'aromatic': [0, 2 * step, step]}, num_attempts=2,
        disorder_cutoff=0, on_generation_failure='raise', verbose=False)
    assert [sequence.count('E') for sequence in fraction_library.values()] == [0, 1, 2]
    assert [sum(sequence.count(aa) for aa in aa_classes['aromatic'])
            for sequence in class_library.values()] == [0, 1, 2]


@pytest.mark.parametrize('fractions', [{'positive': 1}, {'glycine': 1}, {'positive': 0.4, 'negative': 0.6}])
def test_real_class_backend_supports_targets_that_fill_the_entire_sequence(fractions):
    library = generate_class_library(25, class_ranges=fractions, combine=True,
                                     num_attempts=2, disorder_cutoff=0,
                                     on_generation_failure='raise', verbose=False)
    assert len(library) == 1
    sequence = next(iter(library.values()))
    for key, fraction in fractions.items():
        assert sum(sequence.count(aa) for aa in aa_classes[key]) == round(fraction * 25)
