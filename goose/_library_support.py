"""Shared planning and generation helpers for sequence libraries."""

from collections.abc import Mapping
from itertools import product
from math import ceil, floor, isfinite, ulp
from numbers import Integral, Real
import time

from tqdm.auto import tqdm

from goose.backend import goose_tools, parameters
from goose.goose_exceptions import GooseFail, GooseInputError


_EPS = 1e-9
_COMMON_KWARGS = {
    'strict_disorder', 'disorder_cutoff', 'cutoff', 'metapredict_version',
    'max_consecutive_ordered', 'max_total_ordered',
}
FRACTION_KWARGS = _COMMON_KWARGS | {
    'attempts', 'remaining_probabilities', 'max_aa_fractions', 'batch_size',
}
CLASS_KWARGS = _COMMON_KWARGS | {
    'num_attempts', 'remaining_probabilities', 'max_class_fractions',
}
DIMENSION_KWARGS = _COMMON_KWARGS | {
    'attempts', 'allowed_error', 'reduce_pos_charged', 'exclude_aas',
}


class LibraryProgress:
    """One tqdm bar per library, counting completed copies rather than retries."""

    def __init__(self, total, skipped, verbose, description):
        self._enabled = verbose
        self._counts = {'made': 0, 'failed': 0, 'duplicates': 0, 'skipped': skipped}
        self._target = None
        self._last_refresh = 0.0
        self._bar = tqdm(total=total, desc=description, unit='seq',
                         dynamic_ncols=True, disable=not verbose, postfix=self._counts)

    def __enter__(self):
        return self

    def start_copy(self, name):
        self._target = name
        if self._enabled:
            self._bar.set_postfix({**self._counts, 'target': name}, refresh=False)
            # Show which target is starting before expensive generation, while
            # avoiding a terminal/notebook redraw for every fast target.
            now = time.monotonic()
            if now - self._last_refresh >= 0.1:
                self._bar.refresh()
                self._last_refresh = now

    def finish_copy(self, made, failed, duplicates):
        self._counts.update(made=made, failed=failed, duplicates=duplicates)
        if self._enabled:
            self._bar.set_postfix({**self._counts, 'target': self._target}, refresh=False)
        self._bar.update(1)

    def __exit__(self, exc_type, exc_value, traceback):
        if self._enabled and exc_type is None:
            self._bar.set_postfix(self._counts, refresh=False)
        self._bar.close()


def number(value, name, minimum=None, maximum=None):
    """Validate a finite real number, excluding booleans."""
    if isinstance(value, bool) or not isinstance(value, Real) or not isfinite(value):
        raise GooseInputError(f'{name} must be a finite number, got {value!r}.')
    value = float(value)
    if minimum is not None and value < minimum:
        raise GooseInputError(f'{name} must be at least {minimum}, got {value}.')
    if maximum is not None and value > maximum:
        raise GooseInputError(f'{name} must be at most {maximum}, got {value}.')
    return value


def library_length(length):
    if isinstance(length, bool) or not isinstance(length, Integral):
        raise GooseInputError('length must be an integer.')
    length = int(length)
    goose_tools.length_check(length)
    return length


def validate_options(num_copies, require_unique_copies, on_generation_failure,
                     return_details, verbose, kwargs, allowed_kwargs):
    """Reject invalid controls before a potentially expensive generation loop."""
    if isinstance(num_copies, bool) or not isinstance(num_copies, Integral) or num_copies < 1:
        raise GooseInputError('num_copies must be an integer of at least 1.')
    for name, value in (('require_unique_copies', require_unique_copies),
                        ('return_details', return_details), ('verbose', verbose)):
        if not isinstance(value, bool):
            raise GooseInputError(f'{name} must be a boolean.')
    if on_generation_failure not in ('skip', 'raise'):
        raise GooseInputError("on_generation_failure must be 'skip' or 'raise'.")
    goose_tools.check_valid_kwargs(kwargs, allowed_kwargs)
    for name in ('attempts', 'num_attempts', 'batch_size', 'max_consecutive_ordered'):
        value = kwargs.get(name)
        if value is not None:
            minimum = 0 if name == 'max_consecutive_ordered' else 1
            if isinstance(value, bool) or not isinstance(value, Integral) or value < minimum:
                raise GooseInputError(f'{name} must be an integer of at least {minimum}.')
    for name in ('disorder_cutoff', 'cutoff', 'max_total_ordered'):
        if kwargs.get(name) is not None:
            number(kwargs[name], name, 0, 1)
    if 'allowed_error' in kwargs:
        number(kwargs['allowed_error'], 'allowed_error', 0)
    for name in ('strict_disorder', 'reduce_pos_charged'):
        if name in kwargs and not isinstance(kwargs[name], bool):
            raise GooseInputError(f'{name} must be a boolean.')
    if 'metapredict_version' in kwargs:
        version = kwargs['metapredict_version']
        if isinstance(version, bool) or version not in parameters.METAPREDICT_VERSIONS:
            raise GooseInputError(f'metapredict_version must be one of {parameters.METAPREDICT_VERSIONS}.')
    if kwargs.get('remaining_probabilities') is not None:
        if isinstance(kwargs['remaining_probabilities'], dict):
            for aa, probability in kwargs['remaining_probabilities'].items():
                number(probability, f'remaining_probabilities[{aa!r}]', 0, 1)
        goose_tools.handle_custom_probabilities(kwargs['remaining_probabilities'])
    if kwargs.get('exclude_aas') is not None:
        goose_tools.check_basic_parameters(exclude=kwargs['exclude_aas'])


def fraction_limits(defaults, overrides, name):
    limits = dict(defaults)
    if overrides is not None:
        if not isinstance(overrides, Mapping):
            raise GooseInputError(f'{name} must be a dictionary.')
        for key, value in overrides.items():
            if key not in limits:
                raise GooseInputError(f'Unknown key {key!r} in {name}. Options: {list(limits)}.')
            limits[key] = number(value, f'{name}[{key!r}]', 0, 1)
    return limits


def _fraction_counts(length, specification, increment, maximum, name):
    """Expand a range in integer counts, never by rounding repeated fractions."""
    if specification is None:
        start, end = 0.0, maximum
    elif isinstance(specification, Real) and not isinstance(specification, bool):
        start = end = number(specification, name, 0, maximum)
    elif isinstance(specification, (list, tuple)) and 1 <= len(specification) <= 3:
        start = number(specification[0], f'{name} minimum', 0, maximum)
        end = start if len(specification) == 1 else number(specification[1], f'{name} maximum', 0, maximum)
        if len(specification) == 3:
            increment = specification[2]
    else:
        raise GooseInputError(f'{name} must be a number or [minimum, maximum, optional increment].')
    if end < start:
        raise GooseInputError(f'{name} maximum must be at least its minimum.')

    if increment is None:
        # Approximately five percentage points, rounded to the nearest residue.
        count_step = max(1, floor(length * 0.05 + 0.5))
    else:
        increment = number(increment, f'{name} increment', 0, 1)
        if increment * length < 1 - _EPS:
            raise GooseInputError(f'{name} increment must be at least one residue '
                                  f'(1/{length} = {1 / length:g}).')
        # Do not produce spacing finer than requested.
        count_step = max(1, ceil(increment * length - _EPS))

    # Bounds move inward; an exact target must already be representable.
    first = max(0, ceil(start * length - _EPS))
    last = min(length, floor(end * length + _EPS))
    # The epsilon handles multiplication noise, but must not widen a user's
    # actual bounds or exceed the backend's strict maximum-fraction check.
    if start - first / length > 2 * ulp(start):
        first += 1
    if last / length - end > 2 * ulp(end):
        last -= 1
    if last / length > maximum:
        last -= 1
    if first > last:
        raise GooseInputError(f'{name} has no whole-residue target at length {length}. '
                              f'Use multiples of 1/{length}.')
    counts = list(range(first, last + 1, count_step))
    if counts[-1] != last:
        counts.append(last)
    return counts


def format_target(value):
    """Use round-trip precision so distinct numeric targets have distinct names."""
    return repr(float(value))


def enumerate_compositions(length, ranges, increment, selected, limits, combine, selection_name):
    length = library_length(length)
    if not isinstance(combine, bool):
        raise GooseInputError('combine must be a boolean.')
    if ranges is not None and (not isinstance(ranges, Mapping) or not ranges):
        raise GooseInputError('Fraction ranges must be a nonempty dictionary.')
    if combine and ranges is None and selected is None:
        raise GooseInputError(f'For combine=True, specify fraction ranges or {selection_name} '
                              'to select the axes to combine.')
    if selected is None:
        selected = list(limits if ranges is None else ranges)
    if not isinstance(selected, (list, tuple)) or not selected:
        raise GooseInputError(f'{selection_name} must select at least one valid entry.')
    for key in selected:
        if not isinstance(key, str) or key not in limits:
            raise GooseInputError(f'Unknown {selection_name} entry {key!r}. Options: {list(limits)}.')
    if len(set(selected)) != len(selected):
        raise GooseInputError(f'{selection_name} contains duplicate entries.')
    ranges = {} if ranges is None else ranges
    if any(key not in selected for key in ranges):
        raise GooseInputError(f'Every fraction range must also be selected in {selection_name}.')
    if increment is not None:
        number(increment, 'fraction increment', 0, 1)
        if increment * length < 1 - _EPS:
            raise GooseInputError(f'Fraction increment must be at least one residue (1/{length} = {1 / length:g}).')
    axes = {key: _fraction_counts(length, ranges.get(key), increment, limits[key], key)
            for key in selected}

    if combine:
        targets = (dict(zip(selected, counts)) for counts in product(*(axes[key] for key in selected)))
    else:
        targets = ({key: count} for key in selected for count in axes[key])

    possible, skipped = [], []
    for counts in targets:
        combination = {key: count / length for key, count in counts.items()}
        combination['counts'] = counts
        combination['name'] = '_'.join(f'{key}_{format_target(count / length)}' for key, count in counts.items())
        total = sum(counts.values())
        if total > length or (len(counts) == len(limits) and total != length):
            combination['reason'] = ('Specified fractions exceed 1.' if total > length else
                                     'All components are specified but their fractions do not sum to 1.')
            combination['reason_type'] = 'incompatible_fraction_sum'
            skipped.append(combination)
        else:
            possible.append(combination)
    return possible, skipped


def run_library(length, num_copies, possible, skipped, make_sequence,
                require_unique_copies, on_generation_failure, return_details,
                verbose, groups=None, description='Sequence library'):
    """Apply the property-library naming, copy, and failure conventions."""
    library, failed, duplicates = {}, [], []

    def generate(combination):
        sequence = make_sequence(combination)
        if not isinstance(sequence, str) or len(sequence) != length:
            raise GooseFail(f'Generator did not return a sequence of length {length}.')
        if groups is not None:
            for key, count in combination['counts'].items():
                observed = sum(sequence.count(aa) for aa in groups[key])
                if observed != count:
                    raise GooseFail(f'Expected {count} residues for {key}, got {observed}.')
        return sequence

    with LibraryProgress(len(possible) * num_copies, len(skipped), verbose, description) as progress:
        for combination in possible:
            seen = set()
            for copy_number in range(1, num_copies + 1):
                name = f'{combination["name"]}_rep_{copy_number}'
                progress.start_copy(name)
                try:
                    sequence = generate(combination)
                    if require_unique_copies:
                        for _ in range(4):
                            if sequence not in seen:
                                break
                            sequence = generate(combination)
                except GooseInputError:
                    # An invalid request is never a recoverable generation failure.
                    raise
                except Exception as error:
                    if on_generation_failure == 'raise':
                        raise
                    failed.append({**combination, 'name': name,
                                   'reason': f'{type(error).__name__}: {error}'})
                else:
                    if sequence in seen:
                        duplicates.append(name)
                    seen.add(sequence)
                    library[name] = sequence
                progress.finish_copy(len(library), len(failed), len(duplicates))
    if return_details:
        return {'library': library, 'skipped': skipped, 'failed': failed, 'duplicates': duplicates}
    return library
