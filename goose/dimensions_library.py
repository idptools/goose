"""Generate libraries spanning predicted Rg or end-to-end distance."""

from decimal import Decimal
from math import ulp

from goose import create
from goose.backend import parameters
from goose.goose_exceptions import GooseInputError
from goose._library_support import (
    DIMENSION_KWARGS, format_target, library_length, number, run_library,
    validate_options,
)

__all__ = ['generate_dimensions_library', 'enumerate_dimension_combinations']


def enumerate_dimension_combinations(length, rg_or_re='rg', min_rg=None, max_rg=None,
                                     min_re=None, max_re=None, dimension_increment=1.0):
    """Plan dimension targets in angstroms; return ``(possible, skipped)``.

    Omitted bounds use GOOSE's empirical length-dependent limits. The upper
    bound is included even when the interval does not divide the range.
    Entries contain ``rg`` or ``re`` and a ``name``. Invalid bounds raise;
    ``skipped`` is empty. No prediction models are run here.
    """
    length = library_length(length)
    if not isinstance(rg_or_re, str) or rg_or_re.lower() not in ('rg', 're'):
        raise GooseInputError("rg_or_re must be 'rg' or 're'.")
    dimension = rg_or_re.lower()
    if dimension == 'rg':
        if min_re is not None or max_re is not None:
            raise GooseInputError("Re bounds require rg_or_re='re'; only one dimension can be swept.")
        lower, upper = min_rg, max_rg
        possible_min, possible_max = parameters.get_min_rg(length), parameters.get_max_rg(length)
    else:
        if min_rg is not None or max_rg is not None:
            raise GooseInputError("Rg bounds require rg_or_re='rg'; only one dimension can be swept.")
        lower, upper = min_re, max_re
        possible_min, possible_max = parameters.get_min_re(length), parameters.get_max_re(length)
    lower = possible_min if lower is None else number(lower, f'min_{dimension}', possible_min, possible_max)
    upper = possible_max if upper is None else number(upper, f'max_{dimension}', possible_min, possible_max)
    if lower > upper:
        raise GooseInputError(f'min_{dimension} must not exceed max_{dimension}.')
    increment = number(dimension_increment, 'dimension_increment', 0)
    if increment == 0:
        raise GooseInputError('dimension_increment must be greater than zero.')
    if lower < upper and increment < ulp(upper):
        raise GooseInputError('dimension_increment is too small to distinguish numeric targets at these bounds.')

    # Decimal stepping prevents the final value from exceeding the backend's
    # strict empirical bounds due to binary floating-point accumulation.
    first, last, step = (Decimal(str(value)) for value in (lower, upper, increment))
    count = int((last - first) // step)
    values = [float(first + index * step) for index in range(count + 1)]
    if values[-1] != upper:
        values.append(upper)
    possible = [{dimension: value, 'name': f'{dimension}_{format_target(value)}'} for value in values]
    return possible, []


def generate_dimensions_library(length, num_copies=1, rg_or_re='rg', min_rg=None,
                                max_rg=None, min_re=None, max_re=None,
                                dimension_increment=1.0, require_unique_copies=True,
                                on_generation_failure='skip', return_details=False,
                                verbose=True, **kwargs):
    """Generate sequences spanning predicted Rg or end-to-end distance (Re).

    Parameters
    ----------
    length : int
        Length of every sequence.
    num_copies : int, default=1
        Number of independent sequences per target.
    rg_or_re : {'rg', 're'}, default='rg'
        Dimension to sweep. Case insensitive. Only one dimension is targeted.
    min_rg, max_rg : float, optional
        Rg bounds in angstroms, for ``rg_or_re='rg'``. Each omitted bound uses
        ``parameters.get_min_rg(length)`` or ``get_max_rg(length)``.
    min_re, max_re : float, optional
        Re bounds in angstroms, for ``rg_or_re='re'``. Each omitted bound uses
        ``parameters.get_min_re(length)`` or ``get_max_re(length)``.
    dimension_increment : float, default=1.0
        Positive spacing between targets in angstroms. The final upper bound
        is always included, even if the last interval is smaller.
    require_unique_copies : bool, default=True
        Retry duplicate copies up to four times per target; report any that
        remain in ``duplicates`` when returning details.
    on_generation_failure : {'skip', 'raise'}, default='skip'
        Record generation failures and continue, or propagate the exception.
        Invalid inputs always raise.
    return_details : bool, default=False
        Return ``library``, ``skipped``, ``failed``, and ``duplicates`` instead
        of just the name-to-sequence dictionary.
    verbose : bool, default=True
        Show a tqdm progress bar with the current target and generated, failed,
        duplicate, and skipped counts. False suppresses the bar.
    **kwargs
        Forwarded to ``create.seq_by_rg`` or ``create.seq_by_re``. Includes
        ``attempts`` (default 200), ``allowed_error`` (default 0.5 angstroms),
        disorder settings, ``reduce_pos_charged``, and ``exclude_aas``.

    Notes
    -----
    The defaults use GOOSE's empirical generation bounds, not rigorous
    physical extrema or guarantees that every target can be generated.
    Names encode target dimensions; actual predictions are within the
    underlying generator's ``allowed_error`` when generation succeeds.

    Examples
    --------
    >>> library = generate_dimensions_library(100, num_copies=2)
    >>> library = generate_dimensions_library(100, rg_or_re='re',
    ...                                       min_re=30, max_re=70, dimension_increment=5)
    """
    validate_options(num_copies, require_unique_copies, on_generation_failure,
                     return_details, verbose, kwargs, DIMENSION_KWARGS)
    length = library_length(length)
    possible, skipped = enumerate_dimension_combinations(
        length, rg_or_re, min_rg, max_rg, min_re, max_re, dimension_increment)
    dimension = rg_or_re.lower()
    options = {key: value for key, value in kwargs.items() if value is not None}
    if options.get('attempts') is None:
        options['attempts'] = parameters.DEFAULT_ATTEMPTS
    generator = create.seq_by_rg if dimension == 'rg' else create.seq_by_re

    def make_sequence(combination):
        return generator(length, combination[dimension], **options)

    return run_library(length, num_copies, possible, skipped, make_sequence,
                       require_unique_copies, on_generation_failure,
                       return_details, verbose, description=f'{dimension.capitalize()} library')
