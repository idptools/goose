"""Generate libraries spanning individual amino-acid fractions."""

from goose import create
from goose.backend import parameters
from goose._library_support import (
    FRACTION_KWARGS, enumerate_compositions, fraction_limits, library_length,
    run_library, validate_options,
)

__all__ = ['generate_fractions_library', 'enumerate_fraction_combinations']


def enumerate_fraction_combinations(length, fraction_ranges=None, fraction_increment=None,
                                    amino_acids=None, combine=False, max_aa_fractions=None):
    """Plan fraction targets without generating sequences.

    Return ``(possible, skipped)`` lists. Each possible entry contains the
    targeted amino-acid fractions, exact integer ``counts``, and a ``name``.
    See :func:`generate_fractions_library` for range and increment semantics.
    """
    limits = fraction_limits(parameters.MAX_FRACTION_DICT, max_aa_fractions, 'max_aa_fractions')
    if isinstance(amino_acids, str):
        amino_acids = list(amino_acids)
    return enumerate_compositions(length, fraction_ranges, fraction_increment,
                                  amino_acids, limits, combine, 'amino_acids')


def generate_fractions_library(length, num_copies=1, fraction_ranges=None,
                               fraction_increment=None, amino_acids=None, combine=False,
                               require_unique_copies=True, on_generation_failure='skip',
                               return_details=False, verbose=True, **kwargs):
    """Generate disordered sequences spanning amino-acid fractions.

    Parameters
    ----------
    length : int
        Length of every sequence.
    num_copies : int, default=1
        Independent generation calls per target. Unavoidable duplicates are
        retained and reported, including copies of a homopolymer.
    fraction_ranges : dict, optional
        Map amino-acid letters to a fixed fraction, ``[minimum, maximum]``,
        or ``[minimum, maximum, increment]``. When omitted, sweep every amino
        acid independently from zero to its GOOSE maximum. Providing a dict
        selects just its keys unless ``amino_acids`` is also given.
    fraction_increment : float, optional
        Default spacing is approximately 0.05, rounded to the nearest whole
        residue, with a minimum of one residue. At length 25 it is 0.04; at
        length 100 it is 0.05. Explicit increments must be at least 1/length
        and are rounded UP to a whole-residue step. A third range element
        overrides this spacing. The final attainable upper bound is included.
    amino_acids : str or list of str, optional
        Select residues, e.g. ``'DEK'``. Selected residues without a supplied
        range use zero through their GOOSE maximum.
    combine : bool, default=False
        Sweep each selected residue separately. Set True to generate the
        Cartesian product instead; targets summing above one are skipped.
        Explicitly select the axes when combining.
    require_unique_copies : bool, default=True
        Retry duplicate copies up to four times for each target.
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
        Forwarded to ``create.seq_by_fractions``, excluding amino-acid targets
        and ``return_all_sequences``. Includes ``attempts`` (default 200),
        disorder settings, ``remaining_probabilities``, ``batch_size``, and
        ``max_aa_fractions`` overrides, which also change default sweep limits.

    Notes
    -----
    Range bounds move inward to multiples of 1/length. Fixed values must
    already be whole-residue fractions. Unspecified residues remain free;
    they are not constrained to zero. Default maximum fractions are GOOSE's
    empirical generation limits, not guarantees of successful generation.

    Examples
    --------
    >>> library = generate_fractions_library(25, num_copies=2)
    >>> library = generate_fractions_library(100, fraction_ranges={'E': [0, 0.5, 0.1]})
    """
    validate_options(num_copies, require_unique_copies, on_generation_failure,
                     return_details, verbose, kwargs, FRACTION_KWARGS)
    length = library_length(length)
    possible, skipped = enumerate_fraction_combinations(
        length, fraction_ranges, fraction_increment, amino_acids, combine,
        kwargs.get('max_aa_fractions'))
    options = {key: value for key, value in kwargs.items() if value is not None}
    if options.get('attempts') is None:
        options['attempts'] = parameters.DEFAULT_ATTEMPTS

    def make_sequence(combination):
        targets = {aa: combination[aa] for aa in combination['counts']}
        return create.seq_by_fractions(length, **options, **targets)

    return run_library(length, num_copies, possible, skipped, make_sequence,
                       require_unique_copies, on_generation_failure,
                       return_details, verbose,
                       groups={aa: [aa] for aa in parameters.VALID_AMINO_ACIDS},
                       description='Fractions library')
