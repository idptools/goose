"""Generate libraries spanning amino-acid class fractions."""

from goose import create
from goose.backend import parameters
from goose.data.defined_aa_classes import aa_classes
from goose._library_support import (
    CLASS_KWARGS, enumerate_compositions, fraction_limits, library_length,
    run_library, validate_options,
)

__all__ = ['generate_class_library', 'enumerate_class_combinations']


def enumerate_class_combinations(length, class_ranges=None, class_increment=None,
                                 classes=None, combine=False, max_class_fractions=None):
    """Plan class targets without generation; return ``(possible, skipped)``.

    Entries contain targeted class fractions, exact integer ``counts``, and
    a ``name``. See :func:`generate_class_library` for range semantics.
    """
    limits = fraction_limits(parameters.MAX_CLASS_FRACTIONS, max_class_fractions, 'max_class_fractions')
    if isinstance(classes, str):
        classes = [classes]
    return enumerate_compositions(length, class_ranges, class_increment,
                                  classes, limits, combine, 'classes')


def generate_class_library(length, num_copies=1, class_ranges=None,
                           class_increment=None, classes=None, combine=False,
                           require_unique_copies=True, on_generation_failure='skip',
                           return_details=False, verbose=True, **kwargs):
    """Generate disordered sequences spanning amino-acid class fractions.

    Parameters
    ----------
    length : int
        Length of every sequence.
    num_copies : int, default=1
        Number of copies per target; unavoidable duplicates are reported.
    class_ranges : dict, optional
        Map class names to a fixed fraction, ``[minimum, maximum]``, or
        ``[minimum, maximum, increment]``. Without ranges or a class selection,
        sweep all nine GOOSE classes independently from zero to their GOOSE
        maxima. Providing ranges selects just those classes by default.
    class_increment : float, optional
        Approximately 0.05, rounded to the nearest residue, with a minimum of
        one residue. Explicit increments must be at least 1/length and round
        up to whole residues. A third range element overrides this spacing.
        Bounds move inward to attainable fractions; the upper bound is always
        included. A fixed value must already be attainable at this length.
    classes : str or list of str, optional
        Class names: ``aromatic`` (FWY), ``aliphatic`` (AILMV), ``polar``
        (STNQ), ``positive`` (KR), ``negative`` (DE), ``glycine``, ``proline``,
        ``cysteine``, and ``histidine``. Unspecified classes remain free.
    combine : bool, default=False
        Sweep classes individually. True combines the selected axes and skips
        incompatible sums. Select axes explicitly when combining.
    require_unique_copies : bool, default=True
        Retry duplicate copies up to four times for each target.
    on_generation_failure : {'skip', 'raise'}, default='skip'
        Record generation failures or propagate them. Invalid inputs raise.
    return_details : bool, default=False
        Return a dict with ``library``, ``skipped``, ``failed``, and
        ``duplicates``; otherwise return a name-to-sequence dictionary.
    verbose : bool, default=True
        Show a tqdm progress bar with the current target and generated, failed,
        duplicate, and skipped counts. False suppresses the bar.
    **kwargs
        Forwarded to ``create.seq_by_classes``. Includes ``num_attempts``
        (default 200), disorder settings, ``remaining_probabilities``, and
        ``max_class_fractions`` overrides, which also change default limits.

    Notes
    -----
    Targets specify total class content, while residues within each class are
    chosen by the existing class generator. Empirical maxima do not guarantee
    that every target can pass disorder filtering at every length.

    Examples
    --------
    >>> library = generate_class_library(25, num_copies=2)
    >>> library = generate_class_library(100, classes=['aromatic', 'positive'])
    >>> library = generate_class_library(100, class_ranges={'polar': [0.1, 0.5, 0.1]})
    """
    validate_options(num_copies, require_unique_copies, on_generation_failure,
                     return_details, verbose, kwargs, CLASS_KWARGS)
    length = library_length(length)
    possible, skipped = enumerate_class_combinations(
        length, class_ranges, class_increment, classes, combine,
        kwargs.get('max_class_fractions'))
    options = {key: value for key, value in kwargs.items() if value is not None}
    if options.get('num_attempts') is None:
        options['num_attempts'] = parameters.DEFAULT_ATTEMPTS

    def make_sequence(combination):
        targets = {cls: combination[cls] for cls in combination['counts']}
        return create.seq_by_classes(length, **options, **targets)

    return run_library(length, num_copies, possible, skipped, make_sequence,
                       require_unique_copies, on_generation_failure,
                       return_details, verbose, groups=aa_classes,
                       description='Class library')
