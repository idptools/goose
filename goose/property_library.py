##
## property_library.py
##
## User-facing functionality for generating libraries of synthetic IDRs that
## titrate across user-specified ranges of hydropathy, FCR, NCPR, and kappa.
##
## The core idea: the user specifies a range + increment for each property they
## care about, and GOOSE walks the resulting grid of parameter combinations.
## Many combinations in that grid are not possible to make (for example, an
## |NCPR| greater than the FCR, a hydropathy value that cannot be reached at a
## given FCR, or kappa when there are effectively no oppositely charged
## residues to pattern). Those combinations are identified and dropped *before*
## any sequence generation is attempted, and everything that survives is made
## the requested number of times.
##

__all__ = ['generate_property_library', 'enumerate_property_combinations']

import time

from goose import create
from goose.goose_exceptions import GooseInputError
from goose.backend import goose_tools
from goose.backend import parameters
from goose.data.amino_acids import amino_acid_dat


# order the properties are iterated over and the order they appear in names.
_PROPERTY_ORDER = ('hydropathy', 'FCR', 'NCPR', 'kappa')

# what each property is called in the sequence names.
_NAME_PREFIXES = {'hydropathy': 'hydro', 'FCR': 'FCR', 'NCPR': 'NCPR', 'kappa': 'kappa'}

# minimum number of decimal places to use for each property in sequence names.
# FCR, NCPR, and kappa are fractions, so they always get at least one decimal
# place (an NCPR of 0 is named 'NCPR_0.0'). Hydropathy is on a 0 - 9 scale, so
# it gets only the decimal places it needs (a hydropathy of 2 is 'hydro_2').
_MINIMUM_DECIMALS = {'hydropathy': 0, 'FCR': 1, 'NCPR': 1, 'kappa': 1}

# kwargs that are simply handed off to create.sequence() for every sequence.
_SEQUENCE_KWARGS = ['attempts', 'disorder_cutoff', 'exclude', 'strict_disorder',
                    'use_weighted_probabilities', 'custom_probabilities',
                    'metapredict_version', 'max_consecutive_ordered',
                    'max_total_ordered', 'batch_size', 'hydropathy_tolerance',
                    'kappa_tolerance']

# tolerance for float comparisons. Needed because things like 0.3-0.1 != 0.2
# in floating point, and we do a lot of comparisons on values built up by
# repeatedly adding an increment.
_EPS = 1e-9

# how often (in seconds) to report progress while a library is being made.
_SECONDS_BETWEEN_PROGRESS_REPORTS = 30


#-\|/-\|/-\|/-\|/-\|/-\|/-\|/-\|/-\|/-\|/-\|/-\|/-\|/-\|/-\|/-\|/-\|/-\|/-\|/-\|/-\|/-
#-\|/-\|/-\|/-\|/-\|/-\|/-\|/-\|/-\|/             \|/-\|/-\|/-\|/-\|/-\|/-\|/-\|/-\|/-
#-\|/-\|/-\|/-\|/-\|/-\|/-\|/-\|/-\|/  Helper     \|/-\|/-\|/-\|/-\|/-\|/-\|/-\|/-\|/-
#-\|/-\|/-\|/-\|/-\|/-\|/-\|/-\|/-\|/  functions  \|/-\|/-\|/-\|/-\|/-\|/-\|/-\|/-\|/-
#-\|/-\|/-\|/-\|/-\|/-\|/-\|/-\|/-\|/             \|/-\|/-\|/-\|/-\|/-\|/-\|/-\|/-\|/-
#-\|/-\|/-\|/-\|/-\|/-\|/-\|/-\|/-\|/-\|/-\|/-\|/-\|/-\|/-\|/-\|/-\|/-\|/-\|/-\|/-\|/-


def _decimal_places(value):
    '''
    Determine how many decimal places are needed to write a number exactly.
    Used so that the values in a titration - and the names built from them -
    are formatted consistently (an increment of 0.1 gives names like
    'NCPR_0.0' whereas an increment of 1 gives names like 'hydro_2').

    Parameters
    ----------
    value : float or int
        The value to examine.

    Returns
    -------
    int
        Number of decimal places needed to represent value.
    '''
    string_value = repr(float(value))
    # very small values can come back in scientific notation, so reformat those.
    if 'e' in string_value or 'E' in string_value:
        string_value = f'{float(value):.12f}'
    if '.' not in string_value:
        return 0
    return len(string_value.split('.')[1].rstrip('0'))


def _format_value(value, decimals):
    '''
    Format a property value for use in a sequence name.

    Parameters
    ----------
    value : float or int
        The value to format.

    decimals : int
        Number of decimal places to use.

    Returns
    -------
    str
        The formatted value.
    '''
    # get rid of negative zero so we don't end up with names like 'NCPR_-0.0'
    if value == 0:
        value = 0.0
    return f'{value:.{decimals}f}'


def _expand_property_values(property_range, increment, property_name):
    '''
    Turn a user-specified range + increment into the explicit list of values
    that will be titrated over for a single property.

    Parameters
    ----------
    property_range : None, float, int, or list/tuple
        None means the property is not specified at all (GOOSE will not be
        told anything about it). A single number means the property is held
        constant at that value. A two element list/tuple is the
        [start, end] of the range to titrate over. A three element
        list/tuple is [start, end, increment], where the third value
        overrides the increment argument.

    increment : float, int, or None
        The step size to use between start and end. Required if a range
        was specified as a two element list/tuple.

    property_name : str
        Name of the property. Only used for error messages.

    Returns
    -------
    tuple
        (list of values, number of decimal places to use in names). The list
        is [None] if the property was not specified.

    Raises
    ------
    GooseInputError
        If the range or increment is malformed.
    '''
    # property not specified, so nothing to titrate over.
    if property_range is None:
        return [None], 0

    # a single value means hold the property constant.
    if isinstance(property_range, (int, float)) and not isinstance(property_range, bool):
        return [float(property_range)], _decimal_places(property_range)

    if not isinstance(property_range, (list, tuple)):
        raise GooseInputError(f'The range specified for {property_name} must be a single '
                              f'number or a list/tuple, but a {type(property_range)} was passed.')

    # make sure everything in the range is a number.
    for value in property_range:
        if isinstance(value, bool) or not isinstance(value, (int, float)):
            raise GooseInputError(f'The value {value} specified in the range for '
                                  f'{property_name} is not a number.')

    # a single value in a list is the same as holding the property constant.
    if len(property_range) == 1:
        return [float(property_range[0])], _decimal_places(property_range[0])

    if len(property_range) == 2:
        start, end = float(property_range[0]), float(property_range[1])
    elif len(property_range) == 3:
        # a third value in the range overrides the increment argument.
        start, end = float(property_range[0]), float(property_range[1])
        increment = property_range[2]
    else:
        raise GooseInputError(f'The range specified for {property_name} can hold at most '
                              f'3 values ([start, end, increment]) but {len(property_range)} were passed.')

    if end < start:
        raise GooseInputError(f'The range specified for {property_name} has an end value of '
                              f'{end} that is less than the start value of {start}.')

    # if start == end there is only one value, so an increment is not needed.
    if end == start:
        return [start], _decimal_places(start)

    if increment is None:
        raise GooseInputError(f'A range was specified for {property_name} but no increment '
                              f'was provided. Set {property_name}_increment or pass the '
                              f'increment as a third value in the range.')
    if isinstance(increment, bool) or not isinstance(increment, (int, float)):
        raise GooseInputError(f'The increment for {property_name} must be a number, '
                              f'but a {type(increment)} was passed.')
    increment = float(increment)
    if increment <= 0:
        raise GooseInputError(f'The increment for {property_name} must be greater than 0, '
                              f'but {increment} was passed.')

    # number of decimals needed so every value in the titration is exact.
    decimals = max(_decimal_places(start), _decimal_places(end), _decimal_places(increment))

    # build the values. Adding _EPS so that floating point noise doesn't
    # cause us to drop the final value of an evenly divisible range.
    num_steps = int((end - start) / increment + _EPS)
    values = [round(start + (step * increment), decimals) for step in range(num_steps + 1)]

    # if the increment does not divide the range evenly, add the end value so
    # the user always gets the upper bound they asked for.
    if values[-1] < end - _EPS:
        values.append(round(end, decimals))

    return values, decimals


def _hydropathy_extremes(exclude=None):
    '''
    Get the lowest and highest hydropathy values available for positive,
    negative, and uncharged residues. Used to work out the range of mean
    hydropathy values that are compositionally reachable at a given FCR / NCPR.
    Uses the same hydropathy values GOOSE uses everywhere else.

    Parameters
    ----------
    exclude : list or None
        Residues that will not be used during sequence generation. These are
        left out of the calculation.

    Returns
    -------
    dict
        Dictionary mapping 'positive', 'negative', and 'uncharged' to
        (minimum hydropathy, maximum hydropathy) tuples. A residue class that
        has been entirely excluded maps to None.
    '''
    if exclude is None:
        exclude = []
    exclude = [residue.upper() for residue in exclude]

    extremes = {}
    for residue_class, residues in (('positive', ['K', 'R']),
                                    ('negative', ['D', 'E']),
                                    ('uncharged', [aa for aa in parameters.VALID_AMINO_ACIDS
                                                   if aa not in ['K', 'R', 'D', 'E']])):
        values = [amino_acid_dat[aa]['hydrophobicity'] for aa in residues if aa not in exclude]
        extremes[residue_class] = (min(values), max(values)) if values != [] else None
    return extremes


def _hydropathy_limits(fcr, ncpr, extremes):
    '''
    Calculate the minimum and maximum mean hydropathy that a sequence with a
    given FCR and NCPR could possibly have. This is a hard compositional
    limit - a sequence outside of these bounds cannot exist, regardless of how
    many attempts GOOSE makes.

    This is needed because GOOSE's own input checking only bounds the maximum
    FCR for a given hydropathy; it does not catch things like asking for a
    hydropathy of 0.5 at an FCR of 0.2 (impossible, because the 80% of the
    sequence that is uncharged cannot have a hydropathy below 1.0).

    Parameters
    ----------
    fcr : float or None
        The target FCR. None if FCR is not being specified.

    ncpr : float or None
        The target NCPR. None if NCPR is not being specified.

    extremes : dict
        The dict returned by _hydropathy_extremes().

    Returns
    -------
    tuple
        (minimum possible mean hydropathy, maximum possible mean hydropathy).
    '''
    def limits_at(fcr_value, ncpr_value):
        # fraction of the sequence that is positive, negative, and uncharged.
        fraction_positive = (fcr_value + ncpr_value) / 2
        fraction_negative = (fcr_value - ncpr_value) / 2
        fraction_uncharged = 1 - fcr_value
        minimum, maximum = 0.0, 0.0
        for fraction, residue_class in ((fraction_positive, 'positive'),
                                        (fraction_negative, 'negative'),
                                        (fraction_uncharged, 'uncharged')):
            if fraction <= 0:
                continue
            # an entirely excluded residue class means this composition is
            # unreachable. GOOSE catches the charged cases itself, so just
            # skip the class here rather than guessing at a limit.
            if extremes[residue_class] is None:
                continue
            minimum += fraction * extremes[residue_class][0]
            maximum += fraction * extremes[residue_class][1]
        return minimum, maximum

    # work out which FCR / NCPR combinations we need to consider. If a value
    # was not specified, GOOSE is free to choose it, so we take the widest
    # bounds over the values it could choose. The limits are linear in both
    # FCR and NCPR, so checking the corners of the allowed region is enough.
    if fcr is not None and ncpr is not None:
        corners = [(fcr, ncpr)]
    elif fcr is not None:
        # NCPR is free, so it can be anywhere between -FCR and +FCR.
        corners = [(fcr, -fcr), (fcr, fcr)]
    elif ncpr is not None:
        # FCR is free, so it can be anywhere between |NCPR| and 1.
        corners = [(abs(ncpr), ncpr), (1.0, ncpr)]
    else:
        # nothing is pinned down, so the whole scale is available.
        corners = [(0.0, 0.0), (1.0, 1.0), (1.0, -1.0)]

    all_limits = [limits_at(current_fcr, current_ncpr) for current_fcr, current_ncpr in corners]
    return min(limit[0] for limit in all_limits), max(limit[1] for limit in all_limits)


def _charge_counts(length, fcr, ncpr):
    '''
    Work out the number of positive and negative residues implied by an FCR /
    NCPR combination at a given length, and whether those numbers are whole
    residues.

    This matters because GOOSE rounds the number of charged residues to the
    nearest whole residue. Asking for an FCR of 0.15 at a length of 30 (4.5
    charged residues) gives back a sequence whose actual FCR is not 0.15, so a
    library built from those values would be mislabelled.

    Parameters
    ----------
    length : int
        Length of the sequences being made.

    fcr : float or None
        The target FCR.

    ncpr : float or None
        The target NCPR.

    Returns
    -------
    tuple
        (number of positive residues, number of negative residues, exact) where
        exact is True if the requested values land on whole residues. The
        counts are None for values that were not specified.
    '''
    def is_whole(value):
        return abs(value - round(value)) < 1e-6

    if fcr is not None and ncpr is not None:
        num_positive = length * (fcr + ncpr) / 2
        num_negative = length * (fcr - ncpr) / 2
        exact = is_whole(num_positive) and is_whole(num_negative)
        return round(num_positive), round(num_negative), exact
    if fcr is not None:
        return None, None, is_whole(length * fcr)
    if ncpr is not None:
        return None, None, is_whole(length * ncpr)
    return None, None, True


def _kappa_applicable(fcr, ncpr, min_fcr_for_kappa, min_fcr_ncpr_difference):
    '''
    Decide whether kappa can be specified for a given FCR / NCPR combination.

    Kappa describes how positively and negatively charged residues are
    patterned relative to each other, so it is only meaningful when a sequence
    contains both. That makes it undefined when the FCR is 0 and when all of
    the charged residues carry the same sign (FCR == |NCPR|), and it becomes
    increasingly meaningless (and increasingly hard to hit) as the number of
    minority-sign charged residues approaches zero.

    Parameters
    ----------
    fcr : float or None
        The target FCR.

    ncpr : float or None
        The target NCPR.

    min_fcr_for_kappa : float
        Minimum FCR required before kappa will be specified.

    min_fcr_ncpr_difference : float
        Minimum required value of FCR - |NCPR| before kappa will be specified.

    Returns
    -------
    tuple
        (applicable, reason). reason is None when applicable is True.
    '''
    # if FCR was not specified GOOSE chooses it, so there is nothing to check
    # beyond the case where every residue is the same charge.
    if fcr is None:
        if ncpr is not None and abs(abs(ncpr) - 1) < _EPS:
            return False, 'kappa is undefined when |NCPR| is 1 (all charges have the same sign)'
        return True, None

    # kappa is undefined without charged residues.
    if fcr <= _EPS:
        return False, 'kappa is undefined when FCR is 0'

    if fcr < min_fcr_for_kappa - _EPS:
        return False, (f'FCR is below min_FCR_for_kappa '
                       f'({fcr} < {min_fcr_for_kappa})')

    if ncpr is not None:
        difference = fcr - abs(ncpr)
        # all charged residues carry the same sign, so there is no opposite
        # charge patterning to describe.
        if difference <= _EPS:
            return False, 'kappa is undefined when FCR equals |NCPR| (all charges have the same sign)'
        if difference < min_fcr_ncpr_difference - _EPS:
            return False, (f'FCR - |NCPR| is below min_FCR_NCPR_difference_for_kappa '
                           f'({round(difference, 6)} < {min_fcr_ncpr_difference})')

    return True, None


def _goose_input_check(length, property_values, sequence_kwargs):
    '''
    Run a candidate parameter combination through GOOSE's own input checking
    so that anything create.sequence() would reject is caught here instead.

    Parameters
    ----------
    length : int
        Length of the sequences being made.

    property_values : dict
        Dict of FCR, NCPR, hydropathy, and kappa values (any of which can
        be None).

    sequence_kwargs : dict
        The kwargs that will be passed to create.sequence().

    Returns
    -------
    tuple
        (ok, reason). reason is None when ok is True.
    '''
    check_kwargs = dict(sequence_kwargs)
    check_kwargs.update(property_values)
    # fill in the defaults GOOSE expects, then run the same check
    # create.sequence() runs.
    check_kwargs = goose_tools.check_and_correct_props_kwargs(**check_kwargs)
    try:
        goose_tools.check_props_parameters(length, **check_kwargs)
    # GOOSE raises GooseInputError for most bad combinations, but the
    # theoretical max charge calculation raises a bare ValueError for
    # hydropathy values that are out of reach.
    except (GooseInputError, ValueError) as e:
        return False, str(e)
    return True, None


def _count_reasons(skipped):
    '''
    Count how many parameter combinations were skipped for each kind of reason.
    Used for the summary printed when verbose is True, because listing every
    skipped combination is not useful for a large library.

    Parameters
    ----------
    skipped : list
        List of skipped combination dicts.

    Returns
    -------
    dict
        Dict mapping the kind of reason to the number of combinations dropped
        for it, most common first.
    '''
    counts = {}
    for combination in skipped:
        reason_type = combination.get('reason_type', 'unknown')
        counts[reason_type] = counts.get(reason_type, 0) + 1
    return dict(sorted(counts.items(), key=lambda item: item[1], reverse=True))


def _build_name(property_values, decimals):
    '''
    Build the name for a set of property values, e.g.
    'hydro_2_FCR_0.2_NCPR_0.0_kappa_0.1'. Properties that were not specified
    are left out of the name. A property that is part of the library but is not
    applicable to this particular combination (kappa, when there are no
    oppositely charged residues to pattern) is named as 'none'.

    Parameters
    ----------
    property_values : dict
        Dict mapping property names to their values (or None).

    decimals : dict
        Dict mapping property names to the number of decimal places to use.

    Returns
    -------
    str
        The name for this combination of property values.
    '''
    name_pieces = []
    for property_name in _PROPERTY_ORDER:
        if property_name not in property_values:
            continue
        value = property_values[property_name]
        formatted = 'none' if value is None else _format_value(value, decimals[property_name])
        name_pieces.append(f'{_NAME_PREFIXES[property_name]}_{formatted}')
    return '_'.join(name_pieces)


#-\|/-\|/-\|/-\|/-\|/-\|/-\|/-\|/-\|/-\|/-\|/-\|/-\|/-\|/-\|/-\|/-\|/-\|/-\|/-\|/-\|/-
#-\|/-\|/-\|/-\|/-\|/-\|/-\|/-\|/-\|/             \|/-\|/-\|/-\|/-\|/-\|/-\|/-\|/-\|/-
#-\|/-\|/-\|/-\|/-\|/-\|/-\|/-\|/-\|/  Enumerate  \|/-\|/-\|/-\|/-\|/-\|/-\|/-\|/-\|/-
#-\|/-\|/-\|/-\|/-\|/-\|/-\|/-\|/-\|/  possible   \|/-\|/-\|/-\|/-\|/-\|/-\|/-\|/-\|/-
#-\|/-\|/-\|/-\|/-\|/-\|/-\|/-\|/-\|/  parameters \|/-\|/-\|/-\|/-\|/-\|/-\|/-\|/-\|/-
#-\|/-\|/-\|/-\|/-\|/-\|/-\|/-\|/-\|/             \|/-\|/-\|/-\|/-\|/-\|/-\|/-\|/-\|/-
#-\|/-\|/-\|/-\|/-\|/-\|/-\|/-\|/-\|/-\|/-\|/-\|/-\|/-\|/-\|/-\|/-\|/-\|/-\|/-\|/-\|/-


def enumerate_property_combinations(length,
                                    hydropathy_range=None,
                                    hydropathy_increment=None,
                                    FCR_range=None,
                                    FCR_increment=None,
                                    NCPR_range=None,
                                    NCPR_increment=None,
                                    kappa_range=None,
                                    kappa_increment=None,
                                    min_FCR_for_kappa=0.1,
                                    min_FCR_NCPR_difference_for_kappa=0.1,
                                    generate_kappa_free_sequences=True,
                                    require_exact_charge_counts=True,
                                    use_empirical_charge_limit=False,
                                    **kwargs):
    '''
    Work out which combinations of hydropathy, FCR, NCPR, and kappa in a
    requested titration are actually possible to generate. This is the
    parameter selection half of generate_property_library() and does not
    generate any sequences, so it is a cheap way to see how large a library
    will be (and what is being dropped from it) before committing to making it.

    Parameters
    ----------
    length : int
        Length of the sequences that will be made. Needed because whether an
        FCR / NCPR value is possible depends on the length.

    hydropathy_range : float, int, list, tuple, or None
        Range of mean hydropathy values to titrate over. See
        generate_property_library() for the accepted forms.

    hydropathy_increment : float, int, or None
        Increment between hydropathy values.

    FCR_range : float, int, list, tuple, or None
        Range of FCR values to titrate over.

    FCR_increment : float, int, or None
        Increment between FCR values.

    NCPR_range : float, int, list, tuple, or None
        Range of NCPR values to titrate over.

    NCPR_increment : float, int, or None
        Increment between NCPR values.

    kappa_range : float, int, list, tuple, or None
        Range of kappa values to titrate over.

    kappa_increment : float, int, or None
        Increment between kappa values.

    min_FCR_for_kappa : float
        Minimum FCR required for kappa to be specified. Default is 0.1.

    min_FCR_NCPR_difference_for_kappa : float
        Minimum required difference between the FCR and the absolute value of
        the NCPR (FCR - |NCPR|) for kappa to be specified. Default is 0.1.

    generate_kappa_free_sequences : bool
        Whether combinations where kappa cannot be specified should still be
        included without kappa. Default is True.

    require_exact_charge_counts : bool
        Whether to require that the FCR / NCPR values land on whole numbers of
        residues at this length. Default is True.

    use_empirical_charge_limit : bool
        Whether to also apply GOOSE's empirical maximum FCR for a given
        hydropathy. Default is False.

    **kwargs
        Any keyword argument accepted by create.sequence() other than
        return_all_sequences. Only 'exclude', 'hydropathy_tolerance', and
        'kappa_tolerance' change which combinations are possible.

    Returns
    -------
    tuple
        (possible, skipped). possible is a list of dicts, one per usable
        combination, each holding the property values ('hydropathy', 'FCR',
        'NCPR', 'kappa'), a 'name', and the 'charge_residues' the combination
        implies. skipped is a list of dicts holding the property values that
        were dropped and the 'reason' each was dropped.

    Raises
    ------
    GooseInputError
        If any of the ranges, increments, or kappa requirements are invalid.
    '''
    # validate length the same way the rest of GOOSE does.
    goose_tools.length_check(length)
    length = int(length)

    # make sure only kwargs create.sequence() understands were passed. We do
    # not allow return_all_sequences because the library needs one sequence
    # per call.
    goose_tools.check_valid_kwargs(kwargs, _SEQUENCE_KWARGS)
    sequence_kwargs = dict(kwargs)

    # check the kappa requirements.
    for requirement_name, requirement_value in (('min_FCR_for_kappa', min_FCR_for_kappa),
                                                ('min_FCR_NCPR_difference_for_kappa',
                                                 min_FCR_NCPR_difference_for_kappa)):
        if isinstance(requirement_value, bool) or not isinstance(requirement_value, (int, float)):
            raise GooseInputError(f'{requirement_name} must be a number, but a '
                                  f'{type(requirement_value)} was passed.')
        if requirement_value < 0 or requirement_value > 1:
            raise GooseInputError(f'{requirement_name} must be between 0 and 1, but '
                                  f'{requirement_value} was passed.')
    min_FCR_for_kappa = float(min_FCR_for_kappa)
    min_FCR_NCPR_difference_for_kappa = float(min_FCR_NCPR_difference_for_kappa)

    # expand every requested range into the explicit list of values to use.
    values, decimals = {}, {}
    for property_name, property_range, increment in (
            ('hydropathy', hydropathy_range, hydropathy_increment),
            ('FCR', FCR_range, FCR_increment),
            ('NCPR', NCPR_range, NCPR_increment),
            ('kappa', kappa_range, kappa_increment)):
        values[property_name], decimals[property_name] = _expand_property_values(property_range,
                                                                                increment,
                                                                                property_name)
        # never use fewer decimal places in the names than the property calls for.
        decimals[property_name] = max(decimals[property_name], _MINIMUM_DECIMALS[property_name])

    # which properties are actually part of this library. Properties that were
    # not specified are left out of the names entirely.
    specified = [property_name for property_name in _PROPERTY_ORDER
                 if values[property_name] != [None]]
    if specified == []:
        raise GooseInputError('No properties were specified, so there is nothing to titrate '
                              'over. Specify a range for at least one of hydropathy_range, '
                              'FCR_range, NCPR_range, or kappa_range.')

    # hydropathy values available for each class of residue, used for the
    # compositional hydropathy limits below.
    extremes = _hydropathy_extremes(sequence_kwargs.get('exclude'))

    # tolerance GOOSE will allow on hydropathy. A target slightly outside the
    # hard compositional limits is still reachable within tolerance.
    hydropathy_tolerance = sequence_kwargs.get('hydropathy_tolerance')
    if hydropathy_tolerance is None:
        hydropathy_tolerance = parameters.MAXIMUM_HYDRO_ERROR

    possible, skipped = [], []

    # iterate in the same order the properties appear in the names so the
    # returned library reads as a clean titration.
    for hydropathy in values['hydropathy']:
        for fcr in values['FCR']:
            for ncpr in values['NCPR']:

                # ---- checks that do not depend on kappa ----
                charge_reason, charge_reason_type = None, None
                num_positive, num_negative = None, None

                # |NCPR| can never exceed the FCR.
                if fcr is not None and ncpr is not None and abs(ncpr) > fcr + _EPS:
                    charge_reason = f'|NCPR| ({abs(ncpr)}) is greater than the FCR ({fcr})'
                    charge_reason_type = 'NCPR_exceeds_FCR'

                # the FCR / NCPR values have to land on whole residues at this
                # length or the sequence GOOSE makes will not have the values
                # its name claims.
                if charge_reason is None:
                    num_positive, num_negative, exact = _charge_counts(length, fcr, ncpr)
                    if require_exact_charge_counts and not exact:
                        charge_reason = (f'an FCR of {fcr} and NCPR of {ncpr} do not give a whole '
                                         f'number of charged residues at a length of {length}')
                        charge_reason_type = 'charge_values_not_whole_residues'

                # hard compositional limits on mean hydropathy. Nothing GOOSE
                # does can get a sequence outside of these.
                if charge_reason is None and hydropathy is not None:
                    minimum_hydropathy, maximum_hydropathy = _hydropathy_limits(fcr, ncpr, extremes)
                    if hydropathy < minimum_hydropathy - hydropathy_tolerance - _EPS:
                        charge_reason = (f'a hydropathy of {hydropathy} is below the minimum of '
                                         f'{round(minimum_hydropathy, 3)} possible with an FCR of '
                                         f'{fcr} and NCPR of {ncpr}')
                        charge_reason_type = 'hydropathy_too_low_for_charge'
                    elif hydropathy > maximum_hydropathy + hydropathy_tolerance + _EPS:
                        charge_reason = (f'a hydropathy of {hydropathy} is above the maximum of '
                                         f'{round(maximum_hydropathy, 3)} possible with an FCR of '
                                         f'{fcr} and NCPR of {ncpr}')
                        charge_reason_type = 'hydropathy_too_high_for_charge'

                # optionally apply GOOSE's empirical limit on how much charge a
                # sequence can carry at a given hydropathy and still be
                # disordered. Stricter than the checks create.sequence() runs.
                if charge_reason is None and hydropathy is not None and use_empirical_charge_limit:
                    empirical_maximum = parameters.calculate_max_charge_emperical(hydropathy)
                    charge_value = fcr if fcr is not None else (abs(ncpr) if ncpr is not None else None)
                    if charge_value is not None and charge_value > empirical_maximum + _EPS:
                        charge_reason = (f'an FCR of {charge_value} is above the empirical maximum of '
                                         f'{round(empirical_maximum, 3)} for a hydropathy of {hydropathy}')
                        charge_reason_type = 'above_empirical_charge_limit'

                # ---- now walk kappa ----
                # whether kappa can be specified depends only on the charge
                # values, so if it cannot be, the whole kappa dimension
                # collapses to a single sequence made without kappa.
                kappa_ok, kappa_reason = _kappa_applicable(fcr, ncpr, min_FCR_for_kappa,
                                                           min_FCR_NCPR_difference_for_kappa)
                if 'kappa' in specified and not kappa_ok:
                    kappa_values = [None]
                else:
                    kappa_values = values['kappa']

                for kappa in kappa_values:
                    # only include the properties this library is titrating.
                    property_values = {property_name: {'hydropathy': hydropathy, 'FCR': fcr,
                                                       'NCPR': ncpr, 'kappa': kappa}[property_name]
                                       for property_name in specified}
                    combination = dict(property_values)
                    combination['name'] = _build_name(property_values, decimals)

                    # anything the charge / hydropathy checks caught above.
                    if charge_reason is not None:
                        combination['reason'] = charge_reason
                        combination['reason_type'] = charge_reason_type
                        skipped.append(combination)
                        continue

                    # kappa was requested for this library but cannot be used
                    # for this combination.
                    if 'kappa' in specified and not kappa_ok:
                        if not generate_kappa_free_sequences:
                            combination['reason'] = kappa_reason
                            combination['reason_type'] = 'kappa_not_applicable'
                            skipped.append(combination)
                            continue
                        combination['note'] = f'generated without kappa because {kappa_reason}'

                    # finally, run the combination past GOOSE's own input
                    # checking so we do not try to make something
                    # create.sequence() will refuse.
                    ok, reason = _goose_input_check(length,
                                                    {'FCR': fcr, 'NCPR': ncpr,
                                                     'hydropathy': hydropathy,
                                                     'kappa': combination.get('kappa')},
                                                    sequence_kwargs)
                    if not ok:
                        combination['reason'] = reason
                        combination['reason_type'] = 'rejected_by_GOOSE_input_check'
                        skipped.append(combination)
                        continue

                    combination['charge_residues'] = (num_positive, num_negative)
                    possible.append(combination)

    return possible, skipped


#-\|/-\|/-\|/-\|/-\|/-\|/-\|/-\|/-\|/-\|/-\|/-\|/-\|/-\|/-\|/-\|/-\|/-\|/-\|/-\|/-\|/-
#-\|/-\|/-\|/-\|/-\|/-\|/-\|/-\|/-\|/             \|/-\|/-\|/-\|/-\|/-\|/-\|/-\|/-\|/-
#-\|/-\|/-\|/-\|/-\|/-\|/-\|/-\|/-\|/  Generate   \|/-\|/-\|/-\|/-\|/-\|/-\|/-\|/-\|/-
#-\|/-\|/-\|/-\|/-\|/-\|/-\|/-\|/-\|/  the        \|/-\|/-\|/-\|/-\|/-\|/-\|/-\|/-\|/-
#-\|/-\|/-\|/-\|/-\|/-\|/-\|/-\|/-\|/  library    \|/-\|/-\|/-\|/-\|/-\|/-\|/-\|/-\|/-
#-\|/-\|/-\|/-\|/-\|/-\|/-\|/-\|/-\|/             \|/-\|/-\|/-\|/-\|/-\|/-\|/-\|/-\|/-
#-\|/-\|/-\|/-\|/-\|/-\|/-\|/-\|/-\|/-\|/-\|/-\|/-\|/-\|/-\|/-\|/-\|/-\|/-\|/-\|/-\|/-


def generate_property_library(length,
                              num_copies=1,
                              hydropathy_range=None,
                              hydropathy_increment=None,
                              FCR_range=None,
                              FCR_increment=None,
                              NCPR_range=None,
                              NCPR_increment=None,
                              kappa_range=None,
                              kappa_increment=None,
                              min_FCR_for_kappa=0.1,
                              min_FCR_NCPR_difference_for_kappa=0.1,
                              generate_kappa_free_sequences=True,
                              require_exact_charge_counts=True,
                              use_empirical_charge_limit=False,
                              require_unique_copies=True,
                              on_generation_failure='skip',
                              return_details=False,
                              verbose=True,
                              **kwargs):
    '''
    Generate a library of disordered sequences that titrates across
    user-specified ranges of hydropathy, FCR, NCPR, and kappa.

    Every combination of the requested property values is considered. The ones
    that are not possible to make are dropped (see below), and each
    combination that survives is generated num_copies times.

    Combinations are dropped when:
      * |NCPR| is greater than the FCR.
      * The FCR / NCPR values do not give a whole number of charged residues at
        the requested length (turn off with require_exact_charge_counts=False).
      * The hydropathy value is outside the range that is compositionally
        possible at that FCR / NCPR. For example, a hydropathy of 0.5 cannot be
        reached at an FCR of 0.2, because the uncharged 80% of the sequence
        cannot have a hydropathy below 1.0.
      * GOOSE's own input checking rejects the combination (which includes its
        maximum FCR for a given hydropathy).
      * Sequence generation itself fails after all attempts (turn into an
        exception with on_generation_failure='raise').

    Kappa is only specified when it is meaningful. It is skipped when the FCR
    is 0, when all of the charged residues carry the same sign (FCR == |NCPR|),
    when the FCR is below min_FCR_for_kappa, and when FCR - |NCPR| is below
    min_FCR_NCPR_difference_for_kappa. By default those combinations are still
    made, just without kappa specified, and are named with 'kappa_none'.

    Parameters
    ----------
    length : int
        Length of every sequence in the library.

    num_copies : int
        Number of independent sequences to generate for each combination of
        property values. Default is 1.

    hydropathy_range : float, int, list, tuple, or None
        The mean hydropathy values to titrate over. Can be:
          * None - hydropathy is not specified at all (the default).
          * a single number - hydropathy is held constant at that value.
          * [start, end] - titrated from start to end using
            hydropathy_increment.
          * [start, end, increment] - as above, with the increment given as
            the third value.
        If the increment does not divide the range evenly, the end value is
        still included as the final value of the titration.

    hydropathy_increment : float, int, or None
        Increment between hydropathy values. Required if hydropathy_range is
        given as [start, end].

    FCR_range : float, int, list, tuple, or None
        The FCR values to titrate over. Same accepted forms as
        hydropathy_range.

    FCR_increment : float, int, or None
        Increment between FCR values.

    NCPR_range : float, int, list, tuple, or None
        The NCPR values to titrate over. Same accepted forms as
        hydropathy_range.

    NCPR_increment : float, int, or None
        Increment between NCPR values.

    kappa_range : float, int, list, tuple, or None
        The kappa values to titrate over. Same accepted forms as
        hydropathy_range.

    kappa_increment : float, int, or None
        Increment between kappa values.

    min_FCR_for_kappa : float
        The minimum FCR a combination must have for kappa to be specified.
        Default is 0.1.

    min_FCR_NCPR_difference_for_kappa : float
        The minimum difference between the FCR and the absolute value of the
        NCPR (FCR - |NCPR|) a combination must have for kappa to be specified.
        This is the fraction of the sequence available for oppositely charged
        residues, so small values leave very few charges to pattern. Default
        is 0.1.

    generate_kappa_free_sequences : bool
        What to do with combinations where kappa cannot be specified. If True
        (the default) the combination is still made, without kappa, and named
        with 'kappa_none'. If False the combination is dropped.

    require_exact_charge_counts : bool
        If True (the default), combinations whose FCR / NCPR values do not give
        a whole number of positively and negatively charged residues at this
        length are dropped. GOOSE rounds to the nearest residue, so keeping
        these means some sequences will not have the FCR / NCPR values their
        names claim. Note that this depends on both values together: at a
        length of 50, an FCR of 0.1 with an NCPR of 0 needs 2.5 positive and
        2.5 negative residues, so it is dropped, while an FCR of 0.2 with an
        NCPR of 0 (5 and 5) is kept. Picking a length where
        length * increment / 2 is a whole number keeps the most combinations.

    use_empirical_charge_limit : bool
        If True, also drop combinations where the FCR exceeds GOOSE's empirical
        maximum FCR for the requested hydropathy. This accounts for disorder,
        which the check create.sequence() runs does not, but it only prunes
        additional combinations at high hydropathy values (above roughly 5.5).
        Default is False.

    require_unique_copies : bool
        If True (the default), retry when a copy duplicates a sequence already
        made for the same combination of property values.

    on_generation_failure : str
        What to do when GOOSE cannot make a sequence for a combination that
        passed all of the checks. 'skip' (the default) records the failure and
        moves on; 'raise' re-raises the exception.

    return_details : bool
        If True, return a dict of the library plus information about what was
        skipped and what failed instead of just the library. Default is False.

    verbose : bool
        Whether to print progress and a summary of what was skipped or failed.
        Default is True.

    **kwargs
        Any keyword argument accepted by create.sequence() other than
        return_all_sequences. These are passed through for every sequence, so
        this is where to set things like attempts, hydropathy_tolerance,
        kappa_tolerance, exclude, disorder_cutoff, or custom_probabilities.

        One note on attempts: it sets how long GOOSE will keep trying before it
        gives up on a sequence, and the hardest combinations (a high kappa
        combined with a high FCR and a high hydropathy) will use every attempt
        they are given. create.sequence() defaults to 5000 attempts, which can
        mean minutes on a single combination, so this function defaults to
        GOOSE's DEFAULT_ATTEMPTS (currently 200) instead. Lower it to make a
        large library finish sooner, or raise it if combinations you expect to
        be possible are showing up as failures.

    Returns
    -------
    dict
        By default, a dict mapping sequence names to sequences, where each name
        describes the properties the sequence was made with, e.g.
        'hydro_2_FCR_0.2_NCPR_0.0_kappa_0.1_rep_1'. Properties that were not
        part of the library are left out of the names. Every value of a given
        property is written with the same number of decimal places (however
        many its increment needs), so the names sort and parse consistently.

        If return_details is True, returns a dict with the keys:
          * 'library' - the dict of names to sequences described above.
          * 'skipped' - list of dicts for combinations that were not possible,
            each with the property values, the 'reason' it was dropped, and a
            'reason_type' naming the kind of problem.
          * 'failed' - list of dicts for combinations that passed every check
            but that GOOSE could not make, each with the property values and
            the 'reason' generation failed.
          * 'duplicates' - list of names whose sequence duplicates another copy
            of the same combination of property values.

    Raises
    ------
    GooseInputError
        If any of the inputs are invalid.

    GooseFail
        If a sequence cannot be generated and on_generation_failure is 'raise'.

    Examples
    --------
    >>> # titrate hydropathy and FCR, 2 copies of each sequence
    >>> library = generate_property_library(50, num_copies=2,
    ...                                     hydropathy_range=[2, 4], hydropathy_increment=0.5,
    ...                                     FCR_range=[0.1, 0.3], FCR_increment=0.1)
    >>>
    >>> # titrate all four properties
    >>> library = generate_property_library(100, num_copies=3,
    ...                                     hydropathy_range=[2, 4], hydropathy_increment=1,
    ...                                     FCR_range=[0.1, 0.4], FCR_increment=0.1,
    ...                                     NCPR_range=[-0.2, 0.2], NCPR_increment=0.1,
    ...                                     kappa_range=[0.1, 0.9], kappa_increment=0.2,
    ...                                     min_FCR_for_kappa=0.1,
    ...                                     min_FCR_NCPR_difference_for_kappa=0.1)
    '''
    # check num_copies.
    if isinstance(num_copies, bool) or not isinstance(num_copies, int):
        raise GooseInputError(f'num_copies must be an integer, but a {type(num_copies)} was passed.')
    if num_copies < 1:
        raise GooseInputError(f'num_copies must be at least 1, but {num_copies} was passed.')

    # check on_generation_failure.
    if on_generation_failure not in ['skip', 'raise']:
        raise GooseInputError(f"on_generation_failure must be 'skip' or 'raise', but "
                              f"'{on_generation_failure}' was passed.")

    # work out which parameter combinations are actually possible. This also
    # validates length, the ranges, the increments, and the kwargs.
    possible, skipped = enumerate_property_combinations(
        length,
        hydropathy_range=hydropathy_range, hydropathy_increment=hydropathy_increment,
        FCR_range=FCR_range, FCR_increment=FCR_increment,
        NCPR_range=NCPR_range, NCPR_increment=NCPR_increment,
        kappa_range=kappa_range, kappa_increment=kappa_increment,
        min_FCR_for_kappa=min_FCR_for_kappa,
        min_FCR_NCPR_difference_for_kappa=min_FCR_NCPR_difference_for_kappa,
        generate_kappa_free_sequences=generate_kappa_free_sequences,
        require_exact_charge_counts=require_exact_charge_counts,
        use_empirical_charge_limit=use_empirical_charge_limit,
        **kwargs)

    length = int(length)

    if verbose:
        print(f'Generating {len(possible) * num_copies} sequences '
              f'({len(possible)} parameter combinations x {num_copies} '
              f'{"copy" if num_copies == 1 else "copies"}).')
        if skipped != []:
            print(f'{len(skipped)} parameter combinations were not possible and were skipped:')
            for reason_type, count in _count_reasons(skipped).items():
                print(f'  {count} - {reason_type}')

    library, failed, duplicates = {}, [], []
    successful_combinations = 0

    # a library makes a lot of sequences, and a single hard combination can
    # take a long time at create.sequence()'s default of 5000 attempts. Use
    # GOOSE's DEFAULT_ATTEMPTS instead so one difficult combination cannot
    # stall the whole library, unless the user asked for something else.
    if 'attempts' not in kwargs or kwargs['attempts'] is None:
        kwargs['attempts'] = parameters.DEFAULT_ATTEMPTS

    # timing, used to report progress on long runs.
    start_time = time.time()
    last_report_time = start_time

    for combination_number, combination in enumerate(possible, start=1):
        # the property values GOOSE needs, dropping the bookkeeping keys and
        # any property that is not part of this library.
        generation_kwargs = dict(kwargs)
        for property_name in _PROPERTY_ORDER:
            if combination.get(property_name) is not None:
                generation_kwargs[property_name] = combination[property_name]

        sequences_for_combination = []
        for copy_number in range(1, num_copies + 1):
            name = f'{combination["name"]}_rep_{copy_number}'
            try:
                sequence = create.sequence(length, **generation_kwargs)
                # if the copies should be unique, try again a few times when we
                # get a sequence we have already made for this combination.
                if require_unique_copies:
                    for _ in range(4):
                        if sequence not in sequences_for_combination:
                            break
                        sequence = create.sequence(length, **generation_kwargs)
                    if sequence in sequences_for_combination:
                        duplicates.append(name)
            # catching everything here because a combination that passes every
            # check can still fail during generation, and the failure can
            # surface as more than just a GooseFail.
            except Exception as e:
                if on_generation_failure == 'raise':
                    raise
                failure = {property_name: combination[property_name]
                           for property_name in _PROPERTY_ORDER if property_name in combination}
                failure['name'] = name
                failure['reason'] = f'{type(e).__name__}: {e}'
                failed.append(failure)
                continue

            sequences_for_combination.append(sequence)
            library[name] = sequence

        if sequences_for_combination != []:
            successful_combinations += 1

        # report progress on a timer rather than every N combinations. Some
        # combinations (high kappa at high FCR and high hydropathy) take far
        # longer than others, so a timer is what keeps a slow run from looking
        # like a hung one.
        if verbose and (time.time() - last_report_time) > _SECONDS_BETWEEN_PROGRESS_REPORTS:
            print(f'  ...{combination_number} of {len(possible)} parameter combinations, '
                  f'{len(library)} sequences made, {round(time.time() - start_time)}s elapsed.')
            last_report_time = time.time()

    if verbose:
        print(f'\nGenerated {len(library)} sequences across {successful_combinations} '
              f'parameter combinations in {round(time.time() - start_time, 1)}s.')
        if failed != []:
            print(f'{len(failed)} sequences could not be generated despite having valid parameters:')
            for failure in failed[:10]:
                print(f'  {failure["name"]}: {failure["reason"]}')
            if len(failed) > 10:
                print(f'  ...and {len(failed) - 10} more.')
        if duplicates != []:
            print(f'{len(duplicates)} sequences are duplicates of another copy of the '
                  f'same parameter combination.')
        if skipped != [] or failed != []:
            print('Use return_details=True to get the full list of skipped and failed '
                  'parameter combinations.')

    if return_details:
        return {'library': library, 'skipped': skipped, 'failed': failed,
                'duplicates': duplicates}

    return library
