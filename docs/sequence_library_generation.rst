Sequence Libraries
==================

GOOSE generates libraries that sweep physical properties, amino-acid fractions,
amino-acid class fractions, or predicted dimensions. Each generator returns a
dictionary mapping descriptive sequence names to sequences. ``num_copies`` sets
the number of sequences requested per target.

All four generators are available from ``goose.property_library``. The new
generators are also available from their own modules:

.. code-block:: python

   from goose.property_library import (
       generate_property_library,
       generate_fractions_library,
       generate_class_library,
       generate_dimensions_library,
   )

   # Equivalently:
   from goose.fractions_library import generate_fractions_library
   from goose.class_library import generate_class_library
   from goose.dimensions_library import generate_dimensions_library

Amino-acid fractions
--------------------

By default, ``generate_fractions_library`` sweeps all 20 amino acids
**independently**, from zero through the highest whole-residue fraction permitted
by GOOSE's ``MAX_FRACTION_DICT`` at the requested length. Other amino-acid
fractions remain free. For example, the zero-glutamate target contains no E but
does not constrain the fractions of the other 19 amino acids.

.. code-block:: python

   # Sweep every amino acid independently, two copies per target.
   library = generate_fractions_library(25, num_copies=2)

   # Sweep only D, E, and K through their default ranges.
   library = generate_fractions_library(100, amino_acids='DEK')

   # Supply ranges for just the residues of interest.
   library = generate_fractions_library(
       100, num_copies=3,
       fraction_ranges={'E': [0, 0.5, 0.1], 'P': [0.1, 0.4, 0.05]},
   )

``fraction_ranges`` maps residue letters to a fixed value, a
``[minimum, maximum]`` range, or a ``[minimum, maximum, increment]`` range.
Providing ranges selects just those residues, unless ``amino_acids`` also
selects additional residues to sweep through their default ranges.
``max_aa_fractions`` overrides GOOSE's limits and the corresponding default
sweep maxima.

Whole-residue intervals
~~~~~~~~~~~~~~~~~~~~~~~

Fractions are planned in integer residue counts, so every target is exactly
attainable at the specified length. Default intervals are approximately five
percentage points: ``max(1, floor(length * 0.05 + 0.5))`` residues per step.

.. list-table:: Default fraction spacing
   :header-rows: 1

   * - Sequence length
     - Residues per step
     - Fraction increment
   * - 10
     - 1
     - 0.10
   * - 25
     - 1
     - 0.04
   * - 50
     - 3
     - 0.06
   * - 100
     - 5
     - 0.05

Explicit ``fraction_increment`` values must be at least ``1 / length``. Larger
increments round up to whole residues: at length 25, an increment of 0.1 becomes
3 residues, or 0.12. An increment of 0.001 raises ``GooseInputError`` at that
length. A third value in a range sets the increment for that particular axis.

Range minima round up and maxima round down to attainable counts. The highest
attainable value in the range is always included, even if the final step is
shorter. For instance, the default isoleucine maximum of 0.53 becomes 13/25, or
0.52, at length 25. Fixed targets must already correspond to a whole number of
residues; a fixed 0.1 fraction at length 25 raises an input error.

Combinations are opt-in
~~~~~~~~~~~~~~~~~~~~~~~~

``combine=True`` generates the Cartesian product of explicitly selected axes
instead of separate sweeps. Targets whose specified fractions sum above one
are recorded as skipped. If all 20 amino acids are specified, their fractions
must sum to one. Large products can request many sequences; inspect the target
list first using the enumeration helper below.

.. code-block:: python

   library = generate_fractions_library(
       100, fraction_ranges={'E': [0, 0.4, 0.1], 'P': [0, 0.3, 0.1]},
       combine=True,
   )

Amino-acid classes
------------------

``generate_class_library`` uses the same whole-residue intervals and independent
sweeps as the fraction library. By default it sweeps each of GOOSE's nine
classes from zero to its ``MAX_CLASS_FRACTIONS`` limit. Residues within each
targeted class are chosen by the existing class-based sequence generator.

.. list-table:: GOOSE classes
   :header-rows: 1

   * - Class name
     - Residues
   * - ``aromatic``
     - F, W, Y
   * - ``aliphatic``
     - A, I, L, M, V
   * - ``positive``
     - K, R
   * - ``negative``
     - D, E
   * - ``polar``
     - S, T, N, Q
   * - ``glycine``
     - G
   * - ``proline``
     - P
   * - ``cysteine``
     - C
   * - ``histidine``
     - H

.. code-block:: python

   library = generate_class_library(25, num_copies=2)

   library = generate_class_library(
       100, classes=['aromatic', 'positive', 'negative'], class_increment=0.1,
   )

   library = generate_class_library(
       100, class_ranges={'aromatic': [0, 0.3, 0.05], 'polar': [0.1, 0.5, 0.1]},
       combine=True,
   )

``class_ranges`` and ``class_increment`` follow the fraction-library semantics.
``max_class_fractions`` overrides the default class limits. Fully specified class
combinations must sum to one. Unspecified classes remain free.

Rg and end-to-end distance
----------------------------

``generate_dimensions_library`` sweeps either predicted radius of gyration
(``rg``) or predicted end-to-end distance (``re``), in angstroms. It defaults to
Rg and a 1 Å interval. Each omitted minimum or maximum independently uses
GOOSE's empirical bound for the sequence length:
``get_min_rg`` / ``get_max_rg`` or ``get_min_re`` / ``get_max_re`` in
``goose.backend.parameters``. The upper bound is always included.

.. code-block:: python

   # Default length-dependent Rg range, two copies per target, 1 Å spacing.
   library = generate_dimensions_library(100, num_copies=2)

   library = generate_dimensions_library(
       100, min_rg=20, max_rg=35, dimension_increment=2.5,
   )

   library = generate_dimensions_library(
       100, num_copies=3, rg_or_re='re',
       min_re=30, max_re=70, dimension_increment=5,
   )

   # Specify only a minimum; retain the default maximum for this length.
   library = generate_dimensions_library(100, min_rg=20)

Bounds must lie inside GOOSE's permitted range. Set ``rg_or_re='re'`` when
using ``min_re`` or ``max_re``; bounds for both dimensions cannot be combined.
Names encode target dimensions, and the existing generator uses
``allowed_error`` (default 0.5 Å) as its prediction tolerance.

Planning, generation settings, and failures
-------------------------------------------

Planning functions return ``(possible, skipped)`` without generating sequences
or running prediction models. Composition entries include the target
fractions, a name, and exact integer ``counts``; dimension entries contain the
target ``rg`` or ``re`` and a name.

.. code-block:: python

   from goose.property_library import (
       enumerate_fraction_combinations,
       enumerate_class_combinations,
       enumerate_dimension_combinations,
   )

   targets, skipped = enumerate_fraction_combinations(25, amino_acids='DEK')
   print(len(targets))
   print(targets[1])
   # {'D': 0.04, 'counts': {'D': 1}, 'name': 'D_0.04'}

All three new generators forward disorder settings and other supported options
to their corresponding ``create`` functions. Fractions and dimensions use
``attempts``; classes use ``num_attempts``. Their library defaults are 200
attempts per sequence. Lower these for a faster exploratory run, or raise them
for difficult targets. Fraction libraries also accept ``batch_size``.

GOOSE's maximum fractions and dimension bounds are empirical generation limits,
not guarantees of disorder or successful generation at every length. By
default, generation failures are recorded and skipped. Invalid inputs raise
``GooseInputError``. Use ``on_generation_failure='raise'`` to stop at the first
generation failure.

.. code-block:: python

   result = generate_fractions_library(
       25, num_copies=2, amino_acids='DEK',
       attempts=20, return_details=True, verbose=False,
   )
   library = result['library']
   failed = result['failed']
   skipped = result['skipped']
   duplicates = result['duplicates']

Names describe targets and copies, such as ``E_0.2_rep_1``,
``aromatic_0.12_rep_2``, or ``rg_20.0_rep_1``. Composition counts are checked
before a generated sequence is added to the library. With
``require_unique_copies=True`` (default), duplicate copies are retried up to
four times. Unavoidable duplicates, such as multiple copies of a homopolymer,
are retained and reported in ``duplicates``. Uniqueness is checked within each
target, not across separate targets.

All library generators show a ``tqdm`` progress bar by default, using a notebook
widget when available or a terminal bar otherwise. The total is the number of
requested copies after filtering incompatible targets. Each successful or
failed copy advances the bar once; duplicate retries stay within that copy.
The bar includes elapsed time, estimated remaining time, the current sequence
target, and generated, failed, duplicate, and skipped counts. The final bar
retains the counts; full failure and skip reasons are available through
``return_details=True``. Set ``verbose=False`` to suppress the bar.

Physical-property libraries
------------------------------

The existing ``generate_property_library`` sweeps hydropathy, FCR, NCPR, and
kappa. Its interface and behavior are unchanged:

.. code-block:: python

   library = generate_property_library(
       100, num_copies=2,
       hydropathy_range=[2, 4, 0.5],
       FCR_range=[0.1, 0.4, 0.1],
   )

Limiting property-library search time
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

``attempts`` limits the number of candidate **batches** per sequence. It does
not limit the work inside each batch: with default settings, ``attempts=1``
starts ten candidates and can run up to 1000 hydropathy iterations and 20000
kappa iterations. Kappa optimization can be particularly expensive for short
sequences and targets that cannot be reached with their discrete charge counts.

Use ``max_hydropathy_iterations`` and ``max_kappa_iterations`` to set separate
positive integer limits on these inner searches. Their default of ``None``
preserves the existing optimizer limits. A faster exploratory configuration is:

.. code-block:: python

   library = generate_property_library(
       25, num_copies=1,
       hydropathy_range=[1, 4, 0.5],
       FCR_range=[0.2, 0.4, 0.04],
       NCPR_range=[-0.2, 0.2, 0.04],
       kappa_range=[0.1, 0.9, 0.1],
       attempts=1,
       max_hydropathy_iterations=100,
       max_kappa_iterations=100,
       require_exact_charge_counts=True,
       return_details=True,
   )

Lower limits may cause feasible targets to be recorded as failed; they do not
relax the requested tolerances. These are iteration limits, not wall-clock
deadlines. A large Cartesian grid can still take a long time even with short
individual searches. Use ``enumerate_property_combinations`` to count targets
first. Keep ``require_exact_charge_counts=True`` when possible: even if FCR and
NCPR are individually multiples of ``1 / length``, their combination can require
half-residue positive or negative counts. At length 25, FCR=0.04 with NCPR=0
would require half a positive and half a negative residue.
