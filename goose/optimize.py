# before adding in the normalization recalcualton
"""
SequenceOptimizer rewrite. Trying to make it work better
for properties that have variable ranges of values. 
"""

import logging
import math
import random
from collections import deque
from typing import Any, Dict, List, Optional, Tuple, Union
import numpy as np
import sparrow
from tqdm.auto import tqdm

from goose.backend import optimizer_properties
from goose.backend.fast_mutations import (
    cython_make_point_mutation as make_point_mutation,
    cython_make_multi_mutations as make_multi_mutations,
)
from goose.backend.fast_mutations import cython_shuffle_sequence as shuffle_sequence

from goose.backend.optimizer_properties import ProteinProperty
from goose.data import amino_acids


# ==============================================================================
# UTILITY FUNCTIONS
# ==============================================================================

def build_random_sequence(target_length: int) -> str:
    """Build a random amino acid sequence"""
    return ''.join(random.choices(
        list(amino_acids.amino_acid_dat.keys()), 
        k=target_length
    ))


def build_diverse_initial_sequences(target_length: int, num_sequences: int = 5) -> List[str]:
    """Build multiple diverse initial sequences to find the best starting point.
    
    Generates candidates until either ``num_sequences`` unique sequences exist
    or a hard cap of ``5*num_sequences`` attempts is reached. The result is
    de-duplicated but ordering of the first random seed is preserved.
    """
    base_aas = ['G', 'P', 'S', 'T', 'Q', 'N']
    biased_lists = [
        ['A', 'V', 'I', 'L', 'M'],
        ['F', 'W', 'Y'],
        ['D', 'E'],
        ['K', 'R'],
        ['D', 'E', 'K', 'R'],
        ['C'],
        ['H'],
        ['T', 'S', 'N', 'Q', 'P', 'G'],
        ['G', 'P'],
        ['A', 'V', 'I', 'L', 'M', 'E', 'D'],
        ['A', 'V', 'I', 'L', 'M', 'K', 'R'],
        ['F', 'W', 'Y', 'E', 'D'],
        ['F', 'W', 'Y', 'K', 'R'],
    ]
    seen: set = set()
    sequences: List[str] = []
    
    def _add(seq: str) -> None:
        if seq not in seen:
            seen.add(seq)
            sequences.append(seq)
    
    _add(build_random_sequence(target_length))
    max_attempts = max(num_sequences * 5, num_sequences + 10)
    attempts = 0
    while len(sequences) < num_sequences and attempts < max_attempts:
        if len(sequences) % 2 == 0:
            _add(build_random_sequence(target_length))
        else:
            disorder_aas = base_aas + random.choice(biased_lists)
            _add(''.join(random.choices(disorder_aas, k=target_length)))
        attempts += 1
    return sequences




# ==============================================================================
# SEQUENCE OPTIMIZER
# ==============================================================================

class SequenceOptimizer:
    """
    Sequence optimizer that combines the original version [which Jeff crushed] and an attempted rewrite
    that worked better but was over 2,000 lines of code and became unsustainable.
    """
    
    def __init__(
        self,
        target_length: int,
        max_iterations: int = 1000,
        verbose: bool = True, # progress bar and logging
        enable_adaptive_scaling: bool = True,
        enable_shuffling: bool = True,
        shuffle_frequency: int = 50,  # More frequent shuffling for diversity
        stagnation_threshold: int = 25,  # Detect and respond to stagnation sooner
        convergence_tolerance: float = 1e-4,  # Less strict convergence criterion
        convergence_window: int = 20,
        enable_early_convergence: bool = False,
        convergence_patience: int = 20,  # Reduced patience for faster adaptation
        # Error tolerance parameters
        error_tolerance: Optional[float] = 1e-6,  # Stop when total error below this value
        enable_error_tolerance: bool = True,  # Enable error tolerance early stopping
        # Mutation parameters
        num_candidates: int = 5,  # Significantly more candidates for better exploration
        num_starting_candidates: int = 100,  # allow setting different number of starting candidates to explore more seq space.
        min_mutations: int = 1,
        max_mutations: int = 15,  # Allow more aggressive mutations
        mutation_ratio: int = 10,  # More mutations per sequence (length/10 vs length/20)
        # Shuffling parameters
        global_shuffle_probability: float = 0.4,  # More global shuffling for diversity
        shuffle_window_size: int = 15,  # Larger shuffle windows
        # History tracking parameters
        improvement_history_size: int = 20,
        error_history_size: int = 50,
        # Adaptive scaling parameters
        min_trend_samples: int = 5,  # React faster to trends
        max_distance_factor: float = 3.0,  # Allow higher scaling for distant targets
        distance_offset: float = 0.2,  # Lower offset for more sensitive distance calculation
        stagnation_multiplier: float = 1.0,  # Stronger stagnation response
        low_contribution_threshold: float = 0.15,  # More sensitive to low contribution
        high_error_threshold: float = 0.05,  # Lower threshold for high error detection
        boost_factor: float = 2.0,  # Stronger boost for underperforming properties
        scale_momentum: float = 0.5,  # Less momentum for faster adaptation
        scale_learning_rate: float = 0.5,  # Faster learning
        min_scale: float = 0.1,
        max_scale: float = 8.0,  # Allow higher maximum scaling
        # Stagnation parameters
        stagnation_improvement_threshold: float = 0.005,  # More sensitive to stagnation
        improvement_trend_threshold: float = -0.001,
        stagnation_boost_factor: float = 2.0,  # Stronger boost when stagnant
        # --- Multi-objective behavior controls ---
        enforce_raw_monotonicity: bool = False,  # If True, candidates must not
                                                 # increase the unsatisfied raw
                                                 # error sum. Off by default to
                                                 # allow proper multi-objective
                                                 # trade-offs.
        raw_monotonicity_slack: float = 0.0,  # Allowed relative slack when
                                              # enforce_raw_monotonicity is True
                                              # (e.g. 0.05 = 5% leeway).
        scale_freeze_window: int = 0,  # If >0, after stagnation recovery freeze
                                       # adaptive scaling for this many
                                       # iterations so the boost can take
                                       # effect. 0 = use stagnation_threshold.
        max_norm_boost: float = 100.0,  # Cap on the recovery-only norm-boost
                                        # multiplier (separate from clipped
                                        # normalization_factors).
        norm_boost_decay: float = 0.95,  # Per-iteration decay applied to the
                                         # recovery norm-boost (1.0 = no decay).
        recompute_norm_on_emergency: bool = True,  # Blend new normalization
                                                   # factors after emergency
                                                   # diversity injection.
        seed: Optional[int] = None,  # Optional RNG seed for reproducibility.
                                     # Seeds both random and numpy.random.
        max_cache_size: int = 1000,  # Maximum sequences kept in the property
                                     # evaluation cache (LRU eviction). The
                                     # sequence-id map shares this bound.
        # --- Elite pool / population diversity (Round 5, #13) ---
        elite_pool_size: int = 1,  # K elite sequences kept as parents. Default
                                   # 1 reproduces single-best behavior. K>1
                                   # explores neighborhoods of multiple strong
                                   # seeds without adding evaluations.
        parent_selection: str = "weighted",  # "weighted" (rank-weighted),
                                              # "uniform", or "best" (always
                                              # use top of pool). Only matters
                                              # when elite_pool_size > 1.
        aa_fraction_ranges: Optional[Dict[Any, Tuple[float, float]]] = None,
        # Optional amino-acid fraction ranges. Keys may be single residues
        # (e.g. ``'A'``) or residue groups (e.g. ``('W', 'F', 'Y')`` or
        # ``'WFY'``). When provided, candidate generation (initial seeds,
        # point mutations, and multi-mutations) is constrained so each
        # configured per-AA or grouped fraction stays within ``[low, high]``.
        # Residues not listed in a per-AA entry are implicitly bounded
        # ``(0.0, 1.0)``. Shuffling already preserves composition so it is
        # unaffected.
        debugging: bool = False,  # Enable extra logging for debugging
        update_interval: int = 1  # Update progress bar every N iterations
    ):
        """
        Initialize the improved optimizer.
        
        Parameters
        ----------
        target_length : int
            Target sequence length
        max_iterations : int
            Maximum optimization iterations
        verbose : bool
            Enable progress reporting
        enable_adaptive_scaling : bool
            Enable adaptive property scaling
        enable_shuffling : bool
            Enable sequence shuffling for diversity
        shuffle_frequency : int
            How often to shuffle (every N iterations) - DEFAULT: 50 for better diversity
        stagnation_threshold : int
            Iterations before considering property stagnant - DEFAULT: 25 for faster response
        convergence_tolerance : float
            Convergence criterion for early stopping - DEFAULT: 1e-4 (less strict than before)
        convergence_window : int
            Number of recent iterations to check for convergence
        enable_early_convergence : bool
            Enable early stopping when convergence is detected
        convergence_patience : int
            Number of iterations to wait after convergence before stopping - DEFAULT: 20
        error_tolerance : float, optional
            Stop optimization when total error falls below this value (None = disabled)
        enable_error_tolerance : bool
            Enable error tolerance early stopping - DEFAULT: True
        num_candidates : int
            Number of candidate sequences to generate per iteration - DEFAULT: 5 to keep things fast. 
        num_starting_candidates : int
            Number of diverse starting sequences to evaluate - DEFAULT: 100 for better starting point
        min_mutations : int
            Minimum number of mutations for multi-mutation candidates
        max_mutations : int
            Maximum number of mutations for multi-mutation candidates - DEFAULT: 15
        mutation_ratio : int
            Sequence length divisor for calculating maximum mutations - DEFAULT: 10 for more mutations
        global_shuffle_probability : float
            Probability of performing global vs local shuffling - DEFAULT: 0.4 for more diversity
        shuffle_window_size : int
            Size of window for local shuffling - DEFAULT: 15 for larger shuffle regions
        improvement_history_size : int
            Number of recent improvements to track per property
        error_history_size : int
            Number of recent error values to store
        min_trend_samples : int
            Minimum samples required for trend calculation
        max_distance_factor : float
            Maximum scaling factor based on distance from target
        distance_offset : float
            Offset added to raw error for distance calculation
        stagnation_multiplier : float
            Multiplier for stagnation score in adaptive scaling
        low_contribution_threshold : float
            Threshold for identifying low-contributing properties
        high_error_threshold : float
            Threshold for identifying high-error properties
        boost_factor : float
            Factor to boost underperforming properties
        scale_momentum : float
            Momentum factor for scale smoothing (0-1)
        scale_learning_rate : float
            Learning rate for scale updates (0-1)
        min_scale : float
            Minimum allowed property scale
        max_scale : float
            Maximum allowed property scale
        stagnation_improvement_threshold : float
            Minimum improvement threshold to avoid stagnation detection
        improvement_trend_threshold : float
            Threshold for detecting property improvement trends
        stagnation_boost_factor : float
            Factor to boost properties during stagnation recovery
        enforce_raw_monotonicity : bool
            If True, only accept candidate sequences whose unsatisfied raw-error
            sum does not exceed the current best (with optional slack). Off by
            default since it can prevent legitimate multi-objective trade-offs.
        raw_monotonicity_slack : float
            When enforce_raw_monotonicity is True, allow a fractional slack
            (e.g. 0.05 means up to 5%% worse raw error is acceptable).
        scale_freeze_window : int
            After stagnation recovery, freeze adaptive scaling updates for this
            many iterations so recovery boosts have time to take effect. 0
            falls back to ``stagnation_threshold``.
        max_norm_boost : float
            Maximum value for the recovery-only normalization boost multiplier.
        norm_boost_decay : float
            Per-iteration decay applied to recovery norm boosts (1.0 disables).
        recompute_norm_on_emergency : bool
            If True, blend recomputed normalization factors with the existing
            ones after a successful emergency diversity injection.
        seed : int, optional
            Seed for ``random`` and ``numpy.random``. Two runs with the same
            ``seed`` and identical inputs will produce the same trajectory and
            output sequence (assuming external libraries used by properties
            also respect or ignore RNG state).
        max_cache_size : int
            Maximum number of sequences retained in the property-evaluation
            cache. Both the cache and the sequence-id map share this bound
            using LRU eviction.
        elite_pool_size : int
            Number of elite sequences (ranked by unsatisfied raw-error sum,
            which is scale-invariant) kept as parents for mutation. Default 1
            reproduces the original single-best behavior. Larger pools
            (e.g. 3-5) help on multi-objective tasks by maintaining
            structurally diverse strong seeds. No extra property evaluations
            are performed; the pool is updated from the candidates already
            evaluated each iteration.
        parent_selection : str
            Strategy for picking a parent from the elite pool when
            ``elite_pool_size > 1``: ``"weighted"`` (inverse-rank weighted),
            ``"uniform"`` (uniform random), or ``"best"`` (always top of
            pool). Ignored when ``elite_pool_size <= 1``.
        aa_fraction_ranges : dict, optional
            Mapping of amino-acid keys to ``(low, high)`` fractional bounds
            (each in ``[0, 1]``, inclusive). Keys may be:

            * a single-letter string -- per-AA bound, e.g. ``'A': (0.05, 0.15)``
            * a multi-letter string, tuple, set, or frozenset of single
              letters -- bound on the SUM of those AAs' fractions, e.g.
              ``('W', 'F', 'Y'): (0.05, 0.15)`` or ``'WFY': (0.05, 0.15)``.

            When provided, all generated candidates -- initial seeds, point
            mutations, multi-mutations, and emergency-injection seeds -- are
            constrained so each per-AA and per-group fraction stays in its
            configured range. Per-AA bounds for unlisted residues are
            unconstrained. Shuffling preserves composition and is unaffected.
            Bounds are enforced as integer counts derived from
            ``target_length`` (``floor(low * L)`` to ``ceil(high * L)``).
        debugging : bool
            Enable extra logging for debugging
        update_interval : int
            Update progress bar every N iterations
        """
        self.target_length = target_length
        self.max_iterations = max_iterations
        self.verbose = verbose
        self.enable_adaptive_scaling = enable_adaptive_scaling
        self.enable_shuffling = enable_shuffling
        self.shuffle_frequency = shuffle_frequency
        self.stagnation_threshold = stagnation_threshold
        self.convergence_tolerance = convergence_tolerance
        self.convergence_window = convergence_window
        self.enable_early_convergence = enable_early_convergence
        self.convergence_patience = convergence_patience
        
        # Error tolerance parameters
        self.error_tolerance = error_tolerance
        self.enable_error_tolerance = enable_error_tolerance
        
        # Mutation parameters
        self.num_candidates = num_candidates
        self.num_starting_candidates = num_starting_candidates
        self.min_mutations = min_mutations
        self.max_mutations = max_mutations
        self.mutation_ratio = mutation_ratio
        
        # Shuffling parameters
        self.global_shuffle_probability = global_shuffle_probability
        self.shuffle_window_size = shuffle_window_size
        
        # History tracking parameters
        self.improvement_history_size = improvement_history_size
        self.error_history_size = error_history_size
        
        # Adaptive scaling parameters
        self.min_trend_samples = min_trend_samples
        self.max_distance_factor = max_distance_factor
        self.distance_offset = distance_offset
        self.stagnation_multiplier = stagnation_multiplier
        self.low_contribution_threshold = low_contribution_threshold
        self.high_error_threshold = high_error_threshold
        self.boost_factor = boost_factor
        self.scale_momentum = scale_momentum
        self.scale_learning_rate = scale_learning_rate
        self.min_scale = min_scale
        self.max_scale = max_scale
        
        # Stagnation parameters
        self.stagnation_improvement_threshold = stagnation_improvement_threshold
        self.improvement_trend_threshold = improvement_trend_threshold
        self.stagnation_boost_factor = stagnation_boost_factor
        
        # Multi-objective behavior controls
        self.enforce_raw_monotonicity = enforce_raw_monotonicity
        self.raw_monotonicity_slack = raw_monotonicity_slack
        self.scale_freeze_window = scale_freeze_window
        self.max_norm_boost = max_norm_boost
        self.norm_boost_decay = norm_boost_decay
        self.recompute_norm_on_emergency = recompute_norm_on_emergency
        # Seed RNGs for reproducibility (#16). Note that some property classes
        # may call into external libraries with their own RNG state; users
        # seeking strict determinism should also seed those.
        self.seed = seed
        if seed is not None:
            random.seed(seed)
            np.random.seed(seed)
        
        # Cache size bound (#15/#24)
        self.max_cache_size = max(int(max_cache_size), 16)
        
        # Elite pool config (#13). Pool entries are (raw_error_unsatisfied,
        # weighted_error, sequence) tuples. The first field is the sort key
        # and is scale-invariant so adaptive scaling does not invalidate the
        # pool ordering.
        self.elite_pool_size = max(int(elite_pool_size), 1)
        if parent_selection not in ("weighted", "uniform", "best"):
            raise ValueError(
                f"parent_selection must be 'weighted', 'uniform', or 'best'; got {parent_selection!r}"
            )
        self.parent_selection = parent_selection
        self._elite_pool: List[Tuple[float, float, str]] = []
        # Iteration index up to which adaptive scaling is frozen.
        self._scale_freeze_until: int = -1
        # Per-property recovery-only normalization boost multipliers (decayed
        # back toward 1.0 every iteration so they fade once progress resumes).
        self.norm_boost_factors: Dict[str, float] = {}
        
        # Properties and state
        self.properties: List[ProteinProperty] = []
        self.property_scales: Dict[str, float] = {}
        self.property_improvements: Dict[str, deque] = {}
        self.best_error_history: deque = deque(maxlen=self.error_history_size)
        # Raw-error history (sum of raw_error over unsatisfied properties) is
        # used by stagnation/convergence checks because it is invariant to
        # changes in property_scales / normalization_factors during a run.
        self.best_raw_error_history: deque = deque(maxlen=self.error_history_size)
        
        # Dynamic normalization based on initial sequence
        self.initial_property_errors: Dict[str, float] = {}
        self.normalization_factors: Dict[str, float] = {}
        # Whether _calculate_initial_normalization has populated
        # normalization_factors yet. Reset whenever the property set or seed
        # changes so that run() can decide whether to (re)compute.
        self._normalization_initialized: bool = False
        
        # Optimization state
        self.best_sequence: str = ""
        self.best_error: float = float('inf')
        self.iteration: int = 0
        self.convergence_counter: int = 0
        self.starting_sequence: str = None
        
        # Evaluation cache for performance optimization
        # Replace hash-based cache with ID-based cache
        self._raw_property_cache: Dict[int, Dict[str, Dict[str, Any]]] = {}
        self._sequence_id_map: Dict[str, int] = {}  # sequence -> unique ID
        self._next_sequence_id: int = 0
        self._cache_hits: int = 0
        self._cache_misses: int = 0
        
        # Current best sequence property info cache (updated when best_sequence changes)
        self._best_sequence_property_info: Optional[Dict[str, Dict[str, float]]] = None

        # debugging flag
        self.debugging = debugging

        self.update_interval = update_interval

        # store aa list for faster mutation functions (single source of truth:
        # the amino_acids data module)
        self.aa_list = list(amino_acids.amino_acid_dat.keys())
        self.aa_list_without: Dict[str, List[str]] = {}  # Cache filtered AA lists
        for aa in self.aa_list:
            self.aa_list_without[aa] = [x for x in self.aa_list if x != aa]

        # Per-AA fraction range constraints. When ``aa_fraction_ranges`` is
        # set we precompute integer count bounds (min_count, max_count) for
        # every AA in ``self.aa_list`` and any group (multi-AA) constraints,
        # so constrained mutation/initial generation can run in pure-Python
        # without recomputing on each call.
        self.aa_fraction_ranges = aa_fraction_ranges
        self._aa_count_bounds: Optional[Dict[str, Tuple[int, int]]] = None
        self._aa_group_bounds: Optional[List[Tuple[Tuple[str, ...], int, int]]] = None
        self._aa_to_group_indices: Dict[str, List[int]] = {aa: [] for aa in self.aa_list}
        if aa_fraction_ranges is not None:
            self._aa_count_bounds, self._aa_group_bounds = self._compute_aa_count_bounds(
                aa_fraction_ranges, self.target_length
            )
            for gi, (members, _gmin, _gmax) in enumerate(self._aa_group_bounds):
                for m in members:
                    self._aa_to_group_indices[m].append(gi)

        # track property names
        self._property_names: List[str] = []

        # Configure logging
        self._configure_logging()

    def _evict_cache_to_size(self) -> None:
        """Centralized LRU-style eviction for both the sequence-id map and the
        raw property cache. Evicts oldest entries until both are under the
        configured ``max_cache_size``.
        """
        target = self.max_cache_size
        # Evict from the id map first; cascade to the raw cache.
        while len(self._sequence_id_map) > target:
            oldest_seq = next(iter(self._sequence_id_map))
            old_id = self._sequence_id_map.pop(oldest_seq)
            self._raw_property_cache.pop(old_id, None)
        # If for any reason the cache still exceeds target (orphans), trim it.
        while len(self._raw_property_cache) > target:
            oldest_id = next(iter(self._raw_property_cache))
            del self._raw_property_cache[oldest_id]
    
    def _get_sequence_id(self, sequence: str) -> int:
        """
        Get or create a unique integer ID for a sequence.
        
        Uses a simple counter to assign IDs. Much faster than hashing
        because we only need to hash once per unique sequence seen.
        """
        # Check if we've seen this sequence before
        if sequence in self._sequence_id_map:
            return self._sequence_id_map[sequence]
        
        # Assign new ID
        seq_id = self._next_sequence_id
        self._next_sequence_id += 1
        self._sequence_id_map[sequence] = seq_id
        
        # Bound memory using the unified LRU helper.
        if len(self._sequence_id_map) > self.max_cache_size:
            self._evict_cache_to_size()
        
        return seq_id


    def set_initial_sequence(self, sequence: str) -> None:
        """Set the initial sequence and calculate normalization factors"""
        if len(sequence) != self.target_length:
            raise ValueError(f"Initial sequence length {len(sequence)} does not match target length {self.target_length}")
        
        # If per-AA fraction ranges are configured, the starting sequence
        # must satisfy them -- otherwise constrained mutation would be
        # unable to repair it without first violating bounds.
        if self._aa_count_bounds is not None:
            counts: Dict[str, int] = {}
            for c in sequence:
                counts[c] = counts.get(c, 0) + 1
            for aa, (lo, hi) in self._aa_count_bounds.items():
                c = counts.get(aa, 0)
                if c < lo or c > hi:
                    frac = c / len(sequence)
                    raise ValueError(
                        f"Starting sequence violates aa_fraction_ranges for {aa!r}: "
                        f"count={c} (fraction={frac:.3f}), required count in [{lo}, {hi}]."
                    )
            for members, gmin, gmax in (self._aa_group_bounds or []):
                gc = sum(counts.get(a, 0) for a in members)
                if gc < gmin or gc > gmax:
                    frac = gc / len(sequence)
                    raise ValueError(
                        f"Starting sequence violates aa_fraction_ranges group {members}: "
                        f"count={gc} (fraction={frac:.3f}), required count in [{gmin}, {gmax}]."
                    )
        
        self.starting_sequence = sequence
        
        # Calculate initial normalization factors
        self._calculate_initial_normalization(sequence)
        
        # Evaluate initial sequence
        initial_error, _ = self._evaluate_sequence(sequence)
        self.best_error = initial_error
        self.best_error_history.append(initial_error)
        
        if self.debugging:
            self.logger.info(f"Initial sequence set. Initial weighted error: {initial_error:.6f}")
    
    def _configure_logging(self) -> None:
        """Configure a per-instance logger.
        
        Avoids mutating the root logger (#23). When ``verbose`` is on we attach
        a stream handler with our compact format directly to ``self.logger``.
        """
        self.logger = logging.getLogger(f"{__name__}.SequenceOptimizer.{id(self)}")
        # Don't propagate to root so we don't double-print if root is also
        # configured by an embedding application.
        self.logger.propagate = False
        # Idempotent: clear any handlers we previously attached on this logger
        # in case the same instance reconfigures.
        for h in list(self.logger.handlers):
            self.logger.removeHandler(h)
        if self.verbose:
            handler = logging.StreamHandler()
            handler.setFormatter(logging.Formatter('%(message)s'))
            self.logger.addHandler(handler)
            self.logger.setLevel(logging.INFO)
        else:
            self.logger.addHandler(logging.NullHandler())
            self.logger.setLevel(logging.WARNING)
    
    def add_property(self, property_class: type, *args: Any, tolerance: float = 0.0, **kwargs: Any) -> None:
        """
        Add a property to optimize
        
        Parameters
        ----------
        property_class : type
            The property class to instantiate
        *args : Any
            Positional arguments for the property constructor
        tolerance : float, optional
            Error tolerance for this property. If raw_error <= tolerance, 
            the property is considered satisfied and won't contribute to total error.
            Default: 0.0. NOTE: with the default ``tolerance=0.0``, a property
            is only considered "satisfied" when ``|raw_error| < 1e-8`` (a tiny
            numerical epsilon). For floating-point property values this is
            effectively unsatisfiable, so set a small positive tolerance for
            real-valued targets such as hydropathy or Rg.
        **kwargs : Any
            Keyword arguments for the property constructor
        """
        prop = property_class(*args, **kwargs)
        
        # Safety net: if the user supplied ``weight=`` as a kwarg but the
        # subclass __init__ swallowed/ignored it (a common footgun in
        # ``CustomProperty`` subclasses that pass arguments positionally to
        # super().__init__), enforce the user-provided weight on the
        # instance. The ``weight`` setter in ``ProteinProperty`` validates the
        # value, so this is the same as if the constructor had honored it.
        if 'weight' in kwargs and getattr(prop, 'weight', 1.0) != kwargs['weight']:
            if self.verbose:
                self.logger.warning(
                    f"{property_class.__name__}.__init__ did not propagate weight="
                    f"{kwargs['weight']} to the instance (current weight={prop.weight}). "
                    f"Forcing weight={kwargs['weight']}. Check your subclass's "
                    f"super().__init__(...) call -- pass weight as a keyword argument."
                )
            prop.weight = kwargs['weight']
        
        # Set the per-property tolerance
        prop._tolerance = tolerance
        
        # Set sequence length hint for length-dependent calculations
        prop._set_sequence_length_hint(self.target_length)
        
        self.properties.append(prop)
        
        prop_name = prop.tracking_property_name
        if prop_name in self._property_names:
            if prop.multi_target==False:
                raise ValueError(f"Property {prop_name} already added.")
            else:
                # count number of properties with this name
                count = sum(1 for p in self.properties if prop_name in p.__class__.__name__)
                # add a number to the name
                prop_name = f"{prop_name}_{count}"
                prop.tracking_property_name = prop_name
        self._property_names.append(prop_name)
        self.property_scales[prop_name] = 1.0
        self.property_improvements[prop_name] = deque(maxlen=self.improvement_history_size)
        
        # Clear cache when new properties are added since cached values may not include the new property
        self._clear_evaluation_cache()
        # Adding a property invalidates any previously computed normalization.
        self._normalization_initialized = False
        
        if self.verbose:
            tolerance_info = f", tolerance: ±{tolerance}" if tolerance > 0 else ""
            self.logger.info(f"Added property: {prop_name} (target: {prop.target_value}, weight: {prop.weight}{tolerance_info})")
    
    def _assemble_property_info(self, cached_raw_values: Dict[str, Dict[str, Any]]) -> Tuple[float, Dict[str, Dict[str, float]]]:
        """
        Assemble property info from cached raw values with current scaling/normalization.
        Extracted to avoid code duplication between single and batch evaluation.
        """
        property_info = {}
        weighted_errors = []
        # Bind frequently-accessed dicts/methods to locals to avoid repeated
        # attribute lookups in the per-property loop (#25). The previous
        # implementation rebuilt two intermediate dicts per call.
        norm_factors = self.normalization_factors
        scales = self.property_scales
        norm_boosts = self.norm_boost_factors
        
        for prop, prop_name in zip(self.properties, self._property_names):
            
            raw_data = cached_raw_values[prop.tracking_property_name]
            raw_value = raw_data['raw_value']
            raw_error = raw_data['raw_error']
            tolerance = raw_data['tolerance']
            
            within_tolerance = (tolerance == 0.0 and np.abs(raw_error) < 1e-8) or (raw_error <= tolerance)
            
            normalization_factor = norm_factors.get(prop_name, 1.0)
            # Apply optional recovery-only boost on top of the clipped
            # normalization factor (lets repeated stagnation recovery escalate
            # beyond the [0.1, 10.0] clip without polluting the base value).
            boost = norm_boosts.get(prop_name, 1.0)
            effective_norm = normalization_factor * boost
            normalized_error = raw_error * effective_norm
            
            scale = scales.get(prop_name, 1.0)
            weighted_error = normalized_error * prop.weight * scale
            
            effective_weighted_error = 0.0 if within_tolerance else weighted_error
            
            property_info[prop.tracking_property_name] = {
                'raw_value': raw_value,
                'target_value': raw_data['target_value'],
                'raw_error': raw_error,
                'normalized_error': normalized_error,
                'normalization_factor': normalization_factor,
                'norm_boost': boost,
                'effective_normalization': effective_norm,
                'weighted_error': weighted_error,
                'effective_weighted_error': effective_weighted_error,
                'scale': scale,
                'tolerance': tolerance,
                'within_tolerance': within_tolerance,
                'constraint_type': prop.constraint_type,
                'weight': prop.weight
            }
            
            weighted_errors.append(effective_weighted_error)
        
        total_error = sum(weighted_errors)
        return total_error, property_info

    def _evaluate_sequence(self, sequence: Union[str, List[str]], use_cache: bool = True) -> Union[Tuple[float, Dict[str, Dict[str, float]]], List[Tuple[float, Dict[str, Dict[str, float]]]]]:
        """Evaluate sequence(s) and return error and property info."""
        if isinstance(sequence, list):
            return self._evaluate_sequences_batch(sequence, use_cache=use_cache)
        
        raw_cache_key = self._get_sequence_id(sequence)
        
        if use_cache and raw_cache_key in self._raw_property_cache:
            cached_raw_values = self._raw_property_cache[raw_cache_key]
            self._cache_hits += 1
        else:
            protein = sparrow.Protein(sequence)
            cached_raw_values = {}
            
            for prop in self.properties:
                prop_name = prop.tracking_property_name
                raw_value, raw_error = prop.calculate(protein)
                cached_raw_values[prop_name] = {
                    'raw_value': raw_value,
                    'raw_error': raw_error,
                    'target_value': prop.target_value,
                    'tolerance': getattr(prop, '_tolerance', 0.0)
                }
            
            if use_cache:
                self._raw_property_cache[raw_cache_key] = cached_raw_values
                if len(self._raw_property_cache) > self.max_cache_size:
                    self._evict_cache_to_size()
            
            self._cache_misses += 1
        
        # Use shared assembly function
        return self._assemble_property_info(cached_raw_values)

    def _evaluate_sequences_batch(self, sequences: List[str], use_cache: bool = True) -> List[Tuple[float, Dict[str, Dict[str, float]]]]:
        """Evaluate multiple sequences in batch for improved performance."""
        if not sequences:
            return []
        
        sequence_ids = [self._get_sequence_id(seq) for seq in sequences]
        
        # ... existing cache lookup code ...
        cached_results = {}
        uncached_sequences = []
        uncached_indices = []

        for i, seq_id in enumerate(sequence_ids):
            if use_cache and seq_id in self._raw_property_cache:
                cached_results[i] = self._raw_property_cache[seq_id]
                self._cache_hits += 1
            else:
                uncached_sequences.append(sequences[i])
                uncached_indices.append(i)
                self._cache_misses += 1
        
        # ... existing property calculation code ...
        uncached_raw_values = {}
        if uncached_sequences:
            proteins = [sparrow.Protein(seq) for seq in uncached_sequences]
            
            for prop in self.properties:
                prop_name = prop.tracking_property_name
                
                if prop.calculate_in_batch:
                    batch_results = prop.calculate_batch(proteins)
                    for i, (raw_value, raw_error) in enumerate(batch_results):
                        seq_idx = uncached_indices[i]
                        if seq_idx not in uncached_raw_values:
                            uncached_raw_values[seq_idx] = {}
                        uncached_raw_values[seq_idx][prop_name] = {
                            'raw_value': raw_value,
                            'raw_error': raw_error,
                            'target_value': prop.target_value,
                            'tolerance': getattr(prop, '_tolerance', 0.0)
                        }
                else:
                    for i, protein in enumerate(proteins):
                        seq_idx = uncached_indices[i]
                        raw_value, raw_error = prop.calculate(protein)
                        if seq_idx not in uncached_raw_values:
                            uncached_raw_values[seq_idx] = {}
                        uncached_raw_values[seq_idx][prop_name] = {
                            'raw_value': raw_value,
                            'raw_error': raw_error,
                            'target_value': prop.target_value,
                            'tolerance': getattr(prop, '_tolerance', 0.0)
                        }
            
            if use_cache:
                for i, seq_idx in enumerate(uncached_indices):
                    seq_id = sequence_ids[seq_idx]
                    self._raw_property_cache[seq_id] = uncached_raw_values[seq_idx]
                if len(self._raw_property_cache) > self.max_cache_size:
                    self._evict_cache_to_size()
        
        results = []
        for i in range(len(sequences)):
            if i in cached_results:
                cached_raw_values = cached_results[i]
            else:
                cached_raw_values = uncached_raw_values[i]
            
            total_error, property_info = self._assemble_property_info(cached_raw_values)
            results.append((total_error, property_info))
        
        return results
    
    def _get_best_sequence_property_info(self) -> Dict[str, Dict[str, float]]:
        """
        Get property info for the current best sequence (cached).
        Now caches by sequence ID to avoid re-evaluation.
        """
        # Get sequence ID (consistent for same sequence content)
        current_seq_id = self._get_sequence_id(self.best_sequence)
        
        # Check if we've already cached info for this sequence
        if (self._best_sequence_property_info is not None and 
            hasattr(self, '_best_sequence_property_info_cache_id') and
            self._best_sequence_property_info_cache_id == current_seq_id):
            return self._best_sequence_property_info
        
        # Evaluate and cache
        _, self._best_sequence_property_info = self._evaluate_sequence(self.best_sequence)
        self._best_sequence_property_info_cache_id = current_seq_id
        return self._best_sequence_property_info

    def _update_best_sequence(self, new_sequence: str, new_error: float) -> None:
        """Update the best sequence and invalidate cached property info."""
        self.best_sequence = new_sequence
        self.best_error = new_error
        # Invalidate cached property info
        self._best_sequence_property_info = None
        if hasattr(self, '_best_sequence_property_info_cache_id'):
            delattr(self, '_best_sequence_property_info_cache_id')
    
    def _calculate_initial_normalization(self, initial_sequence: str) -> None:
        """
        Calculate dynamic normalization factors based on initial sequence errors.
        This balances properties based on their initial distance from targets.
        """
        protein = sparrow.Protein(initial_sequence)
        initial_errors = []
        
        # Calculate initial errors for all properties
        for prop in self.properties:
            prop_name = prop.tracking_property_name
            raw_value, raw_error = prop.calculate(protein)
            
            # Store initial error (use small minimum to avoid division by zero)
            initial_error = max(raw_error, 1e-6)
            self.initial_property_errors[prop_name] = initial_error
            initial_errors.append(initial_error)
        
        # Calculate normalization factors
        if len(initial_errors) > 1:
            # Use the geometric mean of initial errors as reference scale
            # This provides a balanced reference point across all properties
            reference_error = np.exp(np.mean(np.log(initial_errors)))
            
            for prop in self.properties:
                prop_name = prop.tracking_property_name
                initial_error = self.initial_property_errors[prop_name]
                
                # Normalization factor scales errors to the reference level
                # Properties with larger initial errors get smaller factors (scaled down)
                # Properties with smaller initial errors get larger factors (scaled up)
                normalization_factor = reference_error / initial_error
                
                # Apply reasonable bounds to prevent extreme scaling
                normalization_factor = np.clip(normalization_factor, 0.1, 10.0) 
                
                self.normalization_factors[prop_name] = normalization_factor
                
                if self.debugging:
                    self.logger.info(f"  {prop_name}: initial_error={initial_error:.4f}, "
                                   f"norm_factor={normalization_factor:.3f}")
        else:
            # Single property - no normalization needed
            for prop in self.properties:
                prop_name = prop.tracking_property_name
                self.normalization_factors[prop_name] = 1.0
        
        self._normalization_initialized = True

    def _blend_normalization_with_current_sequence(self, blend: float = 0.5) -> None:
        """Recompute normalization on the current best_sequence and blend it
        with the existing factors.
        
        Used after emergency diversity injection so that normalization adapts
        to the new sequence regime instead of carrying forward potentially
        misleading factors from the prior seed (#12).
        
        Parameters
        ----------
        blend : float
            Weight on the freshly computed factors. ``new = blend*fresh + (1-blend)*old``.
        """
        if not self.properties:
            return
        protein = sparrow.Protein(self.best_sequence)
        fresh_initial_errors: Dict[str, float] = {}
        for prop in self.properties:
            _, raw_error = prop.calculate(protein)
            fresh_initial_errors[prop.tracking_property_name] = max(raw_error, 1e-6)
        if len(fresh_initial_errors) <= 1:
            return
        reference_error = float(np.exp(np.mean(np.log(list(fresh_initial_errors.values())))))
        for prop in self.properties:
            name = prop.tracking_property_name
            fresh_factor = float(np.clip(reference_error / fresh_initial_errors[name], 0.1, 10.0))
            old_factor = self.normalization_factors.get(name, 1.0)
            blended = blend * fresh_factor + (1.0 - blend) * old_factor
            self.normalization_factors[name] = float(np.clip(blended, 0.1, 10.0))
        # Invalidate any cached weighted property info since normalization changed.
        self._best_sequence_property_info = None
        if hasattr(self, '_best_sequence_property_info_cache_id'):
            delattr(self, '_best_sequence_property_info_cache_id')
        if self.debugging:
            self.logger.info(
                "   🔁 Blended normalization with new emergency seed: "
                + ", ".join(f"{n}={self.normalization_factors[n]:.3f}" for n in fresh_initial_errors)
            )
    
    def _update_property_scaling(self, property_info: Dict[str, Dict[str, float]]) -> None:
        """Update adaptive scaling based on property performance"""
        if not self.enable_adaptive_scaling or len(self.properties) <= 1:
            return
        # Honor recovery freeze: skip updates while a stagnation-recovery boost
        # is meant to be in effect (#8).
        if self.iteration < self._scale_freeze_until:
            return
        
        # Calculate relative error contributions using effective_weighted_error so
        # already-satisfied properties (within tolerance) do not deflate the
        # contribution share of unsatisfied properties.
        total_error = sum(
            info.get('effective_weighted_error', info['weighted_error'])
            for info in property_info.values()
        )
        if total_error <= 0:
            return
        
        for prop_name, info in property_info.items():
            # Track improvement history using normalized errors for fair comparison
            self.property_improvements[prop_name].append(info['normalized_error'])
            
            # Calculate stagnation score
            improvements = list(self.property_improvements[prop_name])
            if len(improvements) >= self.min_trend_samples:
                # Closed-form least-squares slope (avoids np.polyfit's SVD
                # call in the per-iteration hot path; ~10-50x faster for
                # short windows).
                _y = np.asarray(improvements, dtype=float)
                _n = _y.shape[0]
                _x = np.arange(_n, dtype=float)
                _xm = (_n - 1) * 0.5
                _ym = float(_y.mean())
                _dx = _x - _xm
                _denom = float((_dx * _dx).sum())
                recent_trend = float((_dx * (_y - _ym)).sum()) / _denom if _denom > 0 else 0.0
                # Normalize slope by the window's mean magnitude so properties
                # at very different scales are comparable (#10).
                window_mean = max(abs(_ym), 1e-8)
                relative_trend = recent_trend / window_mean
                stagnation_score = max(0, -relative_trend)  # Negative slope = improvement
            else:
                stagnation_score = 1.0
            
            # Calculate relative contribution (effective error so satisfied props
            # contribute 0 and unsatisfied props share the full denominator)
            contribution = info.get('effective_weighted_error', info['weighted_error']) / total_error
            
            # Simple adaptive scaling formula
            current_scale = self.property_scales[prop_name]
            
            # Boost properties that are:
            # 1. Far from target (high normalized error)
            # 2. Stagnant (not improving)
            # 3. Contributing little to total error despite high normalized error
            
            # Use normalized error for fair comparison across property types
            distance_factor = min(self.max_distance_factor, info['normalized_error'] + self.distance_offset)
            stagnation_factor = 1.0 + stagnation_score * self.stagnation_multiplier
            
            # If property has low contribution but high normalized error, boost it
            if contribution < self.low_contribution_threshold and info['normalized_error'] > self.high_error_threshold:
                boost_factor = self.boost_factor
            else:
                boost_factor = 1.0
            
            # Calculate new scale with momentum (smooth changes)
            target_scale = distance_factor * stagnation_factor * boost_factor
            new_scale = self.scale_momentum * current_scale + self.scale_learning_rate * target_scale
            
            # Apply bounds
            self.property_scales[prop_name] = np.clip(new_scale, self.min_scale, self.max_scale)

        # Scaling changed -> any cached property info for the current best
        # sequence is now stale (weighted_error / scale fields). Invalidate so
        # the next consumer (stagnation recovery, progress display) recomputes.
        self._best_sequence_property_info = None
        if hasattr(self, '_best_sequence_property_info_cache_id'):
            delattr(self, '_best_sequence_property_info_cache_id')

    # ------------------------------------------------------------------
    # Elite pool helpers (#13)
    # ------------------------------------------------------------------
    @staticmethod
    def _weighted_raw_sum(property_info: Dict[str, Dict[str, Any]],
                          unsatisfied_only: bool = True) -> float:
        """Sum of ``raw_error * weight`` across properties.
        
        This is used wherever we need a scale-invariant priority signal that
        still respects user-provided ``weight`` values. ``raw_error`` is
        independent of the adaptive ``property_scales`` and
        ``normalization_factors`` (so it is stable across a run), but
        weighting by ``prop.weight`` ensures user priorities drive seed
        selection, elite-pool ranking, and stagnation contribution math
        instead of being silently averaged away.
        """
        if unsatisfied_only:
            return sum(
                info['raw_error'] * info.get('weight', 1.0)
                for info in property_info.values()
                if not info.get('within_tolerance', False)
            )
        return sum(
            info['raw_error'] * info.get('weight', 1.0)
            for info in property_info.values()
        )

    def _maybe_add_to_pool(self, sequence: str, raw_unsat: float, weighted: float) -> None:
        """Insert a sequence into the elite pool if it qualifies.
        
        The pool is ranked ascending by ``raw_unsat`` (sum of raw_error over
        unsatisfied properties). This metric is scale-invariant so the pool
        ordering does not become stale when adaptive scaling updates.
        Duplicates are skipped. O(K) per call, K small.
        """
        if self.elite_pool_size <= 1:
            return
        # Skip if the sequence is already present.
        for _, _, s in self._elite_pool:
            if s == sequence:
                return
        if len(self._elite_pool) < self.elite_pool_size:
            self._elite_pool.append((raw_unsat, weighted, sequence))
            self._elite_pool.sort(key=lambda x: x[0])
            return
        # Replace the worst entry if this candidate beats it.
        if raw_unsat < self._elite_pool[-1][0]:
            self._elite_pool[-1] = (raw_unsat, weighted, sequence)
            self._elite_pool.sort(key=lambda x: x[0])
    
    def _select_parent_from_pool(self) -> str:
        """Pick a parent sequence for the next round of mutations.
        
        Always falls back to ``self.best_sequence`` if the pool is empty or
        ``elite_pool_size <= 1`` (the default), preserving the original
        single-best mutation strategy.
        """
        if self.elite_pool_size <= 1 or not self._elite_pool:
            return self.best_sequence
        if self.parent_selection == "best":
            return self._elite_pool[0][2]
        if self.parent_selection == "uniform":
            return random.choice(self._elite_pool)[2]
        # rank-weighted: 1, 1/2, 1/3, ... favors better seeds but still
        # samples weaker ones for diversity.
        n = len(self._elite_pool)
        weights = [1.0 / (i + 1) for i in range(n)]
        return random.choices(
            [entry[2] for entry in self._elite_pool],
            weights=weights,
            k=1,
        )[0]

    # ------------------------------------------------------------------
    # AA fraction-range constraints
    # ------------------------------------------------------------------
    def _compute_aa_count_bounds(
        self,
        aa_fraction_ranges: Dict[Any, Tuple[float, float]],
        target_length: int,
    ) -> Tuple[Dict[str, Tuple[int, int]], List[Tuple[Tuple[str, ...], int, int]]]:
        """Validate ``aa_fraction_ranges`` and convert to integer count bounds.
        
        Keys may be:
          * a single-letter string (per-AA bound), e.g. ``'A': (0.05, 0.15)``
          * a multi-letter string, tuple, frozenset, or set of single letters
            (group bound on the SUM of those AAs), e.g.
            ``('W','F','Y'): (0.05, 0.15)`` or ``'WFY': (0.05, 0.15)``.
        
        AAs not constrained per-letter get ``(0, target_length)``. Min counts
        are ``floor(low*L)``; max counts are ``ceil(high*L)``. Raises if the
        per-AA totals or any group is structurally infeasible.
        
        Returns
        -------
        (per_aa_bounds, group_bounds) where ``group_bounds`` is a list of
        ``(members_tuple, min_count, max_count)``.
        """
        per_aa: Dict[str, Tuple[int, int]] = {}
        group_bounds: List[Tuple[Tuple[str, ...], int, int]] = []
        L = int(target_length)
        valid = set(self.aa_list)
        for key, rng in aa_fraction_ranges.items():
            # Resolve key to an ordered tuple of AAs.
            if isinstance(key, str):
                members = tuple(key)
            elif isinstance(key, (tuple, list, set, frozenset)):
                members = tuple(key)
            else:
                raise ValueError(
                    f"aa_fraction_ranges key {key!r} must be a string, tuple, "
                    f"list, set, or frozenset of single-letter amino acids."
                )
            if not members:
                raise ValueError(f"aa_fraction_ranges key {key!r} is empty.")
            for m in members:
                if m not in valid:
                    raise ValueError(
                        f"aa_fraction_ranges key {key!r} contains invalid "
                        f"amino acid {m!r}; must be one of {sorted(valid)}"
                    )
            if len(set(members)) != len(members):
                raise ValueError(
                    f"aa_fraction_ranges key {key!r} contains duplicate amino acids."
                )
            if not (isinstance(rng, (tuple, list)) and len(rng) == 2):
                raise ValueError(
                    f"aa_fraction_ranges[{key!r}] must be a (low, high) tuple, got {rng!r}"
                )
            low, high = float(rng[0]), float(rng[1])
            if not (0.0 <= low <= high <= 1.0):
                raise ValueError(
                    f"aa_fraction_ranges[{key!r}]=({low}, {high}) must satisfy "
                    f"0 <= low <= high <= 1"
                )
            min_count = int(math.floor(low * L))
            max_count = min(int(math.ceil(high * L)), L)
            if len(members) == 1:
                aa = members[0]
                if aa in per_aa:
                    raise ValueError(
                        f"aa_fraction_ranges has duplicate per-AA entry for {aa!r}"
                    )
                per_aa[aa] = (min_count, max_count)
            else:
                group_bounds.append((members, min_count, max_count))
        for aa in self.aa_list:
            if aa not in per_aa:
                per_aa[aa] = (0, L)
        # Per-AA feasibility.
        total_min = sum(b[0] for b in per_aa.values())
        total_max = sum(b[1] for b in per_aa.values())
        if total_min > L:
            raise ValueError(
                f"aa_fraction_ranges infeasible: sum of minimum counts "
                f"({total_min}) exceeds target_length ({L})."
            )
        if total_max < L:
            raise ValueError(
                f"aa_fraction_ranges infeasible: sum of maximum counts "
                f"({total_max}) is less than target_length ({L})."
            )
        # Group feasibility against per-AA caps/floors.
        for members, gmin, gmax in group_bounds:
            max_from_members = sum(per_aa[a][1] for a in members)
            min_from_members = sum(per_aa[a][0] for a in members)
            if gmin > max_from_members:
                raise ValueError(
                    f"aa_fraction_ranges group {members} infeasible: group min "
                    f"count {gmin} exceeds the sum of member per-AA maximums "
                    f"({max_from_members})."
                )
            if gmax < min_from_members:
                raise ValueError(
                    f"aa_fraction_ranges group {members} infeasible: group max "
                    f"count {gmax} is less than the sum of member per-AA "
                    f"minimums ({min_from_members})."
                )
            if gmin > L:
                raise ValueError(
                    f"aa_fraction_ranges group {members} infeasible: group min "
                    f"count {gmin} exceeds target_length ({L})."
                )
        return per_aa, group_bounds

    def _build_constrained_random_sequence(self, target_length: int) -> str:
        """Build a random sequence whose AA counts respect both per-AA and
        group bounds in ``_aa_count_bounds`` / ``_aa_group_bounds``.
        """
        per_aa = self._aa_count_bounds
        groups = self._aa_group_bounds or []
        aa_to_groups = self._aa_to_group_indices
        counts = {aa: per_aa[aa][0] for aa in self.aa_list}
        group_counts = [sum(counts[a] for a in members) for members, _gmin, _gmax in groups]
        remainder = target_length - sum(counts.values())
        aa_pool = list(self.aa_list)
        while remainder > 0:
            # First satisfy any group still below its minimum -- only members
            # of those groups can take this slot.
            underfill = [
                gi for gi, (_m, gmin, _gmax) in enumerate(groups)
                if group_counts[gi] < gmin
            ]
            if underfill:
                gi = random.choice(underfill)
                members = groups[gi][0]
                pool = list(members)
                random.shuffle(pool)
            else:
                random.shuffle(aa_pool)
                pool = aa_pool
            chosen = None
            for aa in pool:
                if counts[aa] >= per_aa[aa][1]:
                    continue
                # Reject if any containing group is already at its max.
                bad = False
                for gi in aa_to_groups[aa]:
                    if group_counts[gi] >= groups[gi][2]:
                        bad = True
                        break
                if bad:
                    continue
                chosen = aa
                break
            if chosen is None:
                # Should not happen given feasibility checks, but bail out
                # rather than infinite-loop.
                break
            counts[chosen] += 1
            for gi in aa_to_groups[chosen]:
                group_counts[gi] += 1
            remainder -= 1
        chars: List[str] = []
        for aa, c in counts.items():
            if c:
                chars.extend([aa] * c)
        random.shuffle(chars)
        return ''.join(chars)

    def _constrained_point_mutation(self, sequence: str) -> str:
        """Single point mutation that keeps every AA fraction in range."""
        return self._constrained_multi_mutations(sequence, 1)

    def _constrained_multi_mutations(self, sequence: str, n: int) -> str:
        """Apply up to ``n`` point mutations while respecting per-AA and
        group bounds.
        
        Maintains running ``counts`` and per-AA position lists so each
        mutation is O(20 * G) where G is the number of group constraints.
        If at any point no legal swap exists, returns the partially-mutated
        sequence.
        """
        if n <= 0:
            return sequence
        per_aa = self._aa_count_bounds
        groups = self._aa_group_bounds or []
        aa_to_groups = self._aa_to_group_indices
        aa_list = self.aa_list
        seq_list = list(sequence)
        counts: Dict[str, int] = {aa: 0 for aa in aa_list}
        aa_positions: Dict[str, List[int]] = {aa: [] for aa in aa_list}
        for i, c in enumerate(seq_list):
            counts[c] = counts.get(c, 0) + 1
            aa_positions.setdefault(c, []).append(i)
        group_counts = [sum(counts[a] for a in members) for members, _gmin, _gmax in groups]

        def _pair_legal(old_aa: str, new_aa: str) -> bool:
            old_g = aa_to_groups[old_aa]
            new_g = aa_to_groups[new_aa]
            if not old_g and not new_g:
                return True
            old_set = set(old_g)
            new_set = set(new_g)
            for gi in old_set | new_set:
                delta = (1 if gi in new_set else 0) - (1 if gi in old_set else 0)
                if delta == 0:
                    continue
                gc = group_counts[gi] + delta
                gmin = groups[gi][1]
                gmax = groups[gi][2]
                if gc < gmin or gc > gmax:
                    return False
            return True

        for _ in range(n):
            # Per-AA-feasible donors. Group-level legality is checked per pair.
            removable = [
                aa for aa in aa_list
                if counts[aa] > per_aa[aa][0] and aa_positions[aa]
            ]
            if not removable:
                break
            random.shuffle(removable)
            chosen_pair: Optional[Tuple[str, str]] = None
            for old_aa in removable:
                # Per-AA-feasible acceptors.
                acceptors = [
                    a for a in aa_list
                    if a != old_aa and counts[a] < per_aa[a][1]
                    and _pair_legal(old_aa, a)
                ]
                if acceptors:
                    chosen_pair = (old_aa, random.choice(acceptors))
                    break
            if chosen_pair is None:
                break
            old_aa, new_aa = chosen_pair
            positions = aa_positions[old_aa]
            idx = random.randrange(len(positions))
            pos = positions[idx]
            # O(1) removal: swap-pop.
            positions[idx] = positions[-1]
            positions.pop()
            seq_list[pos] = new_aa
            counts[old_aa] -= 1
            counts[new_aa] += 1
            aa_positions[new_aa].append(pos)
            # Update group counts incrementally.
            old_set = set(aa_to_groups[old_aa])
            new_set = set(aa_to_groups[new_aa])
            for gi in old_set - new_set:
                group_counts[gi] -= 1
            for gi in new_set - old_set:
                group_counts[gi] += 1
        return ''.join(seq_list)

    def _generate_candidates(self, sequence: str, iteration: int) -> List[str]:
        """Generate candidate sequences with improved diversity"""
        candidates = []
        
        # Always include current sequence
        candidates.append(sequence)
        
        # Calculate how many candidates for each type
        remaining = self.num_candidates - 1
        
        # Distribute candidates more strategically
        if self.enable_shuffling and iteration % self.shuffle_frequency == 0:
            # When shuffling, allocate more candidates to shuffled variants
            shuffle_count = remaining // 3  # 1/3 for shuffling
            multi_mut_count = remaining // 3  # 1/3 for multi-mutations
            single_mut_count = remaining - shuffle_count - multi_mut_count  # remainder for single
        else:
            # Normal distribution
            shuffle_count = 0
            multi_mut_count = remaining // 2  # 1/2 for multi-mutations
            single_mut_count = remaining - multi_mut_count  # remainder for single
        
        seq_length = len(sequence)

        # Choose mutation primitives. When per-AA fraction ranges are set we
        # need pure-Python mutators that respect those bounds; otherwise we
        # use the much faster Cython implementations.
        if self._aa_count_bounds is not None:
            point_mut = self._constrained_point_mutation
            multi_mut = self._constrained_multi_mutations
        else:
            point_mut = make_point_mutation
            multi_mut = make_multi_mutations

        # Generate single point mutations
        for _ in range(single_mut_count):
            candidates.append(point_mut(sequence))
        
        # Generate multiple mutations with varied mutation counts
        max_muts = max(
            self.min_mutations,
            min(self.max_mutations, seq_length // self.mutation_ratio),
        )
        
        for _ in range(multi_mut_count):
            # Use exponential distribution to favor smaller mutations
            if max_muts > self.min_mutations:
                weights = [2 ** (max_muts - i) for i in range(self.min_mutations, max_muts + 1)]
                num_muts = random.choices(range(self.min_mutations, max_muts + 1), weights=weights)[0]
            else:
                num_muts = self.min_mutations
            
            # Single function call instead of loop over make_point_mutation
            candidates.append(multi_mut(sequence, num_muts))

        # Add shuffled variants if enabled (already efficient)
        if self.enable_shuffling and iteration % self.shuffle_frequency == 0:
            for _ in range(shuffle_count):
                candidates.append(
                    shuffle_sequence(
                    sequence, 
                    window_size=self.shuffle_window_size,
                    global_shuffle_probability=self.global_shuffle_probability
                )
            )
        
        return candidates[:self.num_candidates]
    
    def _compute_convergence_state(self) -> Tuple[bool, float]:
        """
        Compute convergence state without mutating any optimizer state.
        
        Uses the *raw* error history when available (preferred, since raw
        error is invariant to changes in property scaling/normalization);
        falls back to the weighted error history otherwise.
        
        Returns
        -------
        (is_under_tolerance, current_std) :
            is_under_tolerance is True when the recent error history's std is
            below the configured tolerance (history must have at least
            convergence_window entries).
        """
        history = (
            self.best_raw_error_history
            if len(self.best_raw_error_history) >= self.convergence_window
            else self.best_error_history
        )
        if len(history) < self.convergence_window:
            return False, float('inf')
        recent_errors = list(history)[-self.convergence_window:]
        current_std = float(np.std(recent_errors))
        return current_std < self.convergence_tolerance, current_std

    def _check_convergence(self) -> bool:
        """
        Check if optimization has converged based on configurable criteria.
        
        Returns
        -------
        bool
            True if optimization has converged
        """
        if not self.enable_early_convergence:
            return False
        
        is_under_tolerance, _ = self._compute_convergence_state()
        if is_under_tolerance:
            self.convergence_counter += 1
            # Only consider converged after patience period
            return self.convergence_counter >= self.convergence_patience
        else:
            self.convergence_counter = 0  # Reset counter if not converged
            return False
    
    def _check_error_tolerance(self) -> bool:
        """
        Check if optimization should stop due to error tolerance being reached.
        
        For properties with per-property tolerances, uses their individual tolerance.
        For properties without per-property tolerances, uses the global error_tolerance.
        
        Returns
        -------
        bool
            True if error tolerance criterion is met (all properties satisfied)
        """
        if not self.enable_error_tolerance:
            return False
        if not self.properties:
            # Defense-in-depth: an empty property set should not be reported
            # as "tolerance reached".
            return False
        
        # Check if all properties are satisfied (either by per-property or global tolerance)
        property_info = self._get_best_sequence_property_info()
        
        for prop_name, info in property_info.items():
            per_property_tolerance = info.get('tolerance', 0.0)
            
            # If property has per-property tolerance, check against that
            if per_property_tolerance > 0.0:
                if not info.get('within_tolerance', False):
                    return False
            # Otherwise, check against global tolerance (if set)
            elif self.error_tolerance is not None:
                if info['raw_error'] > self.error_tolerance:
                    return False
            # If no tolerance set for this property, it must be exact (can't be satisfied)
            else:
                if info['raw_error'] > 1e-10:  # Allow for floating point precision
                    return False
        
        return True
    
    def _detect_stagnation(self) -> bool:
        """Detect if optimization is stagnant (raw-error-based).
        
        Uses best_raw_error_history rather than the weighted error history
        because adaptive scaling/normalization can shift weighted error
        independently of true progress.
        """
        if len(self.best_raw_error_history) < self.stagnation_threshold:
            return False
        
        initial_error = self.best_raw_error_history[-self.stagnation_threshold]
        final_error = self.best_raw_error_history[-1]

        if initial_error < 1e-10:
            return final_error < 1e-10  # Already at perfect score

        improvement = (initial_error - final_error) / initial_error
        
        return improvement < self.stagnation_improvement_threshold
    
    def _apply_stagnation_recovery(self) -> None:
        """Apply recovery measures when stagnant"""
        if self.debugging:
            self.logger.info("🔄 Applying stagnation recovery...")
        
        # Get current property info to check raw errors and normalization factors
        current_property_info = self._get_best_sequence_property_info()
        
        worst_props = []
        normalization_fixes = []
        
        # Calculate total raw error (weighted) for unsatisfied properties.
        # Using weight here ensures the contribution math respects user
        # priorities -- a property with weight=3000 will dominate the
        # contribution share even if its raw_error is small in absolute terms.
        total_raw_error = self._weighted_raw_sum(current_property_info, unsatisfied_only=True)
        
        property_priorities = []
        
        for prop_name in self.property_scales:
            prop_info = current_property_info.get(prop_name, {})
            raw_error = prop_info.get('raw_error', 0)
            weighted_error = prop_info.get('weighted_error', 0)
            norm_factor = prop_info.get('normalization_factor', 1.0)
            within_tolerance = prop_info.get('within_tolerance', False)
            
            # Skip properties that are already satisfied (within tolerance)
            if within_tolerance:
                if self.debugging:
                    tolerance = prop_info.get('tolerance', 0.0)
                    self.logger.info(f"   ✅ {prop_name} satisfied: error {raw_error:.4f} <= tolerance {tolerance:.4f}")
                continue
            
            # Calculate contribution (weight * raw_error / total_weighted_raw)
            # so user weights drive priority while staying scale-invariant.
            weight = prop_info.get('weight', 1.0)
            contribution = (raw_error * weight) / total_raw_error if total_raw_error > 0 else 0
            
            # Find properties that haven't improved recently (trend-based stagnation)
            improvements = list(self.property_improvements[prop_name])
            trend_stagnant = False
            trend_score = 0.0
            if len(improvements) >= self.min_trend_samples:
                # Closed-form slope (see _update_property_scaling).
                _y = np.asarray(improvements, dtype=float)
                _n = _y.shape[0]
                _x = np.arange(_n, dtype=float)
                _xm = (_n - 1) * 0.5
                _ym = float(_y.mean())
                _dx = _x - _xm
                _denom = float((_dx * _dx).sum())
                trend = float((_dx * (_y - _ym)).sum()) / _denom if _denom > 0 else 0.0
                # Normalize by mean magnitude in the window so the threshold
                # is comparable across properties of different scale (#10).
                window_mean = max(abs(_ym), 1e-8)
                relative_trend = trend / window_mean
                # Use a stricter threshold: only boost if relative trend is
                # positive (getting worse) or very flat.
                if relative_trend >= -0.0001:
                    trend_stagnant = True
                    trend_score = relative_trend  # Higher score = worse trend
            
            # Check for normalization-based stagnation (high raw error + low normalization)
            normalization_stagnant = ((raw_error > 0.3 and norm_factor < 0.2) or 
                                    (raw_error > 0.8))
            
            # Calculate priority score (higher = more need for boost)
            priority_score = 0.0
            if trend_stagnant:
                priority_score += contribution * 2.0 + trend_score * 10.0  # Weight by contribution and trend
            if normalization_stagnant:
                priority_score += raw_error * 5.0  # Weight by raw error magnitude
                
            if trend_stagnant or normalization_stagnant:
                property_priorities.append((prop_name, priority_score, trend_stagnant, normalization_stagnant, raw_error, norm_factor, contribution))
        
        # Sort by priority score (highest first) and only boost high-contribution properties
        property_priorities.sort(key=lambda x: x[1], reverse=True)
        
        # Only boost properties contributing >25% to total error (selective but not too restrictive)
        contribution_threshold = 0.25
        
        boosted_count = 0
        for prop_name, priority_score, trend_stagnant, normalization_stagnant, raw_error, norm_factor, contribution in property_priorities:
            # Only boost if contribution exceeds threshold
            if contribution < contribution_threshold:
                if self.debugging:
                    self.logger.info(f"   ⏭️  Skipping {prop_name}: contribution {contribution*100:.1f}% < {contribution_threshold*100:.0f}% threshold")
                continue
                
            boosted_count += 1
            
            if trend_stagnant:
                worst_props.append(prop_name)
            
            # Check if property needs normalization fix (either explicitly detected or trend stagnant with low norm factor)
            needs_norm_fix = normalization_stagnant or (trend_stagnant and norm_factor < 0.5)
            
            if needs_norm_fix:
                normalization_fixes.append(prop_name)
                if self.debugging:
                    reason = "normalization issue" if normalization_stagnant else "low norm_factor + trend stagnant"
                    self.logger.info(f"   🔧 Boosting {prop_name}: {contribution*100:.1f}% contribution "
                                   f"({reason}, raw_error: {raw_error:.3f}, norm_factor: {norm_factor:.3f})")
        
        # Boost scaling for trend-stagnant properties
        for prop_name in worst_props:
            self.property_scales[prop_name] *= self.stagnation_boost_factor
            self.property_scales[prop_name] = min(self.property_scales[prop_name], self.max_scale)
        
        # Fallback: if no properties meet threshold, boost the top contributor to avoid total stagnation
        if boosted_count == 0 and property_priorities:
            prop_name, priority_score, trend_stagnant, normalization_stagnant, raw_error, norm_factor, contribution = property_priorities[0]
            if self.debugging:
                self.logger.info(f"   🚨 Fallback: No properties >{contribution_threshold*100:.0f}% contribution, boosting top contributor {prop_name} ({contribution*100:.1f}%)")
            
            if trend_stagnant:
                worst_props.append(prop_name)
            
            # Apply same normalization fix logic as main loop
            needs_norm_fix = normalization_stagnant or (trend_stagnant and norm_factor < 0.5)
            if needs_norm_fix:
                normalization_fixes.append(prop_name)
                if self.debugging:
                    reason = "normalization issue" if normalization_stagnant else "low norm_factor + trend stagnant"
                    self.logger.info(f"   🔧 Fallback normalization fix for {prop_name}: {reason} "
                                   f"(raw_error: {raw_error:.3f}, norm_factor: {norm_factor:.3f})")
        
        # Apply normalization fixes for under-normalized properties using the
        # *recovery boost* layer rather than mutating the clipped
        # normalization_factors directly. This lets repeated recoveries keep
        # escalating without hitting the [0.1, 10.0] clip ceiling (#9).
        for prop_name in normalization_fixes:
            # Get current property info to determine appropriate fix
            prop_info = current_property_info.get(prop_name, {})
            raw_error = prop_info.get('raw_error', 0)
            old_boost = self.norm_boost_factors.get(prop_name, 1.0)

            # Be more aggressive with boost based on raw error magnitude
            if raw_error > 0.8:
                new_boost = min(self.max_norm_boost, max(old_boost * 3.0, 5.0))
            elif raw_error > 0.5:
                new_boost = min(self.max_norm_boost, old_boost * 2.5)
            else:
                new_boost = min(self.max_norm_boost, old_boost * 1.5)

            self.norm_boost_factors[prop_name] = new_boost

            if self.debugging:
                base_norm = self.normalization_factors.get(prop_name, 1.0)
                self.logger.info(
                    f"   🔧 Fixed {prop_name} norm boost: {old_boost:.3f} → {new_boost:.3f} "
                    f"(base_norm: {base_norm:.3f}, raw_error: {raw_error:.3f})"
                )

            # Also boost scaling
            self.property_scales[prop_name] *= self.stagnation_boost_factor
            self.property_scales[prop_name] = min(self.property_scales[prop_name], self.max_scale)
        
        all_boosted = worst_props + normalization_fixes
        if self.debugging and all_boosted:
            norm_factors_per_prop=''
            for prop in all_boosted:
                norm_factors_per_prop += f"   - {prop}: norm_factor={self.normalization_factors.get(prop, 1.0):.3f}, scale={self.property_scales.get(prop, 1.0):.3f}\n"
            self.logger.info(f"   Updated normalization factors and scales for boosted properties:\n{norm_factors_per_prop.strip()}")
            self.logger.info(f"   ✅ Boosted {len(all_boosted)} properties: {', '.join(all_boosted)}")
            if normalization_fixes:
                self.logger.info(f"   Applied normalization fixes for: {', '.join(normalization_fixes)}")
        
        # CRITICAL FIX: Recalculate best sequence error with new scaling factors
        # This ensures that candidates in the next iteration are evaluated fairly
        if all_boosted:
            # Freeze adaptive scaling for a window so the recovery boosts have
            # time to take effect before _update_property_scaling smooths them
            # back out (#8).
            freeze_window = self.scale_freeze_window or self.stagnation_threshold
            self._scale_freeze_until = self.iteration + freeze_window
            
            old_error = self.best_error
            # Scaling/normalization changed, but raw property values are
            # independent of those -- the cached raw values are still valid
            # and _assemble_property_info reads the new scales fresh. Reusing
            # the cache here saves N expensive prop.calculate() calls on
            # every stagnation event.
            self.best_error, _ = self._evaluate_sequence(self.best_sequence)
            # Invalidate cached property info since scaling changed
            self._best_sequence_property_info = None
            # Invalidate cached sequence ID
            if hasattr(self, '_best_sequence_property_info_cache_id'):  
                delattr(self, '_best_sequence_property_info_cache_id')

            if self.debugging:
                self.logger.info(f"   🔄 Recalculated best sequence error: {old_error:.6f} → {self.best_error:.6f}")
            
            # Update error history with recalculated value
            if self.best_error_history:
                self.best_error_history[-1] = self.best_error
    
    def _emergency_diversity_injection(self, current_best_raw_error: float) -> float:
        """
        Emergency measure when severely stagnant - try completely different sequences
        Returns the new best raw error (may be updated if better sequence found)
        """
        if self.debugging:
            self.logger.info("🚨 Applying emergency diversity injection...")
        
        # Generate diverse sequences and pick best one
        if self._aa_count_bounds is not None:
            diverse_sequences = []
            seen: set = set()
            attempts = 0
            while len(diverse_sequences) < 10 and attempts < 50:
                seq = self._build_constrained_random_sequence(self.target_length)
                if seq not in seen:
                    seen.add(seq)
                    diverse_sequences.append(seq)
                attempts += 1
        else:
            diverse_sequences = build_diverse_initial_sequences(self.target_length, num_sequences=10)
        
        best_emergency_weighted_error = self.best_error
        best_emergency_raw_error = current_best_raw_error
        best_emergency_seq = self.best_sequence
        
        for seq in diverse_sequences:
            weighted_error, property_info = self._evaluate_sequence(seq)
            # Weighted raw error so user priorities matter (#weight-fix).
            raw_error = self._weighted_raw_sum(property_info, unsatisfied_only=True)
            
            # Acceptance: improves weighted error; raw monotonicity is honored
            # only when the user opted in (consistent with the main loop, #7).
            if weighted_error < best_emergency_weighted_error:
                if self.enforce_raw_monotonicity:
                    raw_limit = current_best_raw_error * (1.0 + self.raw_monotonicity_slack)
                    if raw_error > raw_limit:
                        continue
                best_emergency_weighted_error = weighted_error
                best_emergency_raw_error = raw_error
                best_emergency_seq = seq
        
        # Only switch if we found something better
        if best_emergency_weighted_error < self.best_error:
            old_weighted_error = self.best_error
            old_raw_error = current_best_raw_error
            self._update_best_sequence(best_emergency_seq, best_emergency_weighted_error)
            
            # After jumping to a structurally different sequence, the previous
            # normalization factors may no longer be appropriate. Optionally
            # blend in fresh ones derived from the new seed (#12).
            if self.recompute_norm_on_emergency:
                self._blend_normalization_with_current_sequence()
            
            if self.debugging:
                weighted_improvement = old_weighted_error - best_emergency_weighted_error
                raw_improvement = old_raw_error - best_emergency_raw_error
                self.logger.info(f"   Emergency sequence found! Weighted improvement: {weighted_improvement:.6f}, Raw improvement: {raw_improvement:.6f}")
            
            return best_emergency_raw_error
        else:
            if self.debugging:
                self.logger.info("   No better emergency sequence found")
            return current_best_raw_error
    
    def configure_convergence(
        self,
        tolerance: Optional[float] = None,
        window: Optional[int] = None,
        enable_early_stopping: Optional[bool] = None,
        patience: Optional[int] = None
    ) -> None:
        """
        Configure convergence detection parameters.
        
        Parameters
        ----------
        tolerance : float, optional
            Error standard deviation threshold for convergence detection
        window : int, optional
            Number of recent iterations to analyze for convergence
        enable_early_stopping : bool, optional
            Whether to enable early stopping when convergence is detected
        patience : int, optional
            Number of consecutive convergent iterations required before stopping
        """
        if tolerance is not None:
            self.convergence_tolerance = tolerance
        if window is not None:
            self.convergence_window = window
        if enable_early_stopping is not None:
            self.enable_early_convergence = enable_early_stopping
        if patience is not None:
            self.convergence_patience = patience
        
        # Reset convergence counter when settings change
        self.convergence_counter = 0
        
        if self.verbose:
            self.logger.info(f"Updated convergence settings:")
            self.logger.info(f"  Tolerance: {self.convergence_tolerance}")
            self.logger.info(f"  Window: {self.convergence_window}")
            self.logger.info(f"  Early stopping: {self.enable_early_convergence}")
            self.logger.info(f"  Patience: {self.convergence_patience}")
    
    def configure_error_tolerance(
        self,
        tolerance: Optional[float] = None,
        enable: Optional[bool] = None
    ) -> None:
        """
        Configure error tolerance early stopping parameters.
        
        Parameters
        ----------
        tolerance : float, optional
            Error value below which optimization will stop (None to disable)
        enable : bool, optional
            Whether to enable error tolerance early stopping
        """
        if tolerance is not None:
            self.error_tolerance = tolerance
        if enable is not None:
            self.enable_error_tolerance = enable
        
        if self.verbose:
            if self.error_tolerance is not None and self.enable_error_tolerance:
                self.logger.info(f"Error tolerance stopping enabled: {self.error_tolerance:.6f}")
            else:
                self.logger.info("Error tolerance stopping disabled")
    
    def get_convergence_info(self) -> Dict[str, Any]:
        """
        Get current convergence detection information.
        
        This method is pure: it does not mutate convergence_counter.
        
        Returns
        -------
        Dict[str, Any]
            Dictionary with convergence status and parameters
        """
        is_under_tolerance, current_std = self._compute_convergence_state()
        if current_std > 0 and current_std != float('inf'):
            convergence_ratio = self.convergence_tolerance / current_std
        elif current_std == 0:
            convergence_ratio = float('inf')
        else:
            convergence_ratio = 0.0
        # Predicted convergence (without mutating counter): would the next
        # _check_convergence call return True given current state?
        if self.enable_early_convergence and is_under_tolerance:
            predicted_counter = self.convergence_counter + 1
            is_converged = predicted_counter >= self.convergence_patience
        else:
            is_converged = False
        
        return {
            'is_converged': is_converged,
            'convergence_counter': self.convergence_counter,
            'current_error_std': current_std,
            'tolerance': self.convergence_tolerance,
            'convergence_ratio': convergence_ratio,
            'window_size': self.convergence_window,
            'patience_remaining': max(0, self.convergence_patience - self.convergence_counter),
            'early_stopping_enabled': self.enable_early_convergence
        }
    
    def run(self) -> str:
        """
        Run the optimization.
        
        Returns
        -------
        best_sequence : str
            Optimized sequence
        """
        if not self.properties:
            raise ValueError("No properties defined for optimization")
        
        # Initialize with best of multiple diverse sequences
        if self.starting_sequence is None:
            if self._aa_count_bounds is not None:
                # Honor per-AA fraction ranges by sampling constrained
                # compositions instead of using the unconstrained diverse
                # builder.
                seen: set = set()
                initial_candidates = []
                attempts = 0
                cap = max(self.num_starting_candidates * 5,
                          self.num_starting_candidates + 10)
                while (
                    len(initial_candidates) < self.num_starting_candidates
                    and attempts < cap
                ):
                    seq = self._build_constrained_random_sequence(self.target_length)
                    if seq not in seen:
                        seen.add(seq)
                        initial_candidates.append(seq)
                    attempts += 1
                if not initial_candidates:
                    raise RuntimeError(
                        "Failed to generate any constrained initial candidates; "
                        "check aa_fraction_ranges feasibility."
                    )
            else:
                initial_candidates = build_diverse_initial_sequences(self.target_length, num_sequences=self.num_starting_candidates)
        else:
            initial_candidates = [self.starting_sequence]
        
        # set logging info
        if self.verbose:
            if self.starting_sequence is None:
                self.logger.info(f"Generating and testing {len(initial_candidates)} diverse initial candidate sequences...")
            else:
                self.logger.info(f"Using provided starting sequence for optimization...")
        
        # Batch evaluate initial candidates
        initial_results = self._evaluate_sequence(initial_candidates)

        # Select best initial candidate using a *weight-aware*, scale-invariant
        # metric: sum(raw_error * weight) over all properties. raw_error is
        # independent of the adaptive scales/normalization (so it is stable
        # before normalization is computed) and weighting by prop.weight
        # ensures user priorities drive the seed choice instead of being
        # silently averaged away. (See weight-effectiveness fix.)
        prop_weights = [p.weight for p in self.properties]
        names = self._property_names
        best_initial_idx = 0
        best_initial_raw_sum = float('inf')
        for i, (_, prop_info) in enumerate(initial_results):
            raw_sum = sum(
                prop_info[n]['raw_error'] * w
                for n, w in zip(names, prop_weights)
            )
            if raw_sum < best_initial_raw_sum:
                best_initial_raw_sum = raw_sum
                best_initial_idx = i
        best_initial_seq = initial_candidates[best_initial_idx]
        best_initial_error = initial_results[best_initial_idx][0]
        
        # Seed the elite pool with the top-K of the initial batch (free, no
        # extra evaluations). Pool ranking uses the weight-aware raw sum so
        # user priorities steer pool ordering as well.
        if self.elite_pool_size > 1:
            self._elite_pool = []
            scored = []
            for i, (w_err, prop_info) in enumerate(initial_results):
                raw_sum = sum(
                    prop_info[n]['raw_error'] * w
                    for n, w in zip(names, prop_weights)
                )
                scored.append((raw_sum, w_err, initial_candidates[i]))
            scored.sort(key=lambda x: x[0])
            seen = set()
            for entry in scored:
                if entry[2] in seen:
                    continue
                seen.add(entry[2])
                self._elite_pool.append(entry)
                if len(self._elite_pool) >= self.elite_pool_size:
                    break
        
        # Initialize optimization state
        self.best_sequence = best_initial_seq
        self.best_error = best_initial_error
        self.best_error_history.append(self.best_error)        

        # Calculate dynamic normalization factors based on initial sequence
        if self.verbose:
            self.logger.info(f"Starting optimization (target length: {self.target_length})")
            self.logger.info(f"Starting optimization with initial sequence: {self.best_sequence}")

        # debugging info
        if self.debugging:
            self.logger.info(f"   Initial error: {self.best_error:.4f}")        
            self.logger.info(f"📏 Calculating dynamic normalization factors:")
        
        # calculate initial normalization factors (skip if already computed,
        # e.g. via a prior call to set_initial_sequence)
        if not self._normalization_initialized:
            self._calculate_initial_normalization(self.best_sequence)
        elif self.debugging:
            self.logger.info("   Using pre-computed normalization factors.")
        
        # Recalculate initial error with normalization applied
        self.best_error, initial_property_info = self._evaluate_sequence(self.best_sequence)
        self.best_error_history[-1] = self.best_error  # Update the last entry
        
        if self.debugging:
            self.logger.info(f"   Normalized initial error: {self.best_error:.4f}")
        
        # Progress tracking
        progress_bar = tqdm(
            total=self.max_iterations,
            disable=not self.verbose,
            mininterval=0,
            miniters=self.update_interval,
            leave=True,
        )
        stagnation_counter = 0
        
        # Track best raw error (weight-aware, unsatisfied) so progress
        # tracking and pool ranking honor user priorities. Reuses
        # initial_property_info from the evaluation above (raw values are
        # invariant to scale/normalization, so re-assembly would yield the
        # same raw_error / within_tolerance fields).
        best_raw_error = self._weighted_raw_sum(initial_property_info, unsatisfied_only=True)
        self.best_raw_error_history.append(best_raw_error)
        
        # Main optimization loop
        for _ in range(self.max_iterations):
            self.iteration += 1 # increment iteration count.
            # Pick a parent (best_sequence by default; pool member when
            # elite_pool_size > 1) and generate candidates from it. Always
            # include best_sequence so we never lose the global optimum; the
            # cache makes this essentially free.
            parent = self._select_parent_from_pool()
            candidates = self._generate_candidates(parent, self.iteration)
            if self.elite_pool_size > 1 and self.best_sequence not in candidates:
                candidates.append(self.best_sequence)
            
            # Track improvement in this iteration
            improved_this_iteration = False
            
            # batch calculate
            candidate_results = self._evaluate_sequence(candidates)

            # find best candidate for this iteration
            best_candidate = None
            best_candidate_error = self.best_error
            best_candidate_raw_error = float('inf')
            best_candidate_property_info = None

            # Process batch results to find best candidate
            for i, (candidate, (weighted_error, property_info)) in enumerate(zip(candidates, candidate_results)):
                # Weight-aware raw error over unsatisfied properties. Used
                # for pool ranking and progress tracking; respects priorities
                # while staying invariant to adaptive scaling.
                candidate_raw_error = self._weighted_raw_sum(property_info, unsatisfied_only=True)
                
                # Maintain elite pool from every evaluated candidate (#13).
                # No extra cost: we already have raw_error / weighted_error.
                if self.elite_pool_size > 1:
                    self._maybe_add_to_pool(candidate, candidate_raw_error, weighted_error)
                
                # Skip the first candidate (the parent itself) when it equals
                # the current best -- accepting it would be a no-op.
                if i == 0 and candidate == self.best_sequence:
                    continue
                
                # Acceptance rule. By default we accept any candidate that
                # improves weighted error (true multi-objective behavior). When
                # ``enforce_raw_monotonicity`` is on, additionally require that
                # the unsatisfied raw-error sum does not worsen beyond the
                # configured slack. (#7)
                if weighted_error >= best_candidate_error:
                    continue
                if self.enforce_raw_monotonicity:
                    raw_limit = best_raw_error * (1.0 + self.raw_monotonicity_slack)
                    if candidate_raw_error > raw_limit:
                        continue
                best_candidate = candidate
                best_candidate_error = weighted_error
                best_candidate_raw_error = candidate_raw_error
                best_candidate_property_info = property_info

            # If we found a better candidate, update best sequence
            if best_candidate is not None:
                self._update_best_sequence(best_candidate, best_candidate_error)
                # Reset stagnation counter only when raw error actually improves
                if best_candidate_raw_error < best_raw_error:
                    improved_this_iteration = True
                    stagnation_counter = 0
                    if self.verbose and self.debugging:
                        improvement = best_raw_error - best_candidate_raw_error
                        self.logger.info(f"   \u2705 Raw error improved by {improvement:.6f}")
                
                best_raw_error = best_candidate_raw_error  # Update best raw error
                
                # Update property scaling
                if self.enable_adaptive_scaling:
                    self._update_property_scaling(best_candidate_property_info)
            
            # Track error history
            self.best_error_history.append(self.best_error)
            self.best_raw_error_history.append(best_raw_error)
            
            # Decay any recovery-only normalization boosts back toward 1.0 so
            # they fade as soon as progress resumes (#9). Skip during the
            # adaptive-scaling freeze window so the boost has time to act.
            if (
                self.norm_boost_factors
                and self.norm_boost_decay < 1.0
                and self.iteration >= self._scale_freeze_until
            ):
                changed = False
                for _name in list(self.norm_boost_factors.keys()):
                    new_boost = self.norm_boost_factors[_name] * self.norm_boost_decay
                    if new_boost <= 1.0 + 1e-6:
                        del self.norm_boost_factors[_name]
                    else:
                        self.norm_boost_factors[_name] = new_boost
                    changed = True
                if changed:
                    self._best_sequence_property_info = None
                    if hasattr(self, '_best_sequence_property_info_cache_id'):
                        delattr(self, '_best_sequence_property_info_cache_id')
            
            # Update progress bar every update_interval iterations
            if self.verbose and (
                self.iteration % self.update_interval == 0
                or self.iteration == self.max_iterations
            ):
                # Count satisfied properties for progress display (use cached version)
                current_property_info = self._get_best_sequence_property_info()
                satisfied_count = sum(1 for info in current_property_info.values() if info.get('within_tolerance', False))
                total_properties = len(self.properties)

                postfix = {
                    'raw_error': f'{best_raw_error:.4f}',
                    'satisfied': f'{satisfied_count}/{total_properties}',
                    'stagnation': stagnation_counter
                }

                # Add error tolerance info if enabled
                if self.enable_error_tolerance and self.error_tolerance is not None:
                    postfix['tolerance'] = f'{self.error_tolerance:.4f}'

                # Add convergence info if enabled and enough history
                if self.enable_early_convergence and len(self.best_error_history) >= self.convergence_window:
                    conv_info = self.get_convergence_info()
                    postfix['conv'] = f"{self.convergence_counter}/{self.convergence_patience}"

                progress_bar.set_postfix(postfix, refresh=False)
                # Snap the bar to the current iteration in a single refresh.
                progress_bar.n = self.iteration
                progress_bar.refresh()
            
            # Check for convergence
            if self._check_convergence():
                if self.verbose:
                    conv_info = self.get_convergence_info()
                    self.logger.info(f"\n✅ Converged at iteration {self.iteration}")
                if self.debugging:
                    self.logger.info(f"   Error std: {conv_info['current_error_std']:.6f} < {conv_info['tolerance']:.6f}")
                    self.logger.info(f"   Patience satisfied: {self.convergence_counter}/{self.convergence_patience}")
                break
            
            # Check for error tolerance
            if self._check_error_tolerance():
                if self.verbose:
                    self.logger.info(f"\n🎯 Error tolerance reached at iteration {self.iteration}")
                    self.logger.info(f"   Current error for properties not within tolerance: {self.best_error:.6f} <= {self.error_tolerance:.6f}")
                break
            
            # Handle stagnation
            if not improved_this_iteration:
                stagnation_counter += 1
                
                if stagnation_counter >= self.stagnation_threshold:
                    self._apply_stagnation_recovery()
                    
                    # Recalculate best raw error after stagnation recovery (scaling may have changed)
                    # Use cached property info (already updated by stagnation recovery)
                    updated_property_info = self._get_best_sequence_property_info()
                    best_raw_error = self._weighted_raw_sum(updated_property_info, unsatisfied_only=True)
                    
                    # If severely stagnant, try emergency diversity injection
                    if stagnation_counter >= self.stagnation_threshold * 2:
                        best_raw_error = self._emergency_diversity_injection(best_raw_error)
                    
                    stagnation_counter = 0  # Reset counter after recovery

        # Ensure the bar reflects the final iteration count before closing
        # (covers convergence/tolerance early-exit paths).
        if self.verbose and progress_bar.n < self.iteration:
            progress_bar.n = self.iteration
            progress_bar.refresh()
        progress_bar.close()
        
        # Final evaluation and reporting
        final_property_info = self._get_best_sequence_property_info()
        final_error = self.best_error
        total_raw_error = sum(info['raw_error'] for info in final_property_info.values())
        unsatisfied_raw_error = sum(info['raw_error'] for info in final_property_info.values() 
                                  if not info.get('within_tolerance', False))
        satisfied_count = sum(1 for info in final_property_info.values() if info.get('within_tolerance', False))
        
        if self.debugging:
            cache_stats = self.get_cache_statistics()
            self.logger.info(f"\n\n⚡ Performance:")
            self.logger.info(f"   Cache hit rate: {cache_stats['hit_rate']:.1%}")
            self.logger.info(f"   Evaluations saved: {cache_stats['total_evaluations_saved']}")
            self.logger.info(f"   Total evaluations: {cache_stats['cache_hits'] + cache_stats['cache_misses']}")


        if self.verbose:
            # Show cache performance
            self.logger.info(f"\n\nOptimization complete!")
            self.logger.info(f"   Final raw error for properties not within tolerance: {unsatisfied_raw_error:.6f}")
            self.logger.info(f"   Properties satisfied: {satisfied_count}/{len(self.properties)}")
            self.logger.info(f"   Total iterations: {self.iteration + 1}")
        
        if self.debugging:
            self.logger.info(f"   Final weighted error: {final_error:.6f}")
            self.logger.info(f"   Final raw error (all): {total_raw_error:.6f}")
            
        if self.verbose:
            self.logger.info(f"\n\nFinal property values:")
            
            for prop_name, info in final_property_info.items():
                target = info['target_value']
                actual = info['raw_value']
                raw_error = info['raw_error']
                normalized_error = info['normalized_error']
                norm_factor = info['normalization_factor']
                tolerance = info.get('tolerance', 0.0)
                within_tolerance = info.get('within_tolerance', False)
                initial_error = self.initial_property_errors.get(prop_name, 'N/A')
                
                status_indicator = "✅" if within_tolerance else "❌"
                tolerance_info = f", tolerance: {tolerance:.4f}" if tolerance > 0 else ""
                constraint_type = info.get('constraint_type', 'N/A')
                
                if self.debugging:
                    self.logger.info(
                        f"   {status_indicator} {prop_name}: {actual:.3f} "
                        f"(target: {target:.3f}, raw_error: {raw_error:.4f}), "
                        f"initial_error: {initial_error:.4f}, norm_factor: {norm_factor:.3f}{tolerance_info})"
                    )
                else:
                    self.logger.info(
                        f"   {status_indicator} {prop_name}: {actual:.3f} "
                        f"(target: {constraint_type.__str__()} value of {target:.3f}, raw_error: {raw_error:.4f}), tolerance: {tolerance:.4f}"
                    )
            self.logger.info(f"\nOptimized sequence:\n{self.best_sequence}")
        return self.best_sequence
    
    def _clear_evaluation_cache(self) -> None:
        """Clear the evaluation cache (call when needed)"""
        if hasattr(self, '_raw_property_cache'):
            self._raw_property_cache.clear()
        self._best_sequence_property_info = None
        self._normalization_initialized = False
    
    def get_cache_statistics(self) -> Dict[str, Any]:
        """Get evaluation cache performance statistics"""
        total_requests = self._cache_hits + self._cache_misses
        hit_rate = self._cache_hits / total_requests if total_requests > 0 else 0.0
        
        return {
            'cache_hits': self._cache_hits,
            'cache_misses': self._cache_misses,
            'hit_rate': hit_rate,
            'cache_size': len(getattr(self, '_raw_property_cache', {})),
            'total_evaluations_saved': self._cache_hits
        }
    
    
