# before adding in the normalization recalcualton
"""
SequenceOptimizer rewrite. Trying to make it work better
for properties that have variable ranges of values.
"""

import logging
import math
import random
import time
from collections import defaultdict, deque
from typing import Any, Dict, List, Optional, Tuple, Union

import numpy as np
import sparrow
from tqdm import tqdm

from goose.backend import optimizer_properties
from goose.backend.fast_mutations import (
    cython_make_point_mutation as make_point_mutation,
)
from goose.backend.fast_mutations import cython_shuffle_sequence as shuffle_sequence
from goose.backend.optimizer_properties import ProteinProperty
from goose.data import amino_acids

# ==============================================================================
# UTILITY FUNCTIONS
# ==============================================================================


def build_random_sequence(target_length: int) -> str:
    """Build a random amino acid sequence"""
    return "".join(
        random.choices(list(amino_acids.amino_acid_dat.keys()), k=target_length)
    )


def build_diverse_initial_sequences(
    target_length: int, num_sequences: int = 5
) -> List[str]:
    """Build multiple diverse initial sequences to find the best starting point"""
    sequences = []

    # Pure random
    sequences.append(build_random_sequence(target_length))

    # set base aas. Disorder promoting, not charged, not aromatic.
    base_aas = ["G", "P", "S", "T", "Q", "N"]

    # other lists with varying properties
    biased_lists = [
        ["A", "V", "I", "L", "M"],
        ["F", "W", "Y"],
        ["D", "E"],
        ["K", "R"],
        ["D", "E", "K", "R"],
        ["C"],
        ["H"],
        ["T", "S", "N", "Q", "P", "G"],
        ["G", "P"],
        ["A", "V", "I", "L", "M", "E", "D"],
        ["A", "V", "I", "L", "M", "K", "R"],
        ["F", "W", "Y", "E", "D"],
        ["F", "W", "Y", "K", "R"],
    ]
    # High disorder bias (GP-rich)
    if num_sequences > 1:
        for i in range(num_sequences - 1):
            disorder_aas = base_aas + random.choice(biased_lists)
            seq = "".join(random.choices(disorder_aas, k=target_length))
            sequences.append(seq)
    # Ensure uniqueness
    sequences = list(set(sequences))
    return sequences


# def shuffle_sequence(
#     sequence: str, window_size: int = 10, global_shuffle_probability: float = 0.3
# ) -> str:
#     """Shuffle sequence with window-based approach"""
#     seq_list = list(sequence)

#     # Choose between global and local shuffling
#     if random.random() < global_shuffle_probability:
#         random.shuffle(seq_list)
#     else:  # Local window shuffle
#         if len(sequence) > window_size:
#             start = random.randint(0, len(sequence) - window_size)
#             window = seq_list[start : start + window_size]
#             random.shuffle(window)
#             seq_list[start : start + window_size] = window

#     return "".join(seq_list)


# def make_point_mutation(sequence: str, position: Optional[int] = None) -> str:
#     """Make a single point mutation"""
#     if not sequence:
#         return sequence

#     pos = position if position is not None else random.randint(0, len(sequence) - 1)
#     old_aa = sequence[pos]

#     # Choose new amino acid (different from current)
#     amino_acids_list = list(amino_acids.amino_acid_dat.keys())
#     new_aa = random.choice([aa for aa in amino_acids_list if aa != old_aa])

#     return sequence[:pos] + new_aa + sequence[pos + 1 :]


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
        verbose: bool = True,  # progress bar and logging
        enable_adaptive_scaling: bool = True,
        enable_shuffling: bool = True,
        shuffle_frequency: int = 50,  # More frequent shuffling for diversity
        stagnation_threshold: int = 25,  # Detect and respond to stagnation sooner
        convergence_tolerance: float = 1e-4,  # Less strict convergence criterion
        convergence_window: int = 20,
        enable_early_convergence: bool = False,
        convergence_patience: int = 20,  # Reduced patience for faster adaptation
        # Error tolerance parameters
        error_tolerance: Optional[
            float
        ] = 1e-6,  # Stop when total error below this value
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
        # Analysis parameters
        min_analysis_samples: int = 5,
        improvement_threshold: float = -0.001,
        stability_threshold: float = 0.01,
        debugging: bool = False,  # Enable extra logging for debugging
        update_interval: int = 1,  # Update progress bar every N iterations
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
        min_analysis_samples : int
            Minimum samples required for property analysis
        improvement_threshold : float
            Threshold for determining if a property is improving
        stability_threshold : float
            Variance threshold for determining property stability
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

        # Analysis parameters
        self.min_analysis_samples = min_analysis_samples
        self.improvement_threshold = improvement_threshold
        self.stability_threshold = stability_threshold

        # Properties and state
        self.properties: List[ProteinProperty] = []
        self.property_scales: Dict[str, float] = {}
        self.property_improvements: Dict[str, deque] = {}
        self.best_error_history: deque = deque(maxlen=self.error_history_size)

        # Dynamic normalization based on initial sequence
        self.initial_property_errors: Dict[str, float] = {}
        self.normalization_factors: Dict[str, float] = {}

        # Optimization state
        self.best_sequence: str = ""
        self.best_error: float = float("inf")
        self.iteration: int = 0
        self.convergence_counter: int = 0
        self.starting_sequence: str = None

        # Evaluation cache for performance optimization
        self._raw_property_cache: Dict[str, Dict[str, Dict[str, Any]]] = {}
        self._cache_hits: int = 0
        self._cache_misses: int = 0

        # Current best sequence property info cache (updated when best_sequence changes)
        self._best_sequence_property_info: Optional[Dict[str, Dict[str, float]]] = None

        # debugging flag
        self.debugging = debugging

        self.update_interval = update_interval

        # Configure logging
        self._configure_logging()

    def set_initial_sequence(self, sequence: str) -> None:
        """Set the initial sequence and calculate normalization factors"""
        if len(sequence) != self.target_length:
            raise ValueError(
                f"Initial sequence length {len(sequence)} does not match target length {self.target_length}"
            )

        self.starting_sequence = sequence

        # Calculate initial normalization factors
        self._calculate_initial_normalization(sequence)

        # Evaluate initial sequence
        initial_error, _ = self._evaluate_sequence(sequence)
        self.best_error = initial_error
        self.best_error_history.append(initial_error)

        if self.debugging:
            self.logger.info(
                f"Initial sequence set. Initial weighted error: {initial_error:.6f}"
            )

    def _configure_logging(self) -> None:
        """Configure logging"""
        if self.verbose:
            logging.basicConfig(level=logging.INFO, format="%(message)s")
        self.logger = logging.getLogger(__name__)

    def add_property(
        self, property_class: type, *args: Any, tolerance: float = 0.0, **kwargs: Any
    ) -> None:
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
            Default: 0.0 (no tolerance - property must be exact)
        **kwargs : Any
            Keyword arguments for the property constructor
        """
        prop = property_class(*args, **kwargs)

        # Set the per-property tolerance
        prop._tolerance = tolerance

        # Set sequence length hint for length-dependent calculations
        prop._set_sequence_length_hint(self.target_length)

        self.properties.append(prop)

        prop_name = prop.__class__.__name__
        self.property_scales[prop_name] = 1.0
        self.property_improvements[prop_name] = deque(
            maxlen=self.improvement_history_size
        )

        # Clear cache when new properties are added since cached values may not include the new property
        self._clear_evaluation_cache()

        if self.verbose:
            tolerance_info = f", tolerance: ±{tolerance}" if tolerance > 0 else ""
            self.logger.info(
                f"Added property: {prop_name} (target: {prop.target_value}, weight: {prop.weight}{tolerance_info})"
            )

    def _evaluate_sequences_batch(
        self, sequences: List[str], use_cache: bool = True
    ) -> List[Tuple[float, Dict[str, Dict[str, float]]]]:
        """
        Evaluate multiple sequences in batch for improved performance.

        Parameters
        ----------
        sequences : List[str]
            List of sequences to evaluate
        use_cache : bool, optional
            Whether to use evaluation cache. Default: True

        Returns
        -------
        List[Tuple[float, Dict[str, Dict[str, float]]]]
            List of (total_error, property_info) tuples for each sequence
        """
        if not sequences:
            return []

        # Initialize cache if needed
        if not hasattr(self, "_raw_property_cache"):
            self._raw_property_cache = {}

        # Separate cached and uncached sequences
        cached_results = {}
        uncached_sequences = []
        uncached_indices = []

        for i, seq in enumerate(sequences):
            if use_cache and seq in self._raw_property_cache:
                cached_results[i] = self._raw_property_cache[seq]
                self._cache_hits += 1
            else:
                uncached_sequences.append(seq)
                uncached_indices.append(i)
                self._cache_misses += 1

        # Batch calculate raw property values for uncached sequences
        uncached_raw_values = {}
        if uncached_sequences:
            # Create SPARROW protein objects in batch
            proteins = [sparrow.Protein(seq) for seq in uncached_sequences]

            # Calculate properties for all proteins in batch
            for prop in self.properties:
                prop_name = prop.__class__.__name__

                # Check if property supports batch calculation
                if prop.calculate_in_batch:
                    # Use batch calculation if available
                    batch_results = prop.calculate_batch(proteins)
                    for i, (raw_value, raw_error) in enumerate(batch_results):
                        seq_idx = uncached_indices[i]
                        if seq_idx not in uncached_raw_values:
                            uncached_raw_values[seq_idx] = {}
                        uncached_raw_values[seq_idx][prop_name] = {
                            "raw_value": raw_value,
                            "raw_error": raw_error,
                            "target_value": prop.target_value,
                            "tolerance": getattr(prop, "_tolerance", 0.0),
                        }
                else:
                    # Fall back to individual calculation
                    for i, protein in enumerate(proteins):
                        seq_idx = uncached_indices[i]
                        raw_value, raw_error = prop.calculate(protein)
                        if seq_idx not in uncached_raw_values:
                            uncached_raw_values[seq_idx] = {}
                        uncached_raw_values[seq_idx][prop_name] = {
                            "raw_value": raw_value,
                            "raw_error": raw_error,
                            "target_value": prop.target_value,
                            "tolerance": getattr(prop, "_tolerance", 0.0),
                        }

            # Cache the raw values for uncached sequences
            if use_cache:
                # Limit cache size
                if len(self._raw_property_cache) > 500:
                    # Remove oldest 20% of entries
                    items_to_remove = list(self._raw_property_cache.keys())[:100]
                    for key in items_to_remove:
                        del self._raw_property_cache[key]

                for i, seq in enumerate(uncached_sequences):
                    seq_idx = uncached_indices[i]
                    self._raw_property_cache[seq] = uncached_raw_values[seq_idx]

        # Combine cached and uncached results, then apply scaling/normalization
        results = []
        for i, seq in enumerate(sequences):
            if i in cached_results:
                cached_raw_values = cached_results[i]
            else:
                cached_raw_values = uncached_raw_values[i]

            # Apply current scaling and normalization (same logic as original _evaluate_sequence)
            property_info = {}
            weighted_errors = []

            for prop in self.properties:
                prop_name = prop.__class__.__name__
                constraint_type = getattr(prop, "constraint_type", None)

                # Get cached raw values with safety check
                if prop_name not in cached_raw_values:
                    raise ValueError(
                        f"Property {prop_name} not found in cached values. Available: {list(cached_raw_values.keys())}"
                    )

                raw_data = cached_raw_values[prop_name]
                raw_value = raw_data["raw_value"]
                raw_error = raw_data["raw_error"]
                tolerance = raw_data["tolerance"]

                # Check if property is within tolerance (satisfied)
                if tolerance == 0.0:
                    within_tolerance = np.isclose(raw_error, 0.0, atol=1e-8)
                else:
                    within_tolerance = raw_error <= tolerance

                # Apply dynamic normalization based on initial errors
                normalization_factor = self.normalization_factors.get(prop_name, 1.0)
                normalized_error = raw_error * normalization_factor

                # Apply adaptive scaling
                scale = self.property_scales.get(prop_name, 1.0)
                weighted_error = normalized_error * prop.weight * scale

                # If within tolerance, don't contribute to total error
                effective_weighted_error = 0.0 if within_tolerance else weighted_error

                # Store information
                property_info[prop_name] = {
                    "raw_value": raw_value,
                    "target_value": raw_data["target_value"],
                    "raw_error": raw_error,
                    "normalized_error": normalized_error,
                    "normalization_factor": normalization_factor,
                    "weighted_error": weighted_error,
                    "effective_weighted_error": effective_weighted_error,
                    "scale": scale,
                    "tolerance": tolerance,
                    "within_tolerance": within_tolerance,
                    "satisfied": within_tolerance,
                    "constraint_type": constraint_type,
                }

                weighted_errors.append(effective_weighted_error)

            total_error = sum(weighted_errors)
            results.append((total_error, property_info))

        return results

    def _evaluate_sequence(
        self, sequence: Union[str, List[str]], use_cache: bool = True
    ) -> Union[
        Tuple[float, Dict[str, Dict[str, float]]],
        List[Tuple[float, Dict[str, Dict[str, float]]]],
    ]:
        """
        Evaluate sequence(s) and return error and property info.

        Parameters
        ----------
        sequence : Union[str, List[str]]
            Single sequence or list of sequences to evaluate
        use_cache : bool, optional
            Whether to use evaluation cache. Default: True

        Returns
        -------
        Union[Tuple[float, Dict[str, Dict[str, float]]], List[Tuple[float, Dict[str, Dict[str, float]]]]]
            For single sequence: (total_error, property_info)
            For multiple sequences: List of (total_error, property_info) tuples
        """
        # Handle batch evaluation
        if isinstance(sequence, list):
            return self._evaluate_sequences_batch(sequence, use_cache=use_cache)

        # Handle single sequence evaluation (original logic)
        # For cache efficiency, separate raw property calculations from scaling/weighting
        # Raw property values only depend on sequence, not on scaling factors
        raw_cache_key = sequence

        # Initialize cache if needed
        if not hasattr(self, "_raw_property_cache"):
            self._raw_property_cache = {}

        # Check if we have raw property values cached
        if use_cache and raw_cache_key in self._raw_property_cache:
            cached_raw_values = self._raw_property_cache[raw_cache_key]
            self._cache_hits += 1
        else:
            # Calculate raw property values (expensive part)
            protein = sparrow.Protein(sequence)
            cached_raw_values = {}

            for prop in self.properties:
                prop_name = prop.__class__.__name__
                raw_value, raw_error = prop.calculate(protein)
                cached_raw_values[prop_name] = {
                    "raw_value": raw_value,
                    "raw_error": raw_error,
                    "target_value": prop.target_value,
                    "tolerance": getattr(prop, "_tolerance", 0.0),
                }

            # Cache the raw values
            if use_cache:
                # Limit raw cache size
                if len(self._raw_property_cache) > 500:
                    # Remove oldest 20% of entries
                    items_to_remove = list(self._raw_property_cache.keys())[:100]
                    for key in items_to_remove:
                        del self._raw_property_cache[key]

                self._raw_property_cache[raw_cache_key] = cached_raw_values

            self._cache_misses += 1

        # Now apply current scaling and normalization (fast)
        property_info = {}
        weighted_errors = []

        for prop in self.properties:
            prop_name = prop.__class__.__name__
            constraint_type = getattr(prop, "constraint_type", None)

            # Get cached raw values with safety check
            if prop_name not in cached_raw_values:
                raise ValueError(
                    f"Property {prop_name} not found in cached values. Available: {list(cached_raw_values.keys())}"
                )

            raw_data = cached_raw_values[prop_name]
            raw_value = raw_data["raw_value"]
            raw_error = raw_data["raw_error"]
            tolerance = raw_data["tolerance"]

            # Check if property is within tolerance (satisfied)
            # change to use np.isclose for better numerical stability
            if tolerance == 0.0:
                within_tolerance = np.isclose(raw_error, 0.0, atol=1e-8)
            else:
                within_tolerance = raw_error <= tolerance

            # Apply dynamic normalization based on initial errors
            normalization_factor = self.normalization_factors.get(prop_name, 1.0)
            normalized_error = raw_error * normalization_factor

            # Apply adaptive scaling
            scale = self.property_scales.get(prop_name, 1.0)
            weighted_error = normalized_error * prop.weight * scale

            # If within tolerance, don't contribute to total error
            effective_weighted_error = 0.0 if within_tolerance else weighted_error

            # Store information
            property_info[prop_name] = {
                "raw_value": raw_value,
                "target_value": raw_data["target_value"],
                "raw_error": raw_error,
                "normalized_error": normalized_error,
                "normalization_factor": normalization_factor,
                "weighted_error": weighted_error,
                "effective_weighted_error": effective_weighted_error,
                "scale": scale,
                "tolerance": tolerance,
                "within_tolerance": within_tolerance,
                "satisfied": within_tolerance,
                "constraint_type": constraint_type,
            }

            weighted_errors.append(effective_weighted_error)

        total_error = sum(weighted_errors)
        return total_error, property_info

    def _get_best_sequence_property_info(self) -> Dict[str, Dict[str, float]]:
        """
        Get property info for the current best sequence (cached).

        Returns
        -------
        property_info : Dict[str, Dict[str, float]]
            Cached property information for best sequence
        """
        if self._best_sequence_property_info is None:
            _, self._best_sequence_property_info = self._evaluate_sequence(
                self.best_sequence
            )
        return self._best_sequence_property_info

    def _update_best_sequence(self, new_sequence: str, new_error: float) -> None:
        """
        Update the best sequence and invalidate cached property info.

        Parameters
        ----------
        new_sequence : str
            New best sequence
        new_error : float
            Error for the new best sequence
        """
        self.best_sequence = new_sequence
        self.best_error = new_error
        # Invalidate cached property info since sequence changed
        self._best_sequence_property_info = None

    def _calculate_initial_normalization(self, initial_sequence: str) -> None:
        """
        Calculate dynamic normalization factors based on initial sequence errors.
        This balances properties based on their initial distance from targets.
        """
        protein = sparrow.Protein(initial_sequence)
        initial_errors = []

        # Calculate initial errors for all properties
        for prop in self.properties:
            prop_name = prop.__class__.__name__
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
                prop_name = prop.__class__.__name__
                initial_error = self.initial_property_errors[prop_name]

                # Normalization factor scales errors to the reference level
                # Properties with larger initial errors get smaller factors (scaled down)
                # Properties with smaller initial errors get larger factors (scaled up)
                normalization_factor = reference_error / initial_error

                # Apply reasonable bounds to prevent extreme scaling
                normalization_factor = np.clip(normalization_factor, 0.1, 10.0)

                self.normalization_factors[prop_name] = normalization_factor

                if self.debugging:
                    self.logger.info(
                        f"  {prop_name}: initial_error={initial_error:.4f}, "
                        f"norm_factor={normalization_factor:.3f}"
                    )
        else:
            # Single property - no normalization needed
            for prop in self.properties:
                prop_name = prop.__class__.__name__
                self.normalization_factors[prop_name] = 1.0

    def _update_property_scaling(
        self, property_info: Dict[str, Dict[str, float]]
    ) -> None:
        """Update adaptive scaling based on property performance"""
        if not self.enable_adaptive_scaling or len(self.properties) <= 1:
            return

        # Calculate relative error contributions
        total_error = sum(info["weighted_error"] for info in property_info.values())
        if total_error <= 0:
            return

        for prop_name, info in property_info.items():
            # Track improvement history using normalized errors for fair comparison
            self.property_improvements[prop_name].append(info["normalized_error"])

            # Calculate stagnation score
            improvements = list(self.property_improvements[prop_name])
            if len(improvements) >= self.min_trend_samples:
                recent_trend = np.polyfit(range(len(improvements)), improvements, 1)[0]
                stagnation_score = max(0, -recent_trend)  # Negative slope = improvement
            else:
                stagnation_score = 1.0

            # Calculate relative contribution
            contribution = info["weighted_error"] / total_error

            # Simple adaptive scaling formula
            current_scale = self.property_scales[prop_name]

            # Boost properties that are:
            # 1. Far from target (high normalized error)
            # 2. Stagnant (not improving)
            # 3. Contributing little to total error despite high normalized error

            # Use normalized error for fair comparison across property types
            distance_factor = min(
                self.max_distance_factor,
                info["normalized_error"] + self.distance_offset,
            )
            stagnation_factor = 1.0 + stagnation_score * self.stagnation_multiplier

            # If property has low contribution but high normalized error, boost it
            if (
                contribution < self.low_contribution_threshold
                and info["normalized_error"] > self.high_error_threshold
            ):
                boost_factor = self.boost_factor
            else:
                boost_factor = 1.0

            # Calculate new scale with momentum (smooth changes)
            target_scale = distance_factor * stagnation_factor * boost_factor
            new_scale = (
                self.scale_momentum * current_scale
                + self.scale_learning_rate * target_scale
            )

            # Apply bounds
            self.property_scales[prop_name] = np.clip(
                new_scale, self.min_scale, self.max_scale
            )

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
            single_mut_count = (
                remaining - shuffle_count - multi_mut_count
            )  # remainder for single
        else:
            # Normal distribution
            shuffle_count = 0
            multi_mut_count = remaining // 2  # 1/2 for multi-mutations
            single_mut_count = remaining - multi_mut_count  # remainder for single

        # Generate single point mutations
        for _ in range(single_mut_count):
            candidates.append(make_point_mutation(sequence))

        # Generate multiple mutations with varied mutation counts
        for _ in range(multi_mut_count):
            mutated = sequence
            max_muts = max(
                self.min_mutations,
                min(self.max_mutations, len(sequence) // self.mutation_ratio),
            )

            # Use exponential distribution to favor smaller mutations but allow larger ones
            if max_muts > self.min_mutations:
                # Bias towards smaller mutations but allow larger ones occasionally
                weights = [
                    2 ** (max_muts - i) for i in range(self.min_mutations, max_muts + 1)
                ]
                num_muts = random.choices(
                    range(self.min_mutations, max_muts + 1), weights=weights
                )[0]
            else:
                num_muts = self.min_mutations

            for _ in range(num_muts):
                mutated = make_point_mutation(mutated)
            candidates.append(mutated)

        # Add shuffled variants if enabled
        if self.enable_shuffling and iteration % self.shuffle_frequency == 0:
            for _ in range(shuffle_count):
                candidates.append(
                    shuffle_sequence(
                        sequence,
                        window_size=self.shuffle_window_size,
                        global_shuffle_probability=self.global_shuffle_probability,
                    )
                )

        return candidates[: self.num_candidates]

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

        if len(self.best_error_history) < self.convergence_window:
            return False

        recent_errors = list(self.best_error_history)[-self.convergence_window :]
        error_std = np.std(recent_errors)

        # Check if error standard deviation is below tolerance
        is_converged = error_std < self.convergence_tolerance

        if is_converged:
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

        # Check if all properties are satisfied (either by per-property or global tolerance)
        property_info = self._get_best_sequence_property_info()

        for prop_name, info in property_info.items():
            per_property_tolerance = info.get("tolerance", 0.0)

            # If property has per-property tolerance, check against that
            if per_property_tolerance > 0.0:
                if not info.get("within_tolerance", False):
                    return False
            # Otherwise, check against global tolerance (if set)
            elif self.error_tolerance is not None:
                if info["raw_error"] > self.error_tolerance:
                    return False
            # If no tolerance set for this property, it must be exact (can't be satisfied)
            else:
                if info["raw_error"] > 1e-10:  # Allow for floating point precision
                    return False

        return True

    def _detect_stagnation(self) -> bool:
        """Detect if optimization is stagnant"""
        if len(self.best_error_history) < self.stagnation_threshold:
            return False

        recent_errors = list(self.best_error_history)[-self.stagnation_threshold :]

        # Check if there's been any significant improvement
        initial_error = recent_errors[0]
        final_error = recent_errors[-1]
        improvement = (initial_error - final_error) / max(initial_error, 1e-6)

        return improvement < self.stagnation_improvement_threshold

    def _apply_stagnation_recovery(self) -> None:
        """Apply recovery measures when stagnant"""
        if self.debugging:
            self.logger.info("🔄 Applying stagnation recovery...")

        # Get current property info to check raw errors and normalization factors
        current_property_info = self._get_best_sequence_property_info()

        worst_props = []
        normalization_fixes = []

        # Calculate total raw error for unsatisfied properties only
        # Use raw error instead of weighted error to avoid normalization bias
        total_raw_error = sum(
            info["raw_error"]
            for info in current_property_info.values()
            if not info.get("within_tolerance", False)
        )

        property_priorities = []

        for prop_name in self.property_scales:
            prop_info = current_property_info.get(prop_name, {})
            raw_error = prop_info.get("raw_error", 0)
            weighted_error = prop_info.get("weighted_error", 0)
            norm_factor = prop_info.get("normalization_factor", 1.0)
            within_tolerance = prop_info.get("within_tolerance", False)

            # Skip properties that are already satisfied (within tolerance)
            if within_tolerance:
                if self.debugging:
                    tolerance = prop_info.get("tolerance", 0.0)
                    self.logger.info(
                        f"   ✅ {prop_name} satisfied: error {raw_error:.4f} <= tolerance {tolerance:.4f}"
                    )
                continue

            # Calculate contribution based on raw error (avoids normalization bias)
            contribution = raw_error / total_raw_error if total_raw_error > 0 else 0

            # Find properties that haven't improved recently (trend-based stagnation)
            improvements = list(self.property_improvements[prop_name])
            trend_stagnant = False
            trend_score = 0.0
            if len(improvements) >= self.min_trend_samples:
                trend = np.polyfit(range(len(improvements)), improvements, 1)[0]
                # Use a stricter threshold: only boost if trend is positive (getting worse)
                # or very flat (>= -0.0001 instead of -0.001)
                if trend >= -0.0001:  # Much stricter than -0.001
                    trend_stagnant = True
                    trend_score = trend  # Higher score = worse trend

            # Check for normalization-based stagnation (high raw error + low normalization)
            normalization_stagnant = (raw_error > 0.3 and norm_factor < 0.2) or (
                raw_error > 0.8
            )

            # Calculate priority score (higher = more need for boost)
            priority_score = 0.0
            if trend_stagnant:
                priority_score += (
                    contribution * 2.0 + trend_score * 10.0
                )  # Weight by contribution and trend
            if normalization_stagnant:
                priority_score += raw_error * 5.0  # Weight by raw error magnitude

            if trend_stagnant or normalization_stagnant:
                property_priorities.append(
                    (
                        prop_name,
                        priority_score,
                        trend_stagnant,
                        normalization_stagnant,
                        raw_error,
                        norm_factor,
                        contribution,
                    )
                )

        # Sort by priority score (highest first) and only boost high-contribution properties
        property_priorities.sort(key=lambda x: x[1], reverse=True)

        # Only boost properties contributing >25% to total error (selective but not too restrictive)
        contribution_threshold = 0.25

        boosted_count = 0
        for (
            prop_name,
            priority_score,
            trend_stagnant,
            normalization_stagnant,
            raw_error,
            norm_factor,
            contribution,
        ) in property_priorities:
            # Only boost if contribution exceeds threshold
            if contribution < contribution_threshold:
                if self.debugging:
                    self.logger.info(
                        f"   ⏭️  Skipping {prop_name}: contribution {contribution * 100:.1f}% < {contribution_threshold * 100:.0f}% threshold"
                    )
                continue

            boosted_count += 1

            if trend_stagnant:
                worst_props.append(prop_name)

            # Check if property needs normalization fix (either explicitly detected or trend stagnant with low norm factor)
            needs_norm_fix = normalization_stagnant or (
                trend_stagnant and norm_factor < 0.5
            )

            if needs_norm_fix:
                normalization_fixes.append(prop_name)
                if self.debugging:
                    reason = (
                        "normalization issue"
                        if normalization_stagnant
                        else "low norm_factor + trend stagnant"
                    )
                    self.logger.info(
                        f"   🔧 Boosting {prop_name}: {contribution * 100:.1f}% contribution "
                        f"({reason}, raw_error: {raw_error:.3f}, norm_factor: {norm_factor:.3f})"
                    )

        # Boost scaling for trend-stagnant properties
        for prop_name in worst_props:
            self.property_scales[prop_name] *= self.stagnation_boost_factor
            self.property_scales[prop_name] = min(
                self.property_scales[prop_name], self.max_scale
            )

        # Fallback: if no properties meet threshold, boost the top contributor to avoid total stagnation
        if boosted_count == 0 and property_priorities:
            (
                prop_name,
                priority_score,
                trend_stagnant,
                normalization_stagnant,
                raw_error,
                norm_factor,
                contribution,
            ) = property_priorities[0]
            if self.debugging:
                self.logger.info(
                    f"   🚨 Fallback: No properties >{contribution_threshold * 100:.0f}% contribution, boosting top contributor {prop_name} ({contribution * 100:.1f}%)"
                )

            if trend_stagnant:
                worst_props.append(prop_name)

            # Apply same normalization fix logic as main loop
            needs_norm_fix = normalization_stagnant or (
                trend_stagnant and norm_factor < 0.5
            )
            if needs_norm_fix:
                normalization_fixes.append(prop_name)
                if self.debugging:
                    reason = (
                        "normalization issue"
                        if normalization_stagnant
                        else "low norm_factor + trend stagnant"
                    )
                    self.logger.info(
                        f"   🔧 Fallback normalization fix for {prop_name}: {reason} "
                        f"(raw_error: {raw_error:.3f}, norm_factor: {norm_factor:.3f})"
                    )

        # Apply normalization fixes for under-normalized properties
        for prop_name in normalization_fixes:
            # Get current property info to determine appropriate fix
            prop_info = current_property_info.get(prop_name, {})
            raw_error = prop_info.get("raw_error", 0)
            old_norm = self.normalization_factors.get(prop_name, 1.0)

            # Be more aggressive with normalization fixes based on raw error magnitude
            if raw_error > 0.8:
                # Very high raw error - jump to maximum normalization
                new_norm = 10.0
            elif raw_error > 0.5:
                # High raw error - boost significantly
                new_norm = min(10.0, old_norm * 5.0)
            else:
                # Moderate raw error - boost moderately
                new_norm = min(10.0, old_norm * 3.0)

            self.normalization_factors[prop_name] = new_norm

            if self.debugging:
                self.logger.info(
                    f"   🔧 Fixed {prop_name} normalization: {old_norm:.3f} → {new_norm:.3f} (raw_error: {raw_error:.3f})"
                )

            # Also boost scaling
            self.property_scales[prop_name] *= self.stagnation_boost_factor
            self.property_scales[prop_name] = min(
                self.property_scales[prop_name], self.max_scale
            )

        all_boosted = worst_props + normalization_fixes
        if self.debugging and all_boosted:
            norm_factors_per_prop = ""
            for prop in all_boosted:
                norm_factors_per_prop += f"   - {prop}: norm_factor={self.normalization_factors.get(prop, 1.0):.3f}, scale={self.property_scales.get(prop, 1.0):.3f}\n"
            self.logger.info(
                f"   Updated normalization factors and scales for boosted properties:\n{norm_factors_per_prop.strip()}"
            )
            self.logger.info(
                f"   ✅ Boosted {len(all_boosted)} properties: {', '.join(all_boosted)}"
            )
            if normalization_fixes:
                self.logger.info(
                    f"   Applied normalization fixes for: {', '.join(normalization_fixes)}"
                )

        # CRITICAL FIX: Recalculate best sequence error with new scaling factors
        # This ensures that candidates in the next iteration are evaluated fairly
        if all_boosted:
            old_error = self.best_error
            # Scaling changed, so we need fresh evaluation (don't use cache)
            self.best_error, _ = self._evaluate_sequence(
                self.best_sequence, use_cache=False
            )
            # Invalidate cached property info since scaling changed
            self._best_sequence_property_info = None
            if self.debugging:
                self.logger.info(
                    f"   🔄 Recalculated best sequence error: {old_error:.6f} → {self.best_error:.6f}"
                )

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
        diverse_sequences = build_diverse_initial_sequences(
            self.target_length, num_sequences=10
        )

        best_emergency_weighted_error = self.best_error
        best_emergency_raw_error = current_best_raw_error
        best_emergency_seq = self.best_sequence

        for seq in diverse_sequences:
            weighted_error, property_info = self._evaluate_sequence(seq)
            # Only count raw error from unsatisfied properties
            raw_error = sum(
                info["raw_error"]
                for info in property_info.values()
                if not info.get("within_tolerance", False)
            )

            # Accept if better weighted error and raw error doesn't get worse
            if (
                weighted_error < best_emergency_weighted_error
                and raw_error <= current_best_raw_error
            ):
                best_emergency_weighted_error = weighted_error
                best_emergency_raw_error = raw_error
                best_emergency_seq = seq

        # Only switch if we found something better
        if best_emergency_weighted_error < self.best_error:
            old_weighted_error = self.best_error
            old_raw_error = current_best_raw_error
            self._update_best_sequence(
                best_emergency_seq, best_emergency_weighted_error
            )

            if self.debugging:
                weighted_improvement = (
                    old_weighted_error - best_emergency_weighted_error
                )
                raw_improvement = old_raw_error - best_emergency_raw_error
                self.logger.info(
                    f"   Emergency sequence found! Weighted improvement: {weighted_improvement:.6f}, Raw improvement: {raw_improvement:.6f}"
                )

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
        patience: Optional[int] = None,
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
        self, tolerance: Optional[float] = None, enable: Optional[bool] = None
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
                self.logger.info(
                    f"Error tolerance stopping enabled: {self.error_tolerance:.6f}"
                )
            else:
                self.logger.info("Error tolerance stopping disabled")

    def get_convergence_info(self) -> Dict[str, Any]:
        """
        Get current convergence detection information.

        Returns
        -------
        Dict[str, Any]
            Dictionary with convergence status and parameters
        """
        if len(self.best_error_history) >= self.convergence_window:
            recent_errors = list(self.best_error_history)[-self.convergence_window :]
            current_std = np.std(recent_errors)
            convergence_ratio = (
                self.convergence_tolerance / current_std
                if current_std > 0
                else float("inf")
            )
        else:
            current_std = float("inf")
            convergence_ratio = 0.0

        return {
            "is_converged": self._check_convergence(),
            "convergence_counter": self.convergence_counter,
            "current_error_std": current_std,
            "tolerance": self.convergence_tolerance,
            "convergence_ratio": convergence_ratio,
            "window_size": self.convergence_window,
            "patience_remaining": max(
                0, self.convergence_patience - self.convergence_counter
            ),
            "early_stopping_enabled": self.enable_early_convergence,
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
            initial_candidates = build_diverse_initial_sequences(
                self.target_length, num_sequences=self.num_starting_candidates
            )
        else:
            initial_candidates = [self.starting_sequence]

        # set logging info
        if self.verbose:
            if self.starting_sequence is None:
                self.logger.info(
                    f"Generating and testing {len(initial_candidates)} diverse initial candidate sequences..."
                )
            else:
                self.logger.info(
                    f"Using provided starting sequence for optimization..."
                )

        # Batch evaluate initial candidates
        initial_results = self._evaluate_sequence(initial_candidates)
        best_initial_error = float("inf")
        best_initial_seq = ""

        # Select best initial candidate
        for i, (error, _) in enumerate(initial_results):
            if error < best_initial_error:
                best_initial_error = error
                best_initial_seq = initial_candidates[i]

        # Initialize optimization state
        self.best_sequence = best_initial_seq
        self.best_error = best_initial_error
        self.best_error_history.append(self.best_error)

        # Calculate dynamic normalization factors based on initial sequence
        if self.verbose:
            self.logger.info(
                f"Starting optimization (target length: {self.target_length})"
            )
            self.logger.info(
                f"Starting optimization with initial sequence: {self.best_sequence}"
            )

        # debugging info
        if self.debugging:
            self.logger.info(f"   Initial error: {self.best_error:.4f}")
            self.logger.info(f"📏 Calculating dynamic normalization factors:")

        # calculate initial normalization factors
        self._calculate_initial_normalization(self.best_sequence)

        # Recalculate initial error with normalization applied
        self.best_error, initial_property_info = self._evaluate_sequence(
            self.best_sequence
        )
        self.best_error_history[-1] = self.best_error  # Update the last entry

        if self.debugging:
            self.logger.info(f"   Normalized initial error: {self.best_error:.4f}")

        # Progress tracking
        progress_bar = tqdm(total=self.max_iterations, disable=not self.verbose)
        stagnation_counter = 0

        # Track best raw error to ensure it never increases (only count unsatisfied properties)
        _, initial_property_info = self._evaluate_sequence(self.best_sequence)
        best_raw_error = sum(
            info["raw_error"]
            for info in initial_property_info.values()
            if not info.get("within_tolerance", False)
        )

        # Main optimization loop
        for self.iteration in range(self.max_iterations):
            # Generate candidates
            candidates = self._generate_candidates(self.best_sequence, self.iteration)

            # Track improvement in this iteration
            improved_this_iteration = False

            # batch calculate
            candidate_results = self._evaluate_sequence(candidates)

            # find best candidate for this iteration
            best_candidate = None
            best_candidate_error = self.best_error
            best_candidate_raw_error = float("inf")
            best_candidate_property_info = None

            # Process batch results to find best candidate
            for i, (candidate, (weighted_error, property_info)) in enumerate(
                zip(candidates, candidate_results)
            ):
                # Skip the first candidate (current best sequence) unless it's actually better
                if i == 0 and candidate == self.best_sequence:
                    continue

                # Only count raw error from unsatisfied properties for meaningful progress tracking
                candidate_raw_error = sum(
                    info["raw_error"]
                    for info in property_info.values()
                    if not info.get("within_tolerance", False)
                )

                # Accept candidate if it improves weighted error AND doesn't worsen raw error
                if (
                    weighted_error < best_candidate_error
                    and candidate_raw_error <= best_raw_error
                ):
                    best_candidate = candidate
                    best_candidate_error = weighted_error
                    best_candidate_raw_error = candidate_raw_error
                    best_candidate_property_info = property_info

            # If we found a better candidate, update best sequence
            if best_candidate is not None:
                # Accept candidate if it improves weighted error AND doesn't worsen raw error
                if (
                    best_candidate_error < self.best_error
                    and best_candidate_raw_error <= best_raw_error
                ):
                    self._update_best_sequence(best_candidate, best_candidate_error)
                    # Only reset stagnation counter if raw error actually improves (stricter criteria)
                    if best_candidate_raw_error < best_raw_error:
                        improved_this_iteration = True
                        stagnation_counter = 0
                        if self.verbose:
                            improvement = best_raw_error - best_candidate_raw_error
                            if self.debugging:
                                self.logger.info(
                                    f"   ✅ Raw error improved by {improvement:.6f}"
                                )

                    best_raw_error = best_candidate_raw_error  # Update best raw error

                    # Update property scaling
                    if self.enable_adaptive_scaling:
                        self._update_property_scaling(best_candidate_property_info)

            # Track error history
            self.best_error_history.append(self.best_error)

            # Update progress bar once per iteration
            if self.verbose:
                progress_bar.update(1)

                # Update postfix display every update_interval iterations to avoid excessive overhead
                if self.iteration % self.update_interval == 0:
                    # Count satisfied properties for progress display (use cached version)
                    current_property_info = self._get_best_sequence_property_info()
                    satisfied_count = sum(
                        1
                        for info in current_property_info.values()
                        if info.get("within_tolerance", False)
                    )
                    total_properties = len(self.properties)

                    postfix = {
                        "raw_error": f"{best_raw_error:.4f}",
                        "satisfied": f"{satisfied_count}/{total_properties}",
                        "stagnation": stagnation_counter,
                    }

                    # Add error tolerance info if enabled
                    if self.enable_error_tolerance and self.error_tolerance is not None:
                        postfix["tolerance"] = f"{self.error_tolerance:.4f}"

                    # Add convergence info if enabled and enough history
                    if (
                        self.enable_early_convergence
                        and len(self.best_error_history) >= self.convergence_window
                    ):
                        conv_info = self.get_convergence_info()
                        postfix["conv"] = (
                            f"{self.convergence_counter}/{self.convergence_patience}"
                        )

                    progress_bar.set_postfix(postfix)

            # Check for convergence
            if self._check_convergence():
                if self.verbose:
                    conv_info = self.get_convergence_info()
                    self.logger.info(f"\\n✅ Converged at iteration {self.iteration}")
                if self.debugging:
                    self.logger.info(
                        f"   Error std: {conv_info['current_error_std']:.6f} < {conv_info['tolerance']:.6f}"
                    )
                    self.logger.info(
                        f"   Patience satisfied: {self.convergence_counter}/{self.convergence_patience}"
                    )
                break

            # Check for error tolerance
            if self._check_error_tolerance():
                if self.verbose:
                    self.logger.info(
                        f"\\n🎯 Error tolerance reached at iteration {self.iteration}"
                    )
                    self.logger.info(
                        f"   Current error: {self.best_error:.6f} <= {self.error_tolerance:.6f}"
                    )
                break

            # Handle stagnation
            if not improved_this_iteration:
                stagnation_counter += 1

                if stagnation_counter >= self.stagnation_threshold:
                    self._apply_stagnation_recovery()

                    # Recalculate best raw error after stagnation recovery (scaling may have changed)
                    # Use cached property info (already updated by stagnation recovery)
                    updated_property_info = self._get_best_sequence_property_info()
                    best_raw_error = sum(
                        info["raw_error"]
                        for info in updated_property_info.values()
                        if not info.get("within_tolerance", False)
                    )

                    # If severely stagnant, try emergency diversity injection
                    if stagnation_counter >= self.stagnation_threshold * 2:
                        best_raw_error = self._emergency_diversity_injection(
                            best_raw_error
                        )

                    stagnation_counter = 0  # Reset counter after recovery

        progress_bar.close()

        # Final evaluation and reporting
        final_property_info = self._get_best_sequence_property_info()
        final_error = self.best_error
        total_raw_error = sum(
            info["raw_error"] for info in final_property_info.values()
        )
        unsatisfied_raw_error = sum(
            info["raw_error"]
            for info in final_property_info.values()
            if not info.get("within_tolerance", False)
        )
        satisfied_count = sum(
            1
            for info in final_property_info.values()
            if info.get("within_tolerance", False)
        )

        if self.debugging:
            cache_stats = self.get_cache_statistics()
            self.logger.info(f"\n\n⚡ Performance:")
            self.logger.info(f"   Cache hit rate: {cache_stats['hit_rate']:.1%}")
            self.logger.info(
                f"   Evaluations saved: {cache_stats['total_evaluations_saved']}"
            )
            self.logger.info(
                f"   Total evaluations: {cache_stats['cache_hits'] + cache_stats['cache_misses']}"
            )

        if self.verbose:
            # Show cache performance
            self.logger.info(f"\n\nOptimization complete!")
            self.logger.info(
                f"   Final raw error for properties not within tolerance: {unsatisfied_raw_error:.6f}"
            )
            self.logger.info(
                f"   Properties satisfied: {satisfied_count}/{len(self.properties)}"
            )
            self.logger.info(f"   Total iterations: {self.iteration + 1}")

        if self.debugging:
            self.logger.info(f"   Final weighted error: {final_error:.6f}")
            self.logger.info(f"   Final raw error (all): {total_raw_error:.6f}")

        if self.verbose:
            self.logger.info(f"\n\nFinal property values:")

            for prop_name, info in final_property_info.items():
                target = info["target_value"]
                actual = info["raw_value"]
                raw_error = info["raw_error"]
                normalized_error = info["normalized_error"]
                norm_factor = info["normalization_factor"]
                tolerance = info.get("tolerance", 0.0)
                within_tolerance = info.get("within_tolerance", False)
                initial_error = self.initial_property_errors.get(prop_name, "N/A")

                status_indicator = "✅" if within_tolerance else "❌"
                tolerance_info = (
                    f", tolerance: {tolerance:.4f}" if tolerance > 0 else ""
                )
                constraint_type = info.get("constraint_type", "N/A")

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
        if hasattr(self, "_raw_property_cache"):
            self._raw_property_cache.clear()
        self._best_sequence_property_info = None

    def get_cache_statistics(self) -> Dict[str, Any]:
        """Get evaluation cache performance statistics"""
        total_requests = self._cache_hits + self._cache_misses
        hit_rate = self._cache_hits / total_requests if total_requests > 0 else 0.0

        return {
            "cache_hits": self._cache_hits,
            "cache_misses": self._cache_misses,
            "hit_rate": hit_rate,
            "cache_size": len(getattr(self, "_raw_property_cache", {})),
            "total_evaluations_saved": self._cache_hits,
        }
