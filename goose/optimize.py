import json
import logging
import math
import os
import pickle
import random
from collections import defaultdict
from typing import Any, Dict, List, Optional, Tuple, Union

import sparrow
from tqdm import tqdm

from goose.backend import optimizer_properties
from goose.backend.optimizer_properties import ProteinProperty
from goose.data import amino_acids


# ==============================================================================
# UTILITY FUNCTIONS
# ==============================================================================

def get_data_path(filename: str) -> str:
    """Get path to data file in goose/data directory."""
    root = os.path.abspath(os.path.dirname(__file__))
    return os.path.join(root, 'data', filename)


def get_property_class(class_name: str) -> type:
    """
    Get a property class by name from the optimizer_properties module.
    
    Parameters
    ----------
    class_name : str
        The name of the property class
        
    Returns
    -------
    type
        The property class
        
    Raises
    ------
    ValueError
        If the property class is not found
    """
    try:
        return getattr(optimizer_properties, class_name)
    except AttributeError:
        raise ValueError(f"Property class '{class_name}' not found in optimizer_properties module")


def validate_fixed_ranges(ranges: List[Tuple[int, int]], sequence_length: int) -> List[Tuple[int, int]]:
    """
    Validate and clip fixed residue ranges to sequence bounds.
    
    Parameters
    ----------
    ranges : List[Tuple[int, int]]
        List of (start, end) index pairs (inclusive)
    sequence_length : int
        Length of the sequence
        
    Returns
    -------
    List[Tuple[int, int]]
        Valid ranges clipped to sequence bounds
    """
    if not all(isinstance(r, tuple) and len(r) == 2 for r in ranges):
        raise ValueError("Each range must be a tuple of two integers.")
    
    valid_ranges = []
    for start, end in ranges:
        start = max(0, start)
        end = min(sequence_length - 1, end)
        if start <= end:
            valid_ranges.append((start, end))
    
    return valid_ranges


def ranges_overlap(range1: Tuple[int, int], range2: Tuple[int, int]) -> bool:
    """Check if two ranges overlap."""
    return range1[0] <= range2[1] and range2[0] <= range1[1]


# ==============================================================================
# KMER DICTIONARY CLASS
# ==============================================================================

class KmerDict:
    """
    Dictionary-based storage for k-mers and their associated properties.
    
    This class provides efficient storage and retrieval of k-mer properties
    used in sequence optimization.
    """
    
    def __init__(self, kmer_properties: Optional[Dict[str, dict]] = None):
        """
        Initialize KmerDict.
        
        Parameters
        ----------
        kmer_properties : Dict[str, dict], optional
            Initial k-mer properties to populate the dictionary
        """
        self.kmer_properties: Dict[str, dict] = kmer_properties or {}

    def insert(self, kmer: str, properties: dict) -> None:
        """Insert a k-mer and its properties."""
        self.kmer_properties[kmer] = properties

    def get_properties(self, kmer: str) -> dict:
        """Get properties for a specific k-mer."""
        return self.kmer_properties.get(kmer, {})

    def get_valid_kmers(self, min_k: int, max_k: int) -> List[str]:
        """Get all k-mers with lengths in the specified range."""
        return [
            kmer for kmer in self.kmer_properties 
            if min_k <= len(kmer) <= max_k
        ]

    @classmethod
    def from_dict(cls, kmer_properties: Dict[str, dict]) -> 'KmerDict':
        """Create KmerDict from a dictionary of k-mer properties."""
        return cls(kmer_properties)

    @classmethod
    def from_pickle(cls, filepath: str) -> 'KmerDict':
        """Load KmerDict from pickle file."""
        with open(filepath, "rb") as f:
            kmer_properties = pickle.load(f)
        return cls.from_dict(kmer_properties)

    @classmethod  
    def from_amino_acids(cls) -> 'KmerDict':
        """Create KmerDict using default amino acid properties."""
        return cls.from_dict(amino_acids.amino_acid_dat)

    def __len__(self) -> int:
        """Return number of k-mers stored."""
        return len(self.kmer_properties)

    def __contains__(self, kmer: str) -> bool:
        """Check if k-mer exists in dictionary."""
        return kmer in self.kmer_properties


# ==============================================================================
# SHUFFLE STRATEGIES
# ==============================================================================

class ShuffleStrategy:
    """
    Unified shuffling strategies for sequence optimization.
    
    This class consolidates different shuffling approaches into a single,
    configurable interface, removing the redundancy between shuffle_random_subset
    and get_global_and_window_shuffles.
    """
    
    @staticmethod
    def shuffle_entire(sequence: str) -> str:
        """Shuffle entire sequence randomly."""
        return ''.join(random.sample(sequence, len(sequence)))
    
    @staticmethod
    def shuffle_window(sequence: str, window_size: int, start_pos: Optional[int] = None) -> str:
        """
        Shuffle a specific window within the sequence.
        
        Parameters
        ----------
        sequence : str
            Input sequence
        window_size : int
            Size of window to shuffle
        start_pos : int, optional
            Starting position of window. If None, chosen randomly.
        """
        if not sequence or window_size <= 0:
            return sequence
            
        length = len(sequence)
        if window_size >= length:
            return ShuffleStrategy.shuffle_entire(sequence)
        
        if start_pos is None:
            max_start = length - window_size
            start_pos = random.randint(0, max_start) if max_start >= 0 else 0
        
        start_pos = max(0, min(start_pos, length - window_size))
        end_pos = start_pos + window_size
        
        window = sequence[start_pos:end_pos]
        shuffled_window = ''.join(random.sample(window, len(window)))
        
        return sequence[:start_pos] + shuffled_window + sequence[end_pos:]
    
    @staticmethod
    def shuffle_adaptive(sequence: str, global_prob: float = 0.2) -> str:
        """
        Adaptive shuffling: global shuffle with probability global_prob,
        otherwise shuffle a random window of 10-50% of sequence length.
        """
        if not sequence:
            return sequence
            
        if random.random() < global_prob:
            return ShuffleStrategy.shuffle_entire(sequence)
        
        length = len(sequence)
        min_window = max(1, int(length * 0.1))
        max_window = max(min_window, int(length * 0.5))
        window_size = random.randint(min_window, max_window)
        
        return ShuffleStrategy.shuffle_window(sequence, window_size)
    
    @staticmethod
    def shuffle_global_then_local(sequence: str, window_size: int) -> str:
        """
        First shuffle globally, then apply local window shuffles.
        This combines the functionality of get_global_and_window_shuffles.
        """
        # First global shuffle
        global_shuffled = ShuffleStrategy.shuffle_entire(sequence)
        
        # Then local window shuffles
        result = ""
        for i in range(0, len(global_shuffled), window_size):
            window = global_shuffled[i:i + window_size]
            shuffled_window = ''.join(random.sample(window, len(window)))
            result += shuffled_window
        
        return result
    
    @classmethod
    def generate_shuffled_variants(
        cls, 
        sequence: str, 
        num_shuffles: int, 
        window_size: int = 10,
        strategy: str = "mixed"
    ) -> List[str]:
        """
        Generate multiple shuffled variants using specified strategy.
        
        Parameters
        ----------
        sequence : str
            Input sequence
        num_shuffles : int
            Number of shuffled variants to generate
        window_size : int
            Window size for local shuffles
        strategy : str
            Shuffling strategy: 'adaptive', 'global_local', or 'mixed'
        """
        if num_shuffles <= 0:
            return []
        
        variants = []
        
        for _ in range(num_shuffles):
            if strategy == "adaptive":
                variants.append(cls.shuffle_adaptive(sequence))
            elif strategy == "global_local":
                variants.append(cls.shuffle_global_then_local(sequence, window_size))
            elif strategy == "mixed":
                # Mix different strategies
                if random.random() < 0.5:
                    variants.append(cls.shuffle_adaptive(sequence))
                else:
                    variants.append(cls.shuffle_global_then_local(sequence, window_size))
            else:
                raise ValueError(f"Unknown shuffle strategy: {strategy}")
        
        return variants


# ==============================================================================
# SEQUENCE BUILDER
# ==============================================================================

class SequenceBuilder:
    """Builds initial sequences for optimization."""
    
    @staticmethod
    def build_random_sequence(kmer_dict: KmerDict, target_length: int) -> str:
        """
        Build a random sequence by selecting k-mers from the dictionary.
        
        Parameters
        ----------
        kmer_dict : KmerDict
            K-mer dictionary to sample from
        target_length : int
            Target sequence length
            
        Returns
        -------
        str
            Randomly generated sequence
        """
        sequence = ""
        
        while len(sequence) < target_length:
            remaining_length = target_length - len(sequence)
            valid_kmers = kmer_dict.get_valid_kmers(1, remaining_length)
            
            if not valid_kmers:
                break
                
            new_kmer = random.choice(valid_kmers)
            sequence += new_kmer
        
        return sequence[:target_length]  # Ensure exact length


# ==============================================================================
# MUTATION ENGINE
# ==============================================================================

class MutationEngine:
    """
    Handles sequence mutations through k-mer replacement.
    
    This class consolidates mutation logic and provides a cleaner interface
    for sequence modifications during optimization.
    """
    
    def __init__(self, kmer_dict: KmerDict, min_k: int = 1, max_k: int = 1):
        """
        Initialize mutation engine.
        
        Parameters
        ----------
        kmer_dict : KmerDict
            K-mer dictionary for mutations
        min_k : int
            Minimum k-mer length to consider
        max_k : int  
            Maximum k-mer length to consider
        """
        self.kmer_dict = kmer_dict
        self.min_k = min_k
        self.max_k = max_k
    
    def get_mutable_kmers(self, sequence: str, fixed_ranges: List[Tuple[int, int]]) -> List[str]:
        """
        Get k-mers from sequence that can be mutated (not in fixed ranges).
        
        Parameters
        ----------
        sequence : str
            Current sequence
        fixed_ranges : List[Tuple[int, int]]
            Ranges that cannot be modified
            
        Returns
        -------
        List[str]
            List of k-mers that can be replaced
        """
        mutable_kmers = []
        
        for k in range(self.min_k, self.max_k + 1):
            for i in range(len(sequence) - k + 1):
                kmer = sequence[i:i + k]
                kmer_range = (i, i + k - 1)
                
                # Check if k-mer overlaps with any fixed range
                is_fixed = any(
                    ranges_overlap(kmer_range, fixed_range)
                    for fixed_range in fixed_ranges
                )
                
                if not is_fixed:
                    mutable_kmers.append(kmer)
        
        return mutable_kmers
    
    def filter_candidate_kmers(
        self, 
        directions: Dict[str, Dict[str, Any]]
    ) -> Dict[int, List[str]]:
        """
        Filter k-mers based on their potential to improve properties.
        
        Parameters
        ----------
        directions : Dict[str, Dict[str, Any]]
            Property direction information
            
        Returns
        -------
        Dict[int, List[str]]
            K-mers grouped by length, filtered by improvement potential
        """
        candidate_kmers = defaultdict(list)
        kmer_scores = []
        
        for kmer, prop_dict in self.kmer_dict.kmer_properties.items():
            # Calculate weighted improvement score
            weighted_improvement = 0.0
            total_weight = 0.0
            
            for method_name, info in directions.items():
                if method_name in prop_dict:
                    current_value = info["current_value"]
                    target_value = info["target_value"]
                    weight = info["weight"]
                    scale_range = info["scale_range"]
                    
                    # Calculate improvement potential
                    current_error = abs(current_value - target_value) / scale_range
                    kmer_value = prop_dict[method_name]
                    projected_error = abs(kmer_value - target_value) / scale_range
                    
                    improvement = (current_error - projected_error) * weight
                    weighted_improvement += improvement
                    total_weight += weight
            
            if total_weight > 0:
                normalized_score = weighted_improvement / total_weight
                kmer_scores.append((kmer, normalized_score))
        
        # Sort by improvement score and take top candidates
        kmer_scores.sort(key=lambda x: x[1], reverse=True)
        num_to_keep = max(10, len(kmer_scores) // 2)
        top_kmers = kmer_scores[:num_to_keep]
        
        for kmer, score in top_kmers:
            candidate_kmers[len(kmer)].append(kmer)
        
        # Fallback: if no good k-mers found, allow all k-mers
        if not candidate_kmers:
            for kmer in self.kmer_dict.kmer_properties:
                candidate_kmers[len(kmer)].append(kmer)
        
        return candidate_kmers
    
    def replace_kmer_in_sequence(
        self, 
        sequence: str, 
        old_kmer: str, 
        new_kmer: str,
        fixed_ranges: List[Tuple[int, int]],
        replacement_fraction: float = 0.25
    ) -> str:
        """
        Replace occurrences of a k-mer in the sequence.
        
        Parameters
        ----------
        sequence : str
            Original sequence
        old_kmer : str
            K-mer to replace
        new_kmer : str
            Replacement k-mer
        fixed_ranges : List[Tuple[int, int]]
            Ranges that cannot be modified
        replacement_fraction : float
            Fraction of occurrences to replace
            
        Returns
        -------
        str
            Modified sequence
        """
        if not old_kmer or old_kmer == new_kmer:
            return sequence
        
        # Find all valid replacement positions
        occurrences = []
        start_idx = 0
        kmer_len = len(old_kmer)
        
        while start_idx < len(sequence):
            pos = sequence.find(old_kmer, start_idx)
            if pos == -1:
                break
            
            kmer_range = (pos, pos + kmer_len - 1)
            
            # Check if position overlaps with fixed ranges
            is_fixed = any(
                ranges_overlap(kmer_range, fixed_range)
                for fixed_range in fixed_ranges
            )
            
            if not is_fixed:
                occurrences.append(pos)
            
            start_idx = pos + kmer_len
        
        if not occurrences:
            return sequence
        
        # Select random subset of occurrences to replace
        num_to_replace = max(1, math.ceil(len(occurrences) * replacement_fraction))
        selected_positions = sorted(
            random.sample(occurrences, min(num_to_replace, len(occurrences))),
            reverse=True
        )
        
        # Perform replacements (reverse order to maintain indices)
        new_sequence = sequence
        for pos in selected_positions:
            new_sequence = (
                new_sequence[:pos] + 
                new_kmer + 
                new_sequence[pos + kmer_len:]
            )
        
        return new_sequence
    
    def mutate_sequence(
        self, 
        sequence: str, 
        directions: Dict[str, Dict[str, Any]],
        fixed_ranges: List[Tuple[int, int]]
    ) -> str:
        """
        Perform a single mutation on the sequence.
        
        Parameters
        ----------
        sequence : str
            Current sequence
        directions : Dict[str, Dict[str, Any]]
            Property direction information
        fixed_ranges : List[Tuple[int, int]]
            Fixed ranges that cannot be mutated
            
        Returns
        -------
        str
            Mutated sequence
        """
        # Get candidate k-mers for replacement
        candidate_kmers = self.filter_candidate_kmers(directions)
        
        # Get mutable k-mers from current sequence
        mutable_kmers = self.get_mutable_kmers(sequence, fixed_ranges)

        if not mutable_kmers or not candidate_kmers:
            return sequence
        
        # Select k-mer to replace
        old_kmer = random.choice(mutable_kmers)
        
        # Find suitable replacement
        valid_replacements = candidate_kmers.get(len(old_kmer), [])
        if not valid_replacements:
            return sequence
        
        new_kmer = random.choice(valid_replacements)
        # Perform replacement
        return self.replace_kmer_in_sequence(
            sequence, old_kmer, new_kmer, fixed_ranges
        )


# ==============================================================================
# PROPERTY CALCULATOR
# ==============================================================================

class PropertyCalculator:
    """
    Calculates property errors and manages scaling.
    
    This class provides a clean interface for property calculations
    and error scaling during optimization.
    """
    
    @staticmethod
    def calculate_property_errors(
        protein: sparrow.Protein,
        properties: List[ProteinProperty],
        optimizer: Optional['SequenceOptimizer'] = None
    ) -> Tuple[float, Dict[str, Dict[str, Any]]]:
        """
        Calculate errors for all properties.
        
        Parameters
        ----------
        protein : sparrow.Protein
            Protein object for property calculation
        properties : List[ProteinProperty]
            List of properties to evaluate
        optimizer : SequenceOptimizer, optional
            Optimizer instance for error scaling
            
        Returns
        -------
        Tuple[float, Dict[str, Dict[str, Any]]]
            Combined error and detailed directions information
        """
        errors = []
        directions = {}
        
        for prop in properties:
            property_name = prop.__class__.__name__
            current_value = prop.calculate_raw_value(protein)
            raw_error = prop.calculate_error(current_value)
            
            # Apply scaling if optimizer provided
            if optimizer:
                scaled_error = optimizer._scale_property_error(
                    property_name, raw_error, prop, prop.weight
                )
            else:
                scaled_error = raw_error * prop.weight
            
            errors.append(scaled_error)
            
            # Store comprehensive direction info
            directions[property_name] = {
                'raw_error': raw_error,
                'scaled_error': scaled_error,
                'current_value': current_value,
                'target_value': prop.target_value,
                'weight': prop.weight,
                'is_normalized': prop.is_naturally_normalized(),
                'scale_range': (
                    optimizer.property_scaling_ranges.get(property_name, 1.0) 
                    if optimizer else 1.0
                )
            }
        
        combined_error = sum(errors)
        return combined_error, directions


# ==============================================================================
# MAIN SEQUENCE OPTIMIZER
# ==============================================================================

class SequenceOptimizer:
    """
    Main sequence optimizer class.
    
    This refactored version provides better organization and cleaner interfaces
    while maintaining all original functionality.
    """
    
    # Class-level cache for k-mer dictionaries
    _kmer_dict_cache: Dict[str, KmerDict] = {}
    
    def __init__(
        self, 
        target_length: int,
        kmer_dict_file: Optional[str] = None,
        verbose: bool = True,
        gap_to_report: int = 10,
        num_shuffles: int = 0,
        just_shuffle: bool = False
    ):
        """
        Initialize sequence optimizer.
        
        Parameters
        ----------
        target_length : int
            Target sequence length
        kmer_dict_file : str, optional
            Path to k-mer dictionary file
        verbose : bool
            Enable verbose logging
        gap_to_report : int
            Progress reporting interval
        num_shuffles : int
            Default number of shuffles per iteration
        just_shuffle : bool
            Only perform shuffling (no mutations)
        """
        # Core parameters
        self.target_length = target_length
        self.kmer_dict_file = kmer_dict_file
        self.verbose = verbose
        self.gap_to_report = gap_to_report
        self.num_shuffles = num_shuffles
        self.just_shuffle = just_shuffle
        
        # Optimization parameters (with defaults)
        self.max_iterations = 1000
        self.tolerance = 0.001
        self.window_size = 10
        self.shuffle_interval = 1
        self.initial_sequence: Optional[str] = None
        self.fixed_ranges: List[Tuple[int, int]] = []
        
        # Internal state
        self.properties_dict: Dict[str, ProteinProperty] = {}
        self._property_counter: Dict[str, int] = {}
        self.initial_property_values: Dict[str, float] = {}
        self.property_scaling_ranges: Dict[str, float] = {}
        
        # Initialize components
        self._configure_logger()
        self.kmer_dict = self._load_kmer_dict()
        self.mutation_engine = MutationEngine(self.kmer_dict)
    
    def _configure_logger(self) -> None:
        """Configure logging for the optimizer."""
        if self.verbose and not getattr(self, '_logger_configured', False):
            logging.basicConfig(
                level=logging.INFO,
                format='%(asctime)s - %(levelname)s - %(message)s'
            )
            self._logger_configured = True
        self.logger = logging.getLogger(__name__)
    
    def _load_kmer_dict(self) -> KmerDict:
        """Load k-mer dictionary with caching."""
        cache_key = self.kmer_dict_file or "default"
        
        if cache_key in self._kmer_dict_cache:
            self.logger.info(f"Using cached k-mer dictionary")
            return self._kmer_dict_cache[cache_key]
        
        if self.kmer_dict_file is None:
            # Use default amino acid properties
            kmer_dict = KmerDict.from_amino_acids()
            if self.verbose:
                self.logger.info("Using default amino acid properties")
        else:
            # Handle standard filenames
            if self.kmer_dict_file in ["kmer_properties.pickle", "amino_acids.pkl"]:
                filepath = get_data_path(self.kmer_dict_file)
            else:
                filepath = self.kmer_dict_file
            
            try:
                kmer_dict = KmerDict.from_pickle(filepath)
                if self.verbose:
                    self.logger.info(f"Loaded k-mer dictionary from {filepath}")
            except OSError as e:
                self.logger.error(f"Error loading k-mer dictionary: {e}")
                raise
        
        # Cache the result
        self._kmer_dict_cache[cache_key] = kmer_dict
        return kmer_dict
    
    def _calculate_initial_property_values(self, initial_sequence: str) -> None:
        """Calculate initial property values for scaling."""
        if self.initial_property_values:
            return  # Already calculated
        
        protein = sparrow.Protein(initial_sequence)
        
        for prop in self.properties:
            property_name = prop.__class__.__name__
            initial_value = prop.calculate_raw_value(protein)
            
            # Store initial value and set it in the property
            self.initial_property_values[property_name] = initial_value
            prop.set_initial_value(initial_value)
            
            # Calculate scaling range
            self.property_scaling_ranges[property_name] = prop.get_scaling_range()
    
    def _scale_property_error(
        self, 
        property_name: str, 
        raw_error: float,
        property_instance: ProteinProperty, 
        weight: float
    ) -> float:
        """Scale property error for optimization."""
        if property_instance.is_naturally_normalized():
            normalized_error = raw_error
        else:
            scale_range = self.property_scaling_ranges.get(property_name, 1.0)
            normalized_error = min(raw_error / scale_range, 1.0)

        return normalized_error * weight
    
    # Property management
    def add_property(self, property_class: type, *args: Any, **kwargs: Any) -> None:
        """
        Add a property to optimize.
        
        Parameters
        ----------
        property_class : type
            Property class to instantiate
        *args : Any
            Arguments for property initialization
        **kwargs : Any
            Keyword arguments for property initialization
        """
        new_property = property_class(*args, **kwargs)
        property_name = property_class.__name__
        
        # set seq length hint for length-depoendent props
        if hasattr(new_property, '_set_sequence_length_hint'):
            new_property.set_sequence_length_hint(self.target_length)

        # Handle multi-target properties
        allows_multiple = getattr(new_property, 'multi_target', False)
        constraint = getattr(new_property, 'constraint_type')
        
        if allows_multiple:
            # Generate unique identifier for multi-target properties
            if property_name not in self._property_counter:
                self._property_counter[property_name] = 0
            self._property_counter[property_name] += 1
            
            unique_key = f"{property_name}_{self._property_counter[property_name]}"
            self.properties_dict[unique_key] = new_property
            
            if self.verbose:
                self.logger.info(
                    f"Added property {property_name} targeting {constraint.value} "
                    f"value of {new_property.target_value} as {unique_key}"
                )
        else:
            # Replace existing or add new
            action = "Replaced" if property_name in self.properties_dict else "Added"
            self.properties_dict[property_name] = new_property
            
            if self.verbose:
                self.logger.info(
                    f"{action} property {property_name} targeting {constraint.value} "
                    f"value of {new_property.target_value}"
                )
    
    @property
    def properties(self) -> List[ProteinProperty]:
        """Get list of all properties."""
        return list(self.properties_dict.values())
    
    # Configuration methods
    def set_fixed_ranges(self, ranges: List[Tuple[int, int]]) -> None:
        """Set fixed residue ranges that cannot be mutated."""
        self.fixed_ranges = validate_fixed_ranges(ranges, self.target_length)
    
    def set_optimization_params(
        self,
        max_iterations: Optional[int] = None,
        tolerance: Optional[float] = None,
        window_size: Optional[int] = None,
        num_shuffles: Optional[int] = None,
        shuffle_interval: Optional[int] = None,
        just_shuffle: Optional[bool] = None
    ) -> None:
        """Set optimization parameters."""
        if max_iterations is not None:
            self.max_iterations = max_iterations
        if tolerance is not None:
            self.tolerance = tolerance
        if window_size is not None:
            self.window_size = window_size
        if num_shuffles is not None:
            self.num_shuffles = num_shuffles
        if shuffle_interval is not None:
            self.shuffle_interval = shuffle_interval
        if just_shuffle is not None:
            self.just_shuffle = just_shuffle
    
    def set_initial_sequence(self, sequence: str) -> None:
        """Set initial sequence for optimization."""
        if len(sequence) != self.target_length:
            raise ValueError(
                f"Initial sequence length ({len(sequence)}) must match "
                f"target length ({self.target_length})"
            )
        self.initial_sequence = sequence
    
    # Optimization execution
    def run(self) -> str:
        """
        Run the sequence optimization.
        
        Returns
        -------
        str
            Optimized sequence
        """
        if not self.properties:
            raise ValueError("No properties defined for optimization")
        
        self.logger.info("Starting sequence optimization")
        
        # Get initial sequence
        if self.initial_sequence is not None:
            best_sequence = self.initial_sequence
        else:
            best_sequence = SequenceBuilder.build_random_sequence(
                self.kmer_dict, self.target_length
            )
        
        # Calculate initial property values for scaling
        self._calculate_initial_property_values(best_sequence)
        
        # Initial error calculation
        best_error, directions = PropertyCalculator.calculate_property_errors(
            sparrow.Protein(best_sequence), self.properties, self
        )
        
        # Optimization loop
        best_sequence, final_error = self._optimization_loop(
            best_sequence, best_error, directions
        )
        
        self.logger.info("Sequence optimization completed")
        self._log_results(best_sequence)
        
        return best_sequence
    
    def _optimization_loop(
        self, 
        initial_sequence: str, 
        initial_error: float,
        initial_directions: Dict[str, Dict[str, Any]]
    ) -> Tuple[str, float]:
        """
        Main optimization loop.
        
        Parameters
        ----------
        initial_sequence : str
            Starting sequence
        initial_error : float
            Initial error value
        initial_directions : Dict
            Initial property directions
            
        Returns
        -------
        Tuple[str, float]
            Best sequence and its error
        """
        best_sequence = initial_sequence
        best_error = initial_error
        directions = initial_directions
        
        with tqdm(total=self.max_iterations, desc="Optimizing") as pbar:
            for iteration in range(self.max_iterations):
                # Generate candidate sequences
                candidates = self._generate_candidates(
                    best_sequence, directions, iteration
                )
                
                # Evaluate candidates
                for candidate in candidates:
                    protein = sparrow.Protein(candidate)
                    error, new_directions = PropertyCalculator.calculate_property_errors(
                        protein, self.properties, self
                    )
                    
                    if error < best_error:
                        best_sequence = candidate
                        best_error = error
                        directions = new_directions
                        
                        # Update property states
                        for prop in self.properties:
                            prop_name = prop.__class__.__name__
                            if prop_name in directions:
                                prop._best_value = directions[prop_name]['scaled_error']
                                prop._current_raw_value = directions[prop_name]['current_value']
                
                # Check convergence
                if best_error < self.tolerance:
                    self.logger.info(f"Converged at iteration {iteration}")
                    break
                
                # Update progress
                if iteration % self.gap_to_report == 0:
                    pbar.n = iteration  # Set absolute position
                    pbar.refresh()      # Refresh display
                    pbar.set_description(f"Best Error = {best_error:.5f}")
                    
                    if self.verbose:
                        self._log_iteration_details(iteration, best_error, best_sequence)
        
        return best_sequence, best_error
    
    def _generate_candidates(
        self, 
        sequence: str, 
        directions: Dict[str, Dict[str, Any]], 
        iteration: int
    ) -> List[str]:
        """
        Generate candidate sequences for evaluation.
        
        Parameters
        ----------
        sequence : str
            Current best sequence
        directions : Dict
            Property direction information
        iteration : int
            Current iteration number
            
        Returns
        -------
        List[str]
            List of candidate sequences
        """
        candidates = []
        
        # Determine if shuffling should occur this iteration
        do_shuffle = iteration % self.shuffle_interval == 0
        num_shuffles = self.num_shuffles if do_shuffle else 0
        
        if self.just_shuffle:
            # Only generate shuffled variants
            num_shuffles = max(1, num_shuffles)  # Ensure at least one shuffle
            candidates.extend(
                ShuffleStrategy.generate_shuffled_variants(
                    sequence, num_shuffles, self.window_size
                )
            )
        else:
            # Generate mutated sequence
            mutated = self.mutation_engine.mutate_sequence(
                sequence, directions, self.fixed_ranges
            )
            candidates.append(mutated)
            
            # Add shuffled variants if requested
            if num_shuffles > 0:
                candidates.extend(
                    ShuffleStrategy.generate_shuffled_variants(
                        mutated, num_shuffles, self.window_size
                    )
                )
        
        return candidates
    
    def _log_iteration_details(
        self, 
        iteration: int, 
        best_error: float, 
        best_sequence: str
    ) -> None:
        """Log detailed information about the current iteration."""
        self.logger.info(f"Iteration {iteration}: Best Error = {best_error:.5f}")
        self.logger.info(f"Best Sequence = {best_sequence}")
        
        for prop in self.properties:
            current_value = getattr(prop, '_current_raw_value', None)
            best_value = getattr(prop, '_best_value', None)
            
            current_str = f"{current_value:.3f}" if current_value is not None else "None"
            best_str = f"{best_value:.3f}" if best_value is not None else "None"
            
            self.logger.info(
                f"{prop.__class__.__name__} targeting "
                f"{prop.constraint_type.value} value of {prop.target_value}: "
                f"raw value is {current_str}, scaled error is {best_str}"
            )
    
    def _log_results(self, sequence: str) -> None:
        """Log final optimization results."""
        if self.initial_sequence is not None:
            self.logger.info(f"Initial Sequence: {self.initial_sequence}")
        
        self.logger.info(f"Optimized Sequence: {sequence}")
        
        protein = sparrow.Protein(sequence)
        for prop in self.properties:
            value = prop.calculate_raw_value(protein)
            self.logger.info(
                f"{prop.__class__.__name__}: {value:.2f} "
                f"(Target: {prop.constraint_type.value} {prop.target_value:.2f})"
            )
    
    # Configuration persistence
    def save_configuration(self, filename: str) -> None:
        """Save current configuration to JSON file."""
        properties_config = []
        
        for key, prop in self.properties_dict.items():
            prop_config = {
                "class_name": prop.__class__.__name__,
                "target_value": prop.target_value,
                "weight": prop.weight,
                "init_args": prop.get_init_args()
            }
            
            # Handle multi-target properties
            if hasattr(prop, 'multi_target') and prop.multi_target:
                prop_config.update({
                    "unique_key": key,
                    "multi_target": True
                })
            else:
                prop_config["multi_target"] = False
            
            properties_config.append(prop_config)
        
        config = {
            "target_length": self.target_length,
            "properties": properties_config,
            "fixed_ranges": self.fixed_ranges,
            "max_iterations": self.max_iterations,
            "tolerance": self.tolerance,
            "window_size": self.window_size,
            "num_shuffles": self.num_shuffles,
            "shuffle_interval": self.shuffle_interval,
            "initial_sequence": self.initial_sequence
        }
        
        with open(filename, 'w') as f:
            json.dump(config, f, indent=2)
        
        self.logger.info(f"Configuration saved to {filename}")
    
    @classmethod
    def load_configuration(
        cls, 
        config_file: str, 
        kmer_dict_file: Optional[str] = None
    ) -> 'SequenceOptimizer':
        """Load configuration from JSON file."""
        with open(config_file, 'r') as f:
            config = json.load(f)
        
        # Create optimizer instance
        optimizer = cls(config['target_length'], kmer_dict_file)
        
        # Set parameters
        optimizer.set_fixed_ranges(config['fixed_ranges'])
        optimizer.set_optimization_params(
            max_iterations=config['max_iterations'],
            tolerance=config['tolerance'],
            window_size=config['window_size'],
            num_shuffles=config['num_shuffles'],
            shuffle_interval=config['shuffle_interval']
        )
        
        if config.get('initial_sequence'):
            optimizer.set_initial_sequence(config['initial_sequence'])
        
        # Load properties
        for prop_config in config['properties']:
            try:
                property_class = get_property_class(prop_config['class_name'])
                init_args = prop_config.get('init_args', {})
                new_property = property_class.from_init_args(init_args)
                
                # Handle multi-target properties
                if prop_config.get('multi_target', False) and 'unique_key' in prop_config:
                    optimizer.properties_dict[prop_config['unique_key']] = new_property
                else:
                    optimizer.properties_dict[prop_config['class_name']] = new_property
                    
            except (KeyError, AttributeError, TypeError) as e:
                raise ValueError(
                    f"Error loading property {prop_config['class_name']}: {str(e)}"
                )
        
        optimizer.logger.info(f"Configuration loaded from {config_file}")
        return optimizer
    
    # Utility methods
    def get_defined_properties(self) -> List[Dict[str, Any]]:
        """Get list of defined properties with their details."""
        return [
            {
                "name": prop.__class__.__name__,
                "target_value": prop.target_value,
                "weight": prop.weight
            }
            for prop in self.properties
        ]
    
    def get_optimization_params(self) -> Dict[str, Any]:
        """Get current optimization parameters."""
        return {
            'max_iterations': self.max_iterations,
            'tolerance': self.tolerance,
            'verbose': self.verbose,
            'window_size': self.window_size,
            'num_shuffles': self.num_shuffles,
            'shuffle_interval': self.shuffle_interval,
            'fixed_residue_ranges': self.fixed_ranges,
            'initial_sequence': self.initial_sequence,
            'gap_to_report': self.gap_to_report,
            'just_shuffle': self.just_shuffle
        }
    
    def log_defined_properties(self) -> None:
        """Log all defined properties."""
        properties = self.get_defined_properties()
        if not properties:
            if self.verbose:
                self.logger.info("No properties defined.")
        else:
            if self.verbose:
                self.logger.info("Defined properties:")
                for prop in properties:
                    self.logger.info(
                        f"  - {prop['name']}: target = {prop['target_value']}, "
                        f"weight = {prop['weight']}"
                    )

