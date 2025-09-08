import json
import logging
import math
import os
import pickle
import random
from collections import defaultdict
from typing import Any, Dict, List, Optional, Tuple, Union
import numpy as np
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
    def shuffle_adaptive(sequence: str, global_prob: float = 0.1, conservative: bool = True, 
                        iteration: int = 0, max_iterations: int = 1000) -> str:
        """
        Adaptive shuffling: more aggressive early, more conservative later.
        """
        if not sequence:
            return sequence
        
        # Calculate progress ratio (0.0 = start, 1.0 = end)
        progress = min(1.0, iteration / max_iterations) if max_iterations > 0 else 0.0
        
        # Adjust global shuffle probability based on progress
        if conservative:
            # Start more aggressive, get more conservative
            early_prob = 0.2
            late_prob = 0.05
            adjusted_global_prob = early_prob * (1 - progress) + late_prob * progress
        else:
            adjusted_global_prob = global_prob
            
        if random.random() < adjusted_global_prob:
            return ShuffleStrategy.shuffle_entire(sequence)
        
        length = len(sequence)
        
        # Adjust window sizes based on progress
        if conservative:
            # Early: larger windows, Later: smaller windows
            early_min, early_max = 0.15, 0.4  # 15-40% early
            late_min, late_max = 0.03, 0.1    # 3-10% late
            
            min_frac = early_min * (1 - progress) + late_min * progress
            max_frac = early_max * (1 - progress) + late_max * progress
            
            min_window = max(1, int(length * min_frac))
            max_window = max(min_window, int(length * max_frac))
        else:
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
        strategy: str = "mixed",
        iteration: int = 0,
        max_iterations: int = 1000
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
        iteration : int
            Current iteration number
        max_iterations : int
            Maximum iterations for progress calculation
        """
        if num_shuffles <= 0:
            return []
        
        variants = []
        
        for _ in range(num_shuffles):
            if strategy == "adaptive":
                variants.append(cls.shuffle_adaptive(sequence, iteration=iteration, max_iterations=max_iterations))
            elif strategy == "global_local":
                variants.append(cls.shuffle_global_then_local(sequence, window_size))
            elif strategy == "mixed":
                # Mix different strategies
                if random.random() < 0.5:
                    variants.append(cls.shuffle_adaptive(sequence, iteration=iteration, max_iterations=max_iterations))
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
    
    def __init__(self, kmer_dict: KmerDict, min_k: int = 1, max_k: int = 3, 
                 conservative_mode: bool = True):
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
        conservative_mode : bool
            If True, use more conservative mutation parameters
        """
        self.kmer_dict = kmer_dict
        self.min_k = min_k
        self.max_k = max_k
        self.conservative_mode = conservative_mode
    
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
                    
                    # Calculate improvement potential using raw errors
                    current_error = abs(current_value - target_value)
                    kmer_value = prop_dict[method_name]
                    projected_error = abs(kmer_value - target_value)
                    
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
        replacement_fraction: Optional[float] = None,
        iteration: int = 0
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
        replacement_fraction : float, optional
            Fraction of occurrences to replace. If None, uses adaptive fraction.
        iteration : int
            Current iteration number for adaptive behavior
            
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
        
        # Adaptive replacement fraction based on iteration progress
        if replacement_fraction is None:
            if self.conservative_mode:
                # Early iterations: more aggressive exploration
                if iteration < 300:
                    base_fraction = 0.3 if len(old_kmer) == 1 else 0.2
                elif iteration < 1000:
                    base_fraction = 0.2 if len(old_kmer) == 1 else 0.15
                else:
                    # Later iterations: precise, conservative changes
                    base_fraction = 0.1 if len(old_kmer) == 1 else 0.05
                
                # Reduce fraction if there are many occurrences (avoid too many changes)
                if len(occurrences) > 15:
                    base_fraction *= 0.3
                elif len(occurrences) > 10:
                    base_fraction *= 0.5
                elif len(occurrences) > 5:
                    base_fraction *= 0.7
                
                replacement_fraction = base_fraction
            else:
                # Original behavior
                replacement_fraction = 0.25
        
        # Always replace at least 1, but be more conservative with many occurrences
        num_to_replace = max(1, min(3, math.ceil(len(occurrences) * replacement_fraction)))
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
        fixed_ranges: List[Tuple[int, int]],
        iteration: int = 0
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
        iteration : int
            Current iteration number (for adaptive behavior)
            
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
        
        # Adaptive k-mer selection: start with larger mutations, get more precise later
        if self.conservative_mode:
            # Early iterations: allow larger changes to explore sequence space
            if iteration < 300:
                preferred_length = random.choice([1, 2, 3])  # Explore broadly
            elif iteration < 1000:
                preferred_length = random.choice([1, 2])     # Medium precision
            else:
                # Later iterations: prefer single amino acid changes for fine-tuning
                preferred_length = 1
                
            # Filter mutable k-mers by preferred length
            length_filtered_kmers = [k for k in mutable_kmers if len(k) == preferred_length]
            if length_filtered_kmers:
                old_kmer = random.choice(length_filtered_kmers)
            else:
                old_kmer = random.choice(mutable_kmers)
        else:
            # Original random selection
            old_kmer = random.choice(mutable_kmers)
        
        # Find suitable replacement
        valid_replacements = candidate_kmers.get(len(old_kmer), [])
        if not valid_replacements:
            # Fallback: try different k-mer lengths
            for k_len in sorted(candidate_kmers.keys()):
                if candidate_kmers[k_len]:
                    # Find k-mers of this length in sequence
                    alt_kmers = [k for k in mutable_kmers if len(k) == k_len]
                    if alt_kmers:
                        old_kmer = random.choice(alt_kmers)
                        valid_replacements = candidate_kmers[k_len]
                        break
            
            if not valid_replacements:
                return sequence
        
        new_kmer = random.choice(valid_replacements)
        
        # Perform replacement with adaptive fraction
        return self.replace_kmer_in_sequence(
            sequence, old_kmer, new_kmer, fixed_ranges, iteration=iteration
        )


# ==============================================================================
# PROPERTY CALCULATOR
# ==============================================================================

class PropertyCalculator:
    """
    Calculates property errors.
    
    This class provides a clean interface for property calculations
    during optimization.
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
            Optimizer instance (unused, kept for compatibility)
            
        Returns
        -------
        Tuple[float, Dict[str, Dict[str, Any]]]
            Combined error and detailed directions
        """
        errors = []
        directions = {}
        
        for prop in properties:
            property_name = prop.__class__.__name__
            current_value = prop.calculate_raw_value(protein)
            raw_error = prop.calculate_error(current_value)
            
            # Use raw error with weight only
            weighted_error = raw_error * prop.weight

            # multiply weighted error by the property scaling
            weighted_error = weighted_error * prop._scaling_value
            
            errors.append(weighted_error)
            
            # Update current property value
            prop.current_raw_value = current_value

            # Store direction info
            directions[property_name] = {
                'raw_error': raw_error,
                'weighted_error': weighted_error,
                'current_value': current_value,
                'target_value': prop.target_value,
                'weight': prop.weight
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
        just_shuffle: bool = False,
        adaptive_scaling: bool = True,
        adaptive_scaling_interval: int = 100
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
        adaptive_scaling : bool
            Enable adaptive scaling of property weights
        adaptive_scaling_interval : int
            Number of iterations between adaptive scaling updates
        """
        # Core parameters
        self.target_length = target_length
        self.kmer_dict_file = kmer_dict_file
        self.verbose = verbose
        self.gap_to_report = gap_to_report
        self.num_shuffles = num_shuffles
        self.just_shuffle = just_shuffle
        
        # Adaptive scaling parameters
        self.adaptive_scaling = adaptive_scaling
        self.adaptive_scaling_interval = adaptive_scaling_interval
        
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
        self.iterations_since_improvement = 0
        self.best_error = np.inf
        # to hold generated sequences in case we get stuck during optimization and need to
        # reduce constraints to get out of the stuck state
        self.generated_sequences: List[str] = []
        
        # For adaptive scaling sensitivity analysis
        self._recent_candidate_data: List[Tuple[str, Dict[str, float]]] = []
        self._max_candidate_history = 50  # Keep last 50 candidates for sensitivity analysis

        # Initialize components
        self._configure_logger()
        self.kmer_dict = self._load_kmer_dict()
        self.mutation_engine = MutationEngine(self.kmer_dict, conservative_mode=adaptive_scaling)
    
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
        """Calculate initial property values."""
        if self.initial_property_values:
            return  # Already calculated
        
        protein = sparrow.Protein(initial_sequence)
        
        value_ranges=[]

        for prop in self.properties:
            property_name = prop.__class__.__name__
            initial_value = prop.calculate_raw_value(protein)
            
            # Store initial value and set it in the property
            self.initial_property_values[property_name] = initial_value
            prop.set_initial_value(initial_value)
            value_ranges.append(abs(initial_value-prop.target_value))

        # iterate through properteies and use the value_ranges to set scaling values
        if value_ranges:
            max_range = max(value_ranges)
            for prop in self.properties:
                property_name = prop.__class__.__name__
                initial_value = self.initial_property_values[property_name]
                range_ = abs(initial_value-prop.target_value)
                if range_==0:
                    prop._scaling_value=1.0
                else:
                    prop._scaling_value = max_range/range_
                if self.verbose:
                    print(prop.__class__.__name__, prop._scaling_value, prop.initial_value)

    def _recalculate_adaptive_scaling(self, current_sequence: str, candidate_data: Optional[List] = None) -> None:
        """
        Recalculate scaling values based on sensitivity analysis and proximity to target.
        
        This addresses the issue where highly sensitive properties (like MeanSelfEpsilon) 
        dominate optimization even when they're already close to target, preventing 
        less sensitive properties (like Hydrophobicity) from optimizing.
        
        Parameters
        ----------
        current_sequence : str
            Current best sequence
        candidate_data : List, optional
            List of (candidate_sequence, property_values_dict) from recent evaluations
            Used to estimate property sensitivity to mutations
        """
        if not self.adaptive_scaling or len(self.properties) <= 1:
            return
            
        protein = sparrow.Protein(current_sequence)
        
        # Calculate current state for all properties
        property_info = []
        for prop in self.properties:
            current_value = prop.calculate_raw_value(protein)
            raw_error = prop.calculate_error(current_value)
            
            # Calculate relative distance from target (normalized by target or 1.0 if target is 0)
            target = prop.target_value
            if target != 0:
                relative_distance = abs(current_value - target) / abs(target)
            else:
                relative_distance = abs(current_value - target)
            
            property_info.append((prop, current_value, raw_error, relative_distance))
        
        # Estimate sensitivity from candidate data if available
        sensitivity_factors = self._estimate_property_sensitivities(candidate_data) if candidate_data else {}
        
        # Calculate adaptive scaling
        for prop, current_value, raw_error, relative_distance in property_info:
            prop_name = prop.__class__.__name__
            
            if raw_error == 0:
                # Property is at target - give it very low priority
                prop._scaling_value = 0.05
                continue
            
            # Base scaling starts at 1.0
            base_scaling = 1.0
            
            # 1. Distance-based scaling: properties farther from target get higher priority
            max_relative_distance = max(info[3] for info in property_info)
            if max_relative_distance > 0:
                distance_factor = relative_distance / max_relative_distance
                distance_scaling = 0.5 + (1.5 * distance_factor)  # Scale between 0.5 and 2.0
            else:
                distance_scaling = 1.0
            
            # 2. Proximity penalty: reduce scaling for properties very close to target
            if relative_distance < 0.1:  # Within 10% of target
                proximity_penalty = 0.3  # Reduce scaling significantly
            elif relative_distance < 0.25:  # Within 25% of target  
                proximity_penalty = 0.6
            else:
                proximity_penalty = 1.0  # No penalty
            
            # 3. Sensitivity adjustment: reduce scaling for highly sensitive properties
            sensitivity_factor = sensitivity_factors.get(prop_name, 1.0)
            if sensitivity_factor > 2.0:  # Highly sensitive
                sensitivity_adjustment = 0.4
            elif sensitivity_factor > 1.5:  # Moderately sensitive
                sensitivity_adjustment = 0.7
            else:
                sensitivity_adjustment = 1.0
            
            # Combine all factors
            final_scaling = base_scaling * distance_scaling * proximity_penalty * sensitivity_adjustment
            
            # CRITICAL FIX: Ensure minimum scaling for properties far from target
            # Even sensitive properties need some attention if they're way off target
            if relative_distance > 1.0:  # More than 100% away from target
                min_scaling_for_distant = 0.2
            elif relative_distance > 0.5:  # More than 50% away from target  
                min_scaling_for_distant = 0.1
            else:
                min_scaling_for_distant = 0.05
            
            # Clamp to reasonable bounds, respecting minimum for distant properties
            prop._scaling_value = max(min_scaling_for_distant, min(3.0, final_scaling))
        
        # Recalculate best_error with new scaling values
        self.best_error, _ = PropertyCalculator.calculate_property_errors(
            protein, self.properties, self
        )
        
        if self.verbose:
            self.logger.info("Sensitivity-aware adaptive scaling updated:")
            for prop, current_value, raw_error, relative_distance in property_info:
                prop_name = prop.__class__.__name__
                sensitivity = sensitivity_factors.get(prop_name, "unknown")
                
                # Add status indicators
                status_flags = []
                if relative_distance < 0.1:
                    status_flags.append("NEAR_TARGET")
                if relative_distance > 1.0:
                    status_flags.append("FAR_FROM_TARGET")
                if isinstance(sensitivity, float) and sensitivity > 2.0:
                    status_flags.append("HIGH_SENSITIVITY")
                
                status_str = f" [{', '.join(status_flags)}]" if status_flags else ""
                
                self.logger.info(
                    f"  {prop_name}: scaling = {prop._scaling_value:.3f}, "
                    f"error = {raw_error:.3f}, rel_dist = {relative_distance:.3f}, "
                    f"sensitivity = {sensitivity}{status_str}"
                )
            self.logger.info(f"  New best_error: {self.best_error:.5f}")
            
            # Check for potentially incompatible properties
            incompatible_props = []
            for prop, current_value, raw_error, relative_distance in property_info:
                prop_name = prop.__class__.__name__
                # A property might be incompatible if it's both highly sensitive AND very far from target
                if (isinstance(sensitivity_factors.get(prop_name), float) and 
                    sensitivity_factors.get(prop_name) > 2.5 and 
                    relative_distance > 1.0):
                    incompatible_props.append((prop_name, current_value, prop.target_value))
            
            if incompatible_props:
                self.logger.warning("POTENTIAL INCOMPATIBILITY DETECTED:")
                for prop_name, current_val, target_val in incompatible_props:
                    self.logger.warning(
                        f"  {prop_name}: current={current_val:.2f}, target={target_val:.2f} "
                        f"(highly sensitive + far from target)"
                    )
                self.logger.warning(
                    "Consider: (1) Relaxing targets, (2) Using 'minimum'/'maximum' constraints, "
                    "(3) Checking if target combination is physically possible"
                )
            
    def _estimate_property_sensitivities(self, candidate_data: List) -> Dict[str, float]:
        """
        Estimate how sensitive each property is to sequence mutations.
        
        Parameters
        ----------
        candidate_data : List
            List of (sequence, {property_name: value}) from recent evaluations
            
        Returns
        -------
        Dict[str, float]
            Property name -> sensitivity factor (higher = more sensitive to changes)
        """
        if len(candidate_data) < 2:
            return {}
        
        # Calculate variance in property values across candidates
        property_variances = {}
        property_names = list(candidate_data[0][1].keys()) if candidate_data else []
        
        for prop_name in property_names:
            values = [data[1][prop_name] for data in candidate_data if prop_name in data[1]]
            if len(values) > 1:
                variance = np.var(values)
                property_variances[prop_name] = variance
        
        if not property_variances:
            return {}
        
        # Normalize variances to get relative sensitivity factors
        max_variance = max(property_variances.values())
        if max_variance == 0:
            return {name: 1.0 for name in property_names}
        
        sensitivity_factors = {
            name: variance / max_variance * 3.0 + 0.5  # Scale between 0.5 and 3.5
            for name, variance in property_variances.items()
        }
        
        return sensitivity_factors

    def trigger_adaptive_scaling(self, sequence: str) -> None:
        """
        Manually trigger adaptive scaling recalculation.
        
        Parameters
        ----------
        sequence : str
            Current sequence to evaluate for adaptive scaling
        """
        if not self.adaptive_scaling:
            self.logger.warning("Adaptive scaling is disabled. Enable it first.")
            return
            
        self._recalculate_adaptive_scaling(sequence, self._recent_candidate_data)
    
    def get_current_scaling_info(self) -> Dict[str, float]:
        """
        Get current scaling values for all properties.
        
        Returns
        -------
        Dict[str, float]
            Dictionary mapping property names to their current scaling values
        """
        return {prop.__class__.__name__: prop._scaling_value for prop in self.properties}
        
    def get_property_sensitivity_info(self) -> Dict[str, Dict[str, float]]:
        """
        Get detailed sensitivity and scaling information for all properties.
        
        Returns
        -------
        Dict[str, Dict[str, float]]
            Property name -> dict with scaling, error, distance info
        """
        if not self.properties:
            return {}
            
        # Get current sequence (use initial if no optimization has run)
        if hasattr(self, 'best_sequence'):
            current_sequence = self.best_sequence
        elif self.initial_sequence:
            current_sequence = self.initial_sequence
        else:
            return {}
        
        protein = sparrow.Protein(current_sequence)
        sensitivity_factors = self._estimate_property_sensitivities(self._recent_candidate_data)
        
        info = {}
        for prop in self.properties:
            prop_name = prop.__class__.__name__
            current_value = prop.calculate_raw_value(protein)
            raw_error = prop.calculate_error(current_value)
            
            target = prop.target_value
            if target != 0:
                relative_distance = abs(current_value - target) / abs(target)
            else:
                relative_distance = abs(current_value - target)
            
            info[prop_name] = {
                'current_scaling': prop._scaling_value,
                'raw_error': raw_error,
                'relative_distance': relative_distance,
                'sensitivity_factor': sensitivity_factors.get(prop_name, 'unknown'),
                'current_value': current_value,
                'target_value': target
            }
        
        return info
    
    def boost_distant_properties(self, distance_threshold: float = 1.0) -> List[str]:
        """
        Boost scaling for properties that have drifted very far from their targets.
        This can help recover properties that were over-suppressed by sensitivity scaling.
        
        Parameters
        ----------
        distance_threshold : float
            Relative distance threshold for boosting (default: 1.0 = 100% away from target)
            
        Returns
        -------
        List[str]
            Names of properties that were boosted
        """
        if not self.properties:
            return []
            
        # Get current sequence
        if hasattr(self, 'best_sequence'):
            current_sequence = self.best_sequence
        elif self.initial_sequence:
            current_sequence = self.initial_sequence
        else:
            return []
        
        protein = sparrow.Protein(current_sequence)
        boosted_properties = []
        
        for prop in self.properties:
            prop_name = prop.__class__.__name__
            current_value = prop.calculate_raw_value(protein)
            
            target = prop.target_value
            if target != 0:
                relative_distance = abs(current_value - target) / abs(target)
            else:
                relative_distance = abs(current_value - target)
            
            if relative_distance > distance_threshold:
                # Boost scaling significantly for distant properties
                prop._scaling_value = max(prop._scaling_value, 0.5)
                boosted_properties.append(prop_name)
                
        if self.verbose and boosted_properties:
            self.logger.info(f"Boosted scaling for distant properties: {boosted_properties}")
            
        return boosted_properties
    
    def analyze_optimization_state(self, current_sequence: str) -> Dict[str, Any]:
        """
        Analyze the current optimization state to identify potential issues.
        
        Parameters
        ----------
        current_sequence : str
            Current best sequence
            
        Returns
        -------
        Dict[str, Any]
            Analysis results with recommendations
        """
        protein = sparrow.Protein(current_sequence)
        analysis = {
            'properties': {},
            'scaling_issues': [],
            'recommendations': [],
            'total_error': self.best_error
        }
        
        total_unscaled_error = 0
        dominant_properties = []
        
        for prop in self.properties:
            prop_name = prop.__class__.__name__
            current_value = prop.calculate_raw_value(protein)
            raw_error = prop.calculate_error(current_value)
            weighted_error = raw_error * prop.weight * prop._scaling_value
            
            target = prop.target_value
            if target != 0:
                relative_distance = abs(current_value - target) / abs(target)
            else:
                relative_distance = abs(current_value - target)
            
            total_unscaled_error += raw_error * prop.weight
            
            # Check if this property is dominating the error
            error_contribution = weighted_error / self.best_error if self.best_error > 0 else 0
            if error_contribution > 0.7:
                dominant_properties.append(prop_name)
            
            analysis['properties'][prop_name] = {
                'current_value': current_value,
                'target_value': target,
                'raw_error': raw_error,
                'weighted_error': weighted_error,
                'scaling': prop._scaling_value,
                'weight': prop.weight,
                'relative_distance': relative_distance,
                'error_contribution': error_contribution,
                'iterations_since_improvement': prop.iterations_since_improvement
            }
            
            # Identify scaling issues
            if prop._scaling_value < 0.1 and relative_distance > 0.5:
                analysis['scaling_issues'].append(
                    f"{prop_name}: Very low scaling ({prop._scaling_value:.3f}) but far from target"
                )
            elif prop._scaling_value > 2.0 and relative_distance < 0.1:
                analysis['scaling_issues'].append(
                    f"{prop_name}: High scaling ({prop._scaling_value:.3f}) but already near target"
                )
        
        # Generate recommendations
        if dominant_properties:
            analysis['recommendations'].append(
                f"Properties {dominant_properties} are dominating the error. "
                "Consider reducing their scaling or relaxing their targets."
            )
        
        if analysis['scaling_issues']:
            analysis['recommendations'].append(
                "Scaling issues detected. Consider manual adjustment or disabling adaptive scaling temporarily."
            )
        
        if self.iterations_since_improvement > self.adaptive_scaling_interval * 3:
            analysis['recommendations'].append(
                f"No improvement for {self.iterations_since_improvement} iterations. "
                "Consider increasing mutation rate, changing shuffle parameters, or resetting scaling."
            )
        
        # Compare scaled vs unscaled error
        analysis['total_unscaled_error'] = total_unscaled_error
        analysis['scaling_effect_ratio'] = self.best_error / total_unscaled_error if total_unscaled_error > 0 else 1
        
        return analysis
    
    def reset_to_conservative_scaling(self) -> None:
        """
        Reset all properties to conservative, balanced scaling values.
        Useful when adaptive scaling has gone wrong.
        """
        for prop in self.properties:
            prop._scaling_value = 1.0
            
        if self.verbose:
            self.logger.info("Reset all properties to conservative scaling (1.0)")
    
    def force_balanced_scaling(self, current_sequence: str) -> None:
        """
        Force balanced scaling based only on distance from target, ignoring sensitivity.
        Useful when sensitivity-based scaling is preventing convergence.
        """
        protein = sparrow.Protein(current_sequence)
        
        # Calculate relative distances
        distances = []
        for prop in self.properties:
            current_value = prop.calculate_raw_value(protein)
            raw_error = prop.calculate_error(current_value)
            distances.append(raw_error)
        
        max_error = max(distances) if distances else 1
        
        for i, prop in enumerate(self.properties):
            if distances[i] == 0:
                prop._scaling_value = 0.1  # At target
            else:
                # Scale proportionally to error, between 0.2 and 2.0
                relative_error = distances[i] / max_error
                prop._scaling_value = 0.2 + (1.8 * relative_error)
        
        # Recalculate error
        self.best_error, _ = PropertyCalculator.calculate_property_errors(
            protein, self.properties, self
        )
        
        if self.verbose:
            self.logger.info("Applied balanced scaling based on error magnitude:")
            for prop in self.properties:
                self.logger.info(f"  {prop.__class__.__name__}: scaling = {prop._scaling_value:.3f}")
    
    def print_quick_diagnostic(self, current_sequence: str) -> None:
        """Print a quick diagnostic summary of the current optimization state."""
        analysis = self.analyze_optimization_state(current_sequence)
        
        print("\n=== OPTIMIZATION DIAGNOSTIC ===")
        print(f"Total Error: {analysis['total_error']:.3f}")
        print(f"Iterations since improvement: {self.iterations_since_improvement}")
        print(f"Adaptive scaling: {'ON' if self.adaptive_scaling else 'OFF'}")
        
        print("\nPROPERTY STATUS:")
        for prop_name, info in analysis['properties'].items():
            status = ""
            if info['error_contribution'] > 0.5:
                status += " [DOMINATING]"
            if info['relative_distance'] < 0.1:
                status += " [NEAR_TARGET]"
            if info['relative_distance'] > 1.0:
                status += " [FAR_FROM_TARGET]"
            if info['iterations_since_improvement'] > self.adaptive_scaling_interval:
                status += " [STAGNANT]"
                
            print(f"  {prop_name}: {info['current_value']:.3f} -> {info['target_value']:.3f} "
                  f"(error: {info['raw_error']:.3f}, scaling: {info['scaling']:.3f}){status}")
        
        if analysis['recommendations']:
            print("\nRECOMMENDATIONS:")
            for rec in analysis['recommendations']:
                print(f"  • {rec}")
        
        print("=== END DIAGNOSTIC ===\n")
    
    def adjust_mutation_aggressiveness(self, increase_aggressiveness: bool = False) -> None:
        """
        Adjust mutation engine aggressiveness during optimization.
        
        Parameters
        ----------
        increase_aggressiveness : bool
            If True, make mutations more aggressive. If False, make them more conservative.
        """
        if increase_aggressiveness:
            # Make mutations more aggressive
            self.mutation_engine.conservative_mode = False
            self.mutation_engine.max_k = min(5, self.mutation_engine.max_k + 1)
            if self.verbose:
                self.logger.info("Increased mutation aggressiveness")
        else:
            # Make mutations more conservative
            self.mutation_engine.conservative_mode = True
            self.mutation_engine.max_k = max(1, self.mutation_engine.max_k - 1)
            if self.verbose:
                self.logger.info("Decreased mutation aggressiveness")

    
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
        just_shuffle: Optional[bool] = None,
        adaptive_scaling: Optional[bool] = None,
        adaptive_scaling_interval: Optional[int] = None
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
        if adaptive_scaling is not None:
            self.adaptive_scaling = adaptive_scaling
        if adaptive_scaling_interval is not None:
            self.adaptive_scaling_interval = adaptive_scaling_interval
    
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
        
        # Calculate initial property values
        self._calculate_initial_property_values(best_sequence)
        
        # Initial error calculation
        self.best_error, directions = PropertyCalculator.calculate_property_errors(
            sparrow.Protein(best_sequence), self.properties, self
        )
        
        # Initialize property tracking values with initial results
        for prop in self.properties:
            prop_name = prop.__class__.__name__
            if prop_name in directions:
                prop.best_raw_value = directions[prop_name]['current_value']
                prop.best_raw_error = directions[prop_name]['raw_error'] 
                prop.best_weighted_error = directions[prop_name]['weighted_error']
                prop.iterations_since_improvement = 0
        
        # Optimization loop
        best_sequence, final_error = self._optimization_loop(
            best_sequence, self.best_error, directions
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
        self.best_error = initial_error
        directions = initial_directions
        
        with tqdm(total=self.max_iterations, desc="Optimizing") as pbar:
            for iteration in range(self.max_iterations):
                # Generate candidate sequences
                candidates = self._generate_candidates(
                    best_sequence, directions, iteration
                )

                # Increment iterations since improvement for all properties
                self.iterations_since_improvement += 1
                for prop in self.properties:
                    prop.iterations_since_improvement += 1

                # Evaluate candidates
                for candidate in candidates:
                    protein = sparrow.Protein(candidate)

                    error, new_directions = PropertyCalculator.calculate_property_errors(
                        protein, self.properties, self
                    )
                    
                    # Collect candidate data for sensitivity analysis
                    if self.adaptive_scaling:
                        candidate_property_values = {
                            prop_name: info['current_value'] 
                            for prop_name, info in new_directions.items()
                        }
                        self._recent_candidate_data.append((candidate, candidate_property_values))
                        
                        # Keep only recent candidates
                        if len(self._recent_candidate_data) > self._max_candidate_history:
                            self._recent_candidate_data.pop(0)

                    if error < self.best_error:
                        best_sequence = candidate
                        self.best_error = error
                        directions = new_directions
                        self.iterations_since_improvement = 0

                        # Update property states with best values
                        for prop in self.properties:
                            prop_name = prop.__class__.__name__
                            if prop_name in directions:
                                # Check if this property specifically improved
                                current_error = directions[prop_name]['raw_error']
                                if prop.best_raw_error is None or current_error < prop.best_raw_error:
                                    prop.iterations_since_improvement = 0
                                    
                                # Update best achieved values for this property
                                prop.best_raw_value = directions[prop_name]['current_value']
                                prop.best_raw_error = directions[prop_name]['raw_error']
                                prop.best_weighted_error = directions[prop_name]['weighted_error']

                # Check convergence
                if self.best_error < self.tolerance:
                    self.logger.info(f"Converged at iteration {iteration}")
                    break
                
                # Adaptive scaling check
                if (self.adaptive_scaling and 
                    iteration > 0 and 
                    iteration % self.adaptive_scaling_interval == 0):
                    # Check if we should temporarily disable adaptive scaling if stuck
                    if (self.iterations_since_improvement > self.adaptive_scaling_interval * 4 and
                        hasattr(self, '_adaptive_scaling_failures')):
                        self._adaptive_scaling_failures += 1
                        if self._adaptive_scaling_failures > 3:
                            if self.verbose:
                                self.logger.warning("Adaptive scaling seems to be hindering progress. Temporarily disabling.")
                            self.adaptive_scaling = False
                            self.force_balanced_scaling(best_sequence)
                    else:
                        self._recalculate_adaptive_scaling(best_sequence, self._recent_candidate_data)
                        if not hasattr(self, '_adaptive_scaling_failures'):
                            self._adaptive_scaling_failures = 0
                
                # Auto-adjust mutation aggressiveness if stuck (but be more careful about timing)
                if (iteration > 0 and 
                    iteration % (self.adaptive_scaling_interval * 2) == 0 and
                    self.iterations_since_improvement > self.adaptive_scaling_interval):
                    # Early in optimization: if stuck, try more aggressive exploration
                    if iteration < self.max_iterations * 0.3 and self.iterations_since_improvement > self.adaptive_scaling_interval * 3:
                        self.adjust_mutation_aggressiveness(increase_aggressiveness=True)
                    # Later in optimization: if we became aggressive but still stuck, go back to conservative
                    elif iteration > self.max_iterations * 0.3 and not self.mutation_engine.conservative_mode:
                        self.adjust_mutation_aggressiveness(increase_aggressiveness=False)
                
                # Update progress
                if iteration % self.gap_to_report == 0:
                    pbar.n = iteration  # Set absolute position
                    pbar.refresh()      # Refresh display
                    if iteration > 0:
                        pbar.set_description(f"Best Error = {self.best_error:.5f}")
                    else:
                        pbar.set_description(f"Best Error = N/A")
                    
                    if self.verbose:
                        self._log_iteration_details(iteration, best_sequence)
                
                # if final iteration, make sure that the progress bar is complete
                if iteration == self.max_iterations - 1:
                    pbar.n = self.max_iterations
                    pbar.refresh()
                    pbar.set_description(f"Best Error = {self.best_error:.5f}")
        
        return best_sequence, self.best_error
    
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
            # Use conservative shuffling when adaptive scaling is on
            conservative_shuffling = self.adaptive_scaling
            candidates.extend(
                ShuffleStrategy.generate_shuffled_variants(
                    sequence, num_shuffles, self.window_size, 
                    strategy="adaptive" if conservative_shuffling else "mixed",
                    iteration=iteration, max_iterations=self.max_iterations
                )
            )
        else:
            # Generate mutated sequence
            mutated = self.mutation_engine.mutate_sequence(
                sequence, directions, self.fixed_ranges, iteration
            )
            candidates.append(mutated)
            
            # Add shuffled variants if requested
            if num_shuffles > 0:
                # Use conservative shuffling when adaptive scaling is on
                conservative_shuffling = self.adaptive_scaling
                candidates.extend(
                    ShuffleStrategy.generate_shuffled_variants(
                        mutated, num_shuffles, self.window_size,
                        strategy="adaptive" if conservative_shuffling else "mixed",
                        iteration=iteration, max_iterations=self.max_iterations
                    )
                )
        
        return candidates
    
    def _log_iteration_details(
        self, 
        iteration: int, 
        best_sequence: str
    ) -> None:
        """Log detailed information about the current iteration."""
        self.logger.info(f"Iteration {iteration}:")
        self.logger.info(f"Best Sequence = {best_sequence}")
        
        for prop in self.properties:
            current_raw_value = getattr(prop, 'current_raw_value', None)
            best_raw_value = getattr(prop, 'best_raw_value', None)
            best_weighted_error = getattr(prop, 'best_weighted_error', None)
            
            current_str = f"{current_raw_value:.3f}" if current_raw_value is not None else "None"
            best_value_str = f"{best_raw_value:.3f}" if best_raw_value is not None else "None"
            best_error_str = f"{best_weighted_error:.3f}" if best_weighted_error is not None else "None"
            
            self.logger.info(
                f"{prop.__class__.__name__} targeting "
                f"{prop.constraint_type.value} value of {prop.target_value}: "
                f"current value = {current_str}, best value = {best_value_str}, "
                f"best weighted error = {best_error_str}"
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
            "initial_sequence": self.initial_sequence,
            "adaptive_scaling": self.adaptive_scaling,
            "adaptive_scaling_interval": self.adaptive_scaling_interval
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
            shuffle_interval=config['shuffle_interval'],
            adaptive_scaling=config.get('adaptive_scaling', True),
            adaptive_scaling_interval=config.get('adaptive_scaling_interval', 100)
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
            'just_shuffle': self.just_shuffle,
            'adaptive_scaling': self.adaptive_scaling,
            'adaptive_scaling_interval': self.adaptive_scaling_interval
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

