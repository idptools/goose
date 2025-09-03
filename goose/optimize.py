import random
import math
import os
from collections import defaultdict
from typing import Any, Dict, List, Tuple
from tqdm import tqdm
import json
import logging
import pickle

import sparrow

from goose.data import amino_acids
from goose.backend import optimizer_properties
from goose.backend.optimizer_properties import ProteinProperty


# code that allows access to the data directory
_ROOT = os.path.abspath(os.path.dirname(__file__))
def get_data(path):
    return os.path.join(_ROOT, 'data', path)

def _get_property_class(class_name: str) -> type:
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

# define kmer dict. 
class KmerDict:
    """
    A dictionary-based implementation for storing k-mers and their associated properties.
    """
    def __init__(self) -> None:
        self.kmer_properties: Dict[str, dict] = {}  # Dictionary to store properties for each k-mer

    def insert(self, kmer: str, properties: dict):
        """
        Insert a k-mer and its associated properties into the dictionary.
        Parameters
        ----------
        kmer : str
            The k-mer sequence to be inserted.
        properties : dict
            A dictionary of properties associated with the k-mer.
        """
        self.kmer_properties[kmer] = properties

    def get_kmer_properties(self, kmer: str) -> dict:
        """
        Retrieve the properties associated with a k-mer.
        Parameters
        ----------
        kmer : str
            The k-mer sequence.
        Returns
        -------
        dict
            A dictionary of properties for the k-mer, or an empty dictionary if the k-mer is not found.
        """
        return self.kmer_properties.get(kmer, {})

    def get_valid_kmers(self, k: int, max_length: int) -> List[str]:
        """
        Retrieve a list of valid k-mers of length between k and max_length stored in the dictionary.
        Parameters
        ----------
        k : int
            The minimum length of k-mers to retrieve.
        max_length : int
            The maximum length of k-mers to retrieve.
        Returns
        -------
        List[str]
            A list of valid k-mers with lengths between k and max_length.
        """
        return [kmer for kmer in self.kmer_properties if k <= len(kmer) <= max_length]
    



class SequenceOptimizer:
    '''
    This object minimizes the distance between a generated 
    IDR sequence and a target value stochastically.

    It finds a sequence that matches your property of interest
    via stochastic mutation of a starting reference protein sequence.
    '''
    # Add class-level cache dictionary
    _kmer_dict_cache: Dict[str, KmerDict] = {}

    def __init__(self, target_length: int, kmer_dict_file: str = None, 
                 verbose: bool = True, 
                 gap_to_report: int = 10, 
                 num_shuffles: int = 0, 
                 just_shuffle: bool = False):
        '''Initializes the sequence optimizer
        
        Parameters
        ----------
        target_length : int
            Target sequence length for property matching.
        kmer_dict_file : str, optional
            File location for k-mer bias (default: None).
        verbose : bool
            Enable logging actions/data (default: True).
        gap_to_report : int
            Progress reporting interval (default: 10).
        num_shuffles : int
            Number of shuffles (default: 0).
        just_shuffle : bool
            Only shuffle, do not mutate (default: False).
        '''
        #set some of the values passed in initialization of the seuqence optimizer
        self.target_length = target_length
        self.kmer_dict_file = kmer_dict_file
        self.verbose = verbose
        self.gap_to_report = gap_to_report

        #Initialize the kmer_dict as None
        self.kmer_dict = None 
        self.fixed_ranges: List[Tuple[int, int]] = []
        
        #These are properties that are not passed at initialization
        #These default are set unless the user changes them manually
        self.max_iterations = 1000
        self.tolerance = 0.001
        self.window_size = 10
        self.num_shuffles = num_shuffles
        self.shuffle_interval = 1
        self.initial_sequence = None #initialization of the sequence will occur later
        self.just_shuffle = just_shuffle
        self._configure_logger() #this configures a data logger to ensure that the process is recordered
        self._load_kmer_dict()
        self.properties_dict = {}  # Dictionary to store properties for faster lookup
        self._property_counter = {} # counter for generating unique property ids
        self.initial_property_values = {}
        self.property_scaling_ranges = {}

    def _calculate_initial_property_values(self, initial_sequence: str):
        """Calculate initial property values to establish scaling ranges."""
        if not self.initial_property_values:
            protein = sparrow.Protein(initial_sequence)
            
            for prop in self.properties:
                property_name = prop.__class__.__name__
                initial_value = prop.calculate_raw_value(protein)
                
                # Store initial value and set it in the property
                self.initial_property_values[property_name] = initial_value
                prop.set_initial_value(initial_value)
                
                # Calculate scaling range using property's own method
                self.property_scaling_ranges[property_name] = prop.get_scaling_range()

    def _scale_property_error(self, property_name: str, raw_error: float, 
                             property_instance: ProteinProperty, weight: float) -> float:
        """
        Scale property error to 0-1 range using property's scaling characteristics.
        """
        if property_instance.is_naturally_normalized():
            # Already 0-1 properties: error is naturally 0-1 scaled
            normalized_error = raw_error
        else:
            # Scale based on property's scaling range
            scale_range = self.property_scaling_ranges.get(property_name, 1.0)
            normalized_error = min(raw_error / scale_range, 1.0)  # Cap at 1.0
        
        return normalized_error * weight


    def _configure_logger(self):
        '''Initializes a data logger to store data and capture sequence optimization.'''
        if self.verbose and not getattr(self, '_logger_configured', False):
            logging.basicConfig(level=logging.INFO, format='%(asctime)s - %(levelname)s - %(message)s')
            self._logger_configured = True
        self.logger = logging.getLogger(__name__)

    def _load_kmer_dict(self):
        """
        Load and cache the kmer dictionary from file.
        """
        # Check if the kmer_dict for this file is already cached
        if self.kmer_dict_file in SequenceOptimizer._kmer_dict_cache:
            self.kmer_dict = SequenceOptimizer._kmer_dict_cache[self.kmer_dict_file]
            self.logger.info(f"Using cached kmer dictionary for {self.kmer_dict_file}")
        else:
            # if set to None, just using amino_acids.py. This is defualt because I'm going
            # to remove the use of a .pickle file as default due to the fact that it's slowly
            # becoming less commone due to perceived security risks.
            if self.kmer_dict_file is not None:
                try:
                    if self.kmer_dict_file in ["kmer_properties.pickle", "amino_acids.pkl"]:
                        self.kmer_dict_file = get_data(self.kmer_dict_file)

                    with open(self.kmer_dict_file, "rb") as f:
                        kmer_properties = pickle.load(f)
                    self.kmer_dict = build_kmer_dict(kmer_properties)

                    # Cache the kmer_dict for future use
                    SequenceOptimizer._kmer_dict_cache[self.kmer_dict_file] = self.kmer_dict
                    if self.verbose:
                        self.logger.info(f"Loaded and cached kmer dictionary from {self.kmer_dict_file}")
                except OSError as e:
                    self.logger.error(f"Error loading kmer dictionary: {e}")
                    raise
            else:
                if self.verbose:
                    self.logger.info("Using amino_acids.py for kmer properties")
                self.kmer_dict = build_kmer_dict(amino_acids.amino_acid_dat)

    @classmethod
    def load_configuration(cls, config_file: str, kmer_dict_file: str = None):
        '''This allows you to preload a SequenceOptimizer object from a previous run where you saved the configation in a file.
        
        Parameters
        ----------
        config_file : str
            The path to config file from the previous run that you want to reproduce
        kmer_dict_file : str
            This is the path to the kmer dictionary to used for amino acid frequencies.
            You do not need to load a kmer dictionary as this will go to the default value if you do not specify it.

        Returns
        -------
        SequenceOptimizer
            This is the sequence optimizer with the parameters loaded from the config file and the kmer dictionary that were chosen.
        '''
        with open(config_file, 'r') as f:
            config = json.load(f)

        optimizer = cls(config['target_length'], kmer_dict_file)
        optimizer.set_fixed_ranges(config['fixed_ranges'])
        optimizer.set_optimization_params(
            max_iterations=config['max_iterations'],
            tolerance=config['tolerance'],
            window_size=config['window_size'],
            num_shuffles=config['num_shuffles'],
            shuffle_interval=config['shuffle_interval']
        )
        if config['initial_sequence']:
            optimizer.set_initial_sequence(config['initial_sequence'])

        # Load properties with enhanced multi-target support
        for prop_config in config['properties']:
            try:
                # Get property class using the registry function
                property_class = _get_property_class(prop_config['class_name'])
                
                # Use the stored initialization arguments to recreate the property
                init_args = prop_config.get('init_args', {})
                new_property = property_class.from_init_args(init_args)
                
                # Add the property using the internal method to preserve unique keys for multi-target
                if prop_config.get('multi_target', False) and 'unique_key' in prop_config:
                    optimizer.properties_dict[prop_config['unique_key']] = new_property
                else:
                    optimizer.properties_dict[prop_config['class_name']] = new_property
                
            except (KeyError, AttributeError, TypeError) as e:
                raise ValueError(f"Error loading property {prop_config['class_name']}: {str(e)}")

        optimizer.logger.info(f"Configuration loaded from {config_file}")
        return optimizer

    def add_property(self, property_class: type, *args: Any, **kwargs: Any) -> None:
        '''Updates the optimizer with a new property to try to match on.
        
        This function takes in all the potential arguements and keyword arguments
        you will need passed at initialization for the property class of interest.

        Parameters
        ----------
        property_class : type
            This is the UNinitialized class that you are trying to add a target value for
        *args : Any
            These are any arguments that your class may take in order to initialize it
        **kwargs : Any
            These are any keyword arguments that your class may take to initialize it

        Returns
        -------
        None
            This function returns nothing. It similar adds the property to the stochastic
            minimization function
        '''
        new_property = property_class(*args, **kwargs)
        property_name = property_class.__name__

        # Check if this property type allows multiple targets
        allows_multiple = getattr(new_property, 'multi_target', False)
        constraint = getattr(new_property, 'constraint_type')
        # if allowed to add multiple...
        if allows_multiple:
            # Generate unique identifier for multi-target properties
            if property_name not in self._property_counter:
                self._property_counter[property_name] = 0
            self._property_counter[property_name] += 1
            
            unique_key = f"{property_name}_{self._property_counter[property_name]}"
            self.properties_dict[unique_key] = new_property
            
            if self.verbose:
                self.logger.info(f"Added property constraining {property_name} to {constraint.value} value of {new_property.target_value} as {unique_key}")
        else:
            # Replace if property already exists
            if property_name in self.properties_dict:
                self.properties_dict[property_name] = new_property
                if self.verbose:
                    self.logger.info(f"Replaced existing property {property_name} to {constraint.value} value of {new_property.target_value}")
            else:
                self.properties_dict[property_name] = new_property
                if self.verbose:
                    self.logger.info(f"Added new property {property_name} to {constraint.value} value of {new_property.target_value}")

    def set_fixed_ranges(self, ranges: List[Tuple[int, int]]):
        if not all(isinstance(r, tuple) and len(r) == 2 for r in ranges):
            raise ValueError("Each range must be a tuple of two integers.")
        self.fixed_ranges = ranges

    def set_optimization_params(self, max_iterations: int = None, tolerance: float = None, 
                                window_size: int = None, num_shuffles: int = None, 
                                shuffle_interval: int = None, just_shuffle: bool = None):
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

    def set_initial_sequence(self, sequence: str):
        '''This allows you  to initialize the optimization form a seed sequence that you specify.
        
        Parameters
        ----------
        sequence : str
            This is the seed sequence you want to start from. It MUST be the same length as the 
            target sequence in order to work.
        '''
        assert len(sequence) == self.target_length, f"Initial sequence length must match target length ({self.target_length})"
        self.initial_sequence = sequence

    def run(self) -> str:
        '''This performs the stocastic optimization on the properties you have added to this optimizer object.
        
        Returns
        -------
        str
            This is the sequence that the optimizer found that suits the needs you specified in the 
            protein properties you passed
        '''
        self.logger.info("Starting sequence optimization")
        optimized_sequence = optimize_sequence(
            self.kmer_dict,
            self.target_length,
            list(self.properties_dict.values()),  # Pass the list of properties
            **self.get_optimization_params()
        )
        self.logger.info("Sequence optimization completed")
        self.log_results(optimized_sequence)
        return optimized_sequence

    def get_optimization_params(self) -> dict:
        '''Returns a dictionary of the parameters this object has to find a sequence that matches specifiications.
        
        Returns
        -------
        dict
            This has the object properties as keywords and their values as the value in the dictionary.
        '''
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

    def log_results(self, sequence: str):
        '''Logs each property value for a sequence'''
        if self.initial_sequence is not None:
            self.logger.info("Initial Sequence: " + self.initial_sequence)
        self.logger.info(f"Optimized Sequence: {sequence}")
        protein = sparrow.Protein(sequence)
        for prop in self.properties_dict.values():
            value = prop.calculate_raw_value(protein)
            self.logger.info(f"{prop.__class__.__name__}: {value:.2f} (Target: {prop.constraint_type.value} {prop.target_value:.2f})")
    
    def save_configuration(self, filename: str):
        '''Saves the current configurations for the optimizer to a json file.
        
        Parameters
        ----------
        filename : str
            This is the json filename where you want to store the parameters for this optimization.
        '''
        # Enhanced to handle multi-target properties
        properties_config = []
        for key, prop in self.properties_dict.items():
            # For multi-target properties, store additional metadata
            if hasattr(prop, 'multi_target') and prop.multi_target:
                properties_config.append({
                    "class_name": prop.__class__.__name__,
                    "unique_key": key,
                    "target_value": prop.target_value,
                    "weight": prop.weight,
                    "multi_target": True,
                    # Store initialization args - this would need property-specific handling
                    "init_args": self._get_property_init_args(prop)
                })
            else:
                properties_config.append({
                    "class_name": prop.__class__.__name__,
                    "target_value": prop.target_value, 
                    "weight": prop.weight,
                    "multi_target": False,
                    "init_args": self._get_property_init_args(prop)
                })
        
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

    def _get_property_init_args(self, prop: ProteinProperty) -> dict:
        """
        Extract initialization arguments from a property instance.
        
        Parameters
        ----------
        prop : ProteinProperty
            The property instance to extract arguments from
            
        Returns
        -------
        dict
            Dictionary containing initialization arguments
        """
        return prop.get_init_args()

    def get_defined_properties(self) -> List[Dict[str, Any]]:
        '''This returns a list of the properties you are trying to optimize for during sequence generation.
        
        Returns
        -------
        List[Dict[str,Any]]
            The elements of the list represent each object.
            The keywords in each dictionary are the "name", "target_value", "weight" for each property added.
            The values in each keyword correspond to each of these values.
            No other information is returned.
        '''
        return [
            {
                "name": prop.__class__.__name__,
                "target_value": prop.target_value,
                "weight": prop.weight
            }
            for prop in self.properties_dict.values()
        ]
    
    @property
    def properties(self) -> List[ProteinProperty]:
        """
        Returns the defined properties.
        """
        return list(self.properties_dict.values())

    def log_defined_properties(self) -> None:
        """
        Logs the defined properties.
        """
        defined_properties = self.get_defined_properties()
        if not defined_properties:
            if self.verbose:
                self.logger.info("No properties defined.")
        else:
            if self.verbose:
                self.logger.info("Defined properties:")
                for prop in defined_properties:
                    self.logger.info(f"  - {prop['name']}: target = {prop['target_value']}, weight = {prop['weight']}")


def build_kmer_dict(kmer_properties: Dict[str, dict]) -> KmerDict:
    """
    Build a KmerDict from a dictionary of k-mer properties.

    Parameters
    ----------
    kmer_properties : Dict[str, dict]
        A dictionary mapping k-mers to their associated property dictionaries.

    Returns
    -------
    KmerDict
        The constructed KmerDict instance.
    """
    kmer_dict = KmerDict()
    for kmer, kmer_property in kmer_properties.items():
        kmer_dict.insert(kmer, kmer_property)
    return kmer_dict

def get_fixed_residue_ranges(sequence: str, fixed_ranges: List[Tuple[int, int]]) -> List[Tuple[int, int]]:
    """
    Get the fixed residue ranges in the sequence, ensuring they are within the sequence length.

    Parameters
    ----------
    sequence : str
        The input sequence.
    fixed_ranges : List[Tuple[int, int]]
        A list of tuples representing the start and end indices (inclusive) of the fixed residue ranges.

    Returns
    -------
    List[Tuple[int, int]]
        A list of valid fixed residue ranges within the sequence length.
    """
    valid_ranges = []
    seq_len = len(sequence)

    for start, end in fixed_ranges:
        start = max(0, start)
        end = min(seq_len - 1, end)

        if start <= end:
            valid_ranges.append((start, end))

    return valid_ranges


def get_global_and_window_shuffles(sequence: str, num_shuffles: int, window_size: int) -> List[str]:
    """
    Generate sequences with global and local shuffles.

    Parameters
    ----------
    sequence : str
        The original sequence.
    num_shuffles : int
        The number of global and local shuffled sequences to generate.
    window_size : int
        The size of the window to be shuffled.

    Returns
    -------
    List[str]
        A list containing the sequences with global and local shuffles.
    """
    shuffled_sequences = []

    for _ in range(num_shuffles):
        global_shuffled_sequence = "".join(random.sample(sequence, len(sequence)))

        local_shuffled_sequence = ""
        for i in range(0, len(global_shuffled_sequence), window_size):
            window = global_shuffled_sequence[i:i+window_size]
            shuffled_window = "".join(random.sample(window, len(window)))
            local_shuffled_sequence += shuffled_window

        shuffled_sequences.append(local_shuffled_sequence)

    return shuffled_sequences


def filter_candidate_kmers(sequence: str, kmer_dict: KmerDict, directions: Dict[str, Dict[str, Any]], fixed_residue_ranges: List[Tuple[int, int]]) -> Dict[int, List[str]]:
    """
    Filter candidate k-mers based on their property values, directionality, and fixed residue ranges.

    Parameters
    ----------
    sequence : str
        The current sequence.
    kmer_dict : KmerDict
        The KmerDict instance containing k-mer properties.
    directions : Dict[str, Dict[str, Any]]
        A dictionary containing the direction information for each property.
    fixed_residue_ranges : List[Tuple[int, int]]
        A list of tuples representing the start and end indices (inclusive) of the fixed residue ranges.

    Returns
    -------
    Dict[int, List[str]]
        A dictionary mapping k-mer lengths to lists of valid candidate k-mers.
    """
    current_properties = {method_name: info["current_value"] for method_name, info in directions.items()}
    target_properties = {method_name: info["target_value"] for method_name, info in directions.items()}
    weights = {method_name: info["weight"] for method_name, info in directions.items()}
    
    # Calculate weighted combined improvement for each k-mer
    candidate_kmers = defaultdict(list)
    kmer_scores = []
    
    for kmer, prop_dict in kmer_dict.kmer_properties.items():
        kmer_range = (sequence.find(kmer), sequence.find(kmer) + len(kmer) - 1)
        
        # Check if the k-mer overlaps with any fixed residue range
        if any(start <= kmer_range[0] <= end and start <= kmer_range[1] <= end for start, end in fixed_residue_ranges):
            continue
        
        # Calculate weighted improvement score
        weighted_improvement = 0.0
        total_weight = 0.0
        
        for method_name, info in directions.items():
            if method_name in prop_dict:
                current_value = info["current_value"]
                target_value = info["target_value"]
                weight = info["weight"]
                scale_range = info["scale_range"]
                
                # Calculate current error (normalized)
                current_error = abs(current_value - target_value) / scale_range
                
                # Calculate projected error with this k-mer
                kmer_value = prop_dict[method_name]
                projected_error = abs(kmer_value - target_value) / scale_range
                
                # Improvement is reduction in error (positive = better)
                improvement = (current_error - projected_error) * weight
                weighted_improvement += improvement
                total_weight += weight
        
        if total_weight > 0:
            # Normalize by total weight
            normalized_score = weighted_improvement / total_weight
            kmer_scores.append((kmer, normalized_score))
    
    # Sort k-mers by their improvement score and take top candidates
    kmer_scores.sort(key=lambda x: x[1], reverse=True)
    
    # Take top 50% of k-mers or at least 10 k-mers per length
    num_to_keep = max(10, len(kmer_scores) // 2)
    top_kmers = kmer_scores[:num_to_keep]
    
    for kmer, score in top_kmers:
        candidate_kmers[len(kmer)].append(kmer)
    
    # Fallback: if no good k-mers found, allow all k-mers
    if not candidate_kmers:
        for kmer in kmer_dict.kmer_properties:
            candidate_kmers[len(kmer)].append(kmer)
    
    return candidate_kmers

def select_kmer_to_replace(sequence: str, min_k: int, max_k: int, fixed_residue_ranges: List[Tuple[int, int]]) -> str:
    """
    Select a k-mer from the current sequence to be replaced, excluding k-mers that overlap with fixed residue ranges.

    Parameters
    ----------
    sequence : str
        The current sequence.
    min_k : int
        The minimum length of k-mers to consider.
    max_k : int
        The maximum length of k-mers to consider.
    fixed_residue_ranges : List[Tuple[int, int]]
        A list of tuples representing the start and end indices (inclusive) of the fixed residue ranges.

    Returns
    -------
    str
        The selected k-mer to be replaced.
    """
    current_sequence_kmers = []
    for k in range(min_k, max_k + 1):
        for i in range(len(sequence) - k + 1):
            kmer = sequence[i:i+k]
            kmer_range = (i, i+k-1)

            # Check if any part of the k-mer range intersects with any fixed residue range
            if not any((start <= kmer_range[1] and end >= kmer_range[0]) for start, end in fixed_residue_ranges):
                current_sequence_kmers.append(kmer)


    if current_sequence_kmers:
        kmer_to_replace = random.choice(current_sequence_kmers)
    else:
        # If no k-mers are available for replacement, return None
        kmer_to_replace = None

    return kmer_to_replace


def replace_kmer(sequence: str, kmer_to_replace: str, new_kmer: str, fixed_residue_ranges: List[Tuple[int, int]]) -> str:
    """
    Replace 50% of the occurrences of a k-mer in the given sequence with a new k-mer, respecting fixed residue ranges.

    Parameters
    ----------
    sequence : str
        The original sequence.
    kmer_to_replace : str
        The k-mer to be replaced.
    new_kmer : str
        The new k-mer to replace with.
    fixed_residue_ranges : List[Tuple[int, int]]
        A list of tuples representing the start and end indices (inclusive) of the fixed residue ranges.

    Returns
    -------
    str
        The new sequence with the k-mer replaced, respecting fixed residue ranges.
    """
    if not kmer_to_replace:
        return sequence

    occurrences = []
    start_idx = 0
    kmer_len = len(kmer_to_replace)

    # Find all occurrences of the k-mer in the sequence
    while start_idx < len(sequence):
        kmer_pos = sequence.find(kmer_to_replace, start_idx)
        if kmer_pos == -1:
            break
        kmer_range = (kmer_pos, kmer_pos + kmer_len - 1)

        # Check if any part of the k-mer range intersects with any fixed residue range
        overlaps_fixed_range = any((start <= kmer_range[1] and end >= kmer_range[0]) for start, end in fixed_residue_ranges)

        if not overlaps_fixed_range:
            # Add the index of the occurrence if it doesn't overlap with fixed residue range
            occurrences.append(kmer_pos)

        start_idx = kmer_pos + kmer_len

    if occurrences:
        num_to_replace = math.ceil(len(occurrences)//4+1)
        selected_occurrences = sorted(random.sample(occurrences, num_to_replace), reverse=True)
        
        # Replace the selected occurrences with the new k-mer
        new_sequence = sequence
        for selected_occurrence in selected_occurrences:
            new_sequence = (
                new_sequence[:selected_occurrence] +
                new_kmer +
                new_sequence[selected_occurrence + kmer_len:]
            )
    else:
        # If no suitable occurrences found, return the original sequence
        new_sequence = sequence

    return new_sequence


def shuffle_random_subset(s: str) -> str:
    """
    Shuffle a sequence with:
    - 20% chance to shuffle entire string
    - 80% chance to shuffle a random window of 10-50% of the string
    """
    if not s:
        return s
    length = len(s)
    if random.random() < 0.2:
        return ''.join(random.sample(s, length))
    min_window = max(1, int(length * 0.1))
    max_window = max(min_window, int(length * 0.5))
    window_size = random.randint(min_window, max_window)
    max_start = length - window_size
    if max_start < 0:
        return s
    start_pos = random.randint(0, max_start)
    substring = s[start_pos:start_pos + window_size]
    shuffled_substring = ''.join(random.sample(substring, len(substring)))
    return s[:start_pos] + shuffled_substring + s[start_pos + window_size:]


def mutate_sequence(sequence: str, kmer_dict: KmerDict, 
                    directions: Dict[str, Dict[str, Any]], 
                    properties: List[ProteinProperty], num_shuffles: int = 0, 
                    window_size: int = 10, 
                    fixed_residue_ranges: List[Tuple[int, int]] = [],
                    just_shuffle: bool = False,
                     optimizer: SequenceOptimizer = None) -> Tuple[str, float, Dict[str, Dict[str, Any]]]:
    """
    Mutate a sequence by replacing a k-mer with a new k-mer and generate shuffled sequences.

    Parameters
    ----------
    sequence : str
        The input sequence.
    kmer_dict : KmerDict
        The KmerDict instance containing k-mer properties.
    directions : Dict[str, Dict[str, Any]]
        A dictionary containing the direction information for each property.
    properties : List[ProteinProperty]
        A list of protein properties.
    num_shuffles : int, optional
        The number of shuffles to perform, by default 0.
    window_size : int, optional
        The window size for shuffling, by default 10.
    fixed_residue_ranges : List[Tuple[int, int]], optional
        A list of tuples representing the start and end indices (inclusive) of the fixed residue ranges, by default [].
    just_shuffle : bool, optional
        Whether to just shuffle the sequence, by default False.

    Returns
    -------
    Tuple[str, float, Dict[str, Dict[str, Any]]]
        A tuple containing the new sequence, the combined error, and the updated directions dictionary.
    """
 
 
    # set min and max length for kmer, default is 1 because 
    # default usage is the amino_acids.py file. 
    min_k = 1
    max_k = 1

    # if we aren't just shuffling the sequence, we need to get the candidate kmers
    if not just_shuffle:
        candidate_kmers = filter_candidate_kmers(sequence=sequence, 
                                                 kmer_dict=kmer_dict, 
                                                 directions=directions,
                                                 fixed_residue_ranges = fixed_residue_ranges)

        kmer_to_replace = select_kmer_to_replace(sequence=sequence, 
                                                 min_k=min_k, 
                                                 max_k=max_k, 
                                                 fixed_residue_ranges=fixed_residue_ranges)


        if kmer_to_replace and candidate_kmers:
            valid_kmers = candidate_kmers.get(len(kmer_to_replace), [])
            if valid_kmers:
                new_kmer = random.choice(valid_kmers)
                new_sequence = replace_kmer(sequence, kmer_to_replace, new_kmer, fixed_residue_ranges)     
            else:
                new_sequence = sequence
        else:
            new_sequence = sequence
        
        sequences = [new_sequence] + get_global_and_window_shuffles(new_sequence, num_shuffles, window_size) if num_shuffles > 0 else [new_sequence]

    else:
        if num_shuffles==0:
            num_shuffles=1
        sequences = [shuffle_random_subset(sequence) for _ in range(num_shuffles)]

    min_combined_error = float('inf')
    best_sequence = None
    best_directions = None

    for seq in sequences:
        protein = sparrow.Protein(seq)
        combined_err, directions = calculate_errors(protein, properties, optimizer)

        if combined_err < min_combined_error:
            min_combined_error = combined_err
            best_sequence = seq
            best_directions = directions

    return best_sequence, min_combined_error, best_directions


def optimize_sequence(kmer_dict: KmerDict, target_length: int, properties: List[ProteinProperty],
                      max_iterations: int = 100, tolerance: float = 1e-2, verbose=False,
                      window_size: int = 10, num_shuffles: int = 0, shuffle_interval: int = 50,
                      fixed_residue_ranges: List[Tuple[int, int]] = [], initial_sequence : str = None,
                      gap_to_report: int = 10, just_shuffle : bool =False) -> str:
    """
    Optimize a protein sequence to minimize the combined error between its computed properties and the target property values.

    Parameters
    ----------
    kmer_dict : KmerDict
        The KmerDict instance containing k-mer properties.
    target_length : int
        The target length of the sequence.
    properties : List[ProteinProperty]
        A list of protein properties.
    max_iterations : int, optional
        The maximum number of iterations, by default 100.
    tolerance : float, optional
        The error tolerance for stopping the optimization, by default 1e-2.
    verbose : bool, optional
        Whether to print verbose output, by default False.
    window_size : int, optional
        The window size for shuffling, by default 10.
    num_shuffles : int, optional
        The number of shuffles to perform, by default 100.
    shuffle_interval : int, optional
        The interval at which to perform shuffles, by default 50.
    fixed_residue_ranges : List[Tuple[int, int]], optional
        A list of tuples representing the start and end indices (inclusive) of the fixed residue ranges, by default [].
    initial_sequence : str, optional
        The initial sequence to start the optimization, by default None.
    gap_to_report : int, optional
        The interval at which to report progress, by default 10.
    just_shuffle : bool, optional
        Whether to just shuffle the sequence, by default False.

    Returns
    -------
    str
        The optimized sequence.
    """
    # Create optimizer instance for scaling
    optimizer_tracker = SequenceOptimizer(target_length)
    optimizer_tracker.properties_dict = {prop.__class__.__name__: prop for prop in properties}
    
    # check initial sequence    
    if initial_sequence is not None:
        best_sequence = initial_sequence
    else:
        best_sequence = build_sequence(kmer_dict, target_length)

    # Calculate initial property values for scaling - this is crucial!
    optimizer_tracker._calculate_initial_property_values(best_sequence)    

    best_error, directions = calculate_errors(sparrow.Protein(best_sequence), properties)

    # Start progress bar
    pbar = tqdm(total=max_iterations)

    # Iterate
    for i in range(max_iterations):
        current_num_shuffles = num_shuffles if i % shuffle_interval == 0 else 0
        new_sequence, new_error, directions = mutate_sequence(sequence=best_sequence, 
                                                              kmer_dict=kmer_dict, 
                                                              directions=directions, 
                                                              properties=properties, 
                                                              num_shuffles=current_num_shuffles, 
                                                              window_size=window_size, 
                                                              fixed_residue_ranges=fixed_residue_ranges, 
                                                              just_shuffle=just_shuffle,
                                                              optimizer=optimizer_tracker)
    
        if new_error < best_error:
            best_sequence = new_sequence
            best_error = new_error
            for p in properties:
                p._best_value = directions[p.__class__.__name__]['scaled_error']
                p._current_raw_value=directions[p.__class__.__name__]['current_value']

        if best_error < tolerance:
            break

        # Update progress bar
        if i!=0 and i % gap_to_report == 0:
            pbar.update(gap_to_report)
            pbar.set_description(f"Best Error = {best_error:.5f}")
        
        # make sure we update at last step
        if i == max_iterations - 1:
            pbar.update(gap_to_report)
            pbar.set_description(f"Best Error = {best_error:.5f}")

        # if verbose..
        if verbose and i % gap_to_report == 0:
            print(f"Iteration {i}:\nBest Error = {best_error:.5f}")
            print(f"Best Sequence = {best_sequence}")
            for prop in properties:
                current_value = prop._current_raw_value
                current_str = f"{current_value:.3f}" if current_value is not None else "None"
                best_value = prop._best_value  
                best_str = f"{best_value:.3f}" if best_value is not None else "None"
                print(f"{prop.__class__.__name__} targeting {prop.constraint_type.value} value of {prop.target_value}: raw value is {current_str}, scaled error is {best_str}")

    # Close progress bar
    pbar.close()
    return best_sequence

def calculate_errors(protein: sparrow.Protein, properties: List[ProteinProperty], 
                    optimizer: SequenceOptimizer = None) -> Tuple[float, Dict[str, Dict[str, Any]]]:
    """Enhanced error calculation with property-aware scaling."""
    errors = []
    directions = {}
    
    for prop in properties:
        property_name = prop.__class__.__name__
        current_value = prop.calculate_raw_value(protein)
        raw_error = prop.calculate_error(current_value)  # Use property's error calculation
        
        # Apply scaling if optimizer provided
        if optimizer:
            scaled_error = optimizer._scale_property_error(property_name, raw_error, prop, prop.weight)
        else:
            # Fallback to raw weighted error
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
            'scale_range': optimizer.property_scaling_ranges.get(property_name, 1.0) if optimizer else 1.0
        }
    
    combined_error = sum(errors)
    return combined_error, directions



def build_sequence(kmer_dict: KmerDict, target_length: int) -> str:
    """
    Build a sequence by randomly selecting k-mers from the kmer_dict until the target length is reached.

    Parameters
    ----------
    kmer_dict : KmerDict
        The KmerDict instance containing k-mer properties.
    target_length : int
        The target length of the sequence.

    Returns
    -------
    str
        The constructed sequence.
    """
    sequence = ""

    while len(sequence) < target_length:
        valid_kmers = kmer_dict.get_valid_kmers(1, target_length - len(sequence))

        if not valid_kmers:
            break

        new_kmer = random.choice(valid_kmers)
        sequence += new_kmer

    return sequence
