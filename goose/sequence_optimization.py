import random
import math
from collections import defaultdict
from typing import Any, Dict, List, Tuple
from tqdm import tqdm
import json
import logging
import pickle

import sparrow
from goose import get_data
from goose.properties import (
    FCR,
    NCPR,
    SCD,
    SHD,
    Complexity,
    ComputeIWD,
    EndToEndDistance,
    Hydrophobicity,
    Kappa,
    ProteinProperty,
    RadiusOfGyration,
    TargetAminoAcidFractions,
    epsilon_vector_diff
)


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
    # Add class-level cache dictionary
    _kmer_dict_cache: Dict[str, KmerDict] = {}

    def __init__(self, target_length: int, kmer_dict_file: str = 'amino_acids.pkl', 
                 verbose=False, gap_to_report=10, num_shuffles=0, just_shuffle=False):
        self.target_length = target_length
        self.kmer_dict_file = kmer_dict_file
        self.verbose = verbose
        self.gap_to_report = gap_to_report
        self.kmer_dict = None 
        self.fixed_ranges: List[Tuple[int, int]] = []
        self.max_iterations = 1000
        self.tolerance = 0.001
        self.window_size = 50
        self.num_shuffles = num_shuffles
        self.shuffle_interval = 25
        self.initial_sequence = None
        self.just_shuffle = just_shuffle
        self._configure_logger()
        self._load_kmer_dict()
        self.properties_dict = {}  # Dictionary to store properties for faster lookup

    def _configure_logger(self):
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
            try:
                if self.kmer_dict_file in ["kmer_properties.pickle", "amino_acids.pkl"]:
                    self.kmer_dict_file = get_data(self.kmer_dict_file)

                with open(self.kmer_dict_file, "rb") as f:
                    kmer_properties = pickle.load(f)
                self.kmer_dict = build_kmer_dict(kmer_properties)

                # Cache the kmer_dict for future use
                SequenceOptimizer._kmer_dict_cache[self.kmer_dict_file] = self.kmer_dict
                self.logger.info(f"Loaded and cached kmer dictionary from {self.kmer_dict_file}")
            except OSError as e:
                self.logger.error(f"Error loading kmer dictionary: {e}")
                raise

    @classmethod
    def load_configuration(cls, config_file: str, kmer_dict_file: str = None):
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

        # Load properties
        for prop_name, target_value, weight in config['properties']:
            try:
                # properties being optimized must be in global namespace
                property_class = globals()[prop_name]
            except KeyError:
                raise ValueError(f"property class {prop_name} not found.")

            optimizer.add_property(property_class, target_value, weight)

        optimizer.logger.info(f"Configuration loaded from {config_file}")
        return optimizer

    def add_property(self, property_class: type, *args: Any, **kwargs: Any):
        new_property = property_class(*args, **kwargs)
        property_name = property_class.__name__

        # Replace if property already exists
        if property_name in self.properties_dict:
            self.properties_dict[property_name] = new_property
            self.logger.info(f"Replaced existing property {property_name}")
        else:
            self.properties_dict[property_name] = new_property
            self.logger.info(f"Added new property {property_name}")

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
        assert len(sequence) == self.target_length, f"Initial sequence length must match target length ({self.target_length})"
        self.initial_sequence = sequence

    def run(self) -> str:
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

    def get_optimization_params(self):
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
        if self.initial_sequence is not None:
            self.logger.info("Initial Sequence: " + self.initial_sequence)
        self.logger.info(f"Optimized Sequence: {sequence}")
        protein = sparrow.Protein(sequence)
        for prop in self.properties_dict.values():
            value = prop.calculate(protein)
            self.logger.info(f"{prop.__class__.__name__}: {value:.2f} (Target: {prop.target_value:.2f})")
    
    def save_configuration(self, filename: str):
        config = {
            "target_length": self.target_length,
            "properties": [(prop.__class__.__name__, prop.target_value, prop.weight) for prop in self.properties_dict.values()],
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

    def get_defined_properties(self) -> List[Dict[str, Any]]:
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
            self.logger.info("No properties defined.")
        else:
            self.logger.info("Defined properties:")
            [self.logger.info(f"  - {prop['name']}: target = {prop['target_value']}, weight = {prop['weight']}") for prop in defined_properties]


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



def filter_candidate_kmers(sequence: str, kmer_dict: KmerDict, directions: Dict[str, Dict[str, Any]], min_k: int, max_k: int, fixed_residue_ranges: List[Tuple[int, int]]) -> Dict[int, List[str]]:
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
    min_k : int
        The minimum length of k-mers to consider.
    max_k : int
        The maximum length of k-mers to consider.
    fixed_residue_ranges : List[Tuple[int, int]]
        A list of tuples representing the start and end indices (inclusive) of the fixed residue ranges.

    Returns
    -------
    Dict[int, List[str]]
        A dictionary mapping k-mer lengths to lists of valid candidate k-mers.
    """
    current_properties = {method_name: info["current_value"] for method_name, info in directions.items()}

    # Filter candidate kmers based on their property values, directionality, and fixed residue ranges
    candidate_kmers = defaultdict(list)
    for kmer, prop_dict in kmer_dict.kmer_properties.items():
        valid_kmer = True
        kmer_range = (sequence.find(kmer), sequence.find(kmer) + len(kmer) - 1)

        # Check if the k-mer overlaps with any fixed residue range
        #if any(max(start, kmer_range[0]) <= min(end, kmer_range[1]) for start, end in fixed_residue_ranges):
        if any(start <= kmer_range[0] <= end and start <= kmer_range[1] <= end for start, end in fixed_residue_ranges):
            continue

        for method_name, info in directions.items():
            direction = info["direction"]
            if method_name in prop_dict:
                if (prop_dict[method_name] - current_properties[method_name]) * direction < 0:
                    valid_kmer = False
                    break

        if valid_kmer:
            candidate_kmers[len(kmer)].append(kmer)

    if not candidate_kmers:
        # If no candidate kmers are found in the desired direction, choose from all kmers
        # this maybe bad, but might encourage mixing
        candidate_kmers = None

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
        # Randomly select 25% of the occurrences of the k-mers to flip
        # minimum = 1
        num_to_replace = math.ceil(len(occurrences)//4+1)
        
        selected_occurrences = random.sample(occurrences, num_to_replace)
        
        # Replace the selected occurrences with the new k-mer
        new_sequence = sequence
        for selected_occurrence in selected_occurrences:
            new_sequence = new_sequence[:selected_occurrence] + new_kmer + new_sequence[selected_occurrence + kmer_len:]
    else:
        # If no suitable occurrences found, return the original sequence
        new_sequence = sequence

    return new_sequence

def shuffle_random_subset(s: str) -> str:
    '''
    randomly shuffles a subset of a sequence 
    The frequency that a substring is shuffled is inversely proportional to its length.
    This means shorter shuffles are more common than longer ones. 

    Parameters
    ----------
    s : str
        The input sequence to shuffle.

    Returns
    -------
    str
        The shuffled sequence.

    '''
    length = len(s)
    
    # Define a weighting system where shorter shuffles are more likely.
    # This ensures that shuffling small parts of the string is more frequent than shuffling the entire string.
    # Weights are inversely proportional to the length of the shuffled substring.
    weights = [1 / (i + 1) for i in range(1, length + 1)]
    
    # Randomly select the length of the substring to shuffle, weighted towards shorter lengths.
    shuffle_length = random.choices(range(1, length + 1), weights)[0]
    
    # Randomly pick a starting index for the substring to shuffle
    start_index = random.randint(0, length - shuffle_length)
    
    # Extract the part to shuffle and the rest of the string
    substring = s[start_index:start_index + shuffle_length]
    shuffled_substring = ''.join(random.sample(substring, len(substring)))
    
    # Reconstruct the string with the shuffled part
    shuffled_string = s[:start_index] + shuffled_substring + s[start_index + shuffle_length:]
    
    return shuffled_string



def mutate_sequence(sequence: str, kmer_dict: KmerDict, 
                    target_length: int, directions: Dict[str, Dict[str, Any]], 
                    properties: List[ProteinProperty], num_shuffles: int = 0, 
                    window_size: int = 10, 
                    fixed_residue_ranges: List[Tuple[int, int]] = [],
                    just_shuffle: bool = False) -> Tuple[str, float, Dict[str, Dict[str, Any]]]:
    """
    Mutate a sequence by replacing a k-mer with a new k-mer and generate shuffled sequences.

    Parameters
    ----------
    sequence : str
        The input sequence.
    kmer_dict : KmerDict
        The KmerDict instance containing k-mer properties.
    target_length : int
        The target length of the sequence.
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
 
    min_k = 1
    max_k = 5

    if not just_shuffle:
        candidate_kmers = filter_candidate_kmers(sequence, kmer_dict, directions, min_k, max_k,fixed_residue_ranges)
        kmer_to_replace = select_kmer_to_replace(sequence, min_k, max_k, fixed_residue_ranges)

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
        combined_err, directions = calculate_errors(protein, properties)

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
    if initial_sequence is not None:
        best_sequence = initial_sequence
    else:
        best_sequence = build_sequence(kmer_dict, target_length)
    
    best_error, directions = calculate_errors(sparrow.Protein(best_sequence), properties)

    # Start progress bar
    pbar = tqdm(total=max_iterations)

    # Iterate
    for i in range(max_iterations):
        current_num_shuffles = num_shuffles if i % shuffle_interval == 0 else 0
        new_sequence, new_error, directions = mutate_sequence(best_sequence, kmer_dict, 
                                                              target_length, directions, properties, 
                                                              current_num_shuffles, window_size, 
                                                              fixed_residue_ranges, 
                                                              just_shuffle=just_shuffle)

    
        if new_error < best_error:
            best_sequence = new_sequence
            best_error = new_error

        if best_error < tolerance:
            break

        # Update progress bar
        if i % gap_to_report == 0:
            pbar.update(gap_to_report)
            pbar.set_description(f"Best Error = {best_error:.2f}")

        if verbose and i % gap_to_report == 0:
            print(f"Iteration {i}: Best Error = {best_error}")
            print(f"Iteration {i}: Best Sequence = {best_sequence}")
            for prop_name in directions:
                print(f"Iteration {i}: Target {prop_name} = {directions[prop_name]['target_value']:.3f}")
                print(f"Iteration {i}: Current {prop_name} = {directions[prop_name]['current_value']:.3f}")

    # Close progress bar
    pbar.close()

    return best_sequence


def calculate_errors(protein: 'sparrow.Protein', properties: List[ProteinProperty]) -> Tuple[float, Dict[str, Dict[str, Any]]]:
    """
    Calculate the errors between the computed property values and the target values for a given protein sequence.

    Parameters
    ----------
    protein : sparrow.Protein
        The protein sequence to evaluate.
    properties : List[ProteinProperty]
        A list of protein properties to evaluate.

    Returns
    -------
    Tuple[float, Dict[str, Dict[str, Any]]]
        A tuple containing the combined error and a dictionary of direction information for each property.
    """
    errors = []
    directions = {}

    for prop in properties:
        value = prop.calculate(protein)
        direction = 1 if prop.target_value > value else -1
        direction_info = {"current_value": value, "target_value": prop.target_value, "direction": direction}
        directions[prop.__class__.__name__] = direction_info

        error = abs(value - prop.target_value) * prop.weight
        errors.append(error)

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
