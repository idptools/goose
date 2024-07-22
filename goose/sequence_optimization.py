import json
import logging
import math
import pickle
import random
import types
from collections import defaultdict
from functools import lru_cache
from typing import Any, Dict, List, Tuple

import numpy as np
import sparrow
from IPython import embed

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
)

from tqdm import tqdm, trange

from goose import get_data


class SequenceOptimizer:
    def __init__(self, target_length: int, kmer_dict_file: str = None, verbose=False):
        self.target_length = target_length
        self.kmer_dict_file = kmer_dict_file
        self.verbose = verbose
        self.kmer_dict = None 
        self.properties: List[ProteinProperty] = []
        self.fixed_ranges: List[Tuple[int, int]] = []
        self.max_iterations = 1000
        self.tolerance = 0.001
        self.window_size = 50
        self.num_shuffles = 5000
        self.shuffle_interval = 25
        self.initial_sequence = None

        # Set up logging
        logging.basicConfig(level=logging.INFO, format='%(asctime)s - %(levelname)s - %(message)s')
        self.logger = logging.getLogger(__name__)

        # Load kmer_dict if file is provided
        if self.kmer_dict_file:
            self._load_kmer_dict()
        else:
            self.kmer_dict_file = get_data("kmer_properties.pickle")
            self._load_kmer_dict()

    def _load_kmer_dict(self):
        try:
            with open(self.kmer_dict_file, "rb") as f:
                kmer_properties = pickle.load(f)
            self.kmer_dict = build_kmer_dict(kmer_properties)
            self.logger.info(f"Loaded kmer dictionary from {self.kmer_dict_file}")
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
        optimizer.set_initial_sequence(config['initial_sequence'])

        # Load properties
        for prop_name, target_value, weight in config['properties']:
            # properties being optimized must be in global namespace
            property_class = globals()[prop_name]
            optimizer.add_property(property_class, target_value, weight)

        optimizer.logger.info(f"Configuration loaded from {config_file}")
        return optimizer

    def add_property(self, property_class: type, *args: Any, **kwargs: Any):
        new_property = property_class(*args, **kwargs)
        
        # Check if the property already exists
        for i, prop in enumerate(self.properties):
            if isinstance(prop, property_class):
                # If it exists, replace it
                self.properties[i] = new_property
                self.logger.info(f"Replaced existing property {property_class.__name__}")
                return
        
        # If it doesn't exist, append it
        self.properties.append(new_property)
        self.logger.info(f"Added new property {property_class.__name__}")

    def set_fixed_ranges(self, ranges: List[Tuple[int, int]]):
        self.fixed_ranges = ranges

    def set_optimization_params(self, max_iterations: int = None, tolerance: float = None, 
                                window_size: int = None, num_shuffles: int = None, 
                                shuffle_interval: int = None):
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

    def set_initial_sequence(self, sequence: str):
        assert len(sequence) == self.target_length, f"Initial sequence length must match target length ({self.target_length})"
        self.initial_sequence = sequence

    def run(self) -> str:
        self.logger.info("Starting sequence optimization")
        optimized_sequence = optimize_sequence(
            self.kmer_dict,
            self.target_length,
            self.properties,
            max_iterations=self.max_iterations,
            tolerance=self.tolerance,
            verbose=self.verbose,
            window_size=self.window_size,
            num_shuffles=self.num_shuffles,
            shuffle_interval=self.shuffle_interval,
            fixed_residue_ranges=self.fixed_ranges,
            initial_sequence=self.initial_sequence
        )
        self.logger.info("Sequence optimization completed")
        self.log_results(optimized_sequence)
        return optimized_sequence

    def log_results(self, sequence: str):
        if self.initial_sequence is not None:
            self.logger.info("Initial Sequence: " + self.initial_sequence)
        self.logger.info(f"Optimized Sequence: {sequence}")
        protein = sparrow.Protein(sequence)
        for prop in self.properties:
            value = prop.calculate(protein)
            self.logger.info(f"{prop.__class__.__name__}: {value:.2f} (Target: {prop.target_value:.2f})")

    def save_configuration(self, filename: str):
        config = {
            "target_length": self.target_length,
            "properties": [(prop.__class__.__name__, prop.target_value, prop.weight) for prop in self.properties],
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
            for prop in self.properties
        ]
    
    @property
    def propeties(self):
        defined_properties = self.get_defined_properties()
        if not defined_properties:
            print("No properties defined.")
        else:
            print("Defined properties:")
            for prop in defined_properties:
                print(f"  - {prop['name']}: target = {prop['target_value']}, weight = {prop['weight']}")

class KmerDict:
    """
    A dictionary-based implementation for storing k-mers and their associated properties.
    """
    def __init__(self):
        self.kmer_properties = {}  # Dictionary to store properties for each k-mer

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
        valid_kmers = []
        for kmer in self.kmer_properties:
            if k <= len(kmer) <= max_length:
                valid_kmers.append(kmer)
        return valid_kmers

def build_kmer_dict(kmer_properties: Dict[str, float]) -> KmerDict:
    """
    Build a KmerDict from a dictionary of k-mer properties.

    Parameters
    ----------
    kmer_properties : Dict[str, float]
        A dictionary mapping k-mers to their associated property values.

    Returns
    -------
    KmerDict
        The constructed KmerDict instance.
    """
    kmer_dict = KmerDict()
    for kmer, prop_value in kmer_properties.items():
        kmer_dict.insert(kmer, prop_value)
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
    protein = sparrow.Protein(sequence)
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


def mutate_sequence(sequence: str, kmer_dict: KmerDict, target_length: int, directions: Dict[str, Dict[str, Any]], properties: List[ProteinProperty], num_shuffles: int = 0, window_size: int = 10, fixed_residue_ranges: List[Tuple[int, int]] = []) -> Tuple[str, float, Dict[str, Dict[str, Any]]]:
    """
    Mutate a sequence by replacing a k-mer with a new k-mer and generate shuffled sequences.

    Parameters
    ----------
    ...
    fixed_residue_ranges : List[Tuple[int, int]], optional
        A list of tuples representing the start and end indices (inclusive) of the fixed residue ranges, by default [].

    Returns
    -------
    Tuple[str, float, Dict[str, Dict[str, Any]]]
        A tuple containing the new sequence, the combined error, and the updated directions dictionary.
    """
    min_k = 1
    max_k = 5

    candidate_kmers = filter_candidate_kmers(sequence, kmer_dict, directions, min_k, max_k,fixed_residue_ranges)

    kmer_to_replace = select_kmer_to_replace(sequence, min_k, max_k, fixed_residue_ranges)

    if kmer_to_replace is None:
        # If no k-mer is available for replacement, keep the sequence unchanged
        new_sequence = sequence
    else:
        valid_kmers = candidate_kmers[len(kmer_to_replace)]

        if valid_kmers:
            # Replace a k-mer in the sequence with a new k-mer
            new_kmer = random.choice(valid_kmers)
            new_sequence = replace_kmer(sequence, kmer_to_replace, new_kmer, fixed_residue_ranges)     
        else:
            # If no valid k-mers are found, keep the sequence unchanged
            new_sequence = sequence

    if num_shuffles > 0:
        sequences = [new_sequence] + get_global_and_window_shuffles(new_sequence, num_shuffles, window_size)
    else:
        sequences = [new_sequence]

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

def combined_objective(sequence: str, method_args: List[Tuple[str, Tuple[Any]]], target_values: List[float],
                       weights: List[float] = None, num_shuffles: int = 0, window_size: int = 10) -> Tuple[float, str, Dict[str, Dict[str, Any]]]:
    """
    Compute the combined objective function for a sequence by mutating it and generating shuffled sequences.

    Parameters
    ----------
    sequence : str
        The current sequence.
    method_args : List[Tuple[str, Tuple[Any]]]
        A list of tuples containing the method name and arguments for each property.
    target_values : List[float]
        A list of target values for each property.
    weights : List[float], optional
        A list of weights for each property, by default None (equal weights).
    num_shuffles : int, optional
        The number of shuffled sequences to generate, by default 0.
    window_size : int, optional
        The size of the window to be shuffled, by default 10.

    Returns
    -------
    Tuple[float, str, Dict[str, Dict[str, Any]]]
        A tuple containing the minimum combined error, the best sequence, and the updated directions dictionary.
    """
    if weights is None:
        weights = [1.0 / len(method_args)] * len(method_args)
    else:
        if len(weights) != len(method_args):
            raise ValueError(f"Length of weights should match the number of target methods. Recieved: {len(weights)} weights and {len(method_args)} args")

    protein = sparrow.Protein(sequence)
    combined_err, directions = calculate_errors(protein, method_args, target_values, weights)
    best_sequence, min_combined_error, best_directions = mutate_sequence(sequence, kmer_dict, target_length, directions, num_shuffles, window_size)

    return min_combined_error, best_sequence, best_directions

def optimize_sequence(kmer_dict: KmerDict, target_length: int, properties: List[ProteinProperty],
                      max_iterations: int = 100, tolerance: float = 1e-2, verbose=False,
                      window_size: int = 10, num_shuffles=100, shuffle_interval: int = 50,
                      fixed_residue_ranges: List[Tuple[int, int]] = [], initial_sequence=None) -> str:
    """
    Optimize a protein sequence to minimize the combined error between its computed properties and the target property values.

    Parameters:

    fixed_residue_ranges : List[Tuple[int, int]], optional
        A list of tuples representing the start and end indices (inclusive) of the fixed residue ranges, by default [].
    """
    if initial_sequence is not None:
        best_sequence = initial_sequence
        print(best_sequence)
    else:
        best_sequence = build_sequence(kmer_dict, target_length)
    
    best_error, directions = calculate_errors(sparrow.Protein(best_sequence), properties)

    for i in trange(max_iterations):
        num_shuffles = num_shuffles if i % shuffle_interval == 0 else 0
        new_sequence, new_error, directions = mutate_sequence(best_sequence, kmer_dict, target_length, directions, properties, num_shuffles, window_size, fixed_residue_ranges)

        if new_error < best_error:
            best_sequence = new_sequence
            best_error = new_error

        if best_error < tolerance:
            break

        if verbose:
            if i % 100 == 0:
                print(f"Iteration {i}: Best Error = {best_error}")
                print(f"Iteration {i}: Best Sequence = {best_sequence}")
                for prop_name in directions:
                    print(f"Iteration {i}: Target {prop_name}= {directions[prop_name]['target_value']:.3f}")
                    print(f"Iteration {i}: current {prop_name}= {directions[prop_name]['current_value']:.3f}")

    return best_sequence

def calculate_errors(protein: 'sparrow.Protein', properties: List[ProteinProperty]) -> Tuple[float, Dict[str, Dict[str, Any]]]:
    """
    Calculate the errors between the computed property values and the target values for a given protein sequence.

    Parameters:
    ...
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

        if valid_kmers:
            new_kmer = random.choice(valid_kmers)
            sequence += new_kmer
        else:
            break

    return sequence


if __name__ == "__main__":
    optimizer = SequenceOptimizer(100)
    optimizer.add_property(Hydrophobicity, 0.5, 1.0)
    optimizer.add_property(ComputeIWD,("YFL",), 0.5, 1)
    optimizer.run()
   