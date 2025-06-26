'''
Functionality to create a sequence from the class of amino acids. 
'''
from goose.backend_sequence_generation.seq_by_fractions_vectorized import FractionBasedSequenceGenerator
from typing import Dict, List, Optional, Union
import numpy as np
from goose.data.defined_aa_classes import aa_classes, aa_classes_by_aa, aa_class_indices

def get_individual_amino_acid_fractions(
        length: int,
        aromatic_fraction: float = 0.0,
        aliphatic_fraction: float = 0.0,
        polar_fraction: float = 0.0,
        positive_fraction: float = 0.0,
        negative_fraction: float = 0.0,
        glycine_fraction: float = 0.0,
        proline_fraction: float = 0.0,
        cysteine_fraction: float = 0.0,
        histidine_fraction: float = 0.0) -> Dict[str, float]:
    '''
    Get a dictionary of individual amino acid fractions based on class fractions.
    This takes in length and then the fractions of amino acids by class and then
    returns a dictionary of the individual amino acid fractions for each amino acid
    where the fractions are randomly generated for each amino acid in each class such
    that the sum of the fractions for each class is equal to the class fraction.
    
    IMPORTANT: This function ensures that each amino acid fraction, when multiplied
    by the sequence length, results in an integer number of amino acids. This is
    achieved by first converting class fractions to integer counts, randomly
    distributing these counts within each class, then converting back to fractions.
    
    Parameters:
    -----------
    length : int
        Length of the sequence to be generated.
    aromatic_fraction : float
        Fraction of aromatic amino acids (F, W, Y).
    aliphatic_fraction : float
        Fraction of aliphatic amino acids (A, I, L, M, V).
    polar_fraction : float      
        Fraction of polar amino acids (S, T, N, Q).
    positive_fraction : float
        Fraction of positively charged amino acids (K, R).
    negative_fraction : float
        Fraction of negatively charged amino acids (D, E).
    glycine_fraction : float
        Fraction of glycine (G).
    proline_fraction : float
        Fraction of proline (P).
    cysteine_fraction : float   
        Fraction of cysteine (C).
    histidine_fraction : float
        Fraction of histidine (H).
    Returns:
    --------
    Dict[str, float]
        A dictionary containing the individual amino acid fractions where each
        fraction * length yields an integer number of amino acids.    
    '''
    
    # Initialize the result dictionary with all amino acids set to 0
    aa_fractions = {aa: 0.0 for aa in 'ACDEFGHIKLMNPQRSTVWY'}
    
    # Define the class fractions and their corresponding amino acids
    class_fractions = {
        'aromatic': (aromatic_fraction, aa_classes['aromatic']),
        'aliphatic': (aliphatic_fraction, aa_classes['aliphatic']),
        'polar': (polar_fraction, aa_classes['polar']),
        'positive': (positive_fraction, aa_classes['positive']),
        'negative': (negative_fraction, aa_classes['negative']),
        'glycine': (glycine_fraction, aa_classes['glycine']),
        'proline': (proline_fraction, aa_classes['proline']),
        'cysteine': (cysteine_fraction, aa_classes['cysteine']),
        'histidine': (histidine_fraction, aa_classes['histidine'])
    }
    
    # Process each class - work with integer counts to ensure exact counts
    for class_name, (total_fraction, aa_list) in class_fractions.items():
        if total_fraction > 0.0 and len(aa_list) > 0:
            # Convert fraction to integer count for this class
            total_count = int(round(total_fraction * length))
            
            if len(aa_list) == 1:
                # Single amino acid in class, assign all counts to that amino acid
                aa_fractions[aa_list[0]] = total_count / length
            else:
                # Multiple amino acids in class, randomly distribute the integer counts
                if total_count > 0:
                    # Generate random integer distribution that sums to total_count
                    counts = np.zeros(len(aa_list), dtype=int)
                    
                    # Distribute counts randomly
                    for i in range(total_count):
                        # Randomly select an amino acid to increment
                        selected_aa_idx = np.random.randint(0, len(aa_list))
                        counts[selected_aa_idx] += 1
                    
                    # Convert counts back to fractions
                    for i, aa in enumerate(aa_list):
                        aa_fractions[aa] = counts[i] / length
    
    return aa_fractions

def fill_in_remaining_classes(
        aa_class_fractions: Dict[str, float],
        length: int) -> Dict[str, float]:
    '''
    Fill in the remaining classes of amino acids to make sure that the fractions
    sum to 1.0 but also that the fractions are compatible with the length of sequence such that
    sequence * fraction is an integer.
    Should also ensure that the sum of each fraction times length is equal to the total length.
    The purpose of this function is to fill in fractions that are not specified by the user. 

    Parameters:
    -----------
    aa_class_fractions : Dict[str, float]
        Dictionary of fractions for each class of amino acids. 
        Keys should be: 'aromatic', 'aliphatic', 'polar', 'positive', 'negative', 
        'glycine', 'proline', 'cysteine', 'histidine'
    length : int
        Length of the sequence to be generated.
    Returns:
    --------
    Dict[str, float]
        Updated dictionary of class fractions with remaining classes filled in to sum to 1.0.
    '''
    
    # Make a copy to avoid modifying the original
    filled_fractions = aa_class_fractions.copy()
    
    # Calculate current sum of specified fractions
    current_sum = sum(filled_fractions.values())
    
    # If sum is already 1.0 (within tolerance), return as is
    if abs(current_sum - 1.0) < 1e-10:
        return filled_fractions
    
    # If sum exceeds 1.0, normalize all fractions proportionally
    if current_sum > 1.0:
        for class_name in filled_fractions:
            filled_fractions[class_name] = filled_fractions[class_name] / current_sum
        return filled_fractions
    
    # Calculate remaining fraction to distribute
    remaining_fraction = 1.0 - current_sum
    remaining_count = int(round(remaining_fraction * length))
    
    # Identify classes that are not specified (have 0 fraction)
    all_classes = ['aromatic', 'aliphatic', 'polar', 'positive', 'negative', 
                   'glycine', 'proline', 'cysteine', 'histidine']
    
    unspecified_classes = [cls for cls in all_classes 
                          if filled_fractions.get(cls, 0.0) == 0.0]
    
    # If no unspecified classes, distribute remaining among existing classes proportionally
    if not unspecified_classes:
        if current_sum > 0:
            for class_name in filled_fractions:
                proportion = filled_fractions[class_name] / current_sum
                additional_count = int(round(proportion * remaining_count))
                filled_fractions[class_name] += additional_count / length
    else:
        # Distribute remaining fraction among unspecified classes
        if remaining_count > 0:
            # Use simple random distribution among unspecified classes
            counts_per_class = np.zeros(len(unspecified_classes), dtype=int)
            
            # Randomly distribute the remaining counts
            for i in range(remaining_count):
                selected_class_idx = np.random.randint(0, len(unspecified_classes))
                counts_per_class[selected_class_idx] += 1
            
            # Convert counts back to fractions
            for i, class_name in enumerate(unspecified_classes):
                filled_fractions[class_name] = counts_per_class[i] / length
    
    return filled_fractions

def create_sequence_by_class(
        length: int,
        aromatic_fraction: float = 0.0,
        aliphatic_fraction: float = 0.0,
        polar_fraction: float = 0.0,
        positive_fraction: float = 0.0,
        negative_fraction: float = 0.0,
        glycine_fraction: float = 0.0,
        proline_fraction: float = 0.0,
        cysteine_fraction: float = 0.0,
        histidine_fraction: float = 0.0,
        num_sequences: int = 1,
        convert_to_amino_acids: bool = True) -> Union[str, List[str]]:
    '''
    Create sequence(s) based on the specified amino acid class fractions.
    
    This function first ensures that class fractions sum to 1.0 by filling in 
    unspecified classes, then converts class fractions to individual amino acid
    fractions, and finally generates sequences using the FractionBasedSequenceGenerator.
    
    Parameters:
    -----------
    length : int
        Length of the sequence to be generated.
    aromatic_fraction : float
        Fraction of aromatic amino acids (F, W, Y).
    aliphatic_fraction : float
        Fraction of aliphatic amino acids (A, I, L, M, V).
    polar_fraction : float
        Fraction of polar amino acids (S, T, N, Q).
    positive_fraction : float
        Fraction of positively charged amino acids (K, R).
    negative_fraction : float       
        Fraction of negatively charged amino acids (D, E).
    glycine_fraction : float
        Fraction of glycine (G).
    proline_fraction : float
        Fraction of proline (P).
    cysteine_fraction : float
        Fraction of cysteine (C).
    histidine_fraction : float
        Fraction of histidine (H).
    num_sequences : int
        Number of sequences to generate (default: 1).
    convert_to_amino_acids : bool
        If True, convert the generated sequences to amino acid strings.
    Returns:
    --------
    Union[str, List[str]]
        If num_sequences == 1, returns a single sequence string.
        If num_sequences > 1, returns a list of sequence strings.
    '''
    
    # Validate inputs
    if length <= 0:
        raise ValueError(f"Length must be positive, got {length}")
    if num_sequences <= 0:
        raise ValueError(f"num_sequences must be positive, got {num_sequences}")
    
    # Validate that all fractions are non-negative
    class_fractions = {
        'aromatic': aromatic_fraction,
        'aliphatic': aliphatic_fraction,
        'polar': polar_fraction,
        'positive': positive_fraction,
        'negative': negative_fraction,
        'glycine': glycine_fraction,
        'proline': proline_fraction,
        'cysteine': cysteine_fraction,
        'histidine': histidine_fraction
    }
    
    for class_name, fraction in class_fractions.items():
        if fraction < 0:
            raise ValueError(f"{class_name}_fraction must be non-negative, got {fraction}")
    
    # Check if sum exceeds 1.0
    total_specified = round(sum(class_fractions.values()),8)
    if total_specified > 1.0:
        raise ValueError(f"Sum of class fractions ({total_specified:.4f}) exceeds 1.0")
    
    # Step 1: Fill in remaining classes to ensure fractions sum to 1.0
    filled_class_fractions = fill_in_remaining_classes(class_fractions, length)
    
    # Step 2: Convert class fractions to individual amino acid fractions
    aa_fractions = get_individual_amino_acid_fractions(
        length=length,
        aromatic_fraction=filled_class_fractions['aromatic'],
        aliphatic_fraction=filled_class_fractions['aliphatic'],
        polar_fraction=filled_class_fractions['polar'],
        positive_fraction=filled_class_fractions['positive'],
        negative_fraction=filled_class_fractions['negative'],
        glycine_fraction=filled_class_fractions['glycine'],
        proline_fraction=filled_class_fractions['proline'],
        cysteine_fraction=filled_class_fractions['cysteine'],
        histidine_fraction=filled_class_fractions['histidine']
    )
    
    # Step 3: Generate sequences using FractionBasedSequenceGenerator
    # Normalize fractions to handle potential floating point precision issues
    aa_fractions_sum = sum(aa_fractions.values())
    if abs(aa_fractions_sum - 1.0) > 1e-12:
        for aa in aa_fractions:
            if aa_fractions[aa] > 0:
                aa_fractions[aa] = aa_fractions[aa] / aa_fractions_sum
    
    generator = FractionBasedSequenceGenerator(
        length=length,
        fractions=aa_fractions,
        randomize_unspecified=False,  # We've already filled everything
        normalize=False,  # Don't normalize again since we just did it
        strict=True  # Use strict mode for exact counts
        )
    
    sequences = generator.generate_sequences(num_sequences=num_sequences,
                                             convert_to_amino_acids=convert_to_amino_acids)

    # Return single string if only one sequence requested, otherwise return list
    if num_sequences == 1:
        return sequences[0]
    else:
        return sequences