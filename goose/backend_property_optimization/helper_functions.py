'''
takes in a sequence and makes a specific number of shuffles and returns
the shuffled sequences in a matrix format
'''

import numpy as np
from goose.data import aa_list_probabilities as aa_probs

def shuffle_sequence_return_matrix(seq_array, num_shuffles=1000):
    """
    Shuffle a sequence multiple times and return the results in a matrix format.
    
    Args:
        seq_array (str): The input sequence to shuffle, represented as a np array of ints
        with shape (1, n) where n is the length of the sequence.
        
        num_shuffles (int): The number of times to shuffle the sequence.
        
    Returns:
        np.ndarray: A matrix where each row is a shuffled version of the input sequence.
    """
    seq_length = seq_array.shape[1]
    shuffled_matrix = np.zeros((num_shuffles, seq_length), dtype=seq_array.dtype)
    
    for i in range(num_shuffles):
        # only shuffle if i > 0 to avoid shuffling the original sequence
        if i == 0:
            shuffled_seq = seq_array[0]
        else:
            shuffled_seq = np.random.permutation(seq_array[0])
        shuffled_matrix[i] = shuffled_seq
    
    return shuffled_matrix


def generate_random_sequence(length, probabilities=None):
    """
    Generate a random sequence of a given length based on specified amino acid probabilities.
    
    Args:
        length (int): The length of the sequence to generate.
        probabilities (dict, optional): A dictionary mapping amino acids to their probabilities.
                                    If None, uses default probabilities.
                                    
    Returns:
        str: A random sequence of amino acids.
    """
    if probabilities is None:
        probabilities = aa_probs.IDRProbs
    amino_acids = list(probabilities.keys())
    probabilities = np.array(list(probabilities.values()))
    
    # Normalize probabilities to ensure they sum to 1
    probabilities /= probabilities.sum()
    
    # Generate random sequence
    sequence = np.random.choice(amino_acids, size=length, p=probabilities)
    
    return ''.join(sequence)
