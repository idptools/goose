'''
For calculating properties of a single sequence that is
represented as an array of integers. 
'''

import numpy as np

def sequence_to_array(sequence: str) -> np.ndarray:
    """
    Convert a protein sequence string to a numpy array of integers.
    
    Parameters:
    -----------
    sequence : str
        The input protein sequence.
        
    Returns:
    --------
    np.ndarray
        Numpy array representation of the sequence.
    """
    AA_TO_INT = {
    'A': 0, 'C': 1, 'D': 2, 'E': 3, 'F': 4,
    'G': 5, 'H': 6, 'I': 7, 'K': 8, 'L': 9,
    'M': 10, 'N': 11, 'P': 12, 'Q': 13, 'R': 14,
    'S': 15, 'T': 16, 'V': 17, 'W': 18, 'Y': 19
    }
    
    # Convert sequence to numpy array of integers
    return np.array([AA_TO_INT[aa] for aa in sequence if aa in AA_TO_INT], dtype=int)

def array_to_sequence(seq_array: np.ndarray) -> str:
    """
    Convert a numpy array of integers back to a protein sequence string.
    
    Parameters:
    -----------
    seq_array : np.ndarray
        Numpy array representation of the sequence.
        
    Returns:
    --------
    str
        Protein sequence as a string.
    """
    INT_TO_AA = {0: 'A', 1: 'C', 2: 'D', 3: 'E', 4: 'F', 5: 'G', 
                 6: 'H', 7: 'I', 8: 'K', 9: 'L', 10: 'M', 11: 'N', 
                 12: 'P', 13: 'Q', 14: 'R', 15: 'S', 16: 'T', 
                 17: 'V', 18: 'W', 19: 'Y'}
    # Create a mapping for the integer values to amino acids
    return ''.join([INT_TO_AA[aa] for aa in seq_array])


def calculate_hydropathy_single_sequence(seq_array: np.ndarray) -> float:
    """
    Calculate the hydropathy of a single sequence represented as an array of integers.
    Parameters:
    -----------
    seq_array : np.ndarray
        Array of integers representing the sequence.
    Returns:
    --------
    float
        Hydropathy value of the sequence.
    """
    hydropathy_index = np.array([6.3, 7.0, 1.0, 1.0, 7.3, 4.1, 1.3, 9.0, 
                            0.6, 8.3, 6.4, 1.0, 2.9, 1.0, 0.0, 3.7, 3.8, 8.7, 3.6, 3.2])
    return np.mean(hydropathy_index[seq_array])


def calculate_ncpr_single_sequence(seq_array: np.ndarray) -> float:
    """
    Calculate NCPR for multiple sequences simultaneously.
    
    Parameters:
    -----------
    seq_array : np.ndarray
        Array of sequences represented as integer matrices

    Returns:
    --------
    float
        NCPR value for the sequence
    """
    # Positive: K(8), R(14) = +1, Negative: D(2), E(3) = -1
    positive_mask = (seq_array == 8) | (seq_array == 14)
    negative_mask = (seq_array == 2) | (seq_array == 3)

    net_charge = np.sum(positive_mask) - np.sum(negative_mask)
    return net_charge / seq_array.shape[0]

def calculate_fcr_single_sequence(seq_array: np.ndarray) -> float:
    """
    Calculate FCR for multiple sequences simultaneously.
    
    Parameters:
    -----------
    seq_array : np.ndarray
        Array of sequences represented as integer matrices

    Returns:
    --------
    float
        FCR value for the sequence
    """
    # Charged residues: D(2), E(3), K(8), R(14)
    charged_mask = (
        (seq_array == 2) |  # D
        (seq_array == 3) |  # E
        (seq_array == 8) |  # K
        (seq_array == 14)   # R
    )
    
    return np.mean(charged_mask)

