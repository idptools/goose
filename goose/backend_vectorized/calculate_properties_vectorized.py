'''
Functions for FCR, NCPR, and hydropathy to calculate properties using numpy vectorization.
'''
import numpy as np
from typing import List, Union

# from AA to INT and back
AA_TO_INT = {'A': 0, 'C': 1, 'D': 2, 'E': 3, 'F': 4, 'G': 5, 'H': 6, 'I': 7, 
                'K': 8, 'L': 9, 'M': 10, 'N': 11, 'P': 12, 'Q': 13, 'R': 14, 
                'S': 15, 'T': 16, 'V': 17, 'W': 18, 'Y': 19}
INT_TO_AA = {0: 'A', 1: 'C', 2: 'D', 3: 'E', 4: 'F', 5: 'G', 6: 
             'H', 7: 'I', 8: 'K', 9: 'L', 10: 'M', 11: 'N', 12: 
             'P', 13: 'Q', 14: 'R', 15: 'S', 16: 'T', 17: 'V', 18: 'W', 19: 'Y'}
# constants using the order of amino acids as above. 
HYDROPATHY_SCALE = np.array([6.3, 7.0, 1.0, 1.0, 7.3, 4.1, 1.3, 9.0, 
                            0.6, 8.3, 6.4, 1.0, 2.9, 1.0, 0.0, 3.7, 3.8, 8.7, 3.6, 3.2])
# Constants for FCR and NCPR
FCR_SCALE = np.array([0.0, 0.0, 1.0, 1.0, 0.0, 0.0, 0.0, 0.0, 
                      1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 1.0, 0.0, 
                      0.0, 0.0, 0.0, 0.0])
NCPR_SCALE = np.array([0.0, 0.0, -1.0, -1.0, 0.0, 0.0, 0.0, 0.0, 
                      1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 1.0, 0.0, 
                      0.0, 0.0, 0.0, 0.0])

def sequences_to_matrices(sequences):
    """
    Convert a list of sequences to a list of integer matrices.
    
    Parameters:
    -----------
    sequences : list of str
        List of protein sequences
    
    Returns:
    --------
    list of numpy.ndarray
        List of sequences represented as integer matrices
    """
    seq_matrices = []
    for seq in sequences:
        seq_matrices.append(np.array([AA_TO_INT[aa] for aa in seq]))
    return seq_matrices

def matrices_to_sequences(seq_matrices):
    """
    Convert a list of integer matrices back to sequences.
    
    Parameters:
    -----------
    seq_matrices : list of numpy.ndarray
        List of sequences represented as integer matrices
    
    Returns:
    --------
    list of str
        List of protein sequences
    """
    sequences = []
    for matrix in seq_matrices:
        sequences.append(''.join([INT_TO_AA[aa] for aa in matrix]))
    return sequences


def calculate_fcr_batch(seq_matrices):
    """
    Calculate FCR for multiple sequences simultaneously (fully vectorized).

    Parameters:
    -----------
    seq_matrices : list of numpy.ndarray
        List of sequences represented as integer matrices
        
    Returns:
    --------
    numpy.ndarray
        Array of FCR values for each sequence    
    """
    stacked = np.stack(seq_matrices)  # shape: (num_seqs, seq_len)
    fcr_scores = FCR_SCALE[stacked].mean(axis=1)
    return fcr_scores


def calculate_ncpr_batch(seq_matrices):
    """
    Calculate NCPR for multiple sequences simultaneously.
    
    Parameters:
    -----------
    seq_matrices : list of numpy.ndarray
        List of sequences represented as integer matrices
        
    Returns:
    --------
    numpy.ndarray
        Array of NCPR values for each sequence
    """
    # Calculate mean NCPR for each sequence using vectorized operations
    stacked = np.stack(seq_matrices)
    ncpr_scores = NCPR_SCALE[stacked].mean(axis=1)
    return ncpr_scores


def calculate_hydropathy_batch(seq_matrices):
    """
    Calculate hydropathy scores for multiple sequences simultaneously.
    
    Parameters:
    -----------
    seq_matrices : list of numpy.ndarray
        List of sequences represented as integer matrices
        
    Returns:
    --------
    numpy.ndarray
        Array of hydropathy scores for each sequence
    """
    # Calculate mean hydropathy for each sequence using vectorized operations
    stacked = np.stack(seq_matrices)
    hydropathy_scores = HYDROPATHY_SCALE[stacked].mean(axis=1)
    return hydropathy_scores

