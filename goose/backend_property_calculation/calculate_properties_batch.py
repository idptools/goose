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

# Constants for FCR and NCPR
FCR_SCALE = np.array([0.0, 0.0, 1.0, 1.0, 0.0, 0.0, 0.0, 0.0, 
                      1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 1.0, 0.0, 
                      0.0, 0.0, 0.0, 0.0])
NCPR_SCALE = np.array([0.0, 0.0, -1.0, -1.0, 0.0, 0.0, 0.0, 0.0, 
                      1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 1.0, 0.0, 
                      0.0, 0.0, 0.0, 0.0])

def sequences_to_matrices(sequences):
    """
    Convert a list of sequences to integer matrices (optimized version).
    
    This function is highly optimized for speed:
    - Uses numpy's vectorized operations
    - Pre-allocates output array 
    - Uses direct character-to-integer mapping
    - Minimizes Python loops
    
    Parameters:
    -----------
    sequences : list of str
        List of protein sequences
    
    Returns:
    --------
    numpy.ndarray
        2D array of shape (n_sequences, max_length) with integer amino acid codes
    """
    if isinstance(sequences, np.ndarray):
        sequences = sequences.tolist()

    if sequences==None or len(sequences) == 0:
        return [np.array([])]
    
    # Get dimensions
    n_sequences = len(sequences)
    max_length = max(len(seq) for seq in sequences)
    
    # Pre-allocate output matrix (default to 0, which corresponds to 'A')
    seq_matrix = np.zeros((n_sequences, max_length), dtype=np.int8)
    
    # Create a vectorized lookup array for faster conversion
    # This avoids dictionary lookups in the inner loop
    aa_lookup = np.zeros(256, dtype=np.int8)  # ASCII table size
    for aa, idx in AA_TO_INT.items():
        aa_lookup[ord(aa)] = idx
    
    # Fill the matrix using vectorized operations where possible
    for i, seq in enumerate(sequences):
        if seq:  # Skip empty sequences
            # Convert sequence to numpy array of ASCII codes, then lookup
            ascii_codes = np.frombuffer(seq.encode('ascii'), dtype=np.uint8)
            seq_matrix[i, :len(seq)] = aa_lookup[ascii_codes]
    
    return seq_matrix

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


def calculate_fcr_batch(seq_matrix):
    """
    Calculate FCR for multiple sequences simultaneously (fully vectorized).

    Parameters:
    -----------
    seq_matrix : list of numpy.ndarray
        List of sequences represented as integer matrices
        
    Returns:
    --------
    numpy.ndarray
        Array of FCR values for each sequence    
    """
    # Charged residues: D(2), E(3), K(8), R(14)
    charged_mask = (
        (seq_matrix == 2) |  # D
        (seq_matrix == 3) |  # E
        (seq_matrix == 8) |  # K
        (seq_matrix == 14)   # R
    )
    
    # No need for valid_mask division - just use sequence length
    return np.mean(charged_mask, axis=1)


def calculate_ncpr_batch(seq_matrix):
    """
    Calculate NCPR for multiple sequences simultaneously.
    
    Parameters:
    -----------
    seq_matrix : list of numpy.ndarray
        List of sequences represented as integer matrices
        
    Returns:
    --------
    numpy.ndarray
        Array of NCPR values for each sequence
    """
    # Positive: K(8), R(14) = +1, Negative: D(2), E(3) = -1
    positive_mask = (seq_matrix == 8) | (seq_matrix == 14)
    negative_mask = (seq_matrix == 2) | (seq_matrix == 3)
    
    net_charge = np.sum(positive_mask, axis=1) - np.sum(negative_mask, axis=1)
    return net_charge / seq_matrix.shape[1]


def calculate_hydropathy_batch(seq_matrix):
    """
    Calculate hydropathy scores for multiple sequences simultaneously.
    
    Parameters:
    -----------
    seq_matrix : list of numpy.ndarray
        List of sequences represented as integer matrices
        
    Returns:
    --------
    numpy.ndarray
        Array of hydropathy scores for each sequence
    """
    hydropathy_index = np.array([6.3, 7.0, 1.0, 1.0, 7.3, 4.1, 1.3, 9.0, 
                            0.6, 8.3, 6.4, 1.0, 2.9, 1.0, 0.0, 3.7, 3.8, 8.7, 3.6, 3.2])
    
    # Direct indexing without masking - much faster
    hydropathy_values = hydropathy_index[seq_matrix.astype(int)]
    return np.mean(hydropathy_values, axis=1)


def calculate_aa_fracts_batch(seq_matrix, seq_lengths):
    '''
    Calculate the amino acid fractions for each sequence in the batch.
    Parameters
    ----------
    seq_matrix : np.ndarray
        Matrix of sequences, where each row is a sequence and each column is a position in the sequence.
    seq_lengths : np.ndarray
        Array of lengths of each sequence in the batch.
    Returns
    -------
    list
        List of dictionaries, where each dictionary contains amino acid letters as keys and their fractions as values for each sequence.
    '''
    # Flatten the array and remove NaN values to get only actual amino acids
    all_amino_acids = np.array([0,1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19])
    
    # Mapping from index to amino acid letter
    aa_letters = ['A', 'C', 'D', 'E', 'F', 'G', 'H', 'I', 'K', 'L', 'M', 'N', 'P', 'Q', 'R', 'S', 'T', 'V', 'W', 'Y']

    # Much faster bincount without NaN handling
    seq_length = seq_matrix.shape[1]
    aa_counts = np.apply_along_axis(lambda x: np.bincount(x.astype(int), minlength=20), axis=1, arr=seq_matrix)
    aa_fracts = aa_counts / seq_length
    
    return [dict(zip(aa_letters, fracts)) for fracts in aa_fracts]

def calculate_frac_aro_batch(seq_matrix):
    '''
    Optimized aromatic fraction calculation for fixed-length sequences.
    '''
    aromatic_mask = (seq_matrix == 4) | (seq_matrix == 18) | (seq_matrix == 19)  # F, W, Y
    return np.mean(aromatic_mask, axis=1)

def calculate_frac_polar_batch(seq_matrix):
    '''
    Optimized polar fraction calculation for fixed-length sequences.
    '''
    polar_mask = (
        (seq_matrix == 5) |   # G
        (seq_matrix == 6) |   # H
        (seq_matrix == 11) |  # N
        (seq_matrix == 13) |  # Q
        (seq_matrix == 15) |  # S
        (seq_matrix == 16)    # T
    )
    return np.mean(polar_mask, axis=1)

def calculate_frac_positive_batch(seq_matrix):
    '''
    Optimized positive fraction calculation for fixed-length sequences.
    '''
    positive_mask = (seq_matrix == 8) | (seq_matrix == 14)  # K, R
    return np.mean(positive_mask, axis=1)

def calculate_frac_negative_batch(seq_matrix):
    '''
    Optimized negative fraction calculation for fixed-length sequences.
    '''
    negative_mask = (seq_matrix == 2) | (seq_matrix == 3)  # D, E
    return np.mean(negative_mask, axis=1)

def calculate_frac_aliphatic_batch(seq_matrix):
    '''
    Optimized aliphatic fraction calculation for fixed-length sequences.
    '''
    aliphatic_mask = (
        (seq_matrix == 0) |   # A
        (seq_matrix == 7) |   # I
        (seq_matrix == 9) |   # L
        (seq_matrix == 10) |  # M
        (seq_matrix == 17)    # V
    )
    return np.mean(aliphatic_mask, axis=1)

def calculate_frac_proline_batch(seq_matrix):
    '''
    Optimized proline fraction calculation for fixed-length sequences.
    '''
    proline_mask = (seq_matrix == 12)  # P
    return np.mean(proline_mask, axis=1)

def calculate_complexity_batch(seq_matrix):
    '''
    Optimized complexity calculation for fixed-length sequences.
    
    Parameters
    ----------
    seq_matrix : np.ndarray
        Matrix of sequences with no NaN values.
    
    Returns
    -------
    np.ndarray
        Array of complexity values for each sequence.
    '''
    seq_length = seq_matrix.shape[1]
    
    # Faster bincount without NaN handling
    aa_counts = np.apply_along_axis(
        lambda x: np.bincount(x.astype(int), minlength=20), 
        axis=1, 
        arr=seq_matrix
    )
    
    # Calculate fractions
    aa_fractions = aa_counts / seq_length
    
    # Vectorized complexity calculation
    with np.errstate(divide='ignore', invalid='ignore'):
        log_fractions = np.where(aa_fractions > 0, np.log(aa_fractions) / np.log(20), 0)
    complexities = -np.sum(aa_fractions * log_fractions, axis=1)
    
    return complexities
