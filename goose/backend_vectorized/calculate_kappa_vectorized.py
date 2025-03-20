'''
Numpy vectorized implementation of kappa calculations for multiple sequences
'''

import numpy as np


def vectorized_ternarize(sequences, group1, group2):
    """
    Vectorized function to convert multiple amino acid sequences into ternary arrays.
    
    Parameters
    ----------
    sequences : list or np.ndarray
        List of amino acid sequences (all sequences must be same length)
        
    group1 : list
        List of residues considered as group 1 (assigned value 1)
        
    group2 : list
        List of residues considered as group 2 (assigned value -1)
        
    Returns
    -------
    np.ndarray
        Array of shape (n_sequences, sequence_length) with ternary values
    """
    seq_array = np.array([list(seq) for seq in sequences])
    n_sequences, seq_length = seq_array.shape
    
    # Initialize with zeros
    ternary_array = np.zeros((n_sequences, seq_length), dtype=np.int8)
    
    # Create masks for group1 and group2
    for residue in group1:
        ternary_array[seq_array == residue] = 1
    
    # If group2 is empty, consider all non-group1 residues as group2
    for residue in group2:
        ternary_array[seq_array == residue] = -1
        
    return ternary_array


def vectorized_global_asymmetry(ternary_arrays):
    """
    Compute global compositional asymmetry for multiple sequences.
    
    Parameters
    ----------
    ternary_arrays : np.ndarray
        Array of shape (n_sequences, sequence_length) with ternary values
        
    Returns
    -------
    np.ndarray
        Array of global asymmetry values for each sequence
    """
    n_sequences, seq_length = ternary_arrays.shape
    
    # Count group1 (value 1) and group2 (value -1) occurrences for each sequence
    count_A = np.sum(ternary_arrays == 1, axis=1)
    count_B = np.sum(ternary_arrays == -1, axis=1)
    
    # Calculate fractions
    FA = count_A / seq_length
    FB = count_B / seq_length
    
    # Compute asymmetry, handling division by zero
    denominator = FA + FB
    # Create mask for valid denominators
    valid_mask = denominator > 0
    
    # Initialize result array with zeros
    asymmetry = np.zeros(n_sequences)
    # Calculate only for valid denominators
    asymmetry[valid_mask] = ((FA[valid_mask] - FB[valid_mask])**2) / denominator[valid_mask]
    
    return asymmetry


def vectorized_patterning_asymmetry(ternary_arrays_winged, window_size, overall_asymmetry):
    """
    Compute patterning asymmetry for multiple sequences using sliding windows.
    
    Parameters
    ----------
    ternary_arrays_winged : np.ndarray
        Array of shape (n_sequences, winged_length) with ternary values
        
    window_size : int
        Size of the sliding window
        
    overall_asymmetry : np.ndarray
        Array of global asymmetry values for each sequence
        
    Returns
    -------
    np.ndarray
        Array of patterning asymmetry values for each sequence
    """
    n_sequences, winged_length = ternary_arrays_winged.shape
    n_windows = winged_length - window_size + 1
    
    # Initialize accumulator for squared differences
    delta_sum = np.zeros(n_sequences)
    
    # For each possible window position
    for i in range(n_windows):
        # Extract window for all sequences
        windows = ternary_arrays_winged[:, i:i+window_size]
        
        # Count group1 and group2 occurrences in each window
        window_count_A = np.sum(windows == 1, axis=1)
        window_count_B = np.sum(windows == -1, axis=1)
        
        # Calculate local fractions
        FA_local = window_count_A / window_size
        FB_local = window_count_B / window_size
        
        # Compute window asymmetry, handling division by zero
        denominator = FA_local + FB_local
        window_asymmetry = np.zeros(n_sequences)
        valid_mask = denominator > 0
        window_asymmetry[valid_mask] = ((FA_local[valid_mask] - FB_local[valid_mask])**2) / denominator[valid_mask]
        
        # Accumulate squared differences from overall asymmetry
        delta_sum += (window_asymmetry - overall_asymmetry)**2
    
    # Return average patterning asymmetry
    return delta_sum / n_windows


def vectorized_deltamax_patterning_asymmetry(ternary_arrays, window_size, overall_asymmetry):
    """
    Compute the deltamax patterning asymmetry for multiple sequences.
    
    Parameters
    ----------
    ternary_arrays : np.ndarray
        Array of shape (n_sequences, sequence_length) with ternary values
        
    window_size : int
        Size of the sliding window
        
    overall_asymmetry : np.ndarray
        Array of global asymmetry values for each sequence
        
    Returns
    -------
    np.ndarray
        Array of deltamax values for each sequence
    """
    n_sequences, seq_length = ternary_arrays.shape
    
    # Count occurrences of each value in the arrays
    count_A = np.sum(ternary_arrays == 1, axis=1)  # Group 1
    count_B = np.sum(ternary_arrays == -1, axis=1)  # Group 2
    count_C = np.sum(ternary_arrays == 0, axis=1)  # Neutral
    
    # Build maximally segregated sequences with same composition
    max_segregated_arrays = np.zeros((n_sequences, seq_length), dtype=np.int8)
    
    # For each sequence, build the most segregated version
    for i in range(n_sequences):
        pos = 0
        # First place all group1 residues
        max_segregated_arrays[i, pos:pos+count_A[i]] = 1
        pos += count_A[i]
        
        # Then place all neutral residues
        max_segregated_arrays[i, pos:pos+count_C[i]] = 0
        pos += count_C[i]
        
        # Finally place all group2 residues
        max_segregated_arrays[i, pos:] = -1
    
    # Add wings for consistent window calculation
    winged_length = seq_length + 2*window_size
    max_segregated_winged = np.zeros((n_sequences, winged_length), dtype=np.int8)
    max_segregated_winged[:, window_size:window_size+seq_length] = max_segregated_arrays
    
    # Calculate deltamax using the segregated sequences
    deltamax = vectorized_patterning_asymmetry(max_segregated_winged, window_size, overall_asymmetry)
    
    return deltamax


def vectorized_kappa_x(sequences, group1, group2, window_size, flatten=False, ternarize=True):
    """
    Calculate kappa values for multiple sequences simultaneously using vectorized operations.
    
    Parameters
    ----------
    sequences : list
        List of amino acid sequences (all must be same length)
        
    group1 : list
        List of residues considered as group 1
        
    group2 : list
        List of residues considered as group 2
        
    window_size : int
        Size of the sliding window
        
    flatten : bool
        If True, kappa values above 1 will be set to 1

    ternarize : bool, optional
        If True, sequences will be ternarized before kappa calculation        

    Returns
    -------
    np.ndarray
        Array of kappa values for each sequence
    """
    
    # Validate window size
    if window_size < 2:
        raise ValueError('window_size must be 2 or larger')
    
    if not group1:
        raise ValueError('group1 must contain one or more possible residues')
    
    n_sequences = len(sequences)
    seq_length = len(sequences[0])
    
    # Check if all sequences have the same length
    if not all(len(seq) == seq_length for seq in sequences):
        raise ValueError('All sequences must have the same length')
    
    # Check if window size is valid for the sequences
    if window_size > seq_length:
        return np.full(n_sequences, -1.0)
    
    # Ternarize sequences
    if ternarize:
        ternary_arrays = vectorized_ternarize(sequences, group1, group2)
    else:
        ternary_arrays = sequences
    
    # Create winged arrays for consistent window calculations
    wing_size = window_size
    winged_length = seq_length + 2*wing_size
    ternary_winged = np.zeros((n_sequences, winged_length), dtype=np.int8)
    ternary_winged[:, wing_size:wing_size+seq_length] = ternary_arrays
    
    # Check if any sequence lacks residues from either group
    has_group1 = np.sum(ternary_arrays == 1, axis=1) > 0
    has_group2 = np.sum(ternary_arrays == -1, axis=1) > 0
    valid_sequences = has_group1 & has_group2
    
    # Calculate overall compositional asymmetry
    overall_asymmetry = vectorized_global_asymmetry(ternary_winged)
    
    # Calculate delta (patterning asymmetry)
    delta = vectorized_patterning_asymmetry(ternary_winged, window_size, overall_asymmetry)
    
    # Calculate deltamax
    deltamax = vectorized_deltamax_patterning_asymmetry(ternary_arrays, window_size, overall_asymmetry)
    
    # Calculate kappa
    kappa = np.zeros(n_sequences)
    nonzero_deltamax = deltamax > 0
    valid_calc = valid_sequences & nonzero_deltamax
    kappa[valid_calc] = delta[valid_calc] / deltamax[valid_calc]
    
    # Mark invalid sequences
    kappa[~valid_sequences] = -1.0
    
    # Flatten kappa values if requested
    if flatten:
        kappa[kappa > 1.0] = 1.0
    
    return kappa


def batch_calculate_kappa(sequences, group1=['K', 'R'], group2=['E', 'D'], flatten=True, ternarize=True):
    """
    User-friendly function to calculate kappa values for multiple sequences.
    
    Parameters
    ----------
    sequences : list
        List of amino acid sequences (all must be same length)
        
    group1 : list, optional
        List of residues for group 1, default: positively charged residues
        
    group2 : list, optional
        List of residues for group 2, default: negatively charged residues
        
    window_size : int, optional
        Size of sliding window, default: 5
        
    flatten : bool, optional
        If True, kappa values above 1 will be set to 1
    
    ternarize : bool, optional
        If True, sequences will be ternarized before kappa calculation
        
    Returns
    -------
    np.ndarray
        Array of kappa values for each sequence
    """
    k5=vectorized_kappa_x(sequences, group1, group2, window_size=5, flatten=flatten, ternarize=ternarize)
    k6=vectorized_kappa_x(sequences, group1, group2, window_size=6, flatten=flatten, ternarize=ternarize)
    # get average of k5 and k6
    kappa=(k5+k6)/2
    return kappa

