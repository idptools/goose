from sparrow.protein import Protein
from numpy.lib.stride_tricks import sliding_window_view
import numpy as np
from goose.data.defined_aa_classes import aa_classes, aa_classes_by_aa, aa_class_indices, min_class_hydro, max_class_hydro
from goose.backend_property_calculation.calculate_properties_single_sequence import sequence_to_array, array_to_sequence, calculate_hydropathy_single_sequence, calculate_ncpr_single_sequence, calculate_fcr_single_sequence

def check_properties(
        variant_sequence,
        target_hydropathy=None,
        target_ncpr=None,
        target_fcr=None,
        target_kappa=None,
        hydropathy_tolerance=0.05,
        kappa_tolerance=0.03):
    '''
    Check if properties are what we want them to be. Only
    checks properties that are not None. If all properties
    are None, returns True. If any of the properties are not
    what we want them to be, returns False.
    
    Parameters
    -----------
    variant_sequence : str
        The sequence to check properties for.
        
    target_hydropathy : float or None
        The target hydropathy value. If None, this property is not checked.
        
    target_ncpr : float or None
        The target NCPR value. If None, this property is not checked.
        
    target_fcr : float or None
        The target FCR value. If None, this property is not checked.
        
    target_kappa : float or None
        The target kappa value. If None, this property is not checked.

    hydropathy_tolerance : float
        The tolerance for hydropathy comparison. Default is 0.05.
        
    kappa_tolerance : float
        The tolerance for kappa comparison. Default is 0.03.
        
    Returns
    --------
    bool
        True if all properties match the targets, False otherwise.
    '''
    # if we had changes to FCR or NCPR, we need to allow for a small amount of error
    # as some NCPR / FCR combinations are not possible.
    if target_fcr is not None:
        allowed_ncpr_error = 1/len(variant_sequence)
    else:
        allowed_ncpr_error = 0
    if target_ncpr is not None:
        allowed_fcr_error = 1/len(variant_sequence)
    else:
        allowed_fcr_error = 0

    # Check each property against its target
    if target_hydropathy is not None:
        current_hydropathy = Protein(variant_sequence).hydrophobicity
        if abs(current_hydropathy - target_hydropathy) > hydropathy_tolerance:
            return False
    if target_ncpr is not None:
        current_ncpr = Protein(variant_sequence).NCPR
        if round(abs(current_ncpr - target_ncpr),5) > allowed_ncpr_error:
            return False
    if target_fcr is not None:
        current_fcr = Protein(variant_sequence).FCR
        if round(abs(current_fcr - target_fcr),5) > allowed_fcr_error:
            return False
    if target_kappa is not None:
        current_kappa = Protein(variant_sequence).kappa
        if round(abs(current_kappa - target_kappa),5) > kappa_tolerance:
            return False

    # If all checks passed, return True
    return True


def fraction_net_charge(length, fraction, net_charge):

    """

    Function for figuring out NCPR / FCR values for a sequence based 
    on the sequence length, the wanted fraction, and the wanted
    net_charge. The reason this is a bit convoluted is that the
    function was written to get the fraction and net_charge as
    close as possible, and if the FCR / NCPR / length 
    combination is not possible, to prioritize
    getting the NCPR correct over the FCR. 


    Parameters
    -------------

    length : Int
        How long the sequence is.

    fraction : Float
        The objective FCR as a decimal value.

    net_charge : Float
        The objective NCPR as a decimal value.


    Returns
    ---------

    Dict
        Returns a dictionary holding the NCPR and FCR
        values based on the length, fraction, and net_charge.

    """

    # figure out how many charged residues will be needed
    added_NCPR_residues = round(length * abs(net_charge))
    
    #figure out the remaining FCR fraction
    remaining_FCR_fraction = fraction - abs(net_charge)

    # based on the remamining_FCR_fraction, figure out how many residues to add
    added_FCR_residues = round((length * remaining_FCR_fraction)/2)
    remaining_residues_over_two = round(length - added_NCPR_residues/2)

    # figure ot the number of possible residues to add to the sequence
    remaining_residues = (length - added_NCPR_residues)

    # if the number of added_FCR_residues are less than or equal to 
    # remaining residues over 2, added_residues = added_FCR_residues        
    if added_FCR_residues <= remaining_residues_over_two:
        added_residues = added_FCR_residues 

    # otherwise, added_residues are qualt to remaining_residues_over_two       
    else:
        added_residues = remaining_residues_over_two

    # if the number of residues to add is greater than the remaining number of residues...
    if (added_residues * 2) > remaining_residues:
        added_residues = added_residues - 1

    # otherwise, don't do anything to added_residues       
    else:
        added_residues = added_residues 
    # return the final dict
    return {'NCPR_residues': added_NCPR_residues, 'FCR_residues': added_residues}   


def needed_charged_residues(length, fraction, net_charge):
    '''
    function to figure out how many positive and 
    how many negative charged residues are need
    to satisfy the ncpr and fcr of a 
    given sequence.
    '''
    total_charged = fraction_net_charge(length, fraction, net_charge)

    if net_charge > 0:
        num_positive = total_charged['NCPR_residues'] + round(total_charged['FCR_residues'])
        num_negative = round(total_charged['FCR_residues'])
    elif net_charge < 0:
        num_negative = total_charged['NCPR_residues'] + round(total_charged['FCR_residues'])
        num_positive = round(total_charged['FCR_residues'])
    else:
        num_positive = round(total_charged['FCR_residues'])
        num_negative = round(total_charged['FCR_residues'])

    return {'positive':num_positive, 'negative':num_negative}

def change_all_residues_within_class(input_sequence):
    '''
    function to take in an input sequence and change all residues in
    that sequence to a different residue within the same class. 
    '''
    seq_list = [np.random.choice(aa_classes_by_aa[aa]) for aa in input_sequence]
    return ''.join(seq_list)

def optimize_hydropathy_avoid_original_residues(
        original_sequence,
        variant_sequence,
        target_hydropathy,
        max_iterations=1000,
        tolerance=0.05):
    """
    Optimize hydropathy of a sequence while avoiding original residues.

    Parameters:
    -----------
    original_sequence : str
        The original sequence to avoid.
    variant_sequence : str
        The sequence to optimize.
    target_hydropathy : float
        Target mean hydropathy value to achieve.
    max_iterations : int
        Maximum number of optimization iterations.
    tolerance : float
        Acceptable difference between achieved and target hydropathy.

    Returns:
    --------
    str
        Optimized sequence with hydropathy close to the target.
    """
    
    # Convert sequences to numpy arrays for vectorized operations
    original_seq = sequence_to_array(original_sequence)
    variant_seq = sequence_to_array(variant_sequence)
    
    # set mutable residues. Can't change D,E,K,R,H,G,C,P
    mutable_residues=[0, 4, 7, 8, 9, 10, 11, 13, 14, 15, 16, 17, 18, 19]
    
    current_hydropathy = calculate_hydropathy_single_sequence(variant_seq)
    
    for _ in range(max_iterations):
        if abs(current_hydropathy - target_hydropathy) < tolerance:
            break
        
        # Find positions that can be changed
        changeable_positions = [i for i in range(len(variant_seq)) if variant_seq[i] in mutable_residues]
        
        if not changeable_positions:
            return ''.join(array_to_sequence(variant_seq))
        
        # Randomly select a position to change
        pos_to_change = np.random.choice(changeable_positions)
        
        # Find the best replacement amino acid
        current_aa = variant_seq[pos_to_change]
        best_hydro_error = np.inf
        
        for aa in aa_class_indices[current_aa]:
            if aa not in original_seq:  # Avoid original residues
                new_seq = variant_seq.copy()
                new_seq[pos_to_change] = aa
                new_hydro = calculate_hydropathy_single_sequence(new_seq)

                if abs(new_hydro-target_hydropathy) < best_hydro_error:
                    best_hydro_error = abs(new_hydro - target_hydropathy)
                    variant_seq = new_seq.copy()
                    if best_hydro_error < tolerance:
                        # If we found a suitable amino acid, break early
                        break
    return array_to_sequence(variant_seq)



def calculate_amino_acid_class_fractions(sequence: str) -> dict:
    """
    Calculate the fraction of amino acids in each class for a given sequence.
    
    Parameters:
    -----------
    sequence : str
        Input protein sequence
        
    Returns:
    --------
    dict
        Dictionary with class names as keys and fractions as values
    """
    # Define amino acid classes
    aa_classes = {
        'aliphatic': set(['A', 'I', 'L', 'M', 'V']),
        'aromatic': set(['F', 'W', 'Y']),
        'polar': set(['S', 'T', 'N', 'Q']),
        'positive': set(['D', 'E']),
        'negative': set(['K', 'R']),
        'glycine': set(['G']),
        'cysteine': set(['C']),
        'proline': set(['P']),
        'histidine': set(['H'])
    }
    
    # Convert sequence to uppercase and get length
    sequence = sequence.upper()
    seq_length = len(sequence)
    
    if seq_length == 0:
        return {class_name: 0.0 for class_name in aa_classes.keys()}
    
    # Count amino acids in each class
    class_counts = {class_name: 0 for class_name in aa_classes.keys()}
    
    for aa in sequence:
        for class_name, class_aas in aa_classes.items():
            if aa in class_aas:
                class_counts[class_name] += 1
                break
    
    # Calculate fractions
    class_fractions = {class_name: count / seq_length 
                        for class_name, count in class_counts.items()}
    
    return class_fractions


def binarize_sequence_by_residues(sequence: str,
                       residues: list) -> np.ndarray:
    """
    Convert a sequence of amino acids into a ternary representation.
    
    Parameters:
    -----------
    sequence : str
        Input protein sequence.
        
    residues : list
        List of amino acids to include in the ternary representation.
        Residues in this list will be represented as 1, while all other residues
        will be represented as 0.
        
    Returns:
    --------
    np.ndarray
        Ternary representation of the sequence.
    """
    # Create a binary representation of the sequence
    binary_representation = np.zeros(shape=[len(sequence),], dtype=int)
    
    for i, aa in enumerate(sequence):
        if aa in residues:
            binary_representation[i] = 1
            
    return binary_representation



def get_linear_profile(sequence: str,
                       residues: list,
                       window_size: int,
                       return_binarized_seq=False) -> np.ndarray:
    """
    Get a linear profile of the sequence based on the presence of specified residues.
    The profile is calculated as a sliding window average of the binary representation.

    Parameters:
    -----------
    sequence : str
        Input protein sequence.
        
    residues : list
        List of amino acids to include in the profile.
        
    window_size : int
        Size of the sliding window for averaging.

    return_binarized_seq : bool
        If True, return the binary representation of the sequence instead of the profile.
        
    Returns:
    --------
    np.ndarray
        Linear profile of the sequence with shape (len(sequence), len(residues)).
    """
    # Check if window size is valid
    if window_size <= 0 or window_size > len(sequence):
        raise ValueError("Window size must be a positive integer less than or equal to the sequence length.")
    
    # Convert sequence to binary representation
    binary_array = binarize_sequence_by_residues(sequence, residues)
    
    # Initialize output array
    profile = np.zeros_like(binary_array, dtype=float)
    
    # Calculate sliding window average for each position and residue type
    half_window = window_size // 2
    
    for i in range(len(sequence)):
        # Define window boundaries
        start = max(0, i - half_window)
        end = min(len(sequence), i + half_window + 1)
        
        # Calculate average over the window for each residue type
        profile[i] = np.mean(binary_array[start:end], axis=0)
    
    if return_binarized_seq:
        return profile, binary_array
    
    return profile

def decrease_res_asymmetry(sequence, residues):
    '''
    function to decrease how asymmetrrically specific residues
    are distributed in a sequence. 
    
    parameters
    -----------
    sequence : str
        The input sequence to modify.
    residues : list
        List of residues to consider for modification.
        
    Returns
    --------
    str
        Modified sequence with reduced asymmetry in residue distribution.
    '''
    # Convert to numpy array for easier manipulation
    window_size = 5
    half_window = window_size // 2
    
    # Calculate the linear NCPR for each sequence
    linear_profile, seqs_array = get_linear_profile(sequence, 
                                        residues, 
                                        window_size=window_size,
                                        return_binarized_seq=True)
    
    # get residues in order that are in or not in residues
    target_residues = [a for a in sequence if a in residues]
    other_residues = [a for a in sequence if a not in residues]

    # sort the linear_profile by absolute value
    sorted_indices = np.argsort(linear_profile, axis=0)[::-1]

    # use np.random to get a random index from the first 10% of the sorted indices
    num_to_choose = max(1, len(sorted_indices) // 10)  # at least one position
    chosen_indices = sorted_indices[:num_to_choose]
    random_first_index = np.random.choice(chosen_indices)
    # now get from the last 10%
    last_indices = sorted_indices[-num_to_choose:]
    random_last_index = np.random.choice(last_indices)
    
    # get all positions that are 1 in the sequence
    target_positions = np.where(seqs_array == 1)[0]

    # get a charged position from the window in random_first_index
    first_window_start = max(0, random_first_index - half_window)
    first_window_end = min(len(sequence), random_first_index + half_window + 1)
    first_window_positions = np.arange(first_window_start, first_window_end)
    first_window_target_positions = np.intersect1d(target_positions, first_window_positions)
    # choose a random positon from first_window_target_positions
    if len(first_window_target_positions) > 0:
        random_first_pos = np.random.choice(first_window_target_positions)
    else:
        # If no target positions in window, return the sequence
        return sequence
    
    # get a random position in the random_last_window
    last_window_start = max(0, random_last_index - half_window)
    last_window_end = min(len(sequence), random_last_index + half_window + 1)
    last_window_positions = np.arange(last_window_start, last_window_end)
    random_last_position = np.random.choice(last_window_positions)
    
    # Move the charge from random_first_pos to random_last_position
    if random_first_pos != random_last_position:
        # Store what's currently at the target position
        displaced_value = seqs_array[random_last_position]
        # Move the charge to the target position
        charge_value = seqs_array[random_first_pos]
        seqs_array[random_last_position] = charge_value
        # Put the displaced value where the charge was
        seqs_array[random_first_pos] = displaced_value

    # Convert back to sequence using improved mapping
    # Collect all target and non-target residues from original sequence
    target_residues = [aa for aa in sequence if aa in residues]
    other_residues = [aa for aa in sequence if aa not in residues]
    
    # Build the final sequence based on the modified binary array
    final_sequence = []
    target_idx = 0
    other_idx = 0
    
    for position_value in seqs_array:
        if position_value == 1:
            if target_idx < len(target_residues):
                final_sequence.append(target_residues[target_idx])
                target_idx += 1
        else:
            if other_idx < len(other_residues):
                final_sequence.append(other_residues[other_idx])
                other_idx += 1

    if len(final_sequence) != len(sequence):
        raise ValueError("Final sequence length does not match original sequence length.")
    
    # Join the final sequence into a string  
    return ''.join(final_sequence)

def increase_res_asymmetry(sequence, residues):
    '''
    Function to increase how asymmetrically specific residues
    are distributed in a sequence by moving them toward sequence ends.
    
    Parameters
    -----------
    sequence : str
        The input sequence to modify.
    residues : list
        List of residues to consider for modification.
        
    Returns
    --------
    str
        Modified sequence with increased asymmetry in residue distribution.
    '''
    # Convert to numpy array for easier manipulation
    window_size = 5
    half_window = window_size // 2
    
    # Calculate the linear profile for the sequence
    linear_profile, seqs_array = get_linear_profile(sequence, 
                                        residues, 
                                        window_size=window_size,
                                        return_binarized_seq=True)
    
    # Collect all target and non-target residues from original sequence
    target_residues = [aa for aa in sequence if aa in residues]
    other_residues = [aa for aa in sequence if aa not in residues]
    
    # Calculate current center of mass for target residues
    target_positions = np.where(seqs_array == 1)[0]
    if len(target_positions) == 0:
        return sequence  # No target residues to move
    
    target_center = np.mean(target_positions)
    seq_length = len(sequence)
    
    # Determine which end to push target residues towards
    # If center of mass is already in right half, push further right
    # If center of mass is in left half, push further left
    if target_center > seq_length / 2:
        # Push target residues further right (toward end)
        # Find leftmost target residue to move right
        leftmost_target_pos = np.min(target_positions)
        
        # Find available positions in right half that don't already have target residues
        right_half_start = seq_length // 2
        available_right_positions = []
        for i in range(right_half_start, seq_length):
            if seqs_array[i] == 0 and i > leftmost_target_pos:
                available_right_positions.append(i)
        
        if available_right_positions:
            # Choose random position from available right positions
            target_pos = np.random.choice(available_right_positions)
            
            # Perform the swap
            seqs_array[target_pos] = 1
            seqs_array[leftmost_target_pos] = 0
    else:
        # Push target residues further left (toward beginning)
        # Find rightmost target residue to move left
        rightmost_target_pos = np.max(target_positions)
        
        # Find available positions in left half that don't already have target residues
        left_half_end = seq_length // 2
        available_left_positions = []
        for i in range(left_half_end):
            if seqs_array[i] == 0 and i < rightmost_target_pos:
                available_left_positions.append(i)
        
        if available_left_positions:
            # Choose random position from available left positions
            target_pos = np.random.choice(available_left_positions)
            
            # Perform the swap
            seqs_array[target_pos] = 1
            seqs_array[rightmost_target_pos] = 0
    
    # Convert back to sequence using the same mapping as decrease_res_asymmetry
    final_sequence = []
    target_idx = 0
    other_idx = 0
    
    for position_value in seqs_array:
        if position_value == 1:
            if target_idx < len(target_residues):
                final_sequence.append(target_residues[target_idx])
                target_idx += 1
        else:
            if other_idx < len(other_residues):
                final_sequence.append(other_residues[other_idx])
                other_idx += 1

    if len(final_sequence) != len(sequence):
        raise ValueError("Final sequence length does not match original sequence length.")
    
    return ''.join(final_sequence)


def find_hydro_range_constant_class(
        input_sequence: str):
    """
    Function to find the hydropathy range of a constant class variant.
    
    Parameters:
    -----------
    input_sequence : str
        The input protein sequence.
        
    Returns:
    --------
    tuple
        A tuple containing the minimum and maximum hydropathy values for the constant class variant.
    """
    # convert input sequence to
    min_hydro = sum([min_class_hydro[a] for a in input_sequence])/len(input_sequence)
    max_hydro = sum([max_class_hydro[a] for a in input_sequence])/len(input_sequence)
    return round(min_hydro,5), round(max_hydro,5)

