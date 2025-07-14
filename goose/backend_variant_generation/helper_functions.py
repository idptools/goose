from sparrow.protein import Protein
import numpy as np
from goose.data.defined_aa_classes import aa_classes_by_aa, min_class_hydro, max_class_hydro
import metapredict as meta
from goose.backend import parameters
from goose import goose_exceptions

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

def hydropathy_range(length, fraction, net_charge):
    """
    Determine the mathematically possible hydropathy values based on
    the charge requirements. 
    
    parameters
    ----------
    length : int
        The length of the sequence.
    fraction : float
        The fraction of charged residues in the sequence.
    net_charge : float
        The net charge of the sequence.
    """
    charged_residues = needed_charged_residues(length, fraction, net_charge)
    negative_residues = charged_residues['negative']
    positive_residues = charged_residues['positive']

    # calculate negative residue contribution (D and E are both 1.0 so length = possible value.)
    min_negative_contribution = negative_residues
    max_negative_contribution = negative_residues
    # get min positive contribution (R = 0, so 0 is the min.)
    min_positive_contribution = 0
    # get max positive contribution (K is 0.6 so 0.6 is the max)
    max_positive_contribution = 0.6 * positive_residues

    # figure out non-charged residue number
    non_charged_residues = length - (positive_residues + negative_residues)

    # now get min non-charged contribution (Q and N = 1.0.)
    min_non_charged_contribution = non_charged_residues

    # I = 9.0
    max_non_charged_contribution = non_charged_residues * 9.0

    # calculate the min and max hydropathy values
    min_hydropathy = (min_negative_contribution + min_positive_contribution + min_non_charged_contribution) / length
    max_hydropathy = (max_negative_contribution + max_positive_contribution + max_non_charged_contribution) / length
    return min_hydropathy, max_hydropathy

def change_all_residues_within_class(input_sequence):
    '''
    function to take in an input sequence and change all residues in
    that sequence to a different residue within the same class. 
    '''
    seq_list = [np.random.choice(aa_classes_by_aa[aa]) for aa in input_sequence]
    return ''.join(seq_list)


def check_hydropathy(variant_sequence: str,
                     target_hydropathy: float,
                     hydropathy_tolerance: float) -> bool:
    """
    Check if the hydropathy of a variant sequence is within the specified tolerance
    of the target hydropathy.

    Parameters:
    -----------
    original_hydropathy : float
        The hydropathy of the original sequence.
        
    variant_sequence : str
        The sequence to check.
        
    target_hydropathy : float
        The target hydropathy value.
        
    hydropathy_tolerance : float
        The acceptable deviation from the target hydropathy.
    Returns:
    --------
    bool
        True if the variant sequence's hydropathy is within the tolerance of the target,
        False otherwise.
    """
    return abs(Protein(variant_sequence).hydrophobicity - target_hydropathy) <= hydropathy_tolerance


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
        'negative': set(['D', 'E']),
        'positive': set(['K', 'R']),
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

def decrease_res_asymmetry(sequence, residues, num_attempts=100):
    '''
    function to decrease how asymmetrrically specific residues
    are distributed in a sequence. 
    
    parameters
    -----------
    sequence : str
        The input sequence to modify.
    residues : list
        List of residues to consider for modification.
    num_attempts : int
        Number of attempts to reduce asymmetry. Default is 100.
        
    Returns
    --------
    str
        Modified sequence with reduced asymmetry in residue distribution.
    '''
    # calculate IWD of original protein. 
    starting_IWD = Protein(sequence).compute_iwd(target_residues=residues)

    for _ in range(num_attempts):
        # Convert to numpy array for easier manipulation
        # determine window size based on number of residues in sequence that are specified in residues.
        target_window_size = int(len(sequence)/len([a for a in sequence if a in residues]))+1
        window_size = max(5, target_window_size)
        # make a list of window sizes that are + and - 3 from window_size
        window_sizes = [window_size + i for i in range(-3, 4) if window_size + i > 3]
        if window_sizes == []:
            window_sizes = [5]  # fallback to minimum window size if no valid sizes
        # randomly choose a window size from the list
        window_size = np.random.choice(window_sizes)
        
        # Calculate half window size
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
    
        # Calculate the IWD of the modified sequence
        final_sequence = ''.join(final_sequence)
        final_IWD = Protein(final_sequence).compute_iwd(target_residues=residues)
    
        # If the IWD has not changed, return the original sequence
        if final_IWD < starting_IWD:
            return final_sequence

    # if we failed to reduce asymmetry after all attempts, return None
    return None

def increase_res_asymmetry(sequence, residues, num_attempts=100):
    '''
    Function to increase how asymmetrically specific residues
    are distributed in a sequence by moving them toward sequence ends.
    
    Parameters
    -----------
    sequence : str
        The input sequence to modify.
    residues : list
        List of residues to consider for modification.
    num_attempts : int
        Number of attempts to increase asymmetry. Default is 100.
        
    Returns
    --------
    str
        Modified sequence with increased asymmetry in residue distribution.
    '''
    # calculate IWD of original protein. 
    starting_IWD = Protein(sequence).compute_iwd(target_residues=residues)

    # Convert to numpy array for easier manipulation
    window_size = 5
    half_window = window_size // 2
    
    # binarize sequences
    original_array = binarize_sequence_by_residues(sequence, residues)
    seqs_array = np.copy(original_array)
    
    # Collect all target and non-target residues from original sequence
    target_residues = [aa for aa in sequence if aa in residues]
    other_residues = [aa for aa in sequence if aa not in residues]
    
    # Calculate current center of mass for target residues
    target_positions = np.where(seqs_array == 1)[0]
    if len(target_positions) == 0:
        return sequence  # No target residues to move
    
    target_center = np.mean(target_positions)
    seq_length = len(sequence)
    


    for _ in range(num_attempts):
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
        
        # Calculate the IWD of the modified sequence
        final_sequence = ''.join(final_sequence)
        final_IWD = Protein(final_sequence).compute_iwd(target_residues=residues)
        
        # If the IWD has increased, return the modified sequence
        if final_IWD > starting_IWD:
            return final_sequence
        else:
            # If not, reset seqs_array to original state for next attempt
            seqs_array = np.copy(original_array)

    # if we failed to increase asymmetry after all attempts, return None
    return None



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
    min_hydro = sum([min_class_hydro[a] for a in input_sequence])
    min_hydro = min_hydro / len(input_sequence)
    max_hydro = sum([max_class_hydro[a] for a in input_sequence])
    max_hydro = max_hydro / len(input_sequence)
    if min_hydro < 0:
        min_hydro=0
    if max_hydro > 9:
        max_hydro=9
    return round(min_hydro,5), round(max_hydro,5)



def check_variant_disorder_vectorized(
                              original_sequence,
                              sequences, 
                              strict_disorder=False,
                              disorder_cutoff=parameters.DISORDER_THRESHOLD,
                              metapredict_version=parameters.METAPREDICT_DEFAULT_VERSION,
                              return_best_sequence=False):
    """
    Check disorder of sequences using vectorized operations.
    This function defaults to comparing against a starting variant sequence. 
    By default, the sequence needs to either be above 0.5 or the variant value at that position
    to be disordered. It is the min of these values, so if the variant sequence (original_sequence)
    has a value below 0.5, that will be allowed. However, you can require that the sequence is
    disordered based on a cutoff value.

    Parameters
    ----------
    original_sequence : str
        The original sequence to compare against.
        This is used to determine if the sequence is disordered.
        If the sequence is the same as the original sequence, it is considered disordered.
    sequences : list, str, dict
        list, str, or dict of sequences to check for disorder
    strict_disorder : bool
        If True, all residues must be above the disorder cutoff to be considered disordered.
        default is False
    disorder_cutoff : float
        Cutoff for disorder. Above this value is considered disordered.
    metapredict_version : int
        Version of MetaPredict to use (1, 2, or 3)
        default is 3
    return_best_sequence : bool
        If True, return the sequence with the best disorder score
        default is False
    
    Returns
    -------
    list
        List of sequences that are disordered
    """
    # if sequences are a numpy array, convert to list
    if isinstance(sequences, np.ndarray):
        sequences = sequences.tolist()
    # convert to list if str
    if isinstance(sequences, str):
        sequences = [sequences]
    # check if sequences is a dict
    if isinstance(sequences, dict):
        sequences = list(sequences.values())

    # predict disorder for sequences
    disorder_predictions = meta.predict_disorder(sequences, version=metapredict_version)

    # separate out sequences and predictions into their own lists.
    seqs = [seq[0] for seq in disorder_predictions]
    disorder_scores = np.array([pred[1] for pred in disorder_predictions])

    # if strict_disorder is True, check if all residues are above the cutoff
    if strict_disorder:
        disorder = np.all(disorder_scores > disorder_cutoff, axis=1)

    # if strict_disorder is False, we need to get the disorder of the original sequence
    else:
        # get disorder of original sequence
        original_disorder = meta.predict_disorder([original_sequence], version=metapredict_version)
        original_disorder_scores = np.array(original_disorder[0][1])

        # replace values in original_disorder_scores that are above disorder_cutoff with disorder_cutoff
        original_disorder_scores[original_disorder_scores > disorder_cutoff] = disorder_cutoff
        
        # Create binary mask where True means that the values in the disorder_scores
        # are greater than original_disorder_scores
        mask = disorder_scores > original_disorder_scores

        # now get sequences that have all values above the original disorder scores
        disorder = np.all(mask, axis=1)

    # return sequences that are disordered
    disordered_seqs=[seqs[i] for i in range(len(seqs)) if disorder[i]]

    # if no disorderd sequences and return_best_sequence is True,
    # return the sequence with the highest mean disorder score
    if return_best_sequence and len(disordered_seqs) == 0:
        # return the sequence with the highest mean disorder score
        best_index = np.argmax(np.mean(disorder_scores, axis=1))
        return seqs[best_index]
    
    if disordered_seqs == []:
        # if no disordered sequences, return None
        return None
    
    # return disordered_seqs
    return disordered_seqs