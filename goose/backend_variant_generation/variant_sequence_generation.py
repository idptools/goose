'''
Updated code for variant generation.
This code specifies the changes that are made to each variant sequence. 
If any function fails to generate a valid sequence, it will return None.
This let's us check for 'None' in the variant_generation.py file before checking for disorder. 
'''
import numpy as np
import random

from sparrow.protein import Protein

from goose import parameters
from goose.parameters import get_min_re, get_max_re, get_min_rg, get_max_rg
from goose.backend_sequence_generation.sequence_generation import by_properties
from goose.backend_property_optimization.optimize_kappa import optimize_kappa
from goose.backend_sequence_generation.sequence_by_class import create_sequence_by_class
from goose.backend_property_optimization.optimize_hydropathy import optimize_hydropathy, optimize_hydropathy_within_class, optimize_hydropathy_minimal_changes, optimize_hydropathy_within_class_avoid_original_residues
from goose.backend_variant_generation.helper_functions import needed_charged_residues, change_all_residues_within_class, calculate_amino_acid_class_fractions, decrease_res_asymmetry, increase_res_asymmetry, find_hydro_range_constant_class, check_hydropathy
from goose.backend_property_optimization.optimize_dimensions import optimize_seq_dims, predict_re, predict_rg, batch_predict_re, batch_predict_rg


def shuffle_specific_regions_sequence(input_sequence: str,
                           shuffle_regions: list) -> str:
    """
    Generate a variant of the input sequence by shuffling specified residue indices.
    
    Parameters:
    ----------
    input_sequence : str
        The input protein sequence to be modified.
    shuffle_regions : list of tuples
        List of tuples where each tuple contains start and end indices of regions to shuffle.

    Returns:
    --------
    str or None
        If the variant generation fails, None is returned.
        If successful, a new protein sequence with specified regions shuffled is returned.
    """
    # Shuffle specified regions
    for start, end in shuffle_regions:
        region = list(input_sequence[start:end])
        shuffled_region = ''.join(np.random.permutation(region))
        input_sequence = input_sequence[:start] + shuffled_region + input_sequence[end:]

    return input_sequence

def shuffle_except_specific_regions_sequence(input_sequence: str,
                                    excluded_regions: list) -> str:
    """
    Generate a variant of the input sequence by shuffling residues not specified in excluded_regions.
    Parameters:
    ----------
    input_sequence : str
        The input protein sequence to be modified.
    excluded_regions : list of tuples
        List of tuples where each tuple contains start and end indices of regions to exclude from shuffling.

    Returns:
    --------
    str or None
        If the variant generation fails, None is returned.
        If successful, a new protein sequence with specified regions excluded from shuffling is returned.   
    """
    # Convert input sequence to a list for easier manipulation
    seq_list = list(input_sequence)
    # Get indices of residues to shuffle
    indices_to_shuffle = [i for i in range(len(seq_list)) if not any(start <= i < end for start, end in excluded_regions)]
    # If no residues to shuffle, return original sequence
    if not indices_to_shuffle:
        return input_sequence
    # Shuffle the residues at the target indices
    shuffled_region = np.random.permutation([seq_list[i] for i in indices_to_shuffle])
    # Place shuffled residues back into the sequence
    for idx, new_aa in zip(indices_to_shuffle, shuffled_region):
        seq_list[idx] = new_aa
    return ''.join(seq_list)

def shuffle_except_specific_residues(input_sequence: str,
                             excluded_residues: list) -> str:
    """
    Generate a variant of the input sequence by shuffling residues not specified in exclude_indices.
    
    Parameters:
    ----------
    input_sequence : str
        The input protein sequence to be modified.
    excluded_residues : list
        List of residues to exclude from shuffling.

    Returns:
    --------
    str or None
        If the variant generation fails, None is returned.
        If successful, a new protein sequence with specified residues excluded from shuffling is returned.
    
    Notes:
    -------
    This function can take in class names (e.g. 'charged', 'polar') as well as individual amino acids.
    
    """
    # dict of classes that are possible to choose
    classdict={'charged':['D', 'E', 'K', 'R'], 'polar':['Q', 'N', 'S', 'T'], 'aromatic':
    ['F', 'W', 'Y'], 'aliphatic': ['I', 'V', 'L', 'A', 'M'], 'negative':['D', 'E'], 'positive':['K', 'R']}

    if isinstance(excluded_residues, str):
        if excluded_residues in classdict:
            # if target_aas is a class, get the list of amino acids in that class
            excluded_residues = classdict[excluded_residues]
        else:
            # if target_aas is a single amino acid, convert to list
            excluded_residues = list(excluded_residues)

    # Convert input sequence to a list for easier manipulation
    seq_list = list(input_sequence)
    
    # Get indices of target residues
    indices_to_shuffle = [i for i, aa in enumerate(seq_list) if aa not in excluded_residues]
    # If no residues to shuffle, return original sequence
    if not indices_to_shuffle:
        return input_sequence
    
    # Shuffle the residues at the target indices
    shuffled_region = np.random.permutation([seq_list[i] for i in indices_to_shuffle])
    
    # Place shuffled residues back into the sequence
    for idx, new_aa in zip(indices_to_shuffle, shuffled_region):
        seq_list[idx] = new_aa
    
    return ''.join(seq_list)

def shuffle_specific_residues_sequence(input_sequence: str,
                                      target_residues: list):
    '''
    function that will only shuffle residues specified. All
    other residues will remain in their original positions.
    
    Parameters:
    ----------
    input_sequence : str
        The input protein sequence to be modified.
    target_residues : list
        List of residues to shuffle within their positions.
    
    Returns:
    --------
    str or None
        If the variant generation fails, None is returned.
        If successful, a new protein sequence with specified residues shuffled is returned.

    Notes:
    -------
    This function can take in class names (e.g. 'charged', 'polar') as well as individual amino acids.

    '''
    # dict of classes that are possible to choose
    classdict={'charged':['D', 'E', 'K', 'R'], 'polar':['Q', 'N', 'S', 'T'], 'aromatic':
    ['F', 'W', 'Y'], 'aliphatic': ['I', 'V', 'L', 'A', 'M'], 'negative':['D', 'E'], 'positive':['K', 'R']}

    if isinstance(target_residues, str):
        if target_residues in classdict:
            # if target_aas is a class, get the list of amino acids in that class
            target_residues = classdict[target_residues]
        else:
            # if target_aas is a single amino acid, convert to list
            target_residues = list(target_residues)

    # Convert input sequence to a list for easier manipulation
    seq_list = list(input_sequence)
    
    # Get indices of target residues
    indices_to_shuffle = [i for i, aa in enumerate(seq_list) if aa in target_residues]
    
    # Shuffle the residues at the target indices
    shuffled_region = np.random.permutation([seq_list[i] for i in indices_to_shuffle])
    
    # Place shuffled residues back into the sequence
    for idx, new_aa in zip(indices_to_shuffle, shuffled_region):
        seq_list[idx] = new_aa
    
    return ''.join(seq_list)


def weighted_shuffle_specific_residues_sequence(sequence, target_aas, shuffle_weight):
    '''
    function that will let you perform a weighted shuffle a sequence by 
    specifying residues or classes of residues to shuffle and a weight that
    corresponds to the degree of shuffling that you want to perform. 
    The weight is a number between 0.0-1.0 and corresponds to the probability of moving a residue during shuffling. If you specify target amino acids, only
    those amino acids are included in the shuffling and weighting can still be 
    applied to only those target amino acids.
    parameters
    ----------
    sequence : str
        the amino acid sequence as a string
    target_aas : str or list
        a list of amino acids to target for shuffling
        or a class of amino acids to target for shuffling
        Possible target classes:
            charged : DEKR
            polar : QNST
            aromatic : FYW
            aliphatic : IVLAM
            negative: DE
            positive : KR
            
    shuffle_weight : float
        a weight between 0.0-1.0 representing the probability of 
        moving a residue during shuffling
    '''
    # dict of classes that are possible to choose
    classdict={'charged':['D', 'E', 'K', 'R'], 'polar':['Q', 'N', 'S', 'T'], 'aromatic':
    ['F', 'W', 'Y'], 'aliphatic': ['I', 'V', 'L', 'A', 'M'], 'negative':['D', 'E'], 'positive':['K', 'R']}

    if isinstance(target_aas, str):
        if target_aas in classdict:
            # if target_aas is a class, get the list of amino acids in that class
            target_aas = classdict[target_aas]
        else:
            # if target_aas is a single amino acid, convert to list
            target_aas = list(target_aas)

    # define probabilities of not relocating a residue and relocating a residue, respectively
    shuffle_weights=[1-shuffle_weight, shuffle_weight]

    # get list of target amino acids from the sequence
    target_aa_list = [aa for aa in sequence if aa in target_aas]

    # perform weighted sample to determine residues to shuffle
    target_mask = random.choices([False, True], weights=shuffle_weights, k=len(target_aa_list))
    mask = [target_mask.pop(0) if aa in target_aas else False for aa in sequence]

    # gather positions and identities of residues to shuffle
    orig_scramble_positions = [i for i, val in enumerate(mask) if val == True]
    orig_aas = [sequence[i] for i in orig_scramble_positions]

    # perform Fisher-Yates shuffle only with target_aas marked for shuffling
    for i in range(len(orig_aas) - 1, 0, -1):
        remaining_reposition_sites = orig_scramble_positions[:i]

        # randomly select new position for relocation of the target aa
        new_position = random.choice(remaining_reposition_sites)
        new_position_index = remaining_reposition_sites.index(new_position)

        # swap positions with another target aa marked for shuffling
        orig_aas[i], orig_aas[new_position_index] = orig_aas[new_position_index], orig_aas[i]

    # build final seq. Uses shuffled residues at sites marked for repositioning. Otherwise, uses the original residue at that site
    return ''.join( [orig_aas.pop(0) if i in orig_scramble_positions else aa for i, aa in enumerate(sequence)] )


def targeted_reposition_specific_residues_sequence(input_sequence: str,
                                        target_residues: list) -> str:
    """
    Generate a variant by repositioning specified residues to random new positions
    while maintaining the relative order of non-target residues.
    
    This function differs from targeted shuffle in that target residues can move
    to any position in the sequence, displacing non-target residues, rather than
    just exchanging positions among themselves.
    
    Parameters:
    ----------
    input_sequence : str
        The input protein sequence to be modified.
    target_residues : list or str
        List of amino acids to target for repositioning, or a class name.
        Possible target classes:
            'charged' : ['D', 'E', 'K', 'R']
            'polar' : ['Q', 'N', 'S', 'T'] 
            'aromatic' : ['F', 'W', 'Y']
            'aliphatic' : ['I', 'V', 'L', 'A', 'M']
            'negative' : ['D', 'E']
            'positive' : ['K', 'R']

    Returns:
    --------
    str
        Protein sequence with target residues repositioned to new random locations.
        
    Raises:
    -------
    ValueError
        If target_residues contains invalid amino acids or class names.
    """
    # Define amino acid classes
    amino_acid_classes = {
        'charged': ['D', 'E', 'K', 'R'], 
        'polar': ['Q', 'N', 'S', 'T'], 
        'aromatic': ['F', 'W', 'Y'], 
        'aliphatic': ['I', 'V', 'L', 'A', 'M'], 
        'negative': ['D', 'E'], 
        'positive': ['K', 'R']
    }
    
    # Valid amino acids
    valid_amino_acids = set('ACDEFGHIKLMNPQRSTVWY')
    
    # Process target_residues input
    if isinstance(target_residues, str):
        if target_residues in amino_acid_classes:
            target_residues = amino_acid_classes[target_residues]
        else:
            # Convert single amino acid string to list
            target_residues = list(target_residues.upper())
    
    # Validate target residues
    if not isinstance(target_residues, list):
        raise ValueError("target_residues must be a string (class name or amino acids) or list of amino acids")
    
    for residue in target_residues:
        if residue not in valid_amino_acids:
            raise ValueError(f"Invalid amino acid '{residue}'. Must be one of: {sorted(valid_amino_acids)}")
    
    # Convert to set for faster lookup
    target_residue_set = set(target_residues)
    
    # Separate target and non-target residues while preserving information
    target_residue_list = []
    non_target_residue_list = []
    
    for amino_acid in input_sequence:
        if amino_acid in target_residue_set:
            target_residue_list.append(amino_acid)
        else:
            non_target_residue_list.append(amino_acid)
    
    # Handle edge cases
    if not target_residue_list:
        # No target residues found, return original sequence
        return input_sequence
    
    if not non_target_residue_list:
        # All residues are targets, just shuffle them
        return ''.join(random.sample(target_residue_list, len(target_residue_list)))
    
    # Generate random positions for target residues
    all_positions = list(range(len(input_sequence)))
    new_target_positions = sorted(random.sample(all_positions, k=len(target_residue_list)))
    
    # Build final sequence
    final_sequence = [''] * len(input_sequence)
    target_index = 0
    non_target_index = 0
    
    for position in range(len(input_sequence)):
        if position in new_target_positions:
            final_sequence[position] = target_residue_list[target_index]
            target_index += 1
        else:
            final_sequence[position] = non_target_residue_list[non_target_index]
            non_target_index += 1
    
    return ''.join(final_sequence)



def change_residue_asymmetry_sequence(sequence: str,
                                target_residues: list,
                                num_changes: int = 1,
                                increase_or_decrease='increase') -> str:
    """
    Variant where a class of residues (see below for classes) 
    or a user-specified list of residues is changed to become
    more asymmetrically or less asymmetrically distributed 
    throughout the sequence. Does NOT change sequence composition.
    
    Parameters:
    ----------
    sequence : str
        The input protein sequence to be modified.
    target_residues : list
        List of residues to change to increase or decrease asymmetry.
    num_changes : int, default=1
        Number of residues to change in the sequence.
    increase_or_decrease : str, default='increase'
        Whether to increase or decrease asymmetry. Options are 'increase' or 'decrease'.
    
    Returns:
    --------
    str or None
        If the variant generation fails, None is returned.
        If successful, a new protein sequence with altered residue asymmetry is returned.
    
    Notes:
    -----
    The residue classes are defined as follows:
    - Polar: ['S', 'T', 'N', 'Q']
    - Aliphatic: ['A', 'G', 'V', 'L', 'I', 'M']
    - Aromatic: ['F', 'W', 'Y']
    - Charged: ['D', 'E', 'K', 'R']
    - Cysteine: ['C']
    - Glycine: ['G']
    - Proline: ['P']
    
    This function will change the residues in the sequence to increase or decrease
    the asymmetry of the specified residues. It will not change the overall composition
    of the sequence, but will change the distribution of the specified residues.
    """
    # dict of classes that are possible to choose
    classdict={'charged':['D', 'E', 'K', 'R'], 'polar':['Q', 'N', 'S', 'T'], 'aromatic':
    ['F', 'W', 'Y'], 'aliphatic': ['I', 'V', 'L', 'A', 'M'], 'negative':['D', 'E'], 'positive':['K', 'R']}

    if isinstance(target_residues, str):
        if target_residues in classdict:
            # if target_aas is a class, get the list of amino acids in that class
            target_residues = classdict[target_residues]
        else:
            # if target_aas is a single amino acid, convert to list
            target_residues = list(target_residues)

    for _ in range(num_changes):
        if increase_or_decrease == 'increase':
            new_seq = increase_res_asymmetry(sequence, target_residues, num_attempts=1000)
        elif increase_or_decrease == 'decrease':
            new_seq = decrease_res_asymmetry(sequence, target_residues, num_attempts=1000)
        else:
            raise ValueError("increase_or_decrease must be 'increase' or 'decrease'")
        if new_seq is None:
            return sequence
        sequence = new_seq
    
    return sequence


def constant_residues_and_properties_sequence(sequence: str,
                                  constant_residues: list,
                                  hydropathy_tolerance=parameters.MAXIMUM_HYDRO_ERROR) -> str:
    '''
    Generate a new sequence with the same hydropathy,
    FCR, NCPR, and kappa as the input sequence,
    while keeping the specified residues constant.
    
    Parameters:
    ----------
    sequence : str
        The input protein sequence to be modified.
    constant_residues : list
        List of residues to keep constant in the new sequence.
    hydropathy_tolerance : float, default=0.05
        set by parameters.MAXIMUM_HYDRO_ERROR
        Acceptable difference for hydropathy optimization.
    kappa_tolerance : float, default=0.03
        set by parameters.MAXIMUM_KAPPA_ERROR
        Acceptable difference for kappa optimization.
    
    Returns:
    --------
    str or None
        If the variant generation fails, None is returned.
        If successful, a new protein sequence with specified residues excluded from changing is returned.
    '''
    target_hydropathy = Protein(sequence).hydrophobicity
    # make sequence without specified residues
    modified_sequence = ''.join([aa for aa in sequence if aa not in constant_residues])
    
    # Handle edge case where all residues are constant
    if len(modified_sequence) == 0:
        return sequence  # Return original sequence if nothing to vary
    
    # get properties for the modified sequence
    protein = Protein(modified_sequence)
    original_hydropathy = protein.hydrophobicity
    original_kappa = protein.kappa
    original_FCR = protein.FCR
    original_NCPR = protein.NCPR
    length = len(modified_sequence)

    # Generate a new sequence with the same properties
    for i in range(10):
        new_sequences = by_properties(
            length=length,
            fcr=original_FCR,
            ncpr=original_NCPR,
            hydropathy=original_hydropathy,
            kappa=original_kappa,
            check_disorder=False,
            batch_size=200,
            exclude_residues=constant_residues
        )
        if new_sequences != None:
            break

    if isinstance(new_sequences, str):
        new_sequences = [new_sequences]  # Convert to list for consistency

    if new_sequences is None or len(new_sequences) == 0:
        return None

    all_sequences=[]
    for new_sequence in new_sequences:
        # Create final sequence by inserting constant residues back at their original positions
        final_sequence = list(sequence)  # Start with original sequence as template
        new_seq_index = 0  # Index for the new sequence (without constant residues)
        
        # Iterate through each position in the original sequence
        for i, original_aa in enumerate(sequence):
            if original_aa not in constant_residues:
                # Replace with amino acid from new sequence
                final_sequence[i] = new_sequence[new_seq_index]
                new_seq_index += 1
            else:
                # Keep the constant residue in place
                final_sequence[i] = original_aa
        
        # if constant residues specified as charged, we cannot use kappa optimzation because it will
        # move the charged residues around.
        constant_charged = [aa for aa in constant_residues if aa in ['D', 'E', 'K', 'R']]
        if constant_charged != []:
            charged_residues = ['D', 'E', 'K', 'R']
            positive_residues = ['K', 'R']
            negative_residues = ['D', 'E']
            
            # Identify which charged residue types can be varied (not in constant_residues)
            variable_positive_residues = [aa for aa in positive_residues if aa not in constant_residues]
            variable_negative_residues = [aa for aa in negative_residues if aa not in constant_residues]
            
            # Find positions of variable charged residues in original sequence
            original_variable_positive_positions = [i for i, aa in enumerate(sequence) 
                                                if aa in variable_positive_residues]
            original_variable_negative_positions = [i for i, aa in enumerate(sequence) 
                                                if aa in variable_negative_residues]
            
            # Find variable charged residues in final sequence
            final_variable_positive_residues = [aa for aa in final_sequence 
                                            if aa in variable_positive_residues]
            final_variable_negative_residues = [aa for aa in final_sequence 
                                            if aa in variable_negative_residues]
            
            # If there are variable charged residues, we need to rearrange them
            if original_variable_positive_positions or original_variable_negative_positions:
                # Remove all variable charged residues from their current positions
                for i in range(len(final_sequence)):
                    if (final_sequence[i] in variable_positive_residues or 
                        final_sequence[i] in variable_negative_residues):
                        final_sequence[i] = None  # Mark for replacement
                
                # Shuffle the collected variable charged residues to randomize their assignment
                np.random.shuffle(final_variable_positive_residues)
                np.random.shuffle(final_variable_negative_residues)
                
                # Place variable positive residues back at original variable positive positions
                for i, pos in enumerate(original_variable_positive_positions):
                    if i < len(final_variable_positive_residues):
                        final_sequence[pos] = final_variable_positive_residues[i]
                
                # Place variable negative residues back at original variable negative positions
                for i, pos in enumerate(original_variable_negative_positions):
                    if i < len(final_variable_negative_residues):
                        final_sequence[pos] = final_variable_negative_residues[i]
                
                # Fill any remaining None positions with non-charged amino acids
                non_charged_from_new = [aa for aa in new_sequence if aa not in charged_residues]
                non_charged_index = 0
                
                for i in range(len(final_sequence)):
                    if final_sequence[i] is None:
                        if non_charged_index < len(non_charged_from_new):
                            final_sequence[i] = non_charged_from_new[non_charged_index]
                            non_charged_index += 1
            # make string
            final_sequence = ''.join(final_sequence)
        else:
            final_sequence = ''.join(final_sequence)
        # check hydropathy. Kappa should be fine. 
        if not check_hydropathy(final_sequence, target_hydropathy, hydropathy_tolerance):
            # if not within hydropathy tolerance, skip this sequence
            continue
        all_sequences.append(final_sequence)

    if all_sequences==[]:
        return None

    return all_sequences




def constant_properties_and_class_sequence(sequence: str,
                                        kappa_tolerance: float = parameters.MAXIMUM_KAPPA_ERROR,
                                        hydropathy_tolerance: float = parameters.MAXIMUM_HYDRO_ERROR) -> str:
    '''
    Generate a new sequence where the fractions of amino acids by class
    are constant but their position is not constrained. The returned
    sequence must also have the same hydropathy and kappa as the input sequence.
    
    Parameters:
    ----------
    sequence : str
        The input protein sequence to be modified.  
    kappa_tolerance : float
        Acceptable difference for kappa optimization.
    hydropathy_tolerance : float
        Acceptable difference for hydropathy optimization.    
    
    Returns:
    --------
    str or None
        If the variant generation fails, None is returned.
        If successful, a new protein sequence with the same hydropathy, kappa, and amino acid class fractions is returned.
    
    '''
    # get starting hydrpathy and kappa
    protein = Protein(sequence)
    target_hydropathy = protein.hydrophobicity
    original_kappa = protein.kappa

    # get class fractions
    class_fractions = calculate_amino_acid_class_fractions(sequence)

    # use the create_sequence_by_class function to generate a new sequence
    new_sequence = create_sequence_by_class(
        length=len(sequence),
        aromatic_fraction=class_fractions['aromatic'],
        aliphatic_fraction=class_fractions['aliphatic'],
        polar_fraction=class_fractions['polar'],
        positive_fraction=class_fractions['positive'],
        negative_fraction=class_fractions['negative'],
        glycine_fraction=class_fractions['glycine'],
        proline_fraction=class_fractions['proline'],
        cysteine_fraction=class_fractions['cysteine'],
        histidine_fraction=class_fractions['histidine'],
        num_sequences=5,
        convert_to_amino_acids=True
    )

    # now use within class hydropathy optimization to get the hydropathy correct.
    for _ in range(2):
        new_sequence = optimize_hydropathy_within_class(
            new_sequence,
            target_hydropathy=target_hydropathy,
            max_iterations=1000,
            tolerance=hydropathy_tolerance,
            only_return_within_tolerance=False
        )[0]
        if check_hydropathy(new_sequence, target_hydropathy, hydropathy_tolerance):
            break

    # make sure new_sequence hydropathy is within tolerance
    if not check_hydropathy(new_sequence, target_hydropathy, hydropathy_tolerance):
            # if still not within tolerance, return None
            return None

    # finally, optimize kappa to get it to the original kappa
    new_sequence = optimize_kappa(
        [new_sequence],
        target_kappa=original_kappa,
        num_iterations=1000,
        convert_input_seq_to_matrix=True,
        inputting_matrix=False,
        tolerance=kappa_tolerance,
        only_return_within_tolerance=True,
        return_when_num_hit=1,
        avoid_shuffle=False,
        num_copies=10
    )
    # if we don't have a sequence within kappa tolerance, return None. 
    if new_sequence is None or len(new_sequence) == 0:
        return None
    
    # Get the first sequence from the list
    new_sequence = new_sequence[0]  
    return new_sequence


def constant_properties_and_class_by_order_sequence(input_sequence: str,
                                    hydropathy_tolerance: float = parameters.MAXIMUM_HYDRO_ERROR) -> str:
    """
    Generate a variant of the input sequence that maintains the same number and order of residues by class.
    
    Parameters:
    ----------
    input_sequence : str
        The input protein sequence to be modified.
    hydropathy_tolerance : float, default=0.05
        Acceptable difference for hydropathy optimization.     

    Returns:
    --------
    str or None
        If the variant generation fails, None is returned.
        If successful, a new protein sequence with the same number and order of residues by class is returned.

    """
    # Convert input sequence to a Protein object
    protein = Protein(input_sequence)
    
    # Get the hydropathy, this is the only thing we will need to change back after modifying residues. 
    target_hydropathy = protein.hydrophobicity
    
    # make version where we swap residues within their classes
    variant_sequence = change_all_residues_within_class(input_sequence)

    # now need to try to get hydropathy correct.
    # first try to avoid original residues
    variant_sequence = optimize_hydropathy_within_class_avoid_original_residues(
        original_sequence=input_sequence,
        variant_sequence=variant_sequence,
        target_hydropathy=target_hydropathy,
        max_iterations=1000,
        tolerance=hydropathy_tolerance
    )

    # see if within tolerance, if not use the within class optimization to get it closer.
    if not check_hydropathy(variant_sequence, target_hydropathy, hydropathy_tolerance):
        variant_sequence = optimize_hydropathy_within_class(
            [variant_sequence],
            target_hydropathy=target_hydropathy,
            max_iterations=1000,
            tolerance=hydropathy_tolerance,
            only_return_within_tolerance=True
        )
        # Extract first sequence from the returned list
        if variant_sequence is not None and len(variant_sequence) > 0:
            variant_sequence = variant_sequence[0]
    # if we didn't get anything back from the optimization, return None   
    if variant_sequence is None or (isinstance(variant_sequence, list) and len(variant_sequence) == 0):
        return None
    
    return variant_sequence
    
    
def change_hydropathy_constant_class_sequence(sequence: str,
                                 target_hydropathy: float,
                                 hydropathy_tolerance: float = parameters.MAXIMUM_HYDRO_ERROR):
    '''
    Generate a variant with altered hydropathy while keeping the same order
    and number of residues by class.
    
    Parameters:
    ----------
    sequence : str
        The input protein sequence to be modified.
    target_hydropathy : float
        Target hydropathy value to achieve.
    hydropathy_tolerance : float, default=0.05
        Acceptable difference for hydropathy optimization.
    
    Returns:
    --------
    str
        Protein sequence with altered hydropathy while maintaining class structure. 
    '''
    # check range that is possible
    hydro_range = find_hydro_range_constant_class(sequence)
    if target_hydropathy < hydro_range[0] or target_hydropathy > hydro_range[1]:
        raise ValueError(f"Target hydropathy {target_hydropathy} is outside the possible range {hydro_range} for this sequence.")
    
    # optimize the sequence using the within class hydropathy optimization
    optimized_sequence = optimize_hydropathy_within_class(
        [sequence],
        target_hydropathy=target_hydropathy,
        max_iterations=5000,
        hydropathy_tolerance=hydropathy_tolerance,
        only_return_within_tolerance=True
    )
    if optimized_sequence is None or len(optimized_sequence) == 0:
        return None
    optimized_sequence = optimized_sequence[0]  # Get the first sequence from the list
    
    return optimized_sequence





def change_fcr_minimize_class_changes_sequence(sequence: str,
                               target_FCR: float,
                               target_NCPR: float = None,
                               hydropathy_tolerance: float = parameters.MAXIMUM_HYDRO_ERROR,
                               kappa_tolerance: float = parameters.MAXIMUM_KAPPA_ERROR) -> str:
    '''
    Generate a variant with altered FCR while keeping the same order
    and number of residues by class as much as is possible.
    
    Parameters:
    ----------
    sequence : str
        The input protein sequence to be modified.
    target_FCR : float
        Target FCR value to achieve.
    target_NCPR : float, optional
        Target NCPR value to achieve. If None, will use the original NCPR.
    hydropathy_tolerance : float, default=0.05
        Acceptable difference for hydropathy optimization.
    kappa_tolerance : float, default=0.03
        Acceptable difference for kappa optimization.
    
    Returns:
    --------
    str
        Protein sequence with altered FCR while maintaining class structure. 
    '''
    # set order to change residues in sequence that minimizes changes in major chemistries
    order_to_change = ['Q', 'N', 'S', 'T', 'G', 'A', 'M', 'L', 'V', 'I', 'P', 'C', 'H', 'Y', 'W', 'F']

    # get original NCPR
    if target_NCPR is None:
        starting_NCPR = Protein(sequence).NCPR
    else:
        starting_NCPR = target_NCPR
    starting_hydropathy = Protein(sequence).hydrophobicity
    starting_kappa = Protein(sequence).kappa

    # if target FCR is correct and target ncpr is correct just return the sequence
    if round(target_FCR,8) == Protein(sequence).FCR and round(starting_NCPR,8) == Protein(sequence).NCPR:
        return sequence

    # prioritize FCR.
    if round(target_FCR,8) == 0:
        starting_NCPR=0

    # get num residues needed for fcr and ncpr
    charged_residues = needed_charged_residues(
        length=len(sequence),
        fraction=target_FCR,
        net_charge=starting_NCPR
    )

    # get current number of charged residues
    current_positive = sum(1 for aa in sequence if aa in ['K', 'R'])
    current_negative = sum(1 for aa in sequence if aa in ['D', 'E'])

    # calculate how many residues we need to add or remove
    num_positive_to_add = charged_residues['positive'] - current_positive
    num_negative_to_add = charged_residues['negative'] - current_negative

    # if we need to add positive and remove negative or remove negative and add positive, we will first swap charged residues to try to avoid changing anything else. 
    if num_positive_to_add > 0 and num_negative_to_add < 0:
        # we need to add positive residues and remove negative residues
        # first, swap positive and negative residues
        for _ in range(abs(num_negative_to_add)):
            indices = [i for i, x in enumerate(sequence) if x in ['D', 'E']]
            if indices:
                index_to_change = np.random.choice(indices)
                sequence = sequence[:index_to_change] + np.random.choice(['K', 'R']) + sequence[index_to_change + 1:]
                num_positive_to_add -= 1
                num_negative_to_add += 1

    elif num_negative_to_add > 0 and num_positive_to_add < 0:
        # we need to add negative residues and remove positive residues
        # first, swap positive and negative residues
        for _ in range(abs(num_positive_to_add)):
            indices = [i for i, x in enumerate(sequence) if x in ['K', 'R']]
            if indices:
                index_to_change = np.random.choice(indices)
                sequence = sequence[:index_to_change] + np.random.choice(['D', 'E']) + sequence[index_to_change + 1:]
                num_negative_to_add -= 1
                num_positive_to_add += 1

    # now we can add or remove residues as needed
    # if we need to add or remove residues, we will do so in the order specified
    for _ in range(abs(num_positive_to_add)):
        if num_positive_to_add > 0:
            for aa in order_to_change:
                indices = [i for i, x in enumerate(sequence) if x == aa]
                if indices:
                    index_to_change = np.random.choice(indices)
                    sequence = sequence[:index_to_change] + np.random.choice(['K', 'R']) + sequence[index_to_change + 1:]
                    break   
        else:
            # if we need to remove positive residues, we will remove them in the order specified
            for aa in ['K', 'R']:
                indices = [i for i, x in enumerate(sequence) if x == aa]
                if indices:
                    index_to_change = np.random.choice(indices)
                    sequence = sequence[:index_to_change] + np.random.choice(['Q', 'N', 'S', 'T']) + sequence[index_to_change + 1:]
                    break

    # if we need to add negative residues, we will add them in the order specified
    for _ in range(abs(num_negative_to_add)):
        if num_negative_to_add > 0:
            for aa in order_to_change:
                indices = [i for i, x in enumerate(sequence) if x == aa]
                if indices:
                    index_to_change = np.random.choice(indices)
                    sequence = sequence[:index_to_change] + np.random.choice(['D', 'E']) + sequence[index_to_change + 1:]
                    break
        else:
            # if we need to remove negative residues, we will remove them in the order specified
            for aa in ['D', 'E']:
                indices = [i for i, x in enumerate(sequence) if x == aa]
                if indices:
                    index_to_change = np.random.choice(indices)
                    sequence = sequence[:index_to_change] + np.random.choice(['Q', 'N', 'S', 'T']) + sequence[index_to_change + 1:]
                    break

    # now optimize the hydropathy within class
    sequence = optimize_hydropathy_within_class(
        [sequence],
        target_hydropathy=starting_hydropathy,
        tolerance=hydropathy_tolerance,
        only_return_within_tolerance=False,
        max_iterations=5000
    )
    if sequence is None or len(sequence) == 0:
        return None
    sequence = sequence[0]  # Get the first sequence from the list
    
    # make sure we can have a non negative kappa value. 
    if Protein(sequence).FCR != 0:
        if round(abs(Protein(sequence).NCPR),8) != round(Protein(sequence).FCR,8):
            # finally, optimize kappa to get it to the original kappa
            sequence = optimize_kappa(
                [sequence],
                target_kappa=starting_kappa,
                num_iterations=5000,
                convert_input_seq_to_matrix=True,
                inputting_matrix=False,
                tolerance=kappa_tolerance,
                return_when_num_hit=1,
                avoid_shuffle=True,
                num_copies=10
            )
            if sequence is None or len(sequence) == 0:
                return None
            sequence = sequence[0]

    return sequence


def change_ncpr_constant_class_sequence(sequence: str,
                                target_NCPR: float,
                               hydropathy_tolerance: float = parameters.MAXIMUM_HYDRO_ERROR,
                               kappa_tolerance: float = parameters.MAXIMUM_KAPPA_ERROR) -> str:
    '''
    Generate a variant with altered NCPR while keeping the same order
    and number of residues by class as much as is possible.
    
    Parameters:
    ----------
    sequence : str
        The input protein sequence to be modified.
    target_NCPR : float
        Target NCPR value to achieve.
    hydropathy_tolerance : float, default=0.05
        Acceptable difference for hydropathy optimization.
    kappa_tolerance : float, default=0.03
        Acceptable difference for kappa optimization.
    
    Returns:
    --------
    str
        Protein sequence with altered NCPR while maintaining class structure. 
    '''
    starting_FCR = Protein(sequence).FCR
    seq=change_fcr_minimize_class_changes_sequence(sequence,
                               target_FCR=starting_FCR,
                               target_NCPR=target_NCPR,
                               hydropathy_tolerance=hydropathy_tolerance,
                               kappa_tolerance=kappa_tolerance)
    return seq

def change_kappa_sequence(sequence,
                           target_kappa,
                           kappa_tolerance=parameters.MAXIMUM_KAPPA_ERROR):
    '''
    Generate a sequence with a different kappa value. 

    Parameters:
    ----------
    sequence : str
        The input protein sequence to be modified.
    target_kappa : float
        Target kappa value to achieve.
    kappa_tolerance : float, default=0.03
        Acceptable difference for kappa optimization.
    '''
    # get FCR and NCPR.
    protein = Protein(sequence)
    starting_FCR = protein.FCR
    starting_NCPR = protein.NCPR
    if round(abs(starting_NCPR),8)==round(starting_FCR,8):
        return sequence
    # use kappa optimzer
    sequence = optimize_kappa(
        [sequence],
        target_kappa=target_kappa,
        num_iterations=5000,
        convert_input_seq_to_matrix=True,
        inputting_matrix=False,
        tolerance=kappa_tolerance,
        return_when_num_hit=1,
        avoid_shuffle=True,
        num_copies=10
    )
    if sequence is None or len(sequence) == 0:
        return None
    sequence = sequence[0]
    return sequence

def change_properties_minimize_differences_sequence(input_sequence: str,
                             target_hydropathy: float = None,
                             target_kappa: float = None,
                             target_FCR: float = None,
                             target_NCPR: float = None,
                             max_iterations: int = 1000,
                             hydropathy_tolerance: float = parameters.MAXIMUM_HYDRO_ERROR,
                             kappa_tolerance: float = parameters.MAXIMUM_KAPPA_ERROR) -> str:
    """
    Generate a minimal variant of the input sequence that optimizes hydropathy and kappa.
    
    Parameters:
    ----------
    input_sequence : str
        The input protein sequence to be modified.
    target_hydropathy : float, optional
        Target hydropathy value to achieve.
    target_kappa : float, optional
        Target kappa value to achieve.
    target_FCR : float, optional
        Target FCR value to achieve.
    target_NCPR : float, optional
        Target NCPR value to achieve.
    max_iterations : int, default=1000
        Maximum number of optimization iterations.
    hydropathy_tolerance : float, default=0.07
        set by parameters.MAXIMUM_HYDRO_ERROR
    kappa_tolerance : float, default=0.03
        set by parameters.MAXIMUM_KAPPA_ERROR
        Acceptable difference for kappa optimization.

    Returns:
    --------
    str or None
        If the variant generation fails, None is returned.
        If successful, a new protein sequence with the specified properties is returned.
    """
    # steps will be as follows:
    # 1. Calculate initial properties
    # 2. Modify FCR if specified. 
    # 3. Modify NCPR if specified.
    # 4. Optimize hydropathy if specified.
    # 5. Optimize kappa if specified.
    
    # Step 1: Calculate initial properties
    protein = Protein(input_sequence)
    initial_hydropathy = protein.hydrophobicity
    initial_kappa = protein.kappa
    initial_FCR = protein.FCR
    if target_kappa is None:
        target_kappa = initial_kappa


    # step 2. Modify FCR. We want to change as few residues as possible. 
    # if hydropathy is specified, we can use the charged residues to get 
    # the hydropathy closer to the target. Otherwise we want to change
    # residues that minimize the change in hydropathy.
    if target_FCR is not None:
        # determine the number of residues we need to change to get to target FCR.
        fcr_diff = target_FCR - initial_FCR
        if fcr_diff > 0:
            # if fcr_diff is positive, we need to add more charged residues
            # check what hydropathy change we need to make. 
            hydropathy_diff = target_hydropathy - initial_hydropathy if target_hydropathy is not None else 0

            # choose list of amino acids to change based on hydropathy and chemistry. 
            if hydropathy_diff > 0:
                order_to_change = ['Q', 'N', 'S', 'T', 'G', 'A', 'M', 'L', 'V', 'I', 'P', 'C', 'H', 'Y', 'W', 'F']
            else:
                order_to_change = ['I', 'V', 'L', 'M', 'A','G', 'T', 'S', 'N', 'C', 'P', 'Q','F', 'W', 'Y', 'H']
            # calculate how many residues we need to add
            num_to_add = round(fcr_diff*len(input_sequence))

            # loop over residues to change until we reach the target FCR
            order_to_change = ['Q', 'N', 'S', 'T', 'G', 'A', 'M', 'L', 'V', 'I', 'P', 'C', 'H', 'Y', 'W', 'F']
            for num in range(num_to_add):
                # if num is even, add a 'D' or 'E' else add a 'K' or 'R'
                if num % 2 == 0:
                    # add a negatively charged residue
                    aa_to_add = np.random.choice(['D', 'E'])
                else:
                    # add a positively charged residue
                    aa_to_add = np.random.choice(['K', 'R'])
                # iterate over groups in order_to_change
                for aa in order_to_change:
                    # get occurances of the aa in the sequnece
                    indices = [i for i, x in enumerate(input_sequence) if x == aa]
                    # randomly select one of the indices
                    if indices:
                        index_to_change = np.random.choice(indices)
                        # change the residue at that index
                        input_sequence = input_sequence[:index_to_change] + aa_to_add + input_sequence[index_to_change + 1:]
                        # now that we've added the residue for this round, we need to get out of the loop.
                        break
        elif fcr_diff < 0:
            # if fcr_diff is negative, we need to remove charged residues
            # check what hydropathy change we need to make.
            hydropathy_diff = target_hydropathy - initial_hydropathy if target_hydropathy is not None else 0
            # make list of replacement amino acids based on hydropathy and chemistry
            charge_replacements = ['S', 'T', 'N', 'Q']
            
            # calculate how many residues we need to remove
            num_to_remove = round(abs(fcr_diff)*len(input_sequence))
            # loop over residues to change until we reach the target FCR
            for num in range(num_to_remove):
                # get indices of all charged residues in the sequence
                indices = [i for i, x in enumerate(input_sequence) if x in ['D', 'E', 'K', 'R']]
                # choose a random index to change
                if indices:
                    index_to_change = np.random.choice(indices)
                    # randomly select a replacement amino acid from charge_replacements
                    aa_to_add = np.random.choice(charge_replacements)
                    # change the residue at that index
                    input_sequence = input_sequence[:index_to_change] + aa_to_add + input_sequence[index_to_change + 1:]
        else:
            # if fcr_diff is 0, we don't need to change anything
            pass
                

    # Step 3: Modify NCPR
    if target_NCPR is not None:
        # see if we need to use target FCR or the protein's FCR
        if target_FCR is None:
            use_FCR=initial_FCR
        else:
            use_FCR = target_FCR
        # get needed number of negative and postitive residues.
        needed_residues = needed_charged_residues(len(input_sequence), use_FCR, target_NCPR)
        needed_positive = needed_residues['positive']
        needed_negative = needed_residues['negative']
        # get current number of positive and negative residues
        current_positive = sum(1 for aa in input_sequence if aa in ['K', 'R'])
        current_negative = sum(1 for aa in input_sequence if aa in ['D', 'E'])
        # convert needed number of positive to negative or vice versa
        if needed_positive > current_positive:
            # we need to add more positive residues
            num_to_add = needed_positive - current_positive
            for _ in range(num_to_add):
                # randomly select a position to change
                indices = [i for i, x in enumerate(input_sequence) if x in ['D', 'E']]
                if indices:
                    index_to_change = np.random.choice(indices)
                    input_sequence = input_sequence[:index_to_change] + np.random.choice(['K', 'R']) + input_sequence[index_to_change + 1:]
        elif needed_negative > current_negative:
            # we need to add more negative residues
            num_to_add = needed_negative - current_negative
            for _ in range(num_to_add):
                # randomly select a position to change
                indices = [i for i, x in enumerate(input_sequence) if x in ['K', 'R']]
                if indices:
                    index_to_change = np.random.choice(indices)
                    input_sequence = input_sequence[:index_to_change] + np.random.choice(['D', 'E']) + input_sequence[index_to_change + 1:]

    # Step 4: Optimize hydropathy. First do within class to minimize changes to chemistry. 
    if target_hydropathy is not None:
        input_sequence = optimize_hydropathy_within_class(
            [input_sequence],
            target_hydropathy=target_hydropathy,
            max_iterations=max_iterations,
            tolerance=hydropathy_tolerance,
            only_return_within_tolerance=False
        )[0]
        
        # check hydropathy
        if not check_hydropathy(input_sequence, target_hydropathy, hydropathy_tolerance):
            # This will try to find the minimum number of changes needed to get to the target hydropathy.
            input_sequence = optimize_hydropathy_minimal_changes(
                input_sequence,
                target_hydropathy=target_hydropathy,
                max_iterations=max_iterations,
                tolerance=hydropathy_tolerance,
                preserve_charged=True
            )
        
        # check hydropathy again after minimal changes
        if not check_hydropathy(input_sequence, target_hydropathy, hydropathy_tolerance):
           # use optimize_hydropathy function
            input_sequence = optimize_hydropathy(
                [input_sequence],
                target_hydropathy=target_hydropathy,
                max_iterations=max_iterations,
                tolerance=hydropathy_tolerance,
                preserve_charged=True,
                convert_input_seq_to_matrix=True,
                return_when_num_hit=1,
                avoid_shuffle=True,
                num_copies=10,
                convert_back_to_sequences= True,
                need_to_convert_seqs=True
            )
            if input_sequence is None:
                return None
            input_sequence = input_sequence[0]  # Get the first sequence from the list

    # Step 5: Optimize kappa
    if target_kappa is not None:
        input_sequence = optimize_kappa(
            [input_sequence],
            target_kappa=target_kappa,
            num_iterations=max_iterations,
            tolerance=kappa_tolerance,
            inputting_matrix=False,
            convert_input_seq_to_matrix=True,
            return_when_num_hit=1,
            avoid_shuffle=True,
            num_copies=10
        )
        if input_sequence is None:
            return None
        input_sequence = input_sequence[0]
    
    return input_sequence




def change_any_properties_sequence(sequence,
                                 target_FCR: float = None,
                                 target_NCPR: float = None,
                                 target_hydropathy: float = None,
                                 target_kappa: float = None,
                                 hydropathy_tolerance: float = parameters.MAXIMUM_HYDRO_ERROR,
                                 kappa_tolerance: float = parameters.MAXIMUM_KAPPA_ERROR,
                                 max_iterations: int = 5000) -> str:
    '''
    Generate a variant with altered FCR, NCPR, hydropathy, and kappa while keeping the same order
    and number of residues by class as much as is possible.
    
    Parameters:
    ----------
    sequence : str
        The input protein sequence to be modified.
    target_FCR : float, optional
        Target FCR value to achieve. If None, will not change FCR.
    target_NCPR : float, optional
        Target NCPR value to achieve. If None, will not change NCPR.
    target_hydropathy : float, optional
        Target hydropathy value to achieve. If None, will not change hydropathy.
    target_kappa : float, optional
        Target kappa value to achieve. If None, will not change kappa.
    hydropathy_tolerance : float, default=0.05
        Acceptable difference for hydropathy optimization.
    kappa_tolerance : float, default=0.03
        Acceptable difference for kappa optimization.
    max_iterations : int, default=5000
        Maximum number of iterations for optimization algorithms.
    
    Returns:
    --------
    str
        Protein sequence with altered properties while maintaining class structure. 
    '''
    
    # Step 1: Calculate initial properties
    initial_protein = Protein(sequence)
    if target_FCR is None:
        target_FCR = initial_protein.FCR
    if target_NCPR is None:
        target_NCPR = initial_protein.NCPR

    # use fcr variant function.
    sequence = change_fcr_minimize_class_changes_sequence(sequence,
                                           target_FCR=target_FCR,
                                           target_NCPR=target_NCPR,
                                           hydropathy_tolerance=hydropathy_tolerance,
                                           kappa_tolerance=kappa_tolerance)
    
    # get current hydropathy
    current_hydropathy = Protein(sequence).hydrophobicity
    if target_hydropathy is None:
        target_hydropathy = initial_protein.hydrophobicity
    if abs(current_hydropathy-target_hydropathy) > hydropathy_tolerance:
        # Modify hydropathy using withn class optimization
        sequence = optimize_hydropathy_within_class(
            [sequence],
            target_hydropathy=target_hydropathy,
            max_iterations=max_iterations,
            tolerance=hydropathy_tolerance,
            only_return_within_tolerance=False
        )
    sequence = sequence[0]  # Get the first sequence from the list
    if check_hydropathy(sequence, target_hydropathy, hydropathy_tolerance) is False:
        return None
    
    # now optimize kappa
    # get current kappa
    current_kappa = Protein(sequence).kappa
    if target_kappa is None:
        target_kappa = initial_protein.kappa

    if abs(current_kappa-target_kappa) > kappa_tolerance:
        # make sure we have charged residues otherwise no reason to optimzie kappa
        if Protein(sequence).FCR != 0:
            # make sure sequence has a net charge not equal to fraction charge.
            if round(abs(Protein(sequence).NCPR),8) != round(Protein(sequence).FCR,8):
                # optimize kappa
                # use kappa optimzer
                sequence = optimize_kappa(
                    [sequence],
                    target_kappa=target_kappa,
                    tolerance=kappa_tolerance,
                    return_when_num_hit=1,
                    inputting_matrix=False,
                    convert_input_seq_to_matrix=True,
                    num_copies=10
                )[0]
    # return sequence
    return sequence


def change_dimensions_sequence(sequence, increase_or_decrease, rg_or_re,
                                num_attempts=5, allowed_error=parameters.MAXIMUM_RG_RE_ERROR,
                                reduce_pos_charged=False, exclude_aas=None):
    '''
    Generate a variant of the input sequence with a modified radius of gyration (Rg).

    Parameters
    ----------
    sequence : str
        The input amino acid sequence to modify
    increase_or_decrease : {'increase', 'decrease'}
        Whether to increase or decrease the Rg
    rg_or_re : {'rg', 're'}
        Whether to optimize for radius of gyration (Rg) or radius of elongation (Re)
    num_attempts : int, default=5
        Number of attempts to generate a valid variant
    allowed_error : float, optional
        If specified, the maximum allowed error for the generated Rg or Re compared to the original value
    reduce_pos_charged : bool, default=False
        Whether to reduce positively charged residues (K, R) in the sequence.
    exclude_aas : list, optional
        A list of amino acids to exclude from the optimization process.

    Returns
    -------
    str or list of str
        The generated variant(s) of the input sequence
    '''
    # set number of intervals based on sequence length. 
    if len(sequence)>40:
        num_intervals = 40
    else:
        num_intervals = 20

    if rg_or_re=='re':
        # If we are optimizing for radius of elongation (Re), we need to adjust the intervals accordingly
        dim_intervals = np.linspace(get_min_re(len(sequence)), get_max_re(len(sequence)), num_attempts)
        individual_predictor = predict_re
        batch_predictor = batch_predict_re
    elif rg_or_re=='rg':
        # If we are optimizing for radius of gyration (Rg), we use the standard intervals
        dim_intervals = np.linspace(get_min_rg(len(sequence)), get_max_rg(len(sequence)), num_intervals)
        individual_predictor = predict_rg
        batch_predictor = batch_predict_rg
    else:
        raise ValueError("rg_or_re must be 'rg' or 're'")

    # get the original Rg
    original_dims = individual_predictor(sequence)

    # get intervals greater than original rg
    if increase_or_decrease == 'increase':
        target_intervals = dim_intervals[dim_intervals > original_dims]
    elif increase_or_decrease == 'decrease':
        target_intervals = dim_intervals[dim_intervals < original_dims]
    else:
        raise ValueError("increase_or_decrease must be 'increase' or 'decrease'")
    
    # filter out target_intervals that are within the allowed error range
    target_intervals = target_intervals[np.abs(target_intervals - original_dims) >= allowed_error]

    # Check if target_intervals is empty
    if len(target_intervals) == 0:
        if increase_or_decrease == 'increase':
            closest_interval = [original_dims + (allowed_error+1)]
        elif increase_or_decrease == 'decrease':
            closest_interval = [original_dims - (allowed_error+1)]
    else:
        # now we need to find the closest interval
        closest_interval = target_intervals[np.argmin(np.abs(target_intervals - original_dims))]


    # Initialize variables to avoid scope issues
    optimized_sequence = None
    
    # attempt to make sequence
    for _ in range(num_attempts):
        # now use optimize_seq_dims
        optimized_sequence = optimize_seq_dims(sequence, rg_or_re, closest_interval, 
                                               allowed_error=allowed_error,
                                                return_all_sequences=True,
                                                reduce_pos_charged=reduce_pos_charged,
                                                exclude_aas=exclude_aas)

        if optimized_sequence is None:
            continue  # If optimization failed, skip to the next attempt

        # now get all sequences that are greater than original dimensions if increase or less than if decrease
        if increase_or_decrease == 'increase':
            # predict the dimensions of the optimized sequence
            potential_sequences = batch_predictor(optimized_sequence, return_seq2prediction=True)
            # get sequences greater than the original dimensions
            potential_sequences = {seq: rg for seq, rg in potential_sequences.items() if rg >= original_dims}
        elif increase_or_decrease == 'decrease':
            # predict the dimensions of the optimized sequence
            potential_sequences = batch_predictor(optimized_sequence, return_seq2prediction=True)
            # get sequences less than the original dimensions
            potential_sequences = {seq: rg for seq, rg in potential_sequences.items() if rg <= original_dims}
        
        # if we have no sequences, continue to the next attempt
        if len(potential_sequences) == 0:
            continue

    if len(potential_sequences) > 0:
        return list(potential_sequences.keys())
    return None




