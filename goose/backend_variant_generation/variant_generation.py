'''
Updated code for variant generation. Going to rewrite a lot 
of the functionality to use numpy vectorized operations for efficiency.
'''

from goose.backend_vectorized.optimize_kappa_vectorized import optimize_kappa_vectorized
from goose.backend_vectorized.create_sequence_by_class import create_sequence_by_class
from goose.backend_vectorized.seq_by_probability_vectorized import SequenceGenerator

from goose.backend_variants.optimize_hydropathy_minimal_changes import optimize_hydropathy_minimal_changes
from goose.backend_variants.within_class_hydropathy_optimization import optimize_hydropathy_within_class_vectorized 
from goose.backend_variants.variant_generation_helper_functions import needed_charged_residues, change_all_residues_within_class, optimize_hydropathy_avoid_original_residues, calculate_amino_acid_class_fractions

from sparrow.protein import Protein
import numpy as np


def generate_constant_class_variant(input_sequence: str) -> str:
    """
    Generate a variant of the input sequence that maintains the same number and order of residues by class.
    
    Parameters:
    ----------
    input_sequence : str
        The input protein sequence to be modified.

    Returns:
    --------
    str
        Protein sequence with residues changed within their classes.
    """
    # Convert input sequence to a Protein object
    protein = Protein(input_sequence)
    
    # Get the hydropathy, this is the only thing we will need to change back after modifying residues. 
    original_hydropathy = protein.hydrophobicity
    
    # make version where we swap residues within their classes
    variant_sequence = change_all_residues_within_class(input_sequence)

    # now need to try to get hydropathy correct.
    # first try to avoid original residues
    variant_sequence = optimize_hydropathy_avoid_original_residues(
        original_sequence=input_sequence,
        variant_sequence=variant_sequence,
        target_hydropathy=original_hydropathy,
        max_iterations=1000,
        tolerance=0.05
    )

    # see if within tolerance, if not use the within class optimization to get it closer.
    final_hydropathy = Protein(variant_sequence).hydrophobicity
    if abs(final_hydropathy - original_hydropathy) > 0.05:
        variant_sequence = optimize_hydropathy_within_class_vectorized(
            [variant_sequence],
            target_hydropathy=original_hydropathy,
            max_iterations=1000,
            tolerance=0.05,
            only_return_within_tolerance=True
        )[0]
    
    return variant_sequence
    


def generate_minimal_variant(input_sequence: str,
                             target_hydropathy: float = None,
                             target_kappa: float = None,
                             target_FCR: float = None,
                             target_NCPR: float = None,
                             max_iterations: int = 1000,
                             hydropathy_tolerance: float = 0.05,
                             kappa_tolerance: float = 0.05) -> str:
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
    hydropathy_tolerance : float, default=0.05
        Acceptable difference for hydropathy optimization.
    kappa_tolerance : float, default=0.05
        Acceptable difference for kappa optimization.

    Returns:
    --------
    str
        Optimized protein sequence with minimal changes.
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
    # Step 3: Modify NCPR
    if target_NCPR is not None:
        # calculate NCPR after modifying FCR
        initial_NCPR = protein.NCPR
        # recalculate FCR
        initial_FCR = protein.FCR
        # get needed number of negative and postitive residues.
        needed_residues = needed_charged_residues(len(input_sequence), target_NCPR, initial_NCPR)
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
        input_sequence = optimize_hydropathy_within_class_vectorized(
            [input_sequence],
            target_hydropathy=target_hydropathy,
            max_iterations=max_iterations,
            tolerance=hydropathy_tolerance,
            only_return_within_tolerance=False
        )[0]

    # check if the sequence is within the tolerance of the target hydropathy. If not, take additional steps.
    if target_hydropathy is not None:
        final_hydropathy = Protein(input_sequence).hydrophobicity
        if abs(final_hydropathy - target_hydropathy) > hydropathy_tolerance:
            # This will try to find the minimum number of changes needed to get to the target hydropathy.
            input_sequence = optimize_hydropathy_minimal_changes(
                input_sequence,
                target_hydropathy=target_hydropathy,
                max_iterations=max_iterations,
                tolerance=hydropathy_tolerance,
                preserve_charged=True
            )

    # Step 5: Optimize kappa
    if target_kappa is not None:
        input_sequence = optimize_kappa_vectorized(
            [input_sequence],
            target_kappa=target_kappa,
            max_iterations=max_iterations,
            tolerance=kappa_tolerance
        )[0]

    return input_sequence


def generate_region_shuffle_variant(input_sequence: str,
                           shuffle_regions: list) -> str:
    """
    Generate a variant of the input sequence by shuffling specified regions.
    
    Parameters:
    ----------
    input_sequence : str
        The input protein sequence to be modified.
    shuffle_regions : list of tuples
        List of tuples where each tuple contains start and end indices of regions to shuffle.

    Returns:
    --------
    str
        Protein sequence with specified regions shuffled.
    """
    # Shuffle specified regions
    for start, end in shuffle_regions:
        region = list(input_sequence[start:end])
        shuffled_region = ''.join(np.random.permutation(region))
        input_sequence = input_sequence[:start] + shuffled_region + input_sequence[end:]

    return input_sequence

def generate_excluded_shuffle_variant(input_sequence: str,
                             exclude_indices: list) -> str:
    """
    Generate a variant of the input sequence by shuffling residues not specified in exclude_indices.
    
    Parameters:
    ----------
    input_sequence : str
        The input protein sequence to be modified.
    exclude_indices : list
        List of indices to exclude from shuffling.

    Returns:
    --------
    str
        Protein sequence with residues shuffled except those at excluded indices.
    """
    # Convert input sequence to a list for easier manipulation
    seq_list = list(input_sequence)
    
    # Get indices that are not excluded
    indices_to_shuffle = [i for i in range(len(seq_list)) if i not in exclude_indices]
    
    # Shuffle the residues at the non-excluded indices
    shuffled_region = np.random.permutation([seq_list[i] for i in indices_to_shuffle])
    
    # Place shuffled residues back into the sequence
    for idx, new_aa in zip(indices_to_shuffle, shuffled_region):
        seq_list[idx] = new_aa
    
    return ''.join(seq_list)

def generate_new_seq_constant_class_var(sequence: str,
                                        kappa_tolerance: float,
                                        hydropathy_tolerance: float) -> str:
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
    str
        Protein sequence with residues changed within their classes.
    
    '''
    # get starting hydrpathy and kappa
    protein = Protein(sequence)
    original_hydropathy = protein.hydrophobicity
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
        num_sequences=1
    )

    # now use within class hydropathy optimization to get the hydropathy correct.
    new_sequence = optimize_hydropathy_within_class_vectorized(
        new_sequence,
        target_hydropathy=original_hydropathy,
        max_iterations=1000,
        tolerance=hydropathy_tolerance,
        only_return_within_tolerance=False
    )[0]

    # finally, optimize kappa to get it to the original kappa
    new_sequence = optimize_kappa_vectorized(
        [new_sequence],
        target_kappa=original_kappa,
        max_iterations=1000,
        tolerance=kappa_tolerance
    )[0]

    return new_sequence

def generate_constant_properties_var(sequence: str,
                                     hydropathy_tolerance=0.05,
                                     kappa_tolerance=0.03) -> str:
    '''
    Generate a new sequence with the same hydropathy, 
    FCR, NCPR, and kappa as the input sequence.
    
    Parameters:
    ----------
    sequence : str
        The input protein sequence to be modified.
    
    Returns:
    --------
    str
        A new protein sequence with the same hydropathy, FCR, NCPR, and kappa as the input sequence.
    '''
    # get starting properties
    protein = Protein(sequence)
    original_hydropathy = protein.hydrophobicity
    original_kappa = protein.kappa
    original_FCR = protein.FCR
    original_NCPR = protein.NCPR
    length = len(sequence)

    # use the SequenceGenerator
    # initialize the sequence generators
    seq_gen = SequenceGenerator( 
                 use_weighted_probabilities=None)
    
    # iterate over the number of attempts
    seqs = seq_gen.generate_sequences_vectorized(length, fcr=original_FCR, ncpr=original_NCPR, 
                                                    hydropathy=original_hydropathy, num_sequences=100,
                                                    specific_probabilities=None,
                                                    exclude_residues=None)
    
    # Calculate the hydropathy for the sequences and get the ones within the tolerance
    hydropathies = np.array([Protein(seq).hydrophobicity for seq in seqs])
    within_tolerance = np.abs(hydropathies - original_hydropathy) <= hydropathy_tolerance
    seqs_within_tolerance = seqs[within_tolerance]

    if seqs_within_tolerance.size !=0:
        # get ones within kappa tolerance
        kappas = np.array([Protein(seq).kappa for seq in seqs_within_tolerance])
        within_kappa_tolerance = np.abs(kappas - original_kappa) <= kappa_tolerance
        seqs_within_kappa_tolerance = seqs_within_tolerance[within_kappa_tolerance]
        if seqs_within_kappa_tolerance.size != 0:
            # return the first sequence that meets the criteria
            return seqs_within_kappa_tolerance[0]
   
    # if no sequence was found, return None
    return None

def generate_constant_residue_var(sequence: str,
                                  constant_residues: list,
                                  hydropathy_tolerance=0.05,
                                  kappa_tolerance=0.03) -> str:
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
        Acceptable difference for hydropathy optimization.
    kappa_tolerance : float, default=0.03
        Acceptable difference for kappa optimization.
    
    Returns:
    --------      
    str
        A new protein sequence with the same hydropathy, FCR, NCPR, and kappa as the input sequence,
        while keeping the specified residues constant.
    '''
    # make sequence without specified residues
    modified_sequence = ''.join([aa for aa in sequence if aa not in constant_residues])
    
    # Handle edge case where all residues are constant
    if len(modified_sequence) == 0:
        return sequence  # Return original sequence if nothing to vary
    
    # get constant_residue_indices
    constant_residue_indices = [i for i, aa in enumerate(sequence) if aa in constant_residues]
    # get properties for the modified sequence
    protein = Protein(modified_sequence)
    original_hydropathy = protein.hydrophobicity
    original_kappa = protein.kappa
    original_FCR = protein.FCR
    original_NCPR = protein.NCPR
    length = len(modified_sequence)
    # use the SequenceGenerator to generate a new sequence  
    seq_gen = SequenceGenerator(use_weighted_probabilities=None)
    # iterate over the number of attempts
    seqs = seq_gen.generate_sequences_vectorized(length, fcr=original_FCR, ncpr=original_NCPR, 
                                                    hydropathy=original_hydropathy, num_sequences=100,
                                                    specific_probabilities=None,
                                                    exclude_residues=constant_residues)
    # Calculate the hydropathy for the sequences and get the ones within the tolerance
    hydropathies = np.array([Protein(seq).hydrophobicity for seq in seqs])
    within_tolerance = np.abs(hydropathies - original_hydropathy) <= hydropathy_tolerance
    seqs_within_tolerance = np.array(seqs)[within_tolerance] 

    if len(seqs_within_tolerance) > 0:
        # get ones within kappa tolerance
        kappas = np.array([Protein(seq).kappa for seq in seqs_within_tolerance])
        within_kappa_tolerance = np.abs(kappas - original_kappa) <= kappa_tolerance
        seqs_within_kappa_tolerance = seqs_within_tolerance[within_kappa_tolerance]
        if len(seqs_within_kappa_tolerance) > 0:
            # return the first sequence that meets the criteria
            new_sequence = seqs_within_kappa_tolerance[0]
            
            # Create final sequence by inserting constant residues back at their original positions
            final_sequence = list(sequence)  # Start with original sequence as template
            new_seq_index = 0  # Index for the new sequence (without constant residues)
            
            # Iterate through each position in the original sequence
            for i, original_aa in enumerate(sequence):
                if original_aa not in constant_residues:
                    # Replace with amino acid from new sequence
                    final_sequence[i] = new_sequence[new_seq_index]
                    new_seq_index += 1
                # If it's a constant residue, keep the original amino acid (already in final_sequence)
            
            return ''.join(final_sequence)
        else:
            # No sequences within kappa tolerance, return best hydropathy match
            new_sequence = seqs_within_tolerance[0]
            
            # Create final sequence by inserting constant residues back at their original positions
            final_sequence = list(sequence)  # Start with original sequence as template
            new_seq_index = 0  # Index for the new sequence (without constant residues)
            
            # Iterate through each position in the original sequence
            for i, original_aa in enumerate(sequence):
                if original_aa not in constant_residues:
                    # Replace with amino acid from new sequence
                    final_sequence[i] = new_sequence[new_seq_index]
                    new_seq_index += 1
                # If it's a constant residue, keep the original amino acid (already in final_sequence)
            
            return ''.join(final_sequence)
    else:
        # No sequences within hydropathy tolerance, return original sequence
        return sequence