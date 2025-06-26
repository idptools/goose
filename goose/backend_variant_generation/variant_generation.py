'''
Updated code for variant generation. Going to rewrite a lot 
of the functionality to use numpy vectorized operations for efficiency.
'''
import numpy as np

from sparrow.protein import Protein

from goose.backend_sequence_generation.sequence_generation_vectorized import generate_seq_by_props
from goose.backend_property_calculation.calculate_kappa import kappa
from goose.backend_property_optimization.optimize_kappa import optimize_kappa_vectorized
from goose.backend_sequence_generation.create_sequence_by_class import create_sequence_by_class
from goose.backend_sequence_generation.seq_by_probability_vectorized import SequenceGenerator
from goose import parameters
from goose.backend_property_optimization.optimize_hydropathy_minimal_changes import optimize_hydropathy_minimal_changes
from goose.backend_property_optimization.within_class_hydropathy_optimization import optimize_hydropathy_within_class_vectorized
from goose.backend_property_optimization.optimize_hydropathy import optimize_hydropathy_vectorized 
from goose.backend_variant_generation.helper_functions import needed_charged_residues, change_all_residues_within_class, optimize_hydropathy_avoid_original_residues, calculate_amino_acid_class_fractions




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
        tolerance= parameters.HYDRO_ERROR
    )

    # see if within tolerance, if not use the within class optimization to get it closer.
    final_hydropathy = Protein(variant_sequence).hydrophobicity
    if abs(final_hydropathy - original_hydropathy) > parameters.HYDRO_ERROR:
        variant_sequence = optimize_hydropathy_within_class_vectorized(
            [variant_sequence],
            target_hydropathy=original_hydropathy,
            max_iterations=1000,
            tolerance=parameters.HYDRO_ERROR,
            only_return_within_tolerance=True
        )[0]
    
    return variant_sequence
    


def generate_minimal_variant(input_sequence: str,
                             target_hydropathy: float = None,
                             target_kappa: float = None,
                             target_FCR: float = None,
                             target_NCPR: float = None,
                             max_iterations: int = 1000,
                             hydropathy_tolerance: float = parameters.HYDRO_ERROR,
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
        set by parameters.HYDRO_ERROR
    kappa_tolerance : float, default=0.03
        set by parameters.MAXIMUM_KAPPA_ERROR
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
            num_iterations=max_iterations,
            tolerance=kappa_tolerance,
            inputting_matrix=False,
            convert_input_seq_to_matrix=True,
            return_when_num_hit=1,
            avoid_shuffle=True,
            num_copies=10
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
                                        kappa_tolerance: float = parameters.MAXIMUM_KAPPA_ERROR,
                                        hydropathy_tolerance: float = parameters.HYDRO_ERROR) -> str:
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
        num_sequences=1,
        convert_to_amino_acids=True
    )

    print(new_sequence)

    # now use within class hydropathy optimization to get the hydropathy correct.
    new_sequence = optimize_hydropathy_within_class_vectorized(
        [new_sequence],
        target_hydropathy=original_hydropathy,
        max_iterations=1000,
        tolerance=hydropathy_tolerance,
        only_return_within_tolerance=False
    )[0]

    # finally, optimize kappa to get it to the original kappa
    new_sequence = optimize_kappa_vectorized(
        [new_sequence],
        target_kappa=original_kappa,
        num_iterations=1000,
        convert_input_seq_to_matrix=True,
        inputting_matrix=False,
        tolerance=kappa_tolerance
    )[0]

    return new_sequence

def generate_constant_properties_var(sequence: str,
                                     hydropathy_tolerance=parameters.HYDRO_ERROR,
                                     kappa_tolerance=parameters.MAXIMUM_KAPPA_ERROR) -> str:
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
    return generate_seq_by_props(
        length=length,
        fcr=original_FCR,
        ncpr=original_NCPR,
        hydropathy=original_hydropathy,
        kappa=original_kappa
        )

def generate_constant_residue_var(sequence: str,
                                  constant_residues: list,
                                  hydropathy_tolerance=parameters.HYDRO_ERROR,
                                  kappa_tolerance=parameters.MAXIMUM_KAPPA_ERROR) -> str:
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
        set by parameters.HYDRO_ERROR
        Acceptable difference for hydropathy optimization.
    kappa_tolerance : float, default=0.03
        set by parameters.MAXIMUM_KAPPA_ERROR
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

   # Generate a new sequence with the same properties
    new_sequence = generate_seq_by_props(
        length=length,
        fcr=original_FCR,
        ncpr=original_NCPR,
        hydropathy=original_hydropathy,
        kappa=original_kappa
    )

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
        # If no constant charged residues, we can optimize kappa
        final_sequence = optimize_kappa_vectorized(
            [''.join(final_sequence)],
            target_kappa=original_kappa,
            num_iterations=1000,
            convert_input_seq_to_matrix=True,
            inputting_matrix=False,
            tolerance=kappa_tolerance
        )[0]
    
    return final_sequence
    