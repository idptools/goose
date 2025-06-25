'''
functionality to mutate a sequence that is a numpy array of integers
'''
import numpy as np
from typing import List, Union

AA_TO_INT = {
    'A': 0, 'C': 1, 'D': 2, 'E': 3, 'F': 4,
    'G': 5, 'H': 6, 'I': 7, 'K': 8, 'L': 9,
    'M': 10, 'N': 11, 'P': 12, 'Q': 13, 'R': 14,
    'S': 15, 'T': 16, 'V': 17, 'W': 18, 'Y': 19
}

def generate_mutant_sequences(sequence_array: np.ndarray,
                                             num_mutants: int,
                                             mutations_per_sequence: int = 1,
                                             exclude_residues: List[str] = None,
                                             exclude_positions: List[int] = None,
                                             forbidden_targets: List[str] = None,
                                             add_excluded_residues_to_forbidden_targets=True) -> np.ndarray:
    """
    Generate mutant sequences from a given sequence matrix.
    
    Parameters:
    -----------
    sequence_array : np.ndarray
        Input sequence matrix as a numpy array of integers.
        Should be in shape (seq_length,).
    num_mutants : int
        Number of mutant sequences to generate.
    mutations_per_sequence : int, optional
        Number of mutations to apply to each sequence (default is 1).
    exclude_residues : List[str], optional
        List of residues to exclude from mutation by identity (default is None).
    exclude_positions : List[int], optional
        List of positions (0-indexed) to exclude from mutation (default is None).
    forbidden_targets : List[str], optional
        List of residues that mutations cannot change to (default is None).
    add_excluded_residues_to_forbidden_targets : bool, optional
        If True, adds excluded residues to the forbidden targets (default is True).
        
    Returns:
    --------
    np.ndarray
        Array of mutated sequence matrices with shape (num_mutants, seq_length).
    """
    if exclude_residues is None:
        exclude_residues = []
    if exclude_positions is None:
        exclude_positions = []
    if forbidden_targets is None:
        forbidden_targets = []
    
    if add_excluded_residues_to_forbidden_targets:
        # Add excluded residues to forbidden targets if specified
        forbidden_targets = list(set(forbidden_targets + exclude_residues))

    exclude_indices = np.array([AA_TO_INT[res] for res in exclude_residues if res in AA_TO_INT])
    forbidden_target_indices = np.array([AA_TO_INT[res] for res in forbidden_targets if res in AA_TO_INT])
    seq_length = len(sequence_array)
    
    # Create a mask for positions that can be mutated
    mutation_mask = np.ones(seq_length, dtype=bool)
    
    # Exclude by residue identity
    if len(exclude_indices) > 0:
        residue_mask = np.isin(sequence_array, exclude_indices)
        mutation_mask = mutation_mask & ~residue_mask

    # Exclude by position
    if len(exclude_positions) > 0:
        valid_positions = [pos for pos in exclude_positions if 0 <= pos < seq_length]
        if valid_positions:
            mutation_mask[valid_positions] = False

    # Get mutable positions
    mutable_positions = np.where(mutation_mask)[0]
    
    # If no mutable positions are available after exclusions, handle the case
    if len(mutable_positions) == 0:
        raise ValueError("No mutable positions available after exclusions.")
    # Ensure mutations_per_sequence is less than or equal to the number of mutable positions
    if len(mutable_positions) < mutations_per_sequence:
        # make mutations_per_sequence less than or equal to the number of mutable positions
        mutations_per_sequence = len(mutable_positions)
    
    # Create allowed amino acid indices (excluding forbidden targets)
    all_aa_indices = np.arange(len(AA_TO_INT))
    if len(forbidden_target_indices) > 0:
        allowed_target_indices = all_aa_indices[~np.isin(all_aa_indices, forbidden_target_indices)]
    else:
        allowed_target_indices = all_aa_indices
    
    # Generate all mutant sequences simultaneously
    mutant_sequences = np.tile(sequence_array, (num_mutants, 1))
    
    # Ultra-vectorized approach: generate all mutations at once
    total_mutations = num_mutants * mutations_per_sequence

    # Generate all position selections at once
    all_positions = np.random.choice(
        mutable_positions, 
        size=total_mutations, 
        replace=True  # Allow replacement since we're across different sequences
    ).reshape(num_mutants, mutations_per_sequence)
    
    # Create indices for advanced indexing
    mutant_indices = np.repeat(np.arange(num_mutants), mutations_per_sequence)
    position_indices = all_positions.flatten()

    # Get current amino acids at mutation positions
    current_aas = mutant_sequences[mutant_indices, position_indices]
    
    # Generate new amino acids (excluding current ones and forbidden targets)
    new_aas = np.array([
        np.random.choice(allowed_target_indices[
            (allowed_target_indices != current_aa)
        ])
        for current_aa in current_aas
    ])
    
    # Apply all mutations at once using advanced indexing
    mutant_sequences[mutant_indices, position_indices] = new_aas
    
    return mutant_sequences

