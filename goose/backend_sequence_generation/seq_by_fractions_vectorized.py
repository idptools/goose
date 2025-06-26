'''
Code to generate sequences by specified amino acid fractions using vectorized operations.
This allows for direct control of the amino acid composition in generated sequences.
'''

import numpy as np
from typing import Dict, List, Optional, Union
from dataclasses import dataclass, field
from goose.data import aa_list_probabilities as aa_probs

# Constants
AMINO_ACIDS = "ACDEFGHIKLMNPQRSTVWY"
AMINO_ACIDS_LIST = list(AMINO_ACIDS)
CHARGED_RESIDUES = "DEKR"
POSITIVE_CHARGED = "KR"
NEGATIVE_CHARGED = "DE"
DEFAULT_REMAINING_PROBABILITIES = aa_probs.HumanProbabilitiesByAA

@dataclass
class SequenceParametersByFractions:
    """
    Parameters for sequence generation based on amino acid fractions.
    
    This class handles the specification of amino acid fractions for sequence generation,
    with validation and normalization of the provided fractions.
    
    Attributes:
        length (int): Length of sequences to generate
        fractions (Dict[str, float]): Dictionary mapping amino acids to their fractions
        randomize_unspecified (bool): Whether to randomize unspecified fractions
        normalize (bool): Whether to normalize fractions to sum to 1.0
        strict (bool): Whether to strictly enforce specified fractions
        default_remaining_probabilities (Dict[str, float]): Probability distribution 
            used for unspecified amino acids when randomize_unspecified is False
    """
    length: int
    fractions: Dict[str, float] = field(default_factory=dict)
    randomize_unspecified: bool = False
    normalize: bool = True
    strict: bool = True  # Controls whether specified fractions must be preserved exactly
    default_remaining_probabilities: Dict[str, float] = field(default_factory=lambda: DEFAULT_REMAINING_PROBABILITIES)
    
    def __post_init__(self):
        """Validate and normalize the amino acid fractions."""
        # Make sure we're working with a copy
        self.fractions = self.fractions.copy() if self.fractions else {}
        
        # Ensure all keys are uppercase and valid
        valid_fractions = {}
        for aa, fraction in self.fractions.items():
            aa = aa.upper()
            if aa in AMINO_ACIDS and fraction >= 0:
                valid_fractions[aa] = fraction
        
        self.fractions = valid_fractions
        
        # Validate default_remaining_probabilities
        valid_default_probs = {}
        for aa, prob in self.default_remaining_probabilities.items():
            aa = aa.upper()
            if aa in AMINO_ACIDS and prob >= 0:
                valid_default_probs[aa] = prob
        self.default_remaining_probabilities = valid_default_probs
        
        # Calculate sum of specified fractions
        specified_sum = sum(self.fractions.values())
        
        # Handle empty fractions case - use uniform distribution
        if not self.fractions:
            self.fractions = aa_probs.EvenProbabilities
            return
        
        # If normalize is True, adjust fractions to sum to 1.0
        if self.normalize:
            # If strict mode is on, we never modify specified amino acids
            if self.strict:
                # If sum is too large, we can't enforce strictness
                if specified_sum > 1.0:
                    raise ValueError(f"Specified fractions sum to {specified_sum}, which exceeds 1.0. "
                                     f"Cannot maintain strict fractions. Either set strict=False or "
                                     f"reduce the sum of specified fractions.")
                
                # Find unspecified amino acids
                unspecified = [aa for aa in AMINO_ACIDS if aa not in self.fractions]
                
                # If there are unspecified amino acids, distribute the remaining probability
                if unspecified:
                    remaining = 1.0 - specified_sum
                    
                    # If remaining is effectively 0, avoid division problems
                    if abs(remaining) < 1e-10:
                        for aa in unspecified:
                            self.fractions[aa] = 0.0
                    elif self.randomize_unspecified:
                        # Generate random values for unspecified AAs
                        random_values = np.random.random(len(unspecified))
                        random_values = random_values / np.sum(random_values) * remaining
                        
                        for i, aa in enumerate(unspecified):
                            self.fractions[aa] = random_values[i]
                    else:
                        # Distribute using normalized default_remaining_probabilities
                        # Calculate sum of default_remaining_probabilities for unspecified amino acids
                        unspec_prob_sum = sum(self.default_remaining_probabilities.get(aa, 0) for aa in unspecified)
                        
                        # Distribute with normalization
                        if unspec_prob_sum > 0:
                            for aa in unspecified:
                                self.fractions[aa] = (self.default_remaining_probabilities.get(aa, 0) / unspec_prob_sum) * remaining
                        else:
                            # Fallback to even distribution if no probabilities available
                            per_aa = remaining / len(unspecified)
                            for aa in unspecified:
                                self.fractions[aa] = per_aa
                
                # If all amino acids are specified but sum < 1.0, normalize instead of raising error
                elif specified_sum < 0.999:  # Allow for small floating point errors
                    if self.normalize:
                        for aa in self.fractions:
                            self.fractions[aa] /= specified_sum
                    else:
                        raise ValueError(f"All amino acids are specified but their fractions sum to {specified_sum}, "
                                         f"which is less than 1.0. Please correct the fractions or set normalize=True.")
            else:
                # In non-strict mode, scale everything if sum > 1.0
                if specified_sum > 1.0:
                    for aa in self.fractions:
                        self.fractions[aa] /= specified_sum
                    return
                
                # Handle the case where sum < 1.0 similar to strict mode
                unspecified = [aa for aa in AMINO_ACIDS if aa not in self.fractions]
                if unspecified:
                    remaining = 1.0 - specified_sum
                    
                    if self.randomize_unspecified:
                        # Generate random values for unspecified AAs
                        random_values = np.random.random(len(unspecified))
                        random_values = random_values / np.sum(random_values) * remaining
                        
                        for i, aa in enumerate(unspecified):
                            self.fractions[aa] = random_values[i]
                    else:
                        # Distribute using normalized default_remaining_probabilities
                        # Calculate sum of default_remaining_probabilities for unspecified amino acids
                        unspec_prob_sum = sum(self.default_remaining_probabilities.get(aa, 0) for aa in unspecified)
                        
                        # Distribute with normalization
                        if unspec_prob_sum > 0:
                            for aa in unspecified:
                                self.fractions[aa] = (self.default_remaining_probabilities.get(aa, 0) / unspec_prob_sum) * remaining
                        else:
                            # Fallback to even distribution if no probabilities available
                            per_aa = remaining / len(unspecified)
                            for aa in unspecified:
                                self.fractions[aa] = per_aa


    
    def get_fraction_array(self):
        """Convert fractions dictionary to a numpy array in AMINO_ACIDS order."""
        return np.array([self.fractions.get(aa, 0.0) for aa in AMINO_ACIDS])
    
    def __str__(self):
        """String representation showing fractions and calculated properties."""
        fractions_str = ', '.join(f"{aa}:{self.fractions.get(aa, 0):.3f}" 
                                for aa in sorted(self.fractions.keys()))
        


class FractionBasedSequenceGenerator:
    """
    Class for generating protein sequences with specific amino acid fractions.
    """
    def __init__(self, 
                 length: Optional[int] = None,
                 fractions: Optional[Dict[str, float]] = None,
                 randomize_unspecified: bool = True,
                 normalize: bool = True,
                 strict: bool = True,
                 default_remaining_probabilities: Optional[Dict[str, float]] = None):
        """
        Initialize the sequence generator.
        
        Args:
            length (int, optional): Length of sequences to generate
            fractions (Dict[str, float], optional): Specified amino acid fractions
            randomize_unspecified (bool): Whether to randomize unspecified fractions
            normalize (bool): Whether to normalize fractions to sum to 1.0
            strict (bool): Whether to strictly enforce specified fractions
            default_remaining_probabilities (Dict[str, float], optional): Probability
                distribution used for unspecified amino acids
        """
        self.length = length
        self.fractions = fractions if fractions else {}
        self.randomize_unspecified = randomize_unspecified
        self.normalize = normalize
        self.strict = strict
        if default_remaining_probabilities is None:
            self.default_remaining_probabilities = DEFAULT_REMAINING_PROBABILITIES
        else:
            self.default_remaining_probabilities = default_remaining_probabilities
        
        # Create a mapping from amino acid to index
        self.aa_to_idx = {aa: i for i, aa in enumerate(AMINO_ACIDS)}
        self.idx_to_aa = {i: aa for i, aa in enumerate(AMINO_ACIDS)}
    
    def generate_sequences(self, 
                          num_sequences: int = 1, 
                          length: Optional[int] = None,
                          fractions: Optional[Dict[str, float]] = None,
                          strict: Optional[bool] = None,
                          default_remaining_probabilities: Optional[Dict[str, float]] = None,
                          convert_to_amino_acids: bool = True) -> List[str]:
        """
        Generate protein sequences using the specified amino acid fractions.
        
        Args:
            num_sequences (int): Number of sequences to generate
            length (int, optional): Length of sequences to generate (overrides init value)
            fractions (Dict[str, float], optional): Amino acid fractions (overrides init value)
            strict (bool, optional): Whether to strictly enforce specified fractions
            default_remaining_probabilities (Dict[str, float], optional): Probability
                distribution used for unspecified amino acids
            convert_to_amino_acids (bool): Whether to convert indices back to amino acid strings
            
        Returns:
            List[str]: List of generated protein sequences
        """
        # Validate input parameters
        if num_sequences <= 0:
            raise ValueError(f"num_sequences must be positive, got {num_sequences}")
            
        # Use provided values or fall back to initialized values
        use_length = length if length is not None else self.length
        if use_length is None or use_length <= 0:
            raise ValueError(f"Sequence length must be positive, got {use_length}")
            
        use_fractions = fractions if fractions is not None else self.fractions
        use_strict = strict if strict is not None else self.strict
        use_default_probs = default_remaining_probabilities if default_remaining_probabilities is not None else self.default_remaining_probabilities
        
        # Create parameters object
        params = SequenceParametersByFractions(
            length=use_length,
            fractions=use_fractions,
            randomize_unspecified=self.randomize_unspecified,
            normalize=self.normalize,
            strict=use_strict,
            default_remaining_probabilities=use_default_probs
        )
        
        # Add warning for very short sequences when many amino acids are specified
        if use_strict and use_length < 20:
            specified_aa_count = sum(1 for aa in AMINO_ACIDS if aa in use_fractions and use_fractions[aa] > 0)
            if specified_aa_count > use_length / 2:
                import warnings
                warnings.warn(
                    f"Sequence length ({use_length}) may be too short to accurately represent "
                    f"all {specified_aa_count} specified amino acid fractions. "
                    f"Some amino acids may not appear in the generated sequences."
                )
        
        sequences = []
        
        if use_strict:
            # Deterministic approach for strict mode - exact counts
            for _ in range(num_sequences):
                # More precise algorithm for calculating exact integer counts
                aa_counts = {}
                fractional_counts = {}
                remaining_positions = use_length
                
                # First pass: calculate integer counts for each amino acid
                for aa in AMINO_ACIDS:
                    if aa in params.fractions and params.fractions[aa] > 0:
                        # Calculate exact fractional count
                        frac_count = params.fractions[aa] * use_length
                        # Store integer part
                        aa_counts[aa] = int(frac_count)  # floor, not round
                        # Store fractional part for later allocation
                        fractional_counts[aa] = frac_count - int(frac_count)
                        remaining_positions -= aa_counts[aa]
                
                # Distribute remaining positions based on fractional parts
                if remaining_positions > 0:
                    # Sort amino acids by descending fractional part
                    sorted_aa = sorted(fractional_counts.keys(), 
                                     key=lambda aa: fractional_counts[aa],
                                     reverse=True)
                    
                    # Allocate remaining positions to amino acids with highest fractional parts
                    for i in range(remaining_positions):
                        if i < len(sorted_aa):
                            aa_counts[sorted_aa[i]] += 1
                        else:
                            # If we have more positions than amino acids with fractions,
                            # distribute to unspecified amino acids
                            unspecified = [aa for aa in AMINO_ACIDS if aa not in params.fractions]
                            if unspecified:
                                idx = i % len(unspecified)
                                aa = unspecified[idx]
                                aa_counts[aa] = aa_counts.get(aa, 0) + 1
                            else:
                                # If all amino acids are specified, distribute evenly
                                idx = i % len(AMINO_ACIDS)
                                aa = AMINO_ACIDS[idx]
                                aa_counts[aa] = aa_counts.get(aa, 0) + 1
                
                # Verify we have the right total
                total_count = sum(aa_counts.values())
                if total_count != use_length:
                    # This should not happen, but just in case
                    diff = use_length - total_count
                    if diff > 0:
                        # Add missing positions
                        for i in range(diff):
                            aa = AMINO_ACIDS[i % len(AMINO_ACIDS)]
                            aa_counts[aa] = aa_counts.get(aa, 0) + 1
                    else:
                        # Remove extra positions (least likely to be needed)
                        sorted_counts = sorted(aa_counts.items(), key=lambda x: (x[1], x[0]), reverse=True)
                        for i in range(-diff):
                            aa = sorted_counts[i % len(sorted_counts)][0]
                            if aa_counts[aa] > 0:
                                aa_counts[aa] -= 1
                
                # Create a pool of amino acids based on counts
                aa_pool = []
                for aa, count in aa_counts.items():
                    aa_pool.extend([aa] * count)
                
                # Verify we have the right number of amino acids
                assert len(aa_pool) == use_length, f"Expected {use_length} amino acids, got {len(aa_pool)}"
                
                # Shuffle the pool to randomize order
                np.random.shuffle(aa_pool)
                
                # Convert to string
                if convert_to_amino_acids:
                    sequences.append(''.join(aa_pool))
                else:
                    # Convert to indices if not converting to amino acids
                    sequence_indices = [self.aa_to_idx[aa] for aa in aa_pool]
                    sequences.append(sequence_indices)
        else:
            # Probabilistic approach for non-strict mode
            # Get normalized fractions as array
            fractions_array = params.get_fraction_array()
            
            # Generate sequences using vectorized operations
            for _ in range(num_sequences):
                # Generate sequence indices based on probabilities
                indices = np.random.choice(
                    len(AMINO_ACIDS),
                    size=use_length,
                    p=fractions_array
                )
                
                # Convert indices to amino acids
                if convert_to_amino_acids:
                    sequence = ''.join(AMINO_ACIDS[idx] for idx in indices)
                else:
                    sequence = indices

                sequences.append(sequence)
            
        return sequences

