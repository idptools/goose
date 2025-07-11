'''
Code to generate sequences by specified amino acid fractions using vectorized operations.
This allows for direct control of the amino acid composition in generated sequences.
'''

import numpy as np
from typing import Dict, List, Optional
from dataclasses import dataclass, field
from goose.data import aa_list_probabilities as aa_probs

# Constants
AMINO_ACIDS = "ACDEFGHIKLMNPQRSTVWY"
AMINO_ACIDS_LIST = list(AMINO_ACIDS)
CHARGED_RESIDUES = "DEKR"
POSITIVE_CHARGED = "KR"
NEGATIVE_CHARGED = "DE"
DEFAULT_REMAINING_PROBABILITIES = aa_probs.IDRProbs

@dataclass
class SequenceParametersByFractions:
    """
    Parameters for sequence generation based on amino acid fractions.
    
    This class handles the specification of amino acid fractions for sequence generation,
    with validation and normalization of the provided fractions.
    
    Behavior:
    - Specified fractions are ALWAYS enforced as exact counts (deterministic allocation).
    - Unspecified amino acids are ALWAYS allocated probabilistically based on 
      default_remaining_probabilities.
    - The 'strict' parameter controls other validation behavior but does not affect
      the deterministic nature of specified fractions.
    
    Attributes:
        length (int): Length of sequences to generate
        fractions (Dict[str, float]): Dictionary mapping amino acids to their fractions
        randomize_unspecified (bool): Whether to randomize unspecified fractions using
            uniform random distribution (only applies when strict=False)
        normalize (bool): Whether to normalize fractions to sum to 1.0
        strict (bool): Whether to use strict validation (affects error handling for
            invalid fraction sums, but specified fractions are always deterministic)
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
                
                # In strict mode, we don't pre-assign fractions to unspecified amino acids
                # They will be handled probabilistically during sequence generation
                
                # If all amino acids are specified but sum < 1.0, normalize instead of raising error
                unspecified = [aa for aa in AMINO_ACIDS if aa not in self.fractions]
                if not unspecified and specified_sum < 0.999:  # Allow for small floating point errors
                    if self.normalize:
                        for aa in self.fractions:
                            self.fractions[aa] /= specified_sum
                    else:
                        raise ValueError(f"All amino acids are specified but their fractions sum to {specified_sum}, "
                                         f"which is less than 1.0. Please correct the fractions or set normalize=True.")
            else:
                # In non-strict mode, we still use the same deterministic approach
                # but allow normalization of specified fractions if sum > 1.0
                if specified_sum > 1.0:
                    for aa in self.fractions:
                        self.fractions[aa] /= specified_sum


    
    def get_fraction_array(self):
        """Convert fractions dictionary to a numpy array in AMINO_ACIDS order."""
        return np.array([self.fractions.get(aa, 0.0) for aa in AMINO_ACIDS])
    
    def get_unspecified_probabilities(self):
        """
        Get normalized probabilities for unspecified amino acids.
        
        Returns:
            tuple: (unspecified_aas, normalized_probs) where unspecified_aas is a list
                   of amino acids not specified in fractions, and normalized_probs is
                   the corresponding normalized probability distribution.
        """
        unspecified = [aa for aa in AMINO_ACIDS if aa not in self.fractions or self.fractions[aa] == 0]
        
        if not unspecified:
            return [], []
        
        # Get probabilities for unspecified amino acids
        probs = [self.default_remaining_probabilities.get(aa, 0) for aa in unspecified]
        
        # Normalize probabilities
        prob_sum = sum(probs)
        if prob_sum > 0:
            probs = [p / prob_sum for p in probs]
        else:
            # Fallback to uniform distribution
            probs = [1.0 / len(unspecified)] * len(unspecified)
        
        return unspecified, probs
    
    def __str__(self):
        """String representation showing fractions and calculated properties."""
        fractions_str = ', '.join(f"{aa}:{self.fractions.get(aa, 0):.3f}" 
                                for aa in sorted(self.fractions.keys()))
        


class FractionBasedSequenceGenerator:
    """
    Class for generating protein sequences with specific amino acid fractions.
    
    This generator enforces specified amino acid fractions as exact counts while
    allocating unspecified amino acids probabilistically:
    
    - Specified fractions: Always calculated as exact integer counts (deterministic)
    - Unspecified amino acids: Always allocated probabilistically based on 
      default_remaining_probabilities
    
    For example, with length=100 and fractions={'A':0.1}, you will always get 
    exactly 10 A's in every sequence, while the remaining 90 positions will be
    filled probabilistically based on default_remaining_probabilities.
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
                (these will always be enforced as exact counts)
            randomize_unspecified (bool): Whether to randomize unspecified fractions
                (only affects initialization, not sequence generation behavior)
            normalize (bool): Whether to normalize fractions to sum to 1.0
            strict (bool): Whether to use strict validation for fraction sums
            default_remaining_probabilities (Dict[str, float], optional): Probability
                distribution used for unspecified amino acids (always probabilistic)
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
        
        Specified fractions are always enforced as exact counts, while unspecified
        amino acids are allocated probabilistically based on default_remaining_probabilities.
        
        Args:
            num_sequences (int): Number of sequences to generate
            length (int, optional): Length of sequences to generate (overrides init value)
            fractions (Dict[str, float], optional): Amino acid fractions (overrides init value)
                These will always be enforced as exact counts.
            strict (bool, optional): Whether to use strict validation (does not affect
                the deterministic nature of specified fractions)
            default_remaining_probabilities (Dict[str, float], optional): Probability
                distribution used for unspecified amino acids (always probabilistic)
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
        if use_length < 20:
            specified_aa_count = sum(1 for aa in AMINO_ACIDS if aa in use_fractions and use_fractions[aa] > 0)
            if specified_aa_count > use_length / 2:
                import warnings
                warnings.warn(
                    f"Sequence length ({use_length}) may be too short to accurately represent "
                    f"all {specified_aa_count} specified amino acid fractions. "
                    f"Some amino acids may not appear in the generated sequences."
                )
        
        sequences = []
        
        # Always use deterministic allocation for specified fractions
        # and probabilistic allocation for unspecified amino acids
        for _ in range(num_sequences):
            # Calculate exact counts for specified amino acids
            specified_aa_counts = {}
            fractional_counts = {}
            total_specified_positions = 0
            
            # First pass: calculate integer counts for specified amino acids only
            for aa in AMINO_ACIDS:
                if aa in params.fractions and params.fractions[aa] > 0:
                    # Calculate exact fractional count
                    frac_count = params.fractions[aa] * use_length
                    # Store integer part
                    specified_aa_counts[aa] = int(frac_count)  # floor, not round
                    # Store fractional part for later allocation (only if > 0)
                    fractional_part = frac_count - int(frac_count)
                    if fractional_part > 1e-10:  # Only store if fractional part is significant
                        fractional_counts[aa] = fractional_part
                    total_specified_positions += specified_aa_counts[aa]
            
            # Distribute remaining fractional positions among specified amino acids
            remaining_fractional_positions = use_length - total_specified_positions
            
            # Handle fractional parts for specified amino acids
            if remaining_fractional_positions > 0 and fractional_counts:
                # Sort amino acids by descending fractional part
                sorted_aa = sorted(fractional_counts.keys(), 
                                 key=lambda aa: fractional_counts[aa],
                                 reverse=True)
                
                # Allocate remaining positions to amino acids with highest fractional parts
                positions_to_allocate = min(remaining_fractional_positions, len(sorted_aa))
                for i in range(positions_to_allocate):
                    specified_aa_counts[sorted_aa[i]] += 1
                    total_specified_positions += 1
                
                # Update remaining positions
                remaining_fractional_positions -= positions_to_allocate
            
            # Now handle unspecified amino acids probabilistically
            remaining_positions = use_length - total_specified_positions
            unspecified_aa_counts = {}
            
            if remaining_positions > 0:
                # Get unspecified amino acids and their probabilities
                unspecified, unspec_probs = params.get_unspecified_probabilities()
                
                if unspecified:
                    # Probabilistically allocate remaining positions
                    for _ in range(remaining_positions):
                        # Choose amino acid based on probabilities
                        chosen_idx = np.random.choice(len(unspecified), p=unspec_probs)
                        chosen_aa = unspecified[chosen_idx]
                        unspecified_aa_counts[chosen_aa] = unspecified_aa_counts.get(chosen_aa, 0) + 1
                else:
                    # If no unspecified amino acids, distribute remaining positions
                    # among specified amino acids (this is rare but possible)
                    specified_aas = list(specified_aa_counts.keys())
                    if specified_aas:
                        for i in range(remaining_positions):
                            aa = specified_aas[i % len(specified_aas)]
                            specified_aa_counts[aa] += 1
            
            # Combine all counts
            all_aa_counts = {**specified_aa_counts, **unspecified_aa_counts}
            
            # Create a pool of amino acids based on counts
            aa_pool = []
            for aa, count in all_aa_counts.items():
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
            
        return sequences

