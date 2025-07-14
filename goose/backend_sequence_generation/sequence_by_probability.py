'''
Functionality to use vectorized numpy operations to generate amino acid sequences
using a dictionary of probabilities for each amino acid as a dictioanry where the keys
are the amiono acids and the values are the probabilities for each amino acid as the input
as well as the desired number of sequences and the sequence length. 
'''

import numpy as np
from goose.data import aa_list_probabilities as aa_probs

from typing import List, Dict, Optional
from dataclasses import dataclass

from goose.backend.parameters import calculate_max_charge

# Constants
AMINO_ACIDS = "ACDEFGHIKLMNPQRSTVWY"
HYDROPATHY_SCALE = np.array([6.3, 7.0, 1.0, 1.0, 7.3, 4.1, 1.3, 9.0, 
                            0.6, 8.3, 6.4, 1.0, 2.9, 1.0, 0.0, 3.7, 3.8, 8.7, 3.6, 3.2])
CHARGE_SCALE = np.array([0, 0, -1, -1, 0, 0, 0, 0, 1, 0, 0, 
                        0, 0, 0, 1, 0, 0, 0, 0, 0])
AMINO_ACIDS_ARRAY = np.array(['A', 'C', 'D', 'E', 'F', 'G', 'H', 'I', 'K', 'L', 'M', 'N', 'P', 'Q', 'R', 'S', 'T', 'V', 'W', 'Y'])
DEFAULT_PROBABILITIES = aa_probs.IDRProbs

@dataclass
class SequenceParameters:
    """Parameters for sequence generation with smart defaults."""
    length: int
    target_hydropathy: Optional[float] = None
    target_fcr: Optional[float] = None
    target_ncpr: Optional[float] = None

    def __post_init__(self):
        """Initialize parameter ranges and set defaults."""
        # Define parameter ranges. These are restricted to be chosen
        # to maximize probabilitiy of disordered sequences. 
        self.hydropathy_range = (1, 4.5)
        self.fcr_range = (0, 1)
        self.ncpr_range = (-1, 1)

    @staticmethod
    def _calculate_max_fcr(target_hydropathy: float) -> float:
        """
        Calculate the maximum fraction of charged residues based on target hydropathy.
            This was empirically determined...
        """
        return calculate_max_charge(target_hydropathy)


    @staticmethod
    def _calculate_charged_residues(target_fcr, target_ncpr, objective_length):
        if target_ncpr > target_fcr:
            raise Exception('cannot have objective NCPR greater than objective FCR.')

        positive_charged_residues = round((objective_length * (target_fcr + target_ncpr)) / 2)
        negative_charged_residues = round((objective_length * (target_fcr - target_ncpr)) / 2)

        return positive_charged_residues, negative_charged_residues


    @staticmethod
    def _get_necessary_hydropathy_value_with_charge(sequence_length: float,
                                                    target_hydropathy: float,
                                                    target_fcr: float,
                                                    target_ncpr: float = None) -> float:
        """
        Calculate the hydropathy value needed for the rest of the sequence
        to get the target hydropathy value with the specified charge.
        """
        # calculate the hydropathy value needed for the rest of the sequence
        # to get the target hydropathy value with the specified charge.
        if target_ncpr is None:
            num_charged= round(target_fcr*sequence_length)
            per_charged_residue = 0.65
        else:
            num_positive, num_negative = SequenceParameters._calculate_charged_residues(target_fcr, target_ncpr, sequence_length)
            num_charged = num_positive + num_negative
            if num_charged==0:
                per_charged_residue=0
            else:
                per_charged_residue = ((num_positive*0.3)+(num_negative*1))/(num_positive+num_negative)
        # get charged hydropathy
        charged_hydropathy = per_charged_residue*num_charged
        # get total hydropathy
        total_hydropathy = target_hydropathy*sequence_length
        # get remaining hydropathy
        remaining_hydropathy = total_hydropathy - charged_hydropathy
        if num_charged==sequence_length:
            return 1.0
        objective_hydropathy = remaining_hydropathy/(sequence_length-num_charged)
        return objective_hydropathy
    

class SequenceGenerator:
    """Class for generating protein sequences with specific properties."""
    def __init__(self, 
                 length=None,
                 fcr=None,
                 ncpr=None,
                 hydropathy=None,
                 num_sequences: int = 1,
                 use_weighted_probabilities: bool = False,
                 chosen_probabilities: Dict[str, float] = None,
                 exclude_residues: List[str] = None):
        """Initialize the sequence generator with specified parameters."""
        self.length = length
        self.fcr = fcr
        self.ncpr = ncpr
        self.hydropathy = hydropathy
        self.num_sequences = num_sequences
        self.use_weighted_probabilities = use_weighted_probabilities
        self.chosen_probabilities = chosen_probabilities
        self.exclude_residues = exclude_residues if exclude_residues is not None else []
        # If no probabilities are provided, use the default
        if chosen_probabilities is None:
            self.chosen_probabilities = DEFAULT_PROBABILITIES

        # Initialize the amino acid to integer and integer to amino acid mappings
        self.aa_to_int={'A': 0, 'C': 1, 'D': 2, 'E': 3, 'F': 4, 'G': 5, 'H': 6, 'I': 7, 'K': 8, 'L': 9, 'M': 10, 'N': 11, 'P': 12, 'Q': 13, 'R': 14, 'S': 15, 'T': 16, 'V': 17, 'W': 18, 'Y': 19}
        self.int_to_aa={0: 'A', 1: 'C', 2: 'D', 3: 'E', 4: 'F', 5: 'G', 6: 'H', 7: 'I', 8: 'K', 9: 'L', 10: 'M', 11: 'N', 12: 'P', 13: 'Q', 14: 'R', 15: 'S', 16: 'T', 17: 'V', 18: 'W', 19: 'Y'}
        self.charged_residues_indices = np.where(CHARGE_SCALE != 0)[0]
        self.uncharged_residues_indices = np.where(CHARGE_SCALE == 0)[0]

        # check of the keys in chosen_probabilities are strings, if they are we need to convert them. 
        if isinstance(list(self.chosen_probabilities.keys())[0], str):
            # Convert string keys to integer indices
            self.chosen_probabilities = {self.aa_to_int[aa]: prob for aa, prob in self.chosen_probabilities.items()}
        
        # Pre-compute commonly used values
        self.positive_indices = np.where(CHARGE_SCALE > 0)[0]
        self.negative_indices = np.where(CHARGE_SCALE < 0)[0]
        self.hydropathy_lookup = HYDROPATHY_SCALE[self.uncharged_residues_indices]
        
        # Filter out excluded residues from charge indices
        if self.exclude_residues:
            exclude_indices = [self.aa_to_int[aa] for aa in self.exclude_residues]
            self.positive_indices = np.array([idx for idx in self.positive_indices if idx not in exclude_indices])
            self.negative_indices = np.array([idx for idx in self.negative_indices if idx not in exclude_indices])
        
    def _generate_parameter_batch(self, 
                                  length=None,
                                  fcr=None,
                                  ncpr=None,
                                  hydropathy=None,
                                  num_sequences=None,
                                  exclude_residues=None):
        """Generate a batch of parameters, randomizing unspecified values."""
        # Use specified values if provided
        if length is None:
            length = self.length
        if fcr is None:
            fcr = self.fcr
        if ncpr is None:
            ncpr = self.ncpr
        if hydropathy is None:
            hydropathy = self.hydropathy
        if num_sequences is None:
            num_sequences = self.num_sequences
        if exclude_residues is None:
            exclude_residues = self.exclude_residues
        
        # get excluded residues as numbers. 
        if exclude_residues is not None:
            exclude_residues = [self.aa_to_int[aa] for aa in exclude_residues]
        else:    
            exclude_residues = []
        
        # Initialize lists for parameters
        num_negative = []
        num_positive = []
        sequence_lengths = []  # renamed from length to avoid confusion
        probabilities_list = []  # renamed from probabilities to avoid confusion

        # see if we need to reset values. This lets us get different values for unspecified values
        # every iteration. 
        if hydropathy==None:
            reset_hydro=True
        else:
            reset_hydro=False
        if fcr==None:
            reset_fcr=True
        else:
            reset_fcr=False
        if ncpr==None:
            reset_ncpr=True
        else:
            reset_ncpr=False

        # iterate over the number of sequences and generate the parameters for each sequence
        for _ in range(num_sequences):
            if reset_hydro==True:
                hydropathy=None
            if reset_fcr==True:
                fcr=None
            if reset_ncpr==True:
                ncpr=None
            

            # If hydropathy is specified, calculate other parameters based on it
            if hydropathy is not None:
                # get max FCR, thereotical and empirical
                max_fcr_theoretical = SequenceParameters._calculate_max_fcr(
                        length=length,
                        hydropathy=hydropathy,
                        net_charge=ncpr,
                        empirical=False
                    )
                max_fcr_empirical = SequenceParameters._calculate_max_fcr(
                        length=length,
                        hydropathy=hydropathy,
                        net_charge=ncpr,
                        empirical=True
                    )
                
                # Validate provided FCR against hydropathy constraints
                if fcr is not None and fcr > max_fcr_theoretical:
                    raise ValueError(f"FCR {fcr} exceeds maximum allowed FCR {max_fcr_theoretical:.3f} for hydropathy {hydropathy}")

                # have max FCR already calculated. Just need to get min FCR and NCPR
                if fcr is None:
                    # set max_fcr to be the empirically possible max FCR if user does't specify it. 
                    # this increases our likelihood of generating disordered sequences
                    max_fcr = max_fcr_empirical
                    if ncpr is not None:
                        min_fcr = abs(ncpr)
                        # Check if constraints are satisfiable
                        if min_fcr > max_fcr:
                            # if min_fcr is greater than emperically determined max_fcr, try to use the theoretical max FCR
                            if min_fcr > max_fcr_theoretical:
                                # if min_fcr greater than max_fcr_theoretical, raise an error
                                raise ValueError(f"Cannot satisfy constraints: |NCPR| {abs(ncpr)} > max FCR {max_fcr_theoretical:.3f} for hydropathy {hydropathy}")
                            else:
                                max_fcr = max_fcr_theoretical
                    else:
                        min_fcr = 0.0
                    fcr = np.random.uniform(low=min_fcr, high=max_fcr)
                
                # Validate and set NCPR
                if ncpr is None:
                    ncpr = np.random.uniform(-fcr, fcr)
                else:
                    # Validate that NCPR is within valid range for the FCR
                    if abs(ncpr) > fcr:
                        raise ValueError(f"NCPR {ncpr} is outside valid range [-{fcr}, {fcr}] for FCR {fcr}")
                
                
            # Otherwise, choose parameters based on charge constraints
            else:
                if fcr is None:
                    if ncpr is not None:
                        min_fcr = abs(ncpr) 
                        fcr = np.random.uniform(low=min_fcr, high=1.0)
                    else:
                        min_fcr = 0.0
                        fcr = np.random.beta(1, 2)
                    
                if ncpr is None:
                    ncpr = np.random.uniform(-fcr, fcr)
                else:
                    # Validate that NCPR is within valid range for the FCR
                    if abs(ncpr) > fcr:
                        raise ValueError(f"NCPR {ncpr} is outside valid range [-{fcr}, {fcr}] for FCR {fcr}")
                
                # Calculate hydropathy range based on FCR
                max_hydropathy = (1-fcr) * 9.0
                min_hydropathy = (1-fcr) * 1.0
                hydropathy = np.random.uniform(low=min_hydropathy, high=max_hydropathy)
            
            # Get probabilities for the sequence
            if self.use_weighted_probabilities:
                probs = self._choose_probabilities(
                    target_hydropathy=hydropathy,
                    target_fcr=fcr,
                    target_ncpr=ncpr,
                    length=length
                )
                probs_list = [(probs[aa] if aa not in exclude_residues else 0.0) for aa in range(0,20)]
            else:
                probs_list = [(self.chosen_probabilities[aa] if aa not in exclude_residues else 0.0) for aa in range(0,20)]
            
            # Get charged residues counts
            pos, neg = SequenceParameters._calculate_charged_residues(fcr, ncpr, length)
            
            # Append parameters
            sequence_lengths.append(length)
            num_positive.append(pos)
            num_negative.append(neg)
            probabilities_list.append(probs_list)
       
        # Convert to numpy arrays
        return (np.array(num_positive), 
                np.array(num_negative), 
                np.array(sequence_lengths), 
                np.array(probabilities_list))
    
    def _populate_charged_residues(pos_counts, neg_counts, seq_lengths):
        """
        Randomly populate positions of charged residues in protein sequences using numpy vectorized operations.
        
        Parameters:
        -----------
        pos_counts : array-like
            Array containing the count of positively charged residues for each sequence
        neg_counts : array-like
            Array containing the count of negatively charged residues for each sequence
        seq_lengths : array-like
            Array containing the length of each protein sequence
        
        Returns:
        --------
        list of numpy.ndarray
            List where each element is an array representing a protein sequence with:
            - 1 for positively charged residues
            - -1 for negatively charged residues
            - 0 for other residues
            
        Example:
        --------
        >>> pos_counts = np.array([3, 2])
        >>> neg_counts = np.array([2, 4])
        >>> seq_lengths = np.array([10, 15])
        >>> result = populate_charged_residues(pos_counts, neg_counts, seq_lengths)
        # First sequence has 3 positive and 2 negative charges
        # Second sequence has 2 positive and 4 negative charges
        """
        # Convert inputs to numpy arrays
        pos_counts = np.asarray(pos_counts)
        neg_counts = np.asarray(neg_counts)
        seq_lengths = np.asarray(seq_lengths)
        
        # Validate inputs
        n_sequences = len(seq_lengths)
        if len(pos_counts) != n_sequences or len(neg_counts) != n_sequences:
            raise ValueError("Input arrays must have the same length")
        
        if np.any(pos_counts + neg_counts > seq_lengths):
            raise ValueError("Sum of charged residues cannot exceed sequence length")
        
        result = []
        
        for i in range(n_sequences):
            length = seq_lengths[i]
            pos_count = pos_counts[i]
            neg_count = neg_counts[i]
            
            # Create an array of zeros for this sequence
            sequence = np.zeros(length, dtype=int)
            
            # Total charged residues
            total_charged = pos_count + neg_count
            
            if total_charged > 0:
                # Randomly select positions for charged residues without replacement
                charged_positions = np.random.choice(length, total_charged, replace=False)
                
                # Set positive charges
                sequence[charged_positions[:pos_count]] = 1
                
                # Set negative charges
                sequence[charged_positions[pos_count:pos_count + neg_count]] = -1
            
            result.append(sequence)
        
        return result

    def _generate_uncharged_residues(self, lengths, probabilities):
        """
        Efficiently generate uncharged amino acid indices using vectorized operations.
        
        Parameters:
        -----------
        lengths : array-like
            Array containing the number of uncharged residues needed for each sequence
        probabilities : array-like
            2D array of probabilities for each sequence, with shape (num_sequences, 20)
            
        Returns:
        --------
        list of numpy.ndarray
            List where each element is an array of amino acid indices for uncharged positions
        """
        results = []
        
        for i in range(len(lengths)):
            length = lengths[i]
            probs = probabilities[i]
            
            # Extract and normalize probabilities for uncharged amino acids
            uncharged_probs = np.array([probs[idx] for idx in self.uncharged_residues_indices])
            uncharged_probs = uncharged_probs / np.sum(uncharged_probs)
            
            # Generate uncharged positions using vectorized choice
            uncharged_residues = np.random.choice(
                self.uncharged_residues_indices,
                size=length,
                p=uncharged_probs
            )
            
            results.append(uncharged_residues)
            
        return results

    def _choose_probabilities(self,
                            target_hydropathy,
                            target_fcr,
                            target_ncpr,
                            length):
        """
        Function to choose correction dictionary for probabilities based on the specified properties.
        """
        # get the hydropathy value needed for the rest of the sequence
        # to get the target hydropathy value with the specified charge.
        remaining_hydropathy = SequenceParameters._get_necessary_hydropathy_value_with_charge(
            sequence_length=length,
            target_hydropathy=target_hydropathy,
            target_fcr=target_fcr,
            target_ncpr=target_ncpr
        )
        # get the probabilities 
        dict_val = np.clip(remaining_hydropathy, 1, 9)
        return aa_probs.NeutralHydroDict[float(round(dict_val, 1))]

    def generate_sequences_vectorized(self, 
                                     length=None,
                                     fcr=None,
                                     ncpr=None,
                                     hydropathy=None,
                                     num_sequences=None,
                                     specific_probabilities=None,
                                     exclude_residues=None,
                                     convert_to_amino_acids=False):
        """
        Generate amino acid sequences using fully vectorized numpy operations.
        Uses a fixed sequence length for all generated sequences to enable better vectorization.
        
        Parameters:
        -----------
        length : int, optional
            Length of sequences to generate (same for all sequences)
        fcr : float, optional
            Target fraction of charged residues
        ncpr : float, optional
            Target net charge per residue
        hydropathy : float, optional
            Target hydropathy score
        num_sequences : int, optional
            Number of sequences to generate
        specific_probabilities : dict, optional
            Specific probabilities for amino acids (if not using default)
        exclude_residues : list, optional
            List of residues to exclude from generation
        convert_to_amino_acids : bool, optional
            If True, convert generated indices to amino acid characters
        Returns:
        --------
        list
            List of generated amino acid sequences, all with the same length
        """
        # Ensure length is specified for uniform sequence lengths
        if length is None:
            if self.length is None:
                raise ValueError("Length must be specified for sequence generation.")
            else:
                length = self.length
                
        if num_sequences is None:
            num_sequences = self.num_sequences if self.num_sequences is not None else 1
        
        # if exclude_residues is None, use the instance variable
        if exclude_residues is None:
            exclude_residues = self.exclude_residues

        # set sequence_indices to None
        sequences_array=None

        if specific_probabilities is not None:
            # make sure that all amino acids are in the specific probabilities
            aa_to_int={'A': 0, 'C': 1, 'D': 2, 'E': 3, 'F': 4, 'G': 5, 'H': 6, 'I': 7, 'K': 8, 'L': 9, 'M': 10, 'N': 11, 'P': 12, 'Q': 13, 'R': 14, 'S': 15, 'T': 16, 'V': 17, 'W': 18, 'Y': 19}
            int_to_aa={0: 'A', 1: 'C', 2: 'D', 3: 'E', 4: 'F', 5: 'G', 6: 'H', 7: 'I', 8: 'K', 9: 'L', 10: 'M', 11: 'N', 12: 'P', 13: 'Q', 14: 'R', 15: 'S', 16: 'T', 17: 'V', 18: 'W', 19: 'Y'}
            final_probs={}
            # check if the specific probabilities are in the aa_to_int dictionary
            if list(specific_probabilities.keys())[0] in list(aa_to_int.keys()):
                for aa in aa_to_int:
                    if aa not in exclude_residues:
                        if aa in specific_probabilities.keys():
                            final_probs[aa_to_int[aa]] = specific_probabilities[aa]    
                        # otherwise set to 0
                        else:
                            final_probs[aa_to_int[aa]] = 0.0
                    else:
                        final_probs[aa_to_int[aa]] = 0.0
            else:
                for aa in int_to_aa:
                    if int_to_aa[aa] not in exclude_residues:
                        if aa in specific_probabilities.keys():
                            final_probs[aa] = specific_probabilities[aa]
                        # otherwise set to 0
                        else:
                            final_probs[aa] = 0.0
                    else:
                        final_probs[aa] = 0.0
            # now set specific_probabilities to final_probs
            specific_probabilities = final_probs

            # Always set charged residue probabilities to 0.0 since charged residues
            # are selected directly and not from the probability distribution
            # D=2, E=3, K=8, R=14
            for aa in [2, 3, 8, 14]:
                specific_probabilities[aa] = 0.0
            # now make sure probabilities sum to 1.0
            total = sum(specific_probabilities.values())
            for aa in specific_probabilities.keys():
                specific_probabilities[aa] = specific_probabilities[aa]/total
            
            # get the probabilities
            probs_list = [specific_probabilities[aa] for aa in range(0,20)]
            # get the probabilities for the sequences
            probabilities = np.array([probs_list for _ in range(num_sequences)])

            # if specific probabilities is not None and FCR, NCPR, and hydropathy are None, we want to
            # just generate sequences using the specific probabilities
            # Create 2D array where rows are sequences and columns are positions
            if fcr is None and ncpr is None and hydropathy is None:
                # Generate random indices for amino acids based on probabilities
                sequences_array = np.random.choice(
                    20, 
                    size=(num_sequences, length), 
                    p=probabilities[0]
                )
            
        # if we didn't make the sequences just using the specific probabilities...
        if sequences_array is None:
            # Generate parameters for the batch
            if specific_probabilities is None:
                num_positive, num_negative, _, probabilities = self._generate_parameter_batch(
                    length=length, fcr=fcr, ncpr=ncpr, hydropathy=hydropathy,
                    num_sequences=num_sequences,
                    exclude_residues=exclude_residues
                )
            else:
                num_positive, num_negative, _, _ = self._generate_parameter_batch(
                    length=length, fcr=fcr, ncpr=ncpr, hydropathy=hydropathy,
                    num_sequences=num_sequences,
                    exclude_residues=exclude_residues
                )

            # Create 2D array where rows are sequences and columns are positions
            sequences_array = np.zeros((num_sequences, length), dtype=int)
            
            # Generate charge placements (1 for positive, -1 for negative, 0 for neutral)
            for i in range(num_sequences):
                # Create sequence template
                seq = np.zeros(length, dtype=int)
                
                # Total charged positions
                total_charged = num_positive[i] + num_negative[i]
                
                if total_charged > 0:
                    # Randomly select positions for charges
                    charged_positions = np.random.choice(length, total_charged, replace=False)

                    # Set positive charges
                    if num_positive[i] > 0:
                        seq[charged_positions[:num_positive[i]]] = 1
                    
                    # Set negative charges
                    if num_negative[i] > 0:
                        seq[charged_positions[num_positive[i]:total_charged]] = -1
                
                sequences_array[i] = seq

            # Create local copies of charge indices for this method
            local_positive_indices = self.positive_indices.copy()
            local_negative_indices = self.negative_indices.copy()
            
            # Handle excluded charged residues by creating restricted sets
            if 'E' in exclude_residues:
                local_negative_indices = np.array([2])  # Only D
            if 'D' in exclude_residues:
                local_negative_indices = np.array([3])  # Only E
            if 'R' in exclude_residues:
                local_positive_indices = np.array([8])  # Only K
            if 'K' in exclude_residues:
                local_positive_indices = np.array([14])  # Only R
            
            # Validate that we have available residues for required charges
            if np.any(num_positive > 0) and len(local_positive_indices) == 0:
                raise ValueError("Cannot generate positively charged residues: all positive residues are excluded")
            if np.any(num_negative > 0) and len(local_negative_indices) == 0:
                raise ValueError("Cannot generate negatively charged residues: all negative residues are excluded")

            # Now fill in the actual amino acids based on charge patterns
            for i in range(num_sequences):
                # Get masks for different charge types
                pos_mask = sequences_array[i] == 1
                neg_mask = sequences_array[i] == -1
                neutral_mask = sequences_array[i] == 0

                # Fill positive positions with positive amino acids
                pos_count = np.sum(pos_mask)
                if pos_count > 0:
                    sequences_array[i, pos_mask] = np.random.choice(
                        local_positive_indices, size=pos_count
                    )
                    
                # Fill negative positions with negative amino acids
                neg_count = np.sum(neg_mask)
                if neg_count > 0:
                    sequences_array[i, neg_mask] = np.random.choice(
                        local_negative_indices, size=neg_count
                    )
                    
                # Fill neutral positions using probability distribution
                neutral_count = np.sum(neutral_mask)
                if neutral_count > 0:
                    # Get probabilities for this sequence
                    probs = probabilities[i]
                    
                    # Extract and normalize probabilities for uncharged amino acids
                    uncharged_probs = np.array([probs[idx] for idx in self.uncharged_residues_indices])
                    uncharged_probs = uncharged_probs / np.sum(uncharged_probs)
                    
                    # Choose neutral amino acids based on probabilities
                    neutral_choices = np.random.choice(
                        self.uncharged_residues_indices,
                        size=neutral_count,
                        p=uncharged_probs
                    )
                    
                    sequences_array[i, neutral_mask] = neutral_choices
        
        # Create lookup array for amino acid conversion
        max_idx = max(self.int_to_aa.keys())
        
        # Convert integer arrays to amino acid sequences using vectorized join
        if convert_to_amino_acids:
            sequences = np.array([''.join(self.int_to_aa[aa] for aa in seq) for seq in sequences_array])
        else:
            sequences = sequences_array
        
        return sequences

