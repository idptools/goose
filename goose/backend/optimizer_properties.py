from abc import ABC, abstractmethod
from typing import Any, Tuple, Dict, Callable, Union, List
from enum import Enum
import numpy as np
import statistics as stat
import logging

import sparrow
import metapredict as meta
import finches

from goose.backend.optimizer_tools import MatrixManipulation, VectorManipulation

class ConstraintType(Enum):
    """Enumeration for different constraint types in property optimization."""
    EXACT = "exact"      # Minimize absolute difference from target
    MINIMUM = "minimum"  # Penalize only when below target
    MAXIMUM = "maximum"  # Penalize only when above target

    # add a string representation for easy printing
    def __str__(self):
        return self.value

def _normalize_constraint_type(constraint_type) -> ConstraintType:
    """
    Normalize constraint_type input to ConstraintType enum.
    Follows GOOSE's parameter normalization pattern.
    
    Parameters
    ----------
    constraint_type : str, ConstraintType, or None
        The constraint type specification
        
    Returns
    -------
    ConstraintType
        The normalized constraint type enum
        
    Raises
    ------
    ValueError
        If the constraint type is not recognized
    """
    if constraint_type is None:
        return ConstraintType.EXACT
    
    if isinstance(constraint_type, ConstraintType):
        return constraint_type
    
    if isinstance(constraint_type, str):
        # Normalize string input (case-insensitive)
        constraint_str = constraint_type.lower().strip()
        
        # Handle common aliases following GOOSE's flexible parameter handling
        constraint_mapping = {
            'exact': ConstraintType.EXACT,
            'equal': ConstraintType.EXACT,
            'equals': ConstraintType.EXACT,
            'match': ConstraintType.EXACT,
            
            'minimum': ConstraintType.MINIMUM,
            'min': ConstraintType.MINIMUM,
            'at_least': ConstraintType.MINIMUM,
            'greater_than': ConstraintType.MINIMUM,
            
            'maximum': ConstraintType.MAXIMUM,
            'max': ConstraintType.MAXIMUM,
            'at_most': ConstraintType.MAXIMUM,
            'less_than': ConstraintType.MAXIMUM,
        }
        
        if constraint_str in constraint_mapping:
            return constraint_mapping[constraint_str]
        else:
            valid_options = list(constraint_mapping.keys())
            raise ValueError(f"Invalid constraint_type '{constraint_type}'. Valid options are: {valid_options}")
    
    raise TypeError(f"constraint_type must be a string or ConstraintType enum, got {type(constraint_type)}")

class ProteinProperty(ABC):
    """
    Abstract base class for protein sequence properties to be optimized.
    """
    # the multi_target attribute let's us specify which things can actually be specified towards multiple values. 
    # for example, FCR should only be specified once. However, something like epsilon towards a target 
    # protein could be specified towards multiple targets and therefore we want this to be True for those classes. 
    # default is False
    multi_target = True  # Class-level attribute - override to limit ability to have multiple targets.
    can_be_linear_profile = True  # Class-level attribute - override to disallow linear profiles
    allow_setting_target_by_sequence = True  # Class-level attribute - override to disallow setting target by sequence.
    can_set_target_sequence_and_target_value = False # whether you can set the target sequence and target value at the same time.
    calculate_in_batch = False  # Whether this property can be calculated in batch mode (default False)
    def __init__(self, target_value: float = None, weight: float = 1.0, 
                 constraint_type: ConstraintType = ConstraintType.EXACT,
                 window_size: int = 5, end_mode: str = 'extend',
                 calculate_as_linear_profile: bool = False, target_sequence: str = None,
                 target_profile: np.ndarray = None):
        
        # make sure that if not linear profile, target_value is set. 
        if calculate_as_linear_profile is False:
            if target_value is None and target_sequence is None:
                    raise ValueError("When calculate_as_linear_profile is False, a target_sequence or target_value must be provided.")        
        
        if calculate_as_linear_profile is True:
            if target_sequence is None and target_profile is None:
                raise ValueError("When calculate_as_linear_profile is True, either target_sequence or target_profile must be provided.")
        
        # Validate end_mode
        valid_end_modes = ['extend', 'clip', 'zeros']
        if end_mode not in valid_end_modes:
            raise ValueError(f"Invalid end_mode '{end_mode}'. Valid options are: {valid_end_modes}")
    
        # Validate inputs    
        assert isinstance(weight, (int,float)), f"Weight must be numerical. Received {type(weight)}"
        if self.can_be_linear_profile is False and calculate_as_linear_profile is True:
            raise ValueError(f"{self.__class__.__name__} cannot be calculated as a linear profile.")
        if calculate_as_linear_profile is True and target_profile is None and target_sequence is None:
            raise ValueError("When calculate_as_linear_profile is True, either target_profile or target_sequence must be provided.")
        
        # handle constraint type if user inputs a string instead of Enum
        self.constraint_type = _normalize_constraint_type(constraint_type)
        self._target_value = target_value
        self._weight = weight
        
        # Property tracking attributes with clear names
        self.current_raw_value = None          # Current property value being evaluated
        self.best_raw_value = None             # Property value that gave the best error
        self.best_raw_error = None             # Best raw error achieved
        self.best_weighted_error = None        # Best weighted error achieved
        self.iterations_since_improvement = 0   # Track stagnation
        
        # Initialization values for scaling (though scaling is removed, kept for compatibility)
        self.initial_value = None  # Will be set during optimization setup

        # track scaling value. Start with set to 1. 
        self._scaling_value = 1.0

        # tolerance that can be overridden by the user.
        self._tolerance = 0.0

        # window size for linear profiles
        self.window_size = window_size
        self.end_mode = end_mode  # How to handle ends in linear profiles
        self.calculate_as_linear_profile = calculate_as_linear_profile
        
        # Placeholder for target profile if needed. Allows caching.
        self.target_profile = target_profile  # also allows user to set a target profile if they want to.
        self.target_resized = False # track that target has been resized.
        self.target_sequence = target_sequence # sequence to calculate target profile from if needed.

        # set target value if needed.
        self.set_target_value()

    def set_target_value(self):
        '''
        override target value depending on what is set
        '''
        if self.target_value is not None and self.target_sequence is not None:
            if self.can_set_target_sequence_and_target_value is False:
                raise ValueError("You cannot set both target_value and target_sequence at the same time.")
        if self.target_value is None:
            if self.allow_setting_target_by_sequence is False:    
                raise ValueError("You must provide a target_value.")
            if self.calculate_as_linear_profile:
                self.target_value=0.0
            else:
                if self.target_sequence is None:
                    raise ValueError("You must provide a target_sequence or target_value.")
                else:
                    self.target_value = self.calculate_raw_value(sparrow.Protein(self.target_sequence))
                # else, target_value is already set. Do nothing.
        else:
            assert isinstance(self.target_value, (int,float)), f"target_value must be numerical. Received {type(self.target_value)}"

    
    def set_initial_value(self, value: float) -> None:
        """
        Set the initial value for relative scaling.
        Following GOOSE's parameter validation patterns.
        
        Parameters
        ----------
        value : float
            The initial property value for scaling calculations
            
        Raises
        ------
        TypeError
            If value is not numeric
        ValueError
            If value is NaN or infinite
        """
        # Validate input type following GOOSE patterns
        if not isinstance(value, (int, float)):
            raise TypeError(f"initial_value must be numerical. Received {type(value)}")
        
        # Check for problematic values
        if not np.isfinite(value):
            raise ValueError(f"initial_value must be finite. Received {value}")
        
        self.initial_value = float(value)
    
    def _get_sequence_length_hint(self) -> int:
        """
        Get sequence length hint for length-dependent scaling calculations.
        
        Returns
        -------
        int
            Sequence length hint for calculations
        """
        # Use sequence length set by optimizer (when property is added)
        if hasattr(self, '_sequence_length_hint'):
            return self._sequence_length_hint
        
        # raise error if not set
        raise ValueError("Sequence length hint not set")

    def _set_sequence_length_hint(self, length: int) -> None:
        """Set sequence length hint from optimizer."""
        self._sequence_length_hint = length

    
    @property
    def constraint_type(self) -> ConstraintType:
        return self._constraint_type

    @constraint_type.setter
    def constraint_type(self, value: ConstraintType):
        # Flexible type setter that allows strings or enums as input
        self._constraint_type = _normalize_constraint_type(value)

    def get_init_args(self) -> dict:
        """
        Get the initialization arguments for this property instance.
        Override in subclasses that have additional parameters.
        
        Returns
        -------
        dict
            Dictionary containing the initialization arguments
        """
        vals = {
            "target_value": self.target_value,
            "weight": self.weight,
            "constraint_type": self.constraint_type.value,  # Store as string for JSON compatibility
            "window_size": self.window_size,
            "end_mode": self.end_mode,
            "calculate_as_linear_profile": self.calculate_as_linear_profile,
            "target_sequence": self.target_sequence,
            "target_profile": self.target_profile.tolist() if self.target_profile is not None else None
        }
        return vals

    @classmethod
    def from_init_args(cls, args_dict: dict):
        """
        Create a property instance from initialization arguments.
        
        Parameters
        ----------
        args_dict : dict
            Dictionary containing initialization arguments
            
        Returns
        -------
        ProteinProperty
            New instance of the property class
        """
        return cls(**args_dict)
        

    @property
    def target_value(self) -> float:
        return self._target_value

    @target_value.setter
    def target_value(self, value: float):
        assert isinstance(value, (int,float)), f"target_value must be numerical. Received {type(value)}"
        self._target_value = value

    @property
    def weight(self) -> float:
        return self._weight

    @weight.setter
    def weight(self, value: float):
        assert isinstance(value, (int,float)), f"weight must be numerical. Received {type(value)}"
        self._weight = value

    def calculate_error_batch(self, current_values: Union[List[float], np.ndarray]) -> np.ndarray:
        """
        Calculate errors for a batch of values using vectorized operations.
        Follows GOOSE's vectorized operation patterns for maximum performance.

        Parameters
        ----------
        current_values : Union[List[float], np.ndarray]
            Array or list of current property values to calculate errors for

        Returns
        -------
        np.ndarray
            Array of calculated errors (always positive) corresponding to each input value
            
        Raises
        ------
        ValueError
            If current_values is empty or contains invalid values
        TypeError
            If current_values contains non-numeric values
        """
        # Convert to numpy array for vectorized operations
        current_values = np.asarray(current_values, dtype=float)
        
        # Input validation following GOOSE patterns
        if current_values.size == 0:
            raise ValueError("current_values cannot be empty")
        
        # Check for invalid values (NaN, inf)
        if not np.all(np.isfinite(current_values)):
            invalid_indices = np.where(~np.isfinite(current_values))[0]
            raise ValueError(f"current_values contains invalid values (NaN/inf) at indices: {invalid_indices.tolist()}")
        
        # Vectorized error calculation based on constraint type
        if self.constraint_type == ConstraintType.EXACT:
            # Absolute difference from target (vectorized)
            errors = np.abs(current_values - self.target_value)
            
        elif self.constraint_type == ConstraintType.MINIMUM:
            # Penalty only when below minimum (vectorized)
            # Use np.maximum to handle element-wise max with 0
            errors = np.maximum(0.0, self.target_value - current_values)
            
        elif self.constraint_type == ConstraintType.MAXIMUM:
            # Penalty only when above maximum (vectorized)
            # Use np.maximum to handle element-wise max with 0
            errors = np.maximum(0.0, current_values - self.target_value)
            
        else:
            raise ValueError(f"Unknown constraint type: {self.constraint_type}")
        
        return errors

    def calculate_error(self, current_value: float) -> float:
        """
        Calculate the error between the current value and the target value based on the constraint type.
        Updated to use batch calculation for consistency.

        Parameters
        ----------
        current_value : float
            The current value of the property

        Returns
        -------
        float
            The calculated error (always positive)
        """
        # Use batch calculation for single value to ensure consistency
        # This also provides automatic validation
        batch_result = self.calculate_error_batch([current_value])
        return float(batch_result[0])


    @abstractmethod
    def calculate_raw_value(self, protein: 'sparrow.Protein') -> float:
        """
        Calculate the raw property value for a given protein.
        Subclasses should implement this instead of calculate().

        Parameters
        ----------
        protein : sparrow.Protein
            The protein instance

        Returns
        -------
        float
            The calculated raw property value
        """
        pass

    def calculate_linear_profile_error(self, current_value: np.ndarray) -> float:
        seq_length = self._get_sequence_length_hint()
        if self.constraint_type == ConstraintType.EXACT:
            return np.sum(np.abs(current_value - self.target_profile))/seq_length
        elif self.constraint_type == ConstraintType.MINIMUM:
            # only return the values below the minimum
            return np.sum(np.maximum(0, self.target_profile - current_value))/seq_length
        elif self.constraint_type == ConstraintType.MAXIMUM:
            # only return the values above the maximum
            return np.sum(np.maximum(0, current_value - self.target_profile))/seq_length
        else:
            raise ValueError(f"Unknown constraint type: {self.constraint_type}")
        
    def resize_linear_profile(self) -> np.ndarray:
        if self.target_resized == True:
            return self.target_profile
        else:
            # get seq length using the hint
            seq_length = self._get_sequence_length_hint()
            if self.end_mode == 'clip':
                seq_profile_length = (seq_length - self.window_size) + 1
            else:
                seq_profile_length = seq_length
            if len(self.target_profile) != seq_profile_length:
                # use vector manipulation to resize
                self.target_profile = VectorManipulation.resize_vector(self.target_profile, seq_profile_length)
            # set to True so we don't do this again.
            self.target_resized = True
            return self.target_profile

    def calculate_linear_profile(self, sequence: str) -> np.ndarray:
        """
        Calculate the linear profile for the property across the sequence.
        This is a simple moving window average of the raw property values.

        Parameters
        ----------
        sequence : str
            The protein sequence

        Returns
        -------
        np.ndarray
            The calculated linear profile as a numpy array
        """
        seq_length = len(sequence)
        end = (seq_length-self.window_size)+1
        profile_values=[]
        for i in range(end):
            window_seq = sequence[i:i+self.window_size]
            profile_values.append(self.calculate_raw_value(sparrow.Protein(window_seq)))
        
        # now take care of ends. 
        if self.end_mode == 'extend':
            pad_start = profile_values[0]
            pad_end = profile_values[-1]
        elif self.end_mode == 'zeros':
            pad_start = 0
            pad_end = 0
        else:
            pad_start=None
            pad_end=None
        if pad_start is not None and pad_end is not None:
            start = int(self.window_size/2)
            profile_values = [pad_start]*start + profile_values
            end = seq_length - len(profile_values)  
            profile_values = profile_values + [pad_end]*end
        return np.array(profile_values)

    def get_target_profile(self) -> np.ndarray:
        if self.target_profile is None:
            self.target_profile = self.calculate_linear_profile(self.target_sequence)
            # take care of resizing
            self.resize_linear_profile()
        return self.target_profile

    def calculate_batch(self, proteins: list) -> list:
        """
        Calculate the property for a batch of proteins.
        This is an optional method that can be implemented by subclasses for efficiency.

        Parameters
        ----------
        proteins : list of sparrow.Protein
            List of protein instances

        Returns
        -------
        list of tuples
            Each tuple contains (raw_value, error) for each protein
        """
        results = []
        batch_calculations = self.calculate_raw_value_batch(proteins)
        errors = self.calculate_error_batch(batch_calculations)
        for raw_value, error in zip(batch_calculations, errors):
            results.append((raw_value, error))
        return results

    def calculate(self, protein: 'sparrow.Protein') -> float:
        """
        Calculate the final error value for optimization.
        This applies the constraint type to the raw value.

        Parameters
        ----------
        protein : sparrow.Protein
            The protein instance

        Returns
        -------
        float
            The error value for optimization
        """
        if self.calculate_as_linear_profile:
            self.target_profile = self.get_target_profile()
            raw_profile = self.calculate_linear_profile(protein.sequence)
            error = self.calculate_linear_profile_error(raw_profile)
            average_raw_value = np.sum(raw_profile) / len(raw_profile)
            self._current_raw_value = average_raw_value
            return average_raw_value, error
        else:
            raw_value = self.calculate_raw_value(protein)
            error = self.calculate_error(raw_value)
            self._current_raw_value = raw_value
            return raw_value, error


class CustomProperty(ProteinProperty):
    """
    User-facing custom property class that inherits all functionality from ProteinProperty.
    Users should subclass this and implement calculate_raw_value().
    """
    multi_target = True  # Class-level attribute - override to limit ability to have multiple targets.
    can_be_linear_profile = True  # Class-level attribute - override to disallow linear profiles
    allow_setting_target_by_sequence = True  # Class-level attribute - override to disallow setting target by sequence.
    can_set_target_sequence_and_target_value = True  # whether you can set the target sequence and target value at the same time.
    calculate_in_batch = False  # Whether this property can be calculated in batch mode (default False)
    def __init__(self, target_value: float=None, target_sequence: str = None, weight: float = 1.0,
                 constraint_type: ConstraintType = ConstraintType.EXACT,
                 window_size: int = 5, end_mode: str = 'extend',
                 calculate_as_linear_profile: bool = False, target_profile: np.ndarray = None):

        super().__init__(target_value=target_value, weight=weight, constraint_type=constraint_type,
                         window_size=window_size, end_mode=end_mode,
                         calculate_as_linear_profile=calculate_as_linear_profile,
                         target_sequence=target_sequence, target_profile=target_profile)
        
    def calculate_raw_value(self, protein: 'sparrow.Protein') -> float:
        raise NotImplementedError("CustomProperty requires a user-defined calculation function.")
    

class ComputeIWD(ProteinProperty):
    """
    Compute the Inversed Weighted Distance (IWD) property for the target residues in the sequence.
    """
    can_be_linear_profile = False  # IWD does not make sense as a linear profile
    def __init__(self, residues: Tuple[str, ...], target_value: float=None, 
                 target_sequence: str = None, weight: float = 1.0,
                 constraint_type: ConstraintType = ConstraintType.EXACT):
        # Store the residues parameter before calling super().__init__
        self.residues = residues
        # Pass all other parameters to the parent class
        super().__init__(target_value, weight, constraint_type, target_sequence=target_sequence)
    
    def get_init_args(self) -> dict:
        """Override to include residues parameter"""
        base_args = super().get_init_args()
        base_args['residues'] = self.residues
        return base_args

    def calculate_raw_value(self, protein: 'sparrow.Protein') -> float:
        return protein.compute_iwd(list(self.residues))


class Hydrophobicity(ProteinProperty):
    """
    Calculate the hydrophobicity property.
    """
    def calculate_raw_value(self, protein: 'sparrow.Protein') -> float:
        return protein.hydrophobicity

class FCR(ProteinProperty):
    """
    Calculate the FCR property.
    """
    def calculate_raw_value(self, protein: 'sparrow.Protein') -> float:
        return protein.FCR

class NCPR(ProteinProperty):
    """
    Calculate the NCPR property.
    """
    def calculate_raw_value(self, protein: 'sparrow.Protein') -> float:
        return protein.NCPR


class Kappa(ProteinProperty):
    """
    Calculate the kappa property.
    """
    def calculate_raw_value(self, protein: 'sparrow.Protein') -> float:
        return protein.kappa

class RadiusOfGyration(ProteinProperty):
    """
    Calculate the Radius of gyration
    """
    def calculate_raw_value(self, protein: 'sparrow.Protein') -> float:
        return protein.predictor.radius_of_gyration()

class EndToEndDistance(ProteinProperty):
    """
    Calculate the Radius of gyration
    """
    def calculate_raw_value(self, protein: 'sparrow.Protein') -> float:
        return protein.predictor.end_to_end_distance()


class AminoAcidFractions(ProteinProperty):
    """
    Compute the difference between the target amino acid fractions and the current amino acid fractions.
    """
    can_be_linear_profile = False  # Amino acid fractions do not make sense as a linear profile
    allow_setting_target_by_sequence = False  # Must provide target_fractions explicitly
    def __init__(self, target_fractions: Dict[str, float], weight: float = 1.0,
                 constraint_type: ConstraintType = ConstraintType.EXACT):
        # Store the target_fractions parameter
        self.target_fractions = target_fractions
        # For fractions, target value is not meaningful at class level since we have multiple amino acids
        super().__init__(target_value=0.0, weight=weight, constraint_type=constraint_type)

    def get_init_args(self) -> dict:
        """Override to include target_fractions parameter"""
        base_args = super().get_init_args()
        base_args['target_fractions'] = self.target_fractions
        return base_args

    def calculate_error_for_fraction(self, current_frac: float, target_frac: float) -> float:
        """Calculate error for a single amino acid fraction based on constraint type."""
        if self.constraint_type == ConstraintType.EXACT:
            return abs(target_frac - current_frac)
        elif self.constraint_type == ConstraintType.MINIMUM:
            # Only penalize if below minimum
            return max(0, target_frac - current_frac)
        elif self.constraint_type == ConstraintType.MAXIMUM:
            # Only penalize if above maximum
            return max(0, current_frac - target_frac)
        else:
            raise ValueError(f"Unknown constraint type: {self.constraint_type}")

    def calculate_raw_value(self, protein: 'sparrow.Protein') -> float:
        """
        Calculate the total error across all amino acid fractions.
        This is the final error value, not a raw property value.
        """
        current_fractions = protein.amino_acid_fractions
        total_error = 0.0

        for aa, target_frac in self.target_fractions.items():
            current_frac = current_fractions.get(aa, 0.0)
            total_error += self.calculate_error_for_fraction(current_frac, target_frac)

        return total_error

    def calculate(self, protein: 'sparrow.Protein') -> float:
        """
        Override calculate to handle multi-target amino acid optimization.
        Now supports both single-value and linear profile modes.
        """
        # Use the original amino acid fraction calculation
        final_error = self.calculate_raw_value(protein)
        
        # Update tracking attributes following base class patterns
        self._current_raw_value = final_error
        
        # Return final error (no additional constraint processing needed)
        return final_error, final_error


class SCD(ProteinProperty):
    """
    Returns the default sequence charge decoration (SCD) parameter 
        as defined by Sawle and Ghosh [1]
    
    Reference
    --------
    Sawle, L., & Ghosh, K. (2015). A theoretical method to compute sequence
    dependent configurational properties in charged polymers and proteins.
    The Journal of Chemical Physics, 143(8), 085101.
    """    
    def calculate_raw_value(self, protein: 'sparrow.Protein') -> float:
        return protein.SCD


class SHD(ProteinProperty):
    """
    Returns the default sequence hydropathy decoration (SHD) parameter 
        as defined by Sawle and Ghosh [1]
 
    Reference
    --------
    Sawle, L., & Ghosh, K. (2015). A theoretical method to compute sequence
    dependent configurational properties in charged polymers and proteins.
    The Journal of Chemical Physics, 143(8), 085101.
    """
    def calculate_raw_value(self, protein: 'sparrow.Protein') -> float:
        return protein.SHD

class Complexity(ProteinProperty):
    """
    Calculates the Wootton-Federhen complexity of a sequence (also called
    seg complexity, as this the theory used in the classic SEG algorithm.
    """
    def calculate_raw_value(self, protein: 'sparrow.Protein') -> float:
        return protein.complexity


class FractionDisorder(ProteinProperty):
    """
    Calculates the amount of disorder in
    the current sequence. 
    """
    can_be_linear_profile = False  # Disorder fraction does not make sense as a linear profile
    calculate_in_batch = True  # Disorder can be calculated in batch mode
    def __init__(self, target_value: float=None, weight: float = 1.0, disorder_cutoff: float = 0.5,
                 constraint_type: ConstraintType = ConstraintType.EXACT,
                 target_sequence: str = None):
        # Store the disorder_cutoff parameter
        self.disorder_cutoff = disorder_cutoff
        # Pass all other parameters to the parent class
        super().__init__(target_value=target_value, weight=weight, constraint_type=constraint_type,
                         target_sequence=target_sequence)

    def get_init_args(self) -> dict:
        """Override to include disorder_cutoff parameter"""
        base_args = super().get_init_args()
        base_args['disorder_cutoff'] = self.disorder_cutoff
        return base_args

    def calculate_raw_value(self, protein: 'sparrow.Protein') -> float:
        disorder = meta.predict_disorder(protein.sequence)
        return (disorder > self.disorder_cutoff).sum() / len(disorder)
    
    def calculate_raw_value_batch(self, proteins: list) -> list:
        sequences = [p.sequence for p in proteins]
        disorder_with_seqs = meta.predict_disorder(sequences)
        disorder = np.array([d[1] for d in disorder_with_seqs])
        frac_disorder = (disorder > self.disorder_cutoff).sum(axis=1) / disorder.shape[1]
        results = frac_disorder.tolist()
        return results


class MatchSequenceDisorder(ProteinProperty):
    """
    Matches the disorder of a specified sequence. By default, the 
    specified sequence sets the MINIMUM disorder. 
    You can also try to get the exact disorder.

    """
    # This property already works with profiles inherently, so disable linear profile mode
    can_be_linear_profile = False
    can_set_target_sequence_and_target_value = True
    def __init__(self, target_sequence: str,
                 weight: float = 1.0, constraint_type: ConstraintType = ConstraintType.EXACT):
        # Store the specific parameters for this class
        self.target_sequence = target_sequence
        self.target_disorder = None
        # Pass all other parameters to the parent class
        super().__init__(target_value=0.0, weight=weight, constraint_type=constraint_type,
                         target_sequence=target_sequence)

    def get_init_args(self) -> dict:
        """Override to include target_sequence and exact_match parameters"""
        base_args = super().get_init_args()
        base_args.update({
            'target_sequence': self.target_sequence,
        })
        return base_args

    def set_initial_disorder(self):
        if self.target_disorder is None:
            self.target_disorder = meta.predict_disorder(self.target_sequence)
        return self.target_disorder

    def calculate_raw_value(self, protein: 'sparrow.Protein') -> float:
        target_disorder = self.set_initial_disorder()
        disorder = meta.predict_disorder(protein.sequence)
        if self.constraint_type == ConstraintType.EXACT:
            err = np.sum(np.abs(disorder - target_disorder))
        elif self.constraint_type == ConstraintType.MINIMUM:
            err = np.sum(np.maximum(0, target_disorder - disorder))
        elif self.constraint_type == ConstraintType.MAXIMUM:
            err = np.sum(np.maximum(0, disorder - target_disorder))
        else:
            raise ValueError(f"Unknown constraint type: {self.constraint_type}. Should not have made it this far, contact dev please. ")
        # normalize to length. 
        return err / len(target_disorder)


class MatchingResidues(ProteinProperty):
    '''
    Determines the number of residues that match a target sequence in the current sequence.
    '''
    can_be_linear_profile = False  # Matching residues does not make sense as a linear profile
    allow_setting_target_by_sequence = False  # Must provide target_sequence explicitly
    can_set_target_sequence_and_target_value = True # you can set both target sequence and target
    def __init__(self, target_sequence: str, target_value: float, weight: float = 1.0,
                 constraint_type: ConstraintType = ConstraintType.EXACT,
                 input_fraction=False):
        # Validate required parameters before calling super().__init__
        if target_sequence is None:
            raise ValueError("A target_sequence must be provided for MatchingResidues property.")
        if target_value is None:
            raise ValueError("A target_value must be provided for MatchingResidues property.")
        
        # Pass target_value and all other parameters directly to parent
        super().__init__(target_value=target_value, weight=weight, constraint_type=constraint_type,
                         target_sequence=target_sequence)
        # see if we've verified length
        self.verified_length = False
        if input_fraction:
            if self.target_value > 1.0 or self.target_value < 0.0:
                raise ValueError("If input_fraction is True, target_value must be between 0 and 1.")
            self.target_value = int(round(self.target_value * len(self.target_sequence)))

    def calculate_raw_value(self, protein: 'sparrow.Protein') -> float:
        if self.verified_length is False:
            if self.target_value > len(protein.sequence):
                raise ValueError(f"The target_value is greater than the length of the protein sequence. Please use a value less than or equal to {len(protein.sequence)}.")
            self.verified_length = True
        return sum([1 for i in range(len(protein.sequence)) 
                    if protein.sequence[i] == self.target_sequence[i]])

#=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=    
#=-=-=-=-=-=- EpsilonProperty base class for relatively simple epsilon properties -=-=-=-=-=-=
#=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=

class EpsilonProperty(ProteinProperty):
    """
    Abstract base class for properties that use epsilon calculations.
    Provides common model loading functionality following GOOSE patterns.
    """
    can_be_linear_profile = False
    def __init__(self, target_value: float=None, target_sequence: str = None, weight: float = 1.0, 
                 constraint_type: ConstraintType = ConstraintType.EXACT, 
                 model: str = 'mpipi', preloaded_model=None,
                 window_size: int = 5, end_mode: str = 'extend',
                 calculate_as_linear_profile: bool = False, target_profile: np.ndarray = None):
        # Store epsilon-specific parameters
        self.model = model.lower()
        self._loaded_model = preloaded_model
        self._loaded_imc_object = None
        
        # Handle preloaded model type detection following GOOSE's validation patterns
        if preloaded_model is not None:
            self._detect_preloaded_model_type(preloaded_model)
        
        # Pass all other parameters to the parent class
        super().__init__(target_value=target_value, weight=weight, constraint_type=constraint_type,
                         window_size=window_size, end_mode=end_mode,
                         calculate_as_linear_profile=calculate_as_linear_profile,
                         target_sequence=target_sequence, target_profile=target_profile)

    def get_init_args(self) -> dict:
        """Override to include epsilon-specific parameters"""
        base_args = super().get_init_args()
        base_args.update({
            'model': self.model
        })
        return base_args

    def _detect_preloaded_model_type(self, preloaded_model):
        """
        Detect and set model type from preloaded model instance.
        Follows GOOSE's parameter validation patterns.
        """
        full_model_name = preloaded_model.__class__.__name__
        if full_model_name == 'CALVADOS_frontend':
            self.model = 'calvados'
        elif full_model_name == 'Mpipi_frontend':
            self.model = 'mpipi'
        else:
            # Following GOOSE's error handling patterns
            raise ValueError(f"Unsupported preloaded model type: {full_model_name}")

    def load_model(self, model: str = None):
        """
        Load the epsilon model following GOOSE's lazy loading patterns.
        
        Parameters
        ----------
        model : str, optional
            Model type ('mpipi' or 'calvados'). Uses instance model if None.
            
        Returns
        -------
        Frontend model instance
            The loaded model ready for epsilon calculations
            
        Raises
        ------
        ValueError
            If model type is not supported
        """
        if self._loaded_model is not None:
            return self._loaded_model
            
        # Use instance model if none specified
        if model is None:
            model = self.model
        else:
            model = model.lower()
        
        # Import and load following GOOSE's dynamic import patterns
        if model == 'mpipi':
            self._loaded_model = finches.frontend.mpipi_frontend.Mpipi_frontend()
        elif model == 'calvados':
            self._loaded_model = finches.frontend.calvados_frontend.CALVADOS_frontend()
        else:
            raise ValueError(f"Unsupported model '{model}'. Valid options are: ['mpipi', 'calvados']")
            
        return self._loaded_model

    def load_imc_object(self, model: str = None):
        """
        Load the IMC object for surface calculations.
        
        Parameters
        ----------
        model : str, optional
            Model type. Uses instance model if None.
            
        Returns
        -------
        IMC object
            The loaded IMC object for surface epsilon calculations
        """
        if self._loaded_imc_object is not None:
            return self._loaded_imc_object
            
        # Load the full model first
        loaded_model = self.load_model(model)
        self._loaded_imc_object = loaded_model.IMC_object
        return self._loaded_imc_object


class MeanSelfEpsilon(EpsilonProperty):
    """
    Calculate the self interaction epsilon value of a sequence.
    Note: this simply uses the mean epsilon value. 
    """
    def calculate_raw_value(self, protein: 'sparrow.Protein') -> float:
        loaded_model = self.load_model()
        return loaded_model.epsilon(protein.sequence, protein.sequence)


class MeanEpsilonWithTarget(EpsilonProperty):
    """
    Make a sequence that interacts with a specific sequence with a specific mean epsilon value. 
    """
    can_set_target_sequence_and_target_value = True
    def __init__(self, target_sequence: str, target_value: float=None, weight: float = 1.0, 
                 constraint_type: ConstraintType = ConstraintType.EXACT, model: str = 'mpipi', preloaded_model = None):
        super().__init__(target_value, target_sequence, weight, constraint_type, model, preloaded_model)

    def calculate_raw_value(self, protein: 'sparrow.Protein') -> float:
        loaded_model = self.load_model()
        return loaded_model.epsilon(protein.sequence, self.target_sequence)



class ChemicalFingerprint(EpsilonProperty):
    """
    Uses the chemical foot print from the FINCHES manuscript to generate a sequence with
    a similar chemical fingerprint to the target sequence.

    The chemical fingerprint is calculated by taking the difference between the target
    and the current sequence and summing the differences in the epsilon matrix.
    """
    can_be_linear_profile = False
    can_set_target_sequence_and_target_value = True
    def __init__(self, target_sequence: str, target_value: float = 0.0, 
                 weight: float = 1.0, model: str = 'mpipi', preloaded_model = None, 
                 window_size:int = 15, constraint_type = ConstraintType.EXACT):
        super().__init__(target_value, target_sequence, weight, constraint_type, model, preloaded_model)
        self.target_sequence = target_sequence
        self._target_sequence_fingerprint = None
        self._chemistries = None
        self._matrix_shape=None
        self.window_size = window_size

    def get_init_args(self) -> dict:
        """Override to include target_sequence parameter"""
        base_args = super().get_init_args()
        base_args['target_sequence'] = self.target_sequence
        base_args['window_size'] = self.window_size
        return base_args

    def _get_chemistries(self) -> dict:
        """
        Get the chemistries to be used for the calculation.
        """
        if self._chemistries is None:
            if self.model == 'mpipi':
                self._chemistries = {'c1': 'KRKRKRKRKRKRKRKRKRKR', 'c2': 'RRRRRRRRRRRRRRRRRRRR', 'c3': 'KKKKKKKKKKKKKKKKKKKK', 'c4': 'HRHRHRHRHRHRHRHRHRHR', 'c5': 'HKHKHKHKHKHKHKHKHKHK', 'c6': 'AKAKAKAKAKAKAKAKAKAK', 'c7': 'KQKQKQKQKQKQKQKQKQKQ', 'c8': 'QRQRQRQRQRQRQRQRQRQR', 'c9': 'GRGRGRGRGRGRGRGRGRGR', 'c10': 'MRMRMRMRMRMRMRMRMRMR', 'c11': 'RWRWRWRWRWRWRWRWRWRW', 'c12': 'HHHHHHHHHHHHHHHHHHHH', 'c13': 'HYHYHYHYHYHYHYHYHYHY', 'c14': 'HQHQHQHQHQHQHQHQHQHQ', 'c15': 'HLHLHLHLHLHLHLHLHLHL', 'c16': 'AIAIAIAIAIAIAIAIAIAI', 'c17': 'AGAGAGAGAGAGAGAGAGAG', 'c18': 'ININININININININININ', 'c19': 'GQGQGQGQGQGQGQGQGQGQ', 'c20': 'GSGSGSGSGSGSGSGSGSGS', 'c21': 'DKDKDKDKDKDKDKDKDKDK', 'c22': 'ERERERERERERERERERER', 'c23': 'DHDHDHDHDHDHDHDHDHDH', 'c24': 'QQQQQQQQQQQQQQQQQQQQ', 'c25': 'FGFGFGFGFGFGFGFGFGFG', 'c26': 'QYQYQYQYQYQYQYQYQYQY', 'c27': 'AFAFAFAFAFAFAFAFAFAF', 'c28': 'LWLWLWLWLWLWLWLWLWLW', 'c29': 'FYFYFYFYFYFYFYFYFYFY', 'c30': 'WWWWWWWWWWWWWWWWWWWW', 'c31': 'DWDWDWDWDWDWDWDWDWDW', 'c32': 'EYEYEYEYEYEYEYEYEYEY', 'c33': 'EGEGEGEGEGEGEGEGEGEG', 'c34': 'EPEPEPEPEPEPEPEPEPEP', 'c35': 'AEAEAEAEAEAEAEAEAEAE', 'c36': 'DEDEDEDEDEDEDEDEDEDE'}
            elif self.model == 'calvados':
                self._chemistries = {'c1': 'DEDEDEDEDEDEDEDEDEDE', 'c2': 'EGEGEGEGEGEGEGEGEGEG', 'c3': 'DSDSDSDSDSDSDSDSDSDS', 'c4': 'EVEVEVEVEVEVEVEVEVEV', 'c5': 'DMDMDMDMDMDMDMDMDMDM', 'c6': 'EWEWEWEWEWEWEWEWEWEW', 'c7': 'AVAVAVAVAVAVAVAVAVAV', 'c8': 'APAPAPAPAPAPAPAPAPAP', 'c9': 'AQAQAQAQAQAQAQAQAQAQ', 'c10': 'GQGQGQGQGQGQGQGQGQGQ', 'c11': 'HQHQHQHQHQHQHQHQHQHQ', 'c12': 'QQQQQQQQQQQQQQQQQQQQ', 'c13': 'DRDRDRDRDRDRDRDRDRDR', 'c14': 'DKDKDKDKDKDKDKDKDKDK', 'c15': 'LMLMLMLMLMLMLMLMLMLM', 'c16': 'IMIMIMIMIMIMIMIMIMIM', 'c17': 'LWLWLWLWLWLWLWLWLWLW', 'c18': 'FMFMFMFMFMFMFMFMFMFM', 'c19': 'AWAWAWAWAWAWAWAWAWAW', 'c20': 'IQIQIQIQIQIQIQIQIQIQ', 'c21': 'LQLQLQLQLQLQLQLQLQLQ', 'c22': 'HLHLHLHLHLHLHLHLHLHL', 'c23': 'AMAMAMAMAMAMAMAMAMAM', 'c24': 'FQFQFQFQFQFQFQFQFQFQ', 'c25': 'HWHWHWHWHWHWHWHWHWHW', 'c26': 'FWFWFWFWFWFWFWFWFWFW', 'c27': 'MRMRMRMRMRMRMRMRMRMR', 'c28': 'RWRWRWRWRWRWRWRWRWRW', 'c29': 'KYKYKYKYKYKYKYKYKYKY', 'c30': 'ARARARARARARARARARAR', 'c31': 'KLKLKLKLKLKLKLKLKLKL', 'c32': 'AKAKAKAKAKAKAKAKAKAK', 'c33': 'KSKSKSKSKSKSKSKSKSKS', 'c34': 'KKKKKKKKKKKKKKKKKKKK', 'c35': 'KRKRKRKRKRKRKRKRKRKR', 'c36': 'RRRRRRRRRRRRRRRRRRRR'}
            else:
                raise ValueError(f"Model {self.model} not supported.")
        return self._chemistries

    def _calculate_matrix(self, sequence1: str, sequence2: str, window_size: int) -> np.ndarray:
        """Calculate epsilon matrix for a sequence"""
        loaded_model = self.load_model()
        return loaded_model.intermolecular_idr_matrix(
            sequence1, sequence2,
            disorder_1=False, disorder_2=False, window_size=window_size
        )[0][0]

    def calculate_fingerprint(self, sequence: str):
        """
        Calculate the chemical fingerprint of the target sequence.
        """
        chemistries = self._get_chemistries()

        # now for each chemistry, we calculate the epsilon matrix. 
        sequence_fingerprint = {}
        for chemistry, chemistry_sequence in chemistries.items():
            sequence_fingerprint[chemistry] = self._calculate_matrix(sequence, chemistry_sequence, self.window_size)
        return sequence_fingerprint

    def get_target_fingerprint(self):
        """
        Calculate the chemical fingerprint of the target sequence.
        """
        if self._target_sequence_fingerprint is None:
            self._target_sequence_fingerprint = self.calculate_fingerprint(self.target_sequence)
            # make sure all matrices are the same shape as the current matrix
            for key in self._target_sequence_fingerprint.keys():
                if self._target_sequence_fingerprint[key].shape != self._matrix_shape:
                    self._target_sequence_fingerprint[key] = MatrixManipulation.scale_matrix_to_size(
                        self._target_sequence_fingerprint[key], self._matrix_shape)
                    
        return self._target_sequence_fingerprint
    
    def calculate_raw_value(self, protein: 'sparrow.Protein') -> float:
        current_fingerprint = self.calculate_fingerprint(protein.sequence)
        # get shape of one of the matrices in the current fingerprint.
        if self._matrix_shape is None:
            self._matrix_shape = next(iter(current_fingerprint.values())).shape
        # get the target fingerprint and rescale if necessary. 
        target_fingerprint = self.get_target_fingerprint()
        # calculate the difference between the original and current matrix
        diff = 0
        for key in target_fingerprint.keys():
            diff += np.sum(np.abs(target_fingerprint[key] - current_fingerprint[key]))/len(protein.sequence)
        return diff




#=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=    
#=- Epsilon matrix properties based on adjusting the interaction strength using the matrix. -=
#=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=

class EpsilonMatrixProperty(EpsilonProperty):
    """
    Abstract base class for properties that use epsilon matrix calculations between
    a protein sequence and a target sequence with matrix manipulation capabilities.

    To do a self interaction, set target_sequence to None. This should be done in the
    subclass if you want to do this with self interaction. 

    The only thing you need to set is _manipulate_matrix. This will let you manipulate
    the matrix to (for example) increase the interaction strength.

    """
    can_be_linear_profile = False  # Class-level attribute - override to disallow linear profiles
    can_set_target_sequence_and_target_value = True  # You can set both target sequence and target value
    def __init__(self, sequence: str, target_sequence: str, target_value: float = 0.0,
                 weight: float = 1.0, constraint_type = ConstraintType.EXACT,
                 model: str = 'mpipi', preloaded_model = None,
                 window_size: int = 15, allow_matrix_resizing=True,
                 homotypic_interaction: bool = False):
        super().__init__(target_value, target_sequence, weight, constraint_type, model, preloaded_model)
        self.sequence = sequence
        self.target_sequence = target_sequence
        self.window_size = window_size
        self.target_matrix = None  # Cache for expensive matrix calculations
        self.allow_matrix_resizing = allow_matrix_resizing  # whether to allow matrix resizing. If not allowed,
        self.homotypic_interaction = homotypic_interaction  # whether to do homotypic interaction (same sequence)
        self.window_size_validated=False

    def _verify_window_size(self):
        """Verify and adjust window size if necessary following GOOSE validation patterns."""
        # Get sequence lengths with proper None handling
        generated_length = self._get_sequence_length_hint()
        
        if self.homotypic_interaction:
            min_length = generated_length
        else:
            # Handle sequence lengths with validation
            sequence_length = len(self.sequence) if self.sequence is not None else None
            target_length = len(self.target_sequence) if self.target_sequence is not None else None
            
            # Collect valid lengths for minimum calculation
            valid_lengths = [generated_length]
            if sequence_length is not None:
                valid_lengths.append(sequence_length)
            if target_length is not None:
                valid_lengths.append(target_length)
            
            min_length = min(valid_lengths)
        
        # Validate window size constraints following GOOSE patterns
        if self.window_size < 1:
            raise ValueError("Window size must be at least 1")
        
        if self.window_size >= min_length:
            raise ValueError(f"Window size ({self.window_size}) must be less than minimum sequence length ({min_length})")
        
        # Validate odd window size requirement
        if self.window_size % 2 == 0:
            raise ValueError("Window size must be an odd value")

    def get_init_args(self) -> dict:
        """Override to include matrix-specific parameters"""
        base_args = super().get_init_args()
        base_args.update({
            "sequence": self.sequence,
            "target_sequence": self.target_sequence,
            "window_size": self.window_size,
            "allow_matrix_resizing": self.allow_matrix_resizing,
            "homotypic_interaction": self.homotypic_interaction
        })
        return base_args
    
    def _calculate_matrix(self, sequence1: str, sequence2: str) -> np.ndarray:
        """
        Calculate epsilon matrix between two sequences.
        Uses cached matrix for performance if available.
        
        Parameters
        ----------
        sequence1 : str
            First sequence for matrix calculation
        sequence2 : str
            Second sequence for matrix calculation
            
        Returns
        -------
        np.ndarray
            Epsilon interaction matrix
        """
        loaded_model = self.load_model()
        return loaded_model.intermolecular_idr_matrix(sequence1, sequence2,
                window_size=self.window_size, disorder_1=False, disorder_2=False)[0][0]


    def _manipulate_matrix(self, matrix: np.ndarray) -> np.ndarray:
        """
        Apply matrix manipulations before comparison.
        Override in subclasses to implement specific transformations.
        """
        # Default: return matrix unchanged
        return matrix

    def _match_matrix_size(self) -> np.ndarray:
        """
        Match the size of the matrix to the target sequence length.
        """
        target_matrix_size = self.target_matrix.shape
        generated_seq_length = self._get_sequence_length_hint()

        # Determine the size of the generated matrix based on interaction type
        if self.homotypic_interaction:
            generated_matrix_size = (generated_seq_length - self.window_size + 1, generated_seq_length - self.window_size + 1)
        else:
            generated_matrix_size = (generated_seq_length - self.window_size + 1, len(self.target_sequence) - self.window_size + 1)
        # verify size is same otherwise modify. 
        if target_matrix_size != generated_matrix_size:
            if not self.allow_matrix_resizing:
                raise ValueError("Matrix resizing is not allowed.")
            else:
                # resize target_matrix to match shape of generated_matrix_size
                self.target_matrix = MatrixManipulation.scale_matrix_to_size(self.target_matrix, generated_matrix_size)
        return self.target_matrix

    def _initialize_target_matrix(self) -> np.ndarray:
        """
        Initialize the target matrix for comparison.
        Override in subclasses to implement specific initialization logic.

        Returns
        -------
        np.ndarray
            Initialized target matrix
        """
        if self.target_matrix is None:
            self.target_matrix = self._manipulate_matrix(self._calculate_matrix(
                self.sequence, self.target_sequence))
            # resize so target_matrix matches matrix shape produced from _calculate_matrix for generated sequence
            self.target_matrix = self._match_matrix_size()
        return self.target_matrix
        

    def _calculate_matrix_difference(self, target_matrix: np.ndarray, 
                                   current_matrix: np.ndarray,
                                   normalize_by_length: bool = True) -> float:
        """
        Calculate difference between matrices following GOOSE error patterns.
        
        Parameters
        ----------
        target_matrix : np.ndarray
            Target/reference matrix
        current_matrix : np.ndarray
            Current matrix to compare
        normalize_by_length : bool, optional
            Whether to normalize by sequence length (default: True)
            
        Returns
        -------
        float
            Matrix difference value for optimization
        """

        # Calculate absolute difference
        difference = np.sum(np.abs(target_matrix - current_matrix))
        
        # Normalize by sequence length following GOOSE patterns
        if normalize_by_length:
            sequence_length = min(target_matrix.shape[0], current_matrix.shape[0])
            difference = difference / sequence_length
            
        return difference

    def calculate_raw_value(self, protein: 'sparrow.Protein') -> float:
        """
        Calculate the raw matrix difference value.
        Template method that coordinates matrix calculation and comparison.
        
        Parameters
        ----------
        protein : sparrow.Protein
            The protein sequence being optimized
            
        Returns
        -------
        float
            Raw difference value for optimization
        """
        if self.window_size_validated == False:
            # verify window size
            self._verify_window_size()
            self.window_size_validated = True
        # Get target matrix. this will apply matrix manipulation only on the
        # first calculation. After that, it just returns the target. 
        target_matrix = self._initialize_target_matrix()
        
        # Calculate current interaction matrix
        if self.homotypic_interaction:
            current_matrix = self._calculate_matrix(protein.sequence, protein.sequence)
        else:
            current_matrix = self._calculate_matrix(protein.sequence, 
                                                    self.target_sequence)

        # Calculate and return difference
        return self._calculate_matrix_difference(target_matrix, current_matrix)


class MatchSelfIntermap(EpsilonMatrixProperty):
    """
    Calculates self interaction using matrix representation for closer target matching.
    Uses EpsilonMatrixProperty for consistent matrix handling patterns.
    """
    can_set_target_sequence_and_target_value = True  # You cannot set target sequence for self interaction
    def __init__(self, target_sequence: str, weight: float = 1.0, 
                 model: str = 'mpipi', preloaded_model = None, 
                 inverse: bool = False, window_size: int = 15, 
                 constraint_type = ConstraintType.EXACT, allow_matrix_resizing=True,
                 homotypic_interaction: bool = True):
        # For self-interaction, both sequences are the same initially.
        # set them both to be the same
        super().__init__(target_sequence, target_sequence, 0.0, weight, constraint_type,
                         model, preloaded_model, window_size, allow_matrix_resizing,
                         homotypic_interaction)
        self.inverse = inverse

    def get_init_args(self) -> dict:
        """Override to include inverse parameters"""
        base_args = super().get_init_args()
        base_args.update({
            "inverse": self.inverse
        })
        return base_args

    def _manipulate_matrix(self, matrix: np.ndarray) -> np.ndarray:
        """Apply inverse transformation if requested"""
        if self.inverse:
            return matrix * -1
        else:
            return matrix



class MatchIntermap(EpsilonMatrixProperty):
    """
    Optimize a sequence to match the epsilon matrix that is formed 
    between an original sequence and a target interacting sequence.
    Uses EpsilonMatrixProperty for consistent matrix handling patterns.
    """
    def __init__(self, sequence: str, target_sequence: str, weight: float = 1.0,
                 model: str = 'mpipi', preloaded_model = None, 
                 window_size=15, constraint_type = ConstraintType.EXACT, allow_matrix_resizing=True,
                 homotypic_interaction: bool = False):
        # Initialize with target_sequence as interacting_sequence and interacting_sequence as target
        super().__init__(sequence, target_sequence, 0.0, weight, constraint_type,
                         model, preloaded_model, window_size, allow_matrix_resizing,
                         homotypic_interaction)


class ModifyAttractiveValues(EpsilonMatrixProperty):
    """
    Modify the attractive values. 
    
    Setting attraction_multiplier to a value less than 1 and greater than 0
    will reduce the strength of attractive interactions across the matrix. 

    Setting attraction_multiplier to a value greater than 1 will increase the strength
    of attractive interactions across the matrix. 

    Setting attraction_multiplier to a value less than 0 will invert the attractive
    values to become repulsive. 
        Values less than 0 and greater than -1 will reduce the strength and invert to repulsive. 
        Values less than -1 will flip the repulsive values and increase their strength. 

    """
    def __init__(self, sequence:str, target_sequence: str, multiplier: float,
                 weight: float = 1.0, model: str = 'mpipi', 
                 preloaded_model = None, window_size: int = 15,
                 constraint_type = ConstraintType.EXACT, allow_matrix_resizing=True,
                 homotypic_interaction: bool = False):
        super().__init__(sequence, target_sequence, 0.0, weight, constraint_type,
                         model, preloaded_model, window_size, allow_matrix_resizing,
                         homotypic_interaction)
        self.multiplier = multiplier

    def get_init_args(self) -> dict:
        """Override to include enhancement parameters"""
        base_args = super().get_init_args()
        base_args.update({
            'multiplier': self.multiplier
        })
        return base_args

    def _manipulate_matrix(self, matrix: np.ndarray) -> np.ndarray:
        """Enhance attractive interactions following GOOSE patterns"""
        # Only enhance the target matrix
        return MatrixManipulation.multiply_values_of_matrix(
            matrix, self.multiplier, only_negative=True
        )

class ModifyRepulsiveValues(EpsilonMatrixProperty):
    """
    Increase the repulsive value based on the Epsilon intermap. 

    Setting repulsive_multiplier to a value less than 1 and greater than 0
    will reduce the strength of repulsive interactions across the matrix. 

    Setting repulsive_multiplier to a value greater than 1 will increase the strength
    of repulsive interactions across the matrix. 

    Setting repulsive_multiplier to a value less than 0 will invert the repulsive
    values to become attractive. 
        Values less than 0 and greater than -1 will reduce the strength and invert to attractive. 
        Values less than -1 will flip the repulsive values and increase their strength. 
    """
    def __init__(self, sequence:str, target_sequence: str, 
                 multiplier: float, weight: float = 1.0, model: str = 'mpipi', 
                 preloaded_model = None, window_size: int = 15,
                 constraint_type = ConstraintType.EXACT, allow_matrix_resizing=True,
                 homotypic_interaction: bool = False):
        super().__init__(sequence, target_sequence,
                         0.0, weight, constraint_type, model, preloaded_model,
                         window_size, allow_matrix_resizing, homotypic_interaction)
        self.multiplier = multiplier

    def get_init_args(self) -> dict:
        """Override to include enhancement parameters"""
        base_args = super().get_init_args()
        base_args.update({
            'multiplier': self.multiplier
        })
        return base_args

    def _manipulate_matrix(self, matrix: np.ndarray) -> np.ndarray:
        """Enhance attractive interactions following GOOSE patterns"""
        # Only enhance the target matrix
        return MatrixManipulation.multiply_values_of_matrix(
            matrix, self.multiplier, only_positive=True
        )


class ModifyMatrixValues(EpsilonMatrixProperty):
    """
    combined to allow modulation of attractive and repulsive interactions in a single
    property. Tends to work better than trying to manipulate them individually and is 
    more computationally efficient. 
    """
    def __init__(self, sequence:str, target_sequence: str, 
                 repulsive_multiplier: float, attractive_multiplier: float,
                 weight: float = 1.0, model: str = 'mpipi', 
                 preloaded_model = None, window_size: int = 15,
                 constraint_type = ConstraintType.EXACT, allow_matrix_resizing=True,
                 homotypic_interaction: bool = False):
        super().__init__(sequence, target_sequence, 
                         0.0, weight, constraint_type,
                        model, preloaded_model, window_size, allow_matrix_resizing,
                        homotypic_interaction)
        self.repulsive_multiplier = repulsive_multiplier
        self.attractive_multiplier = attractive_multiplier

    def get_init_args(self) -> dict:
        """Override to include enhancement parameters"""
        base_args = super().get_init_args()
        base_args.update({
            'repulsive_multiplier': self.repulsive_multiplier,
            'attractive_multiplier': self.attractive_multiplier,
        })
        return base_args

    def _manipulate_matrix(self, matrix: np.ndarray) -> np.ndarray:
        """manipulate that matrix. """
        return MatrixManipulation.multiply_values_of_matrix(
            MatrixManipulation.multiply_values_of_matrix(
            matrix, self.attractive_multiplier, only_negative=True),
            self.repulsive_multiplier, only_positive=True)
        
        

class MatchArbitraryMatrix(EpsilonMatrixProperty):
    """
    Make a sequence matched to an arbitrary matrix.
    """
    def __init__(self, arbitrary_matrix: np.ndarray, target_sequence: str = None,
                 weight: float = 1.0, model: str = 'mpipi',
                 preloaded_model=None, window_size: int = 15,
                 constraint_type=ConstraintType.EXACT, allow_matrix_resizing=True,
                 homotypic_interaction: bool = False):
        super().__init__(None, target_sequence, 0.0, weight, constraint_type,
                        model, preloaded_model, window_size, allow_matrix_resizing,
                        homotypic_interaction)
        self.arbitrary_matrix = arbitrary_matrix
        # because this is matching something to a premade matrix, we will handle
        # whether this is homotypic dynamically by checking whether target_sequence is None.
        self.homotypic_interaction = target_sequence is None

    def get_init_args(self) -> dict:
        """Override to include enhancement parameters"""
        base_args = super().get_init_args()
        base_args.update({
            'target_sequence': self.target_sequence,
            'arbitrary_matrix': self.arbitrary_matrix,
        })
        return base_args
    
    def _initialize_target_matrix(self) -> np.ndarray:
        """
        Override because we are inputting the target matrix. 

        Returns
        -------
        np.ndarray
            Initialized target matrix
        """
        if self.target_matrix is None:
            self.target_matrix = self.arbitrary_matrix
            # resize so target_matrix matches matrix shape produced from _calculate_matrix for generated sequence
            self.target_matrix = self._match_matrix_size()
        return self.target_matrix

class FDSurfaceEpsilonProperty(EpsilonProperty):
    """
    Base class for calculations involving FoldedDomain FINCHES objects where we
    are looking at the interactions between a folded domain and an IDR.
    """
    can_set_target_sequence_and_target_value = True  # You can set both target sequence and target value
    def __init__(self, sequence: str,
                 target_value: float = 0, folded_domain: 'finches.folded_domain.FoldedDomain' = None,
                 path_to_folded_domain: str = None, weight: float = 1.0, 
                 model: str = 'mpipi', preloaded_model=None, constraint_type=ConstraintType.EXACT,
                 probe_radius: float = 1.4, surface_thresh: float = 0.10, sasa_mode: str = 'v1', 
                 fd_start: int = None, fd_end: int = None):
        super().__init__(target_value, sequence, weight, constraint_type, model, preloaded_model)
        self.sequence = sequence
        self.folded_domain = folded_domain
        self.path_to_folded_domain = path_to_folded_domain
        self.probe_radius = probe_radius
        self.surface_thresh = surface_thresh
        self.sasa_mode = sasa_mode
        self.fd_start = fd_start
        self.fd_end = fd_end
        self._imc_object = None  # Cache for expensive IMC object loading

    def get_init_args(self) -> dict:
        """Override to include folded_domain parameter"""
        base_args = super().get_init_args()
        base_args['sequence'] = self.sequence
        base_args['folded_domain'] = self.folded_domain
        base_args['path_to_folded_domain'] = self.path_to_folded_domain
        base_args['probe_radius'] = self.probe_radius
        base_args['surface_thresh'] = self.surface_thresh
        base_args['sasa_mode'] = self.sasa_mode
        base_args['fd_start'] = self.fd_start
        base_args['fd_end'] = self.fd_end
        return base_args

    def load_imc_object(self, model: str = None):
        """
        Load the IMC object for surface calculations.
        
        Parameters
        ----------
        model : str, optional
            Model type. Uses instance model if None.
            
        Returns
        -------
        IMC object
            The loaded IMC object for surface epsilon calculations
        """
        if self._imc_object is not None:
            return self._imc_object
            
        # Load the full model first
        loaded_model = self.load_model(model)
        self._imc_object = loaded_model.IMC_object
        return self._imc_object

    def _load_folded_domain(self):
        if self.folded_domain is not None:
            return self.folded_domain
        else:
            self.folded_domain = self._load_folded_domain_from_path()
            return self.folded_domain
    
    def _load_folded_domain_from_path(self):
        if self.path_to_folded_domain is not None:
            # Load the folded domain from the specified path
            self.folded_domain = finches.utils.folded_domain_utils.FoldedDomain(
                                                                    pdbfilename = self.path_to_folded_domain,
                                                                    start=self.fd_start,
                                                                    end=self.fd_end,
                                                                    probe_radius=self.probe_radius,
                                                                    surface_thresh=self.surface_thresh,
                                                                    sasa_mode=self.sasa_mode)
        return self.folded_domain
    
    @abstractmethod
    def calculate_raw_value(self, protein: 'sparrow.Protein') -> float:
        """
        Calculate the raw property value for a given protein.
        Must be implemented by subclasses.
        
        Parameters
        ----------
        protein : sparrow.Protein
            The protein instance
            
        Returns
        -------
        float
            The calculated raw property value
        """
        pass    

class FDMeanSurfaceEpsilon(FDSurfaceEpsilonProperty):
    """
    Class for calculating mean surface epsilon values for folded domain proteins.
    """
    def __init__(self, sequence: str=None, target_value: float = None, folded_domain: 'finches.folded_domain.FoldedDomain' = None,
                 path_to_folded_domain: str = None, weight: float = 1.0,
                 model: str = 'mpipi', preloaded_model=None, constraint_type=ConstraintType.EXACT,
                 probe_radius: float = 1.4, surface_thresh: float = 0.10, sasa_mode: str = 'v1',
                 fd_start: int = None, fd_end: int = None):
        super().__init__(sequence, target_value, folded_domain, path_to_folded_domain, weight,
                         model, preloaded_model, constraint_type, probe_radius, surface_thresh,
                         sasa_mode, fd_start, fd_end)

    def get_init_args(self):
        return {
            "sequence": self.sequence,
            "target_value": self.target_value,
            "folded_domain": self.folded_domain,
            "path_to_folded_domain": self.path_to_folded_domain,
            "weight": self.weight,
            "model": self.model,
            "preloaded_model": self.preloaded_model,
            "constraint_type": self.constraint_type,
            "probe_radius": self.probe_radius,
            "surface_thresh": self.surface_thresh,
            "sasa_mode": self.sasa_mode,
            "fd_start": self.fd_start,
            "fd_end": self.fd_end
        }

    def _initialize_target_epsilon(self):
        if self.sequence is None and self.target_value is None:
            raise ValueError("Either sequence or target_value must be provided.")
        else:
            if self.sequence is not None:
                # load
                fd = self._load_folded_domain()
                # calculate
                self.target_value = fd.calculate_mean_surface_epsilon(self.sequence, self.load_imc_object())
        return self.target_value

    def calculate_raw_value(self, protein):
        fd = self._load_folded_domain()
        return fd.calculate_mean_surface_epsilon(protein.sequence, self.load_imc_object())


class FDSurfaceEpsilon(FDSurfaceEpsilonProperty):
    """
    Class for calculating surface epsilon values for folded domain proteins.
    """
    def __init__(self, sequence: str, attractive_multiplier: float = 1, repulsive_multiplier: float = 1,
                 target_value: float = 0, folded_domain: 'finches.folded_domain.FoldedDomain' = None,
                 path_to_folded_domain: str = None, weight: float = 1.0,
                 model: str = 'mpipi', preloaded_model=None, constraint_type=ConstraintType.EXACT,
                 probe_radius: float = 1.4, surface_thresh: float = 0.10, sasa_mode: str = 'v1',
                 fd_start: int = None, fd_end: int = None):
        super().__init__(sequence, target_value, folded_domain, path_to_folded_domain, weight,
                         model, preloaded_model, constraint_type, probe_radius, surface_thresh,
                         sasa_mode, fd_start, fd_end)
        # store surface epsilon of target to save compute
        self._target_surface_epsilon = None
        self.attractive_multiplier = attractive_multiplier
        self.repulsive_multiplier = repulsive_multiplier

    def get_init_args(self):
        return {
            "sequence": self.sequence,
            "attractive_multiplier": self.attractive_multiplier,
            "repulsive_multiplier": self.repulsive_multiplier,
            "target_value": self.target_value,
            "folded_domain": self.folded_domain,
            "path_to_folded_domain": self.path_to_folded_domain,
            "weight": self.weight,
            "model": self.model,
            "preloaded_model": self.preloaded_model,
            "constraint_type": self.constraint_type,
            "probe_radius": self.probe_radius,
            "surface_thresh": self.surface_thresh,
            "sasa_mode": self.sasa_mode,
            "fd_start": self.fd_start,
            "fd_end": self.fd_end
        }

    def calculate_surface_epsilon(self, sequence: str) -> list:
        """
        Calculate the surface epsilon value for the given sequence.

        Parameters
        ----------
        sequence : str
            The protein sequence to calculate the surface epsilon for.

        Returns
        -------
        list
            List of surface epsilon values.
        """
        imc_object = self.load_imc_object()
        folded_domain = self._load_folded_domain()
        surf_eps_dict = folded_domain.calculate_surface_epsilon(sequence, imc_object)
        return [float(i[2]) for i in list(surf_eps_dict.values())]

    def _initialize_target_epsilon(self):
        if self._target_surface_epsilon is None:
            self._target_surface_epsilon = self.calculate_surface_epsilon(self.sequence)
            self._target_surface_epsilon = VectorManipulation.multiply_attractive_force(self._target_surface_epsilon, self.attractive_multiplier)
            self._target_surface_epsilon = VectorManipulation.multiply_repulsive_force(self._target_surface_epsilon, self.repulsive_multiplier)
        return self._target_surface_epsilon
    
    def calculate_raw_value(self, protein):
        """
        Calculate the raw surface epsilon value for the given protein.

        Parameters
        ----------
        protein : str
            The protein sequence to calculate the surface epsilon for.

        Returns
        -------
        float
            The calculated raw surface epsilon value.
        """
        self._initialize_target_epsilon()
        # calculate value for protein.sequence
        current_surface_epsilon = self.calculate_surface_epsilon(protein.sequence)
        return sum(abs(np.array(current_surface_epsilon)-np.array(self._target_surface_epsilon)))
    
class FDSurfacePatchInteractions(FDSurfaceEpsilonProperty):
    """
    Class for calculating surface patch interactions for folded domain proteins.
    """
    def __init__(self, sequence: str, attractive_multiplier: float = 1, repulsive_multiplier: float = 1,
                 target_value: float = 0, folded_domain: 'finches.folded_domain.FoldedDomain' = None,
                 path_to_folded_domain: str = None, weight: float = 1.0,
                 model: str = 'mpipi', preloaded_model=None, constraint_type=ConstraintType.EXACT,
                 probe_radius: float = 1.4, surface_thresh: float = 0.10, sasa_mode: str = 'v1',
                 fd_start: int = None, fd_end: int = None, idr_tile_size=15, patch_radius=12):
        super().__init__(sequence, target_value, folded_domain, path_to_folded_domain, weight,
                         model, preloaded_model, constraint_type, probe_radius, surface_thresh,
                         sasa_mode, fd_start, fd_end)
        # store surface epsilon of target to save compute
        self._target_surface_epsilon = None
        self.attractive_multiplier = attractive_multiplier
        self.repulsive_multiplier = repulsive_multiplier
        self.idr_tile_size = idr_tile_size
        self.patch_radius = patch_radius

    def get_init_args(self):
        return {
            "sequence": self.sequence,
            "attractive_multiplier": self.attractive_multiplier,
            "repulsive_multiplier": self.repulsive_multiplier,
            "target_value": self.target_value,
            "folded_domain": self.folded_domain,
            "path_to_folded_domain": self.path_to_folded_domain,
            "weight": self.weight,
            "model": self.model,
            "preloaded_model": self.preloaded_model,
            "constraint_type": self.constraint_type,
            "probe_radius": self.probe_radius,
            "surface_thresh": self.surface_thresh,
            "sasa_mode": self.sasa_mode,
            "fd_start": self.fd_start,
            "fd_end": self.fd_end,
            "idr_tile_size": self.idr_tile_size,
            "patch_radius": self.patch_radius
        }

    def calculate_surface_epsilon(self, sequence: str) -> list:
        """
        Calculate the surface epsilon value for the given sequence.

        Parameters
        ----------
        sequence : str
            The protein sequence to calculate the surface epsilon for.

        Returns
        -------
        list
            List of surface epsilon values.
        """
        imc_object = self.load_imc_object()
        folded_domain = self._load_folded_domain()
        surf_eps=folded_domain.calculate_idr_surface_patch_interactions(sequence, imc_object,
                                                                               idr_tile_size=self.idr_tile_size,
                                                                              patch_radius=self.patch_radius)[1]
        # return the matrix of values only after changing nan to 0 
        return np.nan_to_num(surf_eps)


    def _initialize_target_epsilon(self) -> np.ndarray:
        """Initialize and cache the target surface epsilon matrix."""
        if self._target_surface_epsilon is None:
            # Calculate target matrix
            target_matrix = self.calculate_surface_epsilon(self.sequence)
            
            # Apply multipliers using matrix manipulation
            if self.attractive_multiplier != 1.0:
                target_matrix = MatrixManipulation.multiply_values_of_matrix(
                    target_matrix, self.attractive_multiplier, only_negative=True
                )
            
            if self.repulsive_multiplier != 1.0:
                target_matrix = MatrixManipulation.multiply_values_of_matrix(
                    target_matrix, self.repulsive_multiplier, only_positive=True
                )
            
            # Resize if needed based on sequence length hint
            generated_sequence_length = self._get_sequence_length_hint()
            expected_idr_tiles = generated_sequence_length - self.idr_tile_size + 1
            
            if target_matrix.shape[1] != expected_idr_tiles:
                new_shape = (target_matrix.shape[0], expected_idr_tiles)
                target_matrix = MatrixManipulation.scale_matrix_to_size(target_matrix, new_shape)
            
            self._target_surface_epsilon = target_matrix
            
        return self._target_surface_epsilon
    
    def calculate_raw_value(self, protein):
        """
        Calculate the raw surface epsilon value for the given protein.

        Parameters
        ----------
        protein : str
            The protein sequence to calculate the surface epsilon for.

        Returns
        -------
        float
            The calculated raw surface epsilon value.
        """
            # Initialize target matrix
        target_matrix = self._initialize_target_epsilon()
        
        # Calculate current matrix
        current_matrix = self.calculate_surface_epsilon(protein.sequence)
        
        # Ensure both matrices are numpy arrays
        target_matrix = np.asarray(target_matrix)
        current_matrix = np.asarray(current_matrix)
        
        # Calculate absolute difference
        difference = np.sum(np.abs(target_matrix - current_matrix))
        
        # Normalize by sequence length following GOOSE patterns
        sequence_length = min(target_matrix.shape[1], current_matrix.shape[1])
        if sequence_length > 0:
            difference = difference / sequence_length
        return float(difference)