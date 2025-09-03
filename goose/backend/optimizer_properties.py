from abc import ABC, abstractmethod
from typing import Any, Tuple, Dict, Callable
from enum import Enum
import inspect
import numpy as np

import sparrow
import metapredict as meta
import finches
from finches.utils.folded_domain_utils import FoldedDomain

from goose.backend.optimizer_tools import MatrixManipulation

class ConstraintType(Enum):
    """Enumeration for different constraint types in property optimization."""
    EXACT = "exact"      # Minimize absolute difference from target
    MINIMUM = "minimum"  # Penalize only when below target
    MAXIMUM = "maximum"  # Penalize only when above target

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
    multi_target = True  # Class-level attribute - override to limit ability to have multiple targest. 

    # currently testing this to get rid of the problem whereby properties that have larger
    # values are overweighted due to their having larger absolute error... Hope it worsk. 
    # Class-level normalization attributes - override in subclasses
    IS_NATURALLY_NORMALIZED = False  # True if property naturally ranges 0-1

    def __init__(self, target_value: float, weight: float = 1.0, constraint_type: ConstraintType = ConstraintType.EXACT):
        assert isinstance(target_value, (int,float)), f"target_value must be numerical. Received {type(target_value)}"
        assert isinstance(weight, (int,float)), f"Weight must be numerical. Received {type(weight)}"

        # handle constraint type if user inputs a string instead of Enum
        self.constraint_type = _normalize_constraint_type(constraint_type)
        self._target_value = target_value
        self._weight = weight
        self._current_raw_value=None
        self._best_value=None
        # this should update at initialization to allow for relative error values 
        # for non 0 to 1 properties. 
        self.initial_value = None  # Will be set during optimization setup
        

    ## note: the is_naturally_normalized, set_initial_value and get_scaling_range, 
    # are all things to be checked carefully if there are problems. I'm still not confident
    # that relative error value normalization is actually the best way to do this. 
    def is_naturally_normalized(self) -> bool:
        """Check if this property is naturally 0-1 scaled."""
        return self.IS_NATURALLY_NORMALIZED
    
    def set_initial_value(self, value: float):
        """Set the initial value for relative scaling"""
        self.initial_value = value
    
    def get_scaling_range(self) -> float:
        """
        Get the scaling range based on initial vs target values.
        Used for relative error normalization.
        """
        if self.is_naturally_normalized():
            return 1.0  # Already normalized, max possible error is 1.0
        
        if self.initial_value is not None:
            # Use actual initial-to-target range for scaling
            return max(abs(self.target_value - self.initial_value), 0.01)
        
        # Fallback if no initial value (shouldn't happen in normal usage)
        return 1.0    
    
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
        return {
            "target_value": self.target_value,
            "weight": self.weight,
            "constraint_type": self.constraint_type.value  # Store as string for JSON compatibility
        }

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

    def calculate_error(self, current_value: float) -> float:
        """
        Calculate the error between the current value and the target value based on the constraint type.

        Parameters:
        current_value (float): The current value of the property.

        Returns:
        float: The calculated error (always positive).
        """
        if self.constraint_type == ConstraintType.EXACT:
            return abs(current_value - self.target_value)
        elif self.constraint_type == ConstraintType.MINIMUM:
            if current_value >= self.target_value:
                return 0.0  # No penalty if above minimum
            else:
                return self.target_value - current_value  # Penalty for being below minimum
        elif self.constraint_type == ConstraintType.MAXIMUM:
            if current_value <= self.target_value:
                return 0.0  # No penalty if below maximum
            else:
                return current_value - self.target_value  # Penalty for being above maximum
        else:
            raise ValueError(f"Unknown constraint type: {self.constraint_type}")

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
        
        raw_value = self.calculate_raw_value(protein)
        self._current_raw_value = raw_value
        if self._best_value is None:
            self._best_value = raw_value
        return self.calculate_error(raw_value)



class ComputeIWD(ProteinProperty):
    """
    Compute the Inversed Weighted Distance (IWD) property for the target residues in the sequence.
    """
    IS_NATURALLY_NORMALIZED = False  
    
    def __init__(self, residues: Tuple[str, ...], target_value: float, 
                 weight: float = 1.0, constraint_type: ConstraintType = ConstraintType.EXACT):
        super().__init__(target_value, weight, constraint_type)
        self.residues = residues

    # note: This only needs to be overridden in subclasses that have additional parameters
    def get_init_args(self) -> dict:
        """Override to include residues parameter"""
        return {
            "residues": self.residues,
            "target_value": self.target_value,
            "weight": self.weight,
            "constraint_type": self.constraint_type.value  # Store as string for JSON compatibility
        }    

    def calculate_raw_value(self, protein: 'sparrow.Protein') -> float:
        return protein.compute_iwd(list(self.residues))


class Hydrophobicity(ProteinProperty):
    """
    Calculate the hydrophobicity property.
    """
    IS_NATURALLY_NORMALIZED = False  # Hydrophobicity scale is typically 0-6.6
    
    def __init__(self, target_value: float, weight: float = 1.0, constraint_type: ConstraintType = ConstraintType.EXACT):
        super().__init__(target_value, weight, constraint_type)

    def calculate_raw_value(self, protein: 'sparrow.Protein') -> float:
        return protein.hydrophobicity

class FCR(ProteinProperty):
    """
    Calculate the FCR property.
    """
    IS_NATURALLY_NORMALIZED = True  # FCR ranges 0-1
    
    def __init__(self, target_value: float, weight: float = 1.0, constraint_type: ConstraintType = ConstraintType.EXACT):
        super().__init__(target_value, weight, constraint_type)

    def calculate_raw_value(self, protein: 'sparrow.Protein') -> float:
        return protein.FCR

class NCPR(ProteinProperty):
    """
    Calculate the NCPR property.
    """
    IS_NATURALLY_NORMALIZED = True  # NCPR ranges -1 to 1, but effectively normalized
    
    def __init__(self, target_value: float, weight: float = 1.0, constraint_type: ConstraintType = ConstraintType.EXACT):
        super().__init__(target_value, weight, constraint_type)

    def calculate_raw_value(self, protein: 'sparrow.Protein') -> float:
        return protein.NCPR


class Kappa(ProteinProperty):
    """
    Calculate the kappa property.
    """
    IS_NATURALLY_NORMALIZED = True  # Kappa ranges 0-1
    
    def __init__(self, target_value: float, weight: float = 1.0, constraint_type: ConstraintType = ConstraintType.EXACT):
        super().__init__(target_value, weight, constraint_type)

    def calculate_raw_value(self, protein: 'sparrow.Protein') -> float:
        return protein.kappa

class RadiusOfGyration(ProteinProperty):
    """
    Calculate the Radius of gyration
    """
    IS_NATURALLY_NORMALIZED = False  # Rg values depend on sequence length
    
    def __init__(self, target_value: float, weight: float = 1.0, 
                 constraint_type: ConstraintType = ConstraintType.EXACT):
        super().__init__(target_value, weight, constraint_type)

    def calculate_raw_value(self, protein: 'sparrow.Protein') -> float:
        return protein.predictor.radius_of_gyration()

class EndToEndDistance(ProteinProperty):
    """
    Calculate the Radius of gyration
    """
    IS_NATURALLY_NORMALIZED = False  # End-to-end distance depends on sequence length
    
    def __init__(self, target_value: float, weight: float = 1.0, 
                 constraint_type: ConstraintType = ConstraintType.EXACT):
        super().__init__(target_value, weight, constraint_type)

    def calculate_raw_value(self, protein: 'sparrow.Protein') -> float:
        return protein.predictor.end_to_end_distance()


class AminoAcidFractions(ProteinProperty):
    """
    Compute the difference between the target amino acid fractions and the current amino acid fractions.
    """
    IS_NATURALLY_NORMALIZED = False  # Fractions range 0-1
    
    def __init__(self, target_fractions: Dict[str, float], weight: float = 1.0, 
                 constraint_type: ConstraintType = ConstraintType.EXACT):
        # For fractions, target value is not meaningful at class level since we have multiple amino acids
        super().__init__(0.0, weight, constraint_type)
        self.target_fractions = target_fractions

    def get_init_args(self) -> dict:
        """Override to include target_fractions parameter"""
        return {
            "target_fractions": self.target_fractions,
            "weight": self.weight,
            "constraint_type": self.constraint_type.value
        }

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
        Following GOOSE patterns for complex properties.
        """
        # Calculate the final error value
        final_error = self.calculate_raw_value(protein)
        
        # Update tracking attributes following base class patterns
        self._current_raw_value = final_error
        if self._best_value is None:
            self._best_value = final_error
        
        # Return final error (no additional constraint processing needed)
        return final_error


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
    IS_NATURALLY_NORMALIZED = False  # SCD values are arbitrary scale
    
    def __init__(self, target_value: float, weight: float = 1.0, 
                 constraint_type: ConstraintType = ConstraintType.EXACT):
        super().__init__(target_value, weight, constraint_type)

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
    IS_NATURALLY_NORMALIZED = False  # SHD values are arbitrary scale
    
    def __init__(self, target_value: float, weight: float = 1.0, 
                 constraint_type: ConstraintType = ConstraintType.EXACT):
        super().__init__(target_value, weight, constraint_type)

    def calculate_raw_value(self, protein: 'sparrow.Protein') -> float:
        return protein.SHD

class Complexity(ProteinProperty):
    """
    Calculates the Wootton-Federhen complexity of a sequence (also called
    seg complexity, as this the theory used in the classic SEG algorithm.
    """
    IS_NATURALLY_NORMALIZED = True  # Complexity values are arbitrary scale
    
    def __init__(self, target_value: float, weight: float = 1.0, 
                 constraint_type: ConstraintType = ConstraintType.EXACT):
        super().__init__(target_value, weight, constraint_type)

    def calculate_raw_value(self, protein: 'sparrow.Protein') -> float:
        return protein.complexity


class FractionDisorder(ProteinProperty):
    """
    Calculates the amount of disorder in
    the current sequence. 
    """
    IS_NATURALLY_NORMALIZED = True  # Fraction ranges 0-1
    
    def __init__(self, target_value: float, weight: float = 1.0, disorder_cutoff = 0.5,
                 constraint_type: ConstraintType = ConstraintType.EXACT):
        super().__init__(target_value, weight, constraint_type)
        self.disorder_cutoff = disorder_cutoff

    def get_init_args(self) -> dict:
        """Override to include disorder_cutoff parameter"""
        return {
            "target_value": self.target_value,
            "weight": self.weight,
            "disorder_cutoff": self.disorder_cutoff,
            "constraint_type": self.constraint_type.value
        }

    def calculate_raw_value(self, protein: 'sparrow.Protein') -> float:
        disorder=meta.predict_disorder(protein.sequence)
        return (disorder>self.disorder_cutoff).sum()/len(disorder)

        
    

class MatchSequenceDisorder(ProteinProperty):
    """
    Matches the disorder of a specified sequence. By default, the 
    specified sequence sets the MINIMUM disorder. 
    You can also try to get the exact disorder.

    """
    IS_NATURALLY_NORMALIZED = True  # Normalized by sequence length, ranges 0-1
    
    def __init__(self, target_sequence: str, exact_match: bool = False, target_value : float = 0,
                 weight: float = 1.0, constraint_type: ConstraintType = ConstraintType.EXACT):
        super().__init__(target_value, weight, constraint_type)
        self.target_sequence = target_sequence
        self.exact_match = exact_match
        self.target_disorder = None

    def get_init_args(self) -> dict:
        """Override to include target_sequence and exact_match parameters"""
        return {
            "target_sequence": self.target_sequence,
            "exact_match": self.exact_match,
            "target_value": self.target_value,
            "weight": self.weight,
            "constraint_type": self.constraint_type.value
        }

    def set_initial_disorder(self):
        if self.target_disorder is None:
            self.target_disorder = meta.predict_disorder(self.target_sequence)
        return self.target_disorder

    def calculate_raw_value(self, protein: 'sparrow.Protein') -> float:
        target_disorder=self.set_initial_disorder()
        disorder=meta.predict_disorder(protein.sequence)
        if self.exact_match==True:
            err = np.sum(np.abs(disorder-target_disorder))
        else:
            err = np.sum(disorder<target_disorder)
        # normalize to length. 
        return err/len(target_disorder)


class MatchingResidues(ProteinProperty):
    '''
    Determines the number of residues that match a target sequence in the current sequence.
    '''
    IS_NATURALLY_NORMALIZED = False  # Count of residues, depends on sequence length
    
    def __init__(self, target_sequence: str, target_value: float, weight: float = 1.0,
                 constraint_type: ConstraintType = ConstraintType.EXACT):
        super().__init__(target_value, weight, constraint_type)
        self.target_sequence = target_sequence

    def get_init_args(self):
        """Override to include target_sequence parameter"""
        return {
            "target_sequence": self.target_sequence,
            "target_value": self.target_value,
            "weight": self.weight,
            "constraint_type": self.constraint_type.value
        }

    def calculate_raw_value(self, protein: 'sparrow.Protein') -> float:
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
    IS_NATURALLY_NORMALIZED = False  # Epsilon values are arbitrary scale
    
    def __init__(self, target_value: float, weight: float = 1.0, 
                 constraint_type = ConstraintType.EXACT, model: str = 'mpipi', 
                 preloaded_model = None):
        super().__init__(target_value, weight, constraint_type)
        self.model = model.lower()
        self._loaded_model = preloaded_model
        self._loaded_imc_object = None
        
        # Handle preloaded model type detection following GOOSE's validation patterns
        if preloaded_model is not None:
            self._detect_preloaded_model_type(preloaded_model)

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

    def get_init_args(self) -> dict:
        """
        Get initialization arguments including model parameter.
        Override in subclasses with additional parameters.
        """
        base_args = super().get_init_args()
        base_args['model'] = self.model
        return base_args


class MeanSelfEpsilon(EpsilonProperty):
    """
    Calculate the self interaction epsilon value of a sequence.
    Note: this simply uses the mean epsilon value. 
    """
    IS_NATURALLY_NORMALIZED = False  # Epsilon values are arbitrary scale
    
    def __init__(self, target_value: float, weight: float = 1.0,
                 model: str = 'mpipi', preloaded_model = None, 
                 constraint_type = ConstraintType.EXACT):
        super().__init__(target_value, weight, constraint_type, model, preloaded_model)

    def calculate_raw_value(self, protein: 'sparrow.Protein') -> float:
        loaded_model = self.load_model()
        return loaded_model.epsilon(protein.sequence, protein.sequence)


class MatchSelfIntermap(EpsilonProperty):
    """
    Calculates self interaction using matrix representation for closer target matching.
    """
    IS_NATURALLY_NORMALIZED = False  # Matrix differences are arbitrary scale
    
    def __init__(self, target_sequence: str, weight: float = 1.0, 
                 model: str = 'mpipi', preloaded_model = None, 
                 inverse: bool = False, window_size=15, constraint_type = ConstraintType.EXACT):
        # Always target 0.0 since we minimize matrix differences
        super().__init__(0.0, weight, constraint_type, model, preloaded_model)
        self.target_sequence = target_sequence
        self.inverse = inverse
        self.window_size = window_size
        self._original_epsilon_matrix = None

    def get_init_args(self) -> dict:
        """Override to include target_sequence and inverse parameters"""
        base_args = super().get_init_args()
        base_args.update({
            "target_sequence": self.target_sequence,
            "inverse": self.inverse,
            "window_size": self.window_size
        })
        return base_args

    def _calculate_matrix(self, sequence: str, window_size: int) -> np.ndarray:
        """Calculate epsilon matrix for a sequence"""
        loaded_model = self.load_model()
        return loaded_model.intermolecular_idr_matrix(
            sequence, sequence, window_size=window_size,
            disorder_1=False, disorder_2=False
        )[0][0]

    def _get_original_epsilon_matrix(self) -> np.ndarray:
        """Get cached original epsilon matrix"""
        if self._original_epsilon_matrix is None:
            matrix = self._calculate_matrix(self.target_sequence, self.window_size)
            if self.inverse:
                self._original_epsilon_matrix = matrix * -1
            else:
                self._original_epsilon_matrix = matrix    
        return self._original_epsilon_matrix
    
    def calculate_raw_value(self, protein: 'sparrow.Protein') -> float:
        original_matrix = self._get_original_epsilon_matrix()
        current_matrix = self._calculate_matrix(protein.sequence, self.window_size)
        # Normalize by sequence length 
        return np.sum(np.abs(original_matrix - current_matrix)) / len(protein.sequence)


class MeanEpsilonWithTarget(EpsilonProperty):
    """
    Make a sequence that interacts with a specific sequence with a specific mean epsilon value. 
    """
    IS_NATURALLY_NORMALIZED = False  # Epsilon values are arbitrary scale
    
    def __init__(self, target_value: float, target_sequence: str, weight: float = 1.0,
                 model: str = 'mpipi', preloaded_model = None, 
                 constraint_type = ConstraintType.EXACT):
        super().__init__(target_value, weight, constraint_type, model, preloaded_model)
        self.target_sequence = target_sequence

    def get_init_args(self) -> dict:
        """Override to include target_sequence parameter"""
        base_args = super().get_init_args()
        base_args['target_sequence'] = self.target_sequence
        return base_args

    def calculate_raw_value(self, protein: 'sparrow.Protein') -> float:
        loaded_model = self.load_model()
        return loaded_model.epsilon(protein.sequence, self.target_sequence)


class MatchIntermap(EpsilonProperty):
    """
    Optimize a sequence to match the epsilon matrix that is formed 
    between an original sequence and a target interacting sequence.
    Note: the original sequence must be the same length as the sequence being optimized. 
    """
    IS_NATURALLY_NORMALIZED = False  # Matrix differences are arbitrary scale
    
    def __init__(self, target_sequence: str, interacting_sequence: str, weight: float = 1.0,
                 model: str = 'mpipi', preloaded_model = None, 
                 window_size=15, constraint_type = ConstraintType.EXACT):
        # target_value always 0.0
        super().__init__(0.0, weight, constraint_type, model, preloaded_model)
        self.target_sequence = target_sequence
        self.interacting_sequence = interacting_sequence
        self.window_size = window_size
        self._target_matrix = None

    def get_init_args(self) -> dict:
        """Override to include sequence parameters"""
        base_args = super().get_init_args()
        base_args.update({
            "target_sequence": self.target_sequence,
            "interacting_sequence": self.interacting_sequence,
            "window_size": self.window_size
        })
        return base_args

    def _calculate_matrix(self, sequence1: str, sequence2: str, window_size: int) -> np.ndarray:
        """Calculate epsilon matrix for a sequence"""
        loaded_model = self.load_model()
        return loaded_model.intermolecular_idr_matrix(sequence1, sequence2,
            disorder_1=False, disorder_2=False, window_size=window_size )[0][0]

    def _get_target_matrix(self) -> float:
        """Get cached target matrix"""
        if self._target_matrix is None:
            self._target_matrix = self._calculate_matrix(self.target_sequence, self.interacting_sequence, self.window_size)
        return self._target_matrix
    
    def calculate_raw_value(self, protein):
        target_matrix = self._get_target_matrix()
        current_matrix = self._calculate_matrix(protein.sequence, self.interacting_sequence, self.window_size)
        # make sure the matrices are the same size otherwise scale.
        if current_matrix.shape != target_matrix.shape:
            # scale target to be same dims as current
            target_matrix = MatrixManipulation.scale_matrix_to_size(target_matrix, current_matrix.shape)
            # store scaled target so we don't need to recalculate
            self._target_matrix = target_matrix
        return np.sum(np.abs(target_matrix - current_matrix)) / len(protein.sequence)



class ChemicalFingerprint(EpsilonProperty):
    """
    Uses the chemical foot print from the FINCHES manuscript to generate a sequence with
    a similar chemical fingerprint to the target sequence.

    The chemical fingerprint is calculated by taking the difference between the target
    and the current sequence and summing the differences in the epsilon matrix.
    """
    IS_NATURALLY_NORMALIZED = False  # Fingerprint differences are arbitrary scale
    
    def __init__(self, target_sequence: str, target_value: float = 0.0, 
                 weight: float = 1.0, model: str = 'mpipi', preloaded_model = None, 
                 window_size:int = 15, constraint_type = ConstraintType.EXACT):
        super().__init__(target_value, weight, constraint_type, model, preloaded_model)
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
    the matrix to (for example) increase hte interaction strength. 

    """
    IS_NATURALLY_NORMALIZED = False  # Matrix manipulations are arbitrary scale
    
    def __init__(self, interacting_sequence:str, target_interacting_sequence: str, 
                 target_value: float = 0.0, 
                 weight: float = 1.0, constraint_type = ConstraintType.EXACT,
                 model: str = 'mpipi', preloaded_model = None,
                 window_size: int = 15):
        super().__init__(target_value, weight, constraint_type, model, preloaded_model)
        self.interacting_sequence = interacting_sequence
        self.target_interacting_sequence = target_interacting_sequence
        self.window_size = window_size
        self.target_matrix = None  # Cache for expensive matrix calculations
        
    def get_init_args(self) -> dict:
        """Override to include matrix-specific parameters"""
        base_args = super().get_init_args()
        base_args.update({
            "interacting_sequence": self.interacting_sequence,
            "target_interacting_sequence": self.target_interacting_sequence,
            "window_size": self.window_size
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
                self.interacting_sequence, self.target_interacting_sequence))
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
        # if set to None, we are doing a self interaction. 
        if self.interacting_sequence is None:
            self.interacting_sequence = protein.sequence
        
        # Get target matrix. this will apply matrix manipulation only on the
        # first calculation. After that, it just returns the target. 
        target_matrix = self._initialize_target_matrix()
        
        # Calculate current interaction matrix
        current_matrix = self._calculate_matrix(protein.sequence, 
                                                self.target_interacting_sequence)
        
        # Calculate and return difference
        return self._calculate_matrix_difference(target_matrix, current_matrix)


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
    IS_NATURALLY_NORMALIZED = False  # Modified matrix values are arbitrary scale
    
    def __init__(self, interacting_sequence:str, target_interacting_sequence: str, 
                 multiplier: float,
                 weight: float = 1.0, model: str = 'mpipi', 
                 preloaded_model = None, window_size: int = 15,
                 constraint_type = ConstraintType.EXACT):
        super().__init__(interacting_sequence,target_interacting_sequence, 
                         0.0, weight, constraint_type,
                        model, preloaded_model, window_size)
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
    IS_NATURALLY_NORMALIZED = False  # Modified matrix values are arbitrary scale
    
    def __init__(self, interacting_sequence:str, target_interacting_sequence: str, 
                 multiplier: float,
                 weight: float = 1.0, model: str = 'mpipi', 
                 preloaded_model = None, window_size: int = 15,
                 constraint_type = ConstraintType.EXACT):
        super().__init__(interacting_sequence,target_interacting_sequence, 
                         0.0, weight, constraint_type,
                        model, preloaded_model, window_size)
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
    IS_NATURALLY_NORMALIZED = False  # Modified matrix values are arbitrary scale
    
    def __init__(self, interacting_sequence:str, target_interacting_sequence: str, 
                 repulsive_multiplier: float, attractive_multiplier,
                 weight: float = 1.0, model: str = 'mpipi', 
                 preloaded_model = None, window_size: int = 15,
                 constraint_type = ConstraintType.EXACT):
        super().__init__(interacting_sequence,target_interacting_sequence, 
                         0.0, weight, constraint_type,
                        model, preloaded_model, window_size)
        self.repulsive_multiplier = repulsive_multiplier
        self.attractive_multiplier = attractive_multiplier

    def get_init_args(self) -> dict:
        """Override to include enhancement parameters"""
        base_args = super().get_init_args()
        base_args.update({
            'repulsive_multiplier': self.multiplier,
            'attractive_multiplier': self. attractive_multiplier,
        })
        return base_args

    def _manipulate_matrix(self, matrix: np.ndarray) -> np.ndarray:
        """manipulate that matrix. """
        return MatrixManipulation.multiply_values_of_matrix(
            MatrixManipulation.multiply_values_of_matrix(
            matrix, self.attractive_multiplier, only_negative=True),
            self.repulsive_multiplier, only_positive=True)
        
        


'''
=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-

NOTE TO SELF
NOTE TO SELF
NOTE TO SELF
NOTE TO SELF

Still need to update all code below!!

NOTE TO SELF
NOTE TO SELF
NOTE TO SELF
NOTE TO SELF

=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
'''

class FDSurfaceInteractionByRepulsiveAndAttractiveValues(ProteinProperty):
    """
    Try to get a specific surface repulsion and attraction

    Example usage:
    # imports
    import finches
    from finches.utils.folded_domain_utils import FoldedDomain

    import numpy as np

    import goose
    from sparrow.protein import Protein as pr

    # preload the FD to speed things up.
    path_to_pdb=f'/Users/thisUser/Desktop/test_FD.pdb'
    preloaded_fd = FoldedDomain(path_to_pdb)

    attractive_target=-5
    repulsive_target=5

    optimizer = goose.SequenceOptimizer(target_length=50, verbose=True, gap_to_report=100)
    optimizer.add_property(goose.FDSurfaceInteractionByValue, weight=1.0,
                            repulsive_target=repulsive_target, 
                            attractive_target=attractive_target,
                            model='mpipi', preloaded_fd=preloaded_fd)

    optimizer.set_optimization_params(max_iterations=5000, tolerance=2)
    best_seq=optimizer.run()
    print(best_seq)

    # double check
    loaded_IMC_object = finches.frontend.mpipi_frontend.Mpipi_frontend().IMC_object
    current_epsilon = preloaded_fd.calculate_surface_epsilon(best_seq, loaded_IMC_object)
    current_repulsive = np.sum([a[2] for a in current_epsilon.values() if a[2]>0])
    current_attractive = np.sum([a[2] for a in current_epsilon.values() if a[2]<0])
    print(current_repulsive)
    print(current_attractive)
    print(np.abs(repulsive_target-current_repulsive) + np.abs(attractive_target-current_attractive))


    """
    IS_NATURALLY_NORMALIZED = False  # Epsilon surface values are arbitrary scale
    
    def __init__(self,  
                 repulsive_target : float,
                 attractive_target : float,
                 weight: float = 1.0,
                 target_value: float = 0,
                 model = 'mpipi',
                 path_to_pdb: str = None,
                 probe_radius: float = 1.4,
                 surface_thresh: float = 0.10,
                 sasa_mode: str = 'v1',
                 fd_start : int = None,
                 fd_end : int = None,
                 preloaded_fd = None):
        super().__init__(target_value, weight)
        self.repulsive_target = repulsive_target
        self.attractive_target = attractive_target
        self.model = model
        self.path_to_pdb = path_to_pdb
        self.probe_radius = probe_radius
        self.sasa_mode = sasa_mode
        self.surface_thresh = surface_thresh
        self.fd_start = fd_start
        self.fd_end = fd_end

        # stuff to set to None and then update when this class is initiated    
        self.loaded_IMC_object = None
        if preloaded_fd != None:
            self.folded_domain = preloaded_fd
        else:
            self.folded_domain = None
        self.target_epsilon = None
    
    def load_IMC_object(self, model: str = None):
        """
        Load the IMC_object to be used for the calculation.

        Parameters:
        model (str): The model to be used for the calculation.
        """
        if model is None:
            model = self.model.lower()
        
        if self.loaded_IMC_object==None:
            if model == 'mpipi':
                self.loaded_IMC_object = finches.frontend.mpipi_frontend.Mpipi_frontend().IMC_object
            elif model == 'calvados':
                self.loaded_IMC_object = finches.frontend.calvados_frontend.CALVADOS_frontend().IMC_object
            else:
                raise ValueError(f"Model {model} not supported.")
        return self.loaded_IMC_object

    def load_folded_domain(self, path_to_pdb: str = None, start=None, end=None,
                           probe_radius: float = None, surface_thresh: float = None,
                           sasa_mode: str = None):
        """
        Load the folded domain to be used for the calculation.

        Parameters:
        path_to_pdb (str): The path to the pdb file for the folded domain.
        start (int): The start residue of the folded domain.
        end (int): The end residue of the folded domain.
        probe_radius (float): The probe radius to use for the calculation.
        surface_thresh (float): The surface threshold to use for the calculation.
        sasa_mode (str): The SASA mode to use for the calculation.
        """
        # only run if we haven't already loaded the folded domain
        if self.folded_domain is None:
            # get everything from the class if not provided
            if path_to_pdb is None:
                path_to_pdb = self.path_to_pdb
            if start is None:
                start = self.fd_start
            if end is None:
                end = self.fd_end
            if probe_radius is None:
                probe_radius = self.probe_radius
            if surface_thresh is None:
                surface_thresh = self.surface_thresh
            if sasa_mode is None:
                sasa_mode = self.sasa_mode
            
            # get folded domain
            self.folded_domain = FoldedDomain(path_to_pdb,
                                                        start=start,
                                                        end=end,
                                                        probe_radius=probe_radius,
                                                        surface_thresh=surface_thresh,
                                                        sasa_mode=sasa_mode)
        
        return self.folded_domain
    
    def calculate(self, protein: 'sparrow.Protein') -> float:
        # make sure we have the original epsilon vectors
        self.loadedloaded_IMC_object_model = self.load_IMC_object()
        self.folded_domain = self.load_folded_domain()
        current_epsilon = self.folded_domain.calculate_surface_epsilon(protein.sequence, self.loaded_IMC_object)
        current_repulsive = np.sum([a[2] for a in current_epsilon.values() if a[2]>0])
        current_attractive = np.sum([a[2] for a in current_epsilon.values() if a[2]<0])
        # get the diff
        return np.abs(self.repulsive_target-current_repulsive) + np.abs(self.attractive_target-current_attractive)
        

# make new class that makes sequences by net. 
class FDSurfaceInteractionByMean(ProteinProperty):
    """
    Try to get a specific Net attraction. 

    """
    IS_NATURALLY_NORMALIZED = False  # Mean surface epsilon values are arbitrary scale
    
    def __init__(self,  
                 weight: float = 1.0,
                 target_value: float = 0,
                 model = 'mpipi',
                 path_to_pdb: str = None,
                 probe_radius: float = 1.4,
                 surface_thresh: float = 0.10,
                 sasa_mode: str = 'v1',
                 fd_start : int = None,
                 fd_end : int = None,
                 preloaded_fd = None):
        super().__init__(target_value, weight)
        self.model = model
        self.path_to_pdb = path_to_pdb
        self.probe_radius = probe_radius
        self.sasa_mode = sasa_mode
        self.surface_thresh = surface_thresh
        self.fd_start = fd_start
        self.fd_end = fd_end

        # stuff to set to None and then update when this class is initiated    
        self.loaded_IMC_object = None
        if preloaded_fd != None:
            self.folded_domain = preloaded_fd
        else:
            self.folded_domain = None
        self.target_epsilon = None
    
    def load_IMC_object(self, model: str = None):
        """
        Load the IMC_object to be used for the calculation.

        Parameters:
        model (str): The model to be used for the calculation.
        """
        if model is None:
            model = self.model.lower()
        
        if self.loaded_IMC_object==None:
            if model == 'mpipi':
                self.loaded_IMC_object = finches.frontend.mpipi_frontend.Mpipi_frontend().IMC_object
            elif model == 'calvados':
                self.loaded_IMC_object = finches.frontend.calvados_frontend.CALVADOS_frontend().IMC_object
            else:
                raise ValueError(f"Model {model} not supported.")
        return self.loaded_IMC_object

    def load_folded_domain(self, path_to_pdb: str = None, start=None, end=None,
                           probe_radius: float = None, surface_thresh: float = None,
                           sasa_mode: str = None):
        """
        Load the folded domain to be used for the calculation.

        Parameters:
        path_to_pdb (str): The path to the pdb file for the folded domain.
        start (int): The start residue of the folded domain.
        end (int): The end residue of the folded domain.
        probe_radius (float): The probe radius to use for the calculation.
        surface_thresh (float): The surface threshold to use for the calculation.
        sasa_mode (str): The SASA mode to use for the calculation.
        """
        # only run if we haven't already loaded the folded domain
        if self.folded_domain is None:
            # get everything from the class if not provided
            if path_to_pdb is None:
                path_to_pdb = self.path_to_pdb
            if start is None:
                start = self.fd_start
            if end is None:
                end = self.fd_end
            if probe_radius is None:
                probe_radius = self.probe_radius
            if surface_thresh is None:
                surface_thresh = self.surface_thresh
            if sasa_mode is None:
                sasa_mode = self.sasa_mode
            
            # get folded domain
            self.folded_domain = FoldedDomain(path_to_pdb,
                                                        start=start,
                                                        end=end,
                                                        probe_radius=probe_radius,
                                                        surface_thresh=surface_thresh,
                                                        sasa_mode=sasa_mode)
        
        return self.folded_domain
    
    def calculate(self, protein: 'sparrow.Protein') -> float:
        # make sure we have the original epsilon vectors
        self.loaded_IMC_object = self.load_IMC_object()
        self.folded_domain = self.load_folded_domain()
        cur_net = self.folded_domain.calculate_mean_surface_epsilon(protein.sequence, self.loaded_IMC_object)
        # get the diff
        return cur_net


