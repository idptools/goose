from abc import ABC, abstractmethod
from typing import Any, Tuple, Dict
import sparrow
import numpy as np
import finches
import metapredict as meta

class ProteinProperty(ABC):
    """
    Abstract base class for protein sequence properties to be optimized.
    """
    def __init__(self, target_value: float, weight: float = 1.0):
        assert isinstance(target_value, (int,float)), f"target_value must be numerical. Received {type(target_value)}"
        assert isinstance(weight, (int,float)), f"Weight must be numerical. Received {type(weight)}"
        self._target_value = target_value
        self._weight = weight

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

    @abstractmethod
    def calculate(self, protein: 'sparrow.Protein') -> float:
        """
        Calculate the property value for a given protein.

        Parameters:
        protein (sparrow.Protein): The protein instance.

        Returns:
        float: The calculated property value.
        """
        pass


class ComputeIWD(ProteinProperty):
    """
    Compute the Inversed Weighted Distance (IWD) property for the target residues in the sequence.
    """
    def __init__(self, residues: Tuple[str, ...], target_value: float, weight: float = 1.0):
        super().__init__(target_value, weight)
        self.residues = residues

    def calculate(self, protein: 'sparrow.Protein') -> float:
        return protein.compute_iwd(*self.residues)

class Hydrophobicity(ProteinProperty):
    """
    Calculate the hydrophobicity property.
    """
    def __init__(self, target_value: float, weight: float = 1.0):
        super().__init__(target_value, weight)

    def calculate(self, protein: 'sparrow.Protein') -> float:
        return protein.hydrophobicity

class FCR(ProteinProperty):
    """
    Calculate the FCR property.
    """
    def __init__(self, target_value: float, weight: float = 1.0):
        super().__init__(target_value, weight)

    def calculate(self, protein: 'sparrow.Protein') -> float:
        return protein.FCR

class NCPR(ProteinProperty):
    """
    Calculate the NCPR property.
    """
    def __init__(self, target_value: float, weight: float = 1.0):
        super().__init__(target_value, weight)

    def calculate(self, protein: 'sparrow.Protein') -> float:
        return protein.NCPR

class Hydrophobicity(ProteinProperty):
    """
    Calculate the hydrophobicity property.
    """
    def __init__(self, target_value: float, weight: float = 1.0):
        super().__init__(target_value, weight)

    def calculate(self, protein: 'sparrow.Protein') -> float:
        return protein.hydrophobicity

class Kappa(ProteinProperty):
    """
    Calculate the kappa property.
    """
    def __init__(self, target_value: float, weight: float = 1.0):
        super().__init__(target_value, weight)

    def calculate(self, protein: 'sparrow.Protein') -> float:
        return protein.kappa

class RadiusOfGyration(ProteinProperty):
    """
    Calculate the Radius of gyration
    """
    def __init__(self, target_value: float, weight: float = 1.0):
        super().__init__(target_value, weight)

    def calculate(self, protein: 'sparrow.Protein') -> float:
        return protein.predictor.radius_of_gyration()

class EndToEndDistance(ProteinProperty):
    """
    Calculate the Radius of gyration
    """
    def __init__(self, target_value: float, weight: float = 1.0):
        super().__init__(target_value, weight)

    def calculate(self, protein: 'sparrow.Protein') -> float:
        return protein.predictor.end_to_end_distance()


class TargetAminoAcidFractions(ProteinProperty):
    """
    Compute the difference between the target amino acid fractions and the current amino acid fractions.
    """
    def __init__(self, target_fractions: Dict[str, float], weight: float = 1.0):
        """
        Parameters:
        target_fractions (Dict[str, float]): A dictionary where keys are single-letter amino acid codes,
            and values are the target fractions for each amino acid.
        weight (float): The weight of this property in the combined objective function.
        """
        super().__init__(0.0, weight)  # Target value is set to 0.0 since we want to minimize the difference
        self.target_fractions = target_fractions

    def calculate(self, protein: 'sparrow.Protein') -> float:
        """
        Calculate the difference between the target amino acid fractions and the current amino acid fractions.

        Parameters:
        protein (sparrow.Protein): The protein instance.

        Returns:
        float: The difference between the target and current amino acid fractions.
        """
        current_fractions = protein.amino_acid_fractions

        diff_sum = 0.0
        for aa, target_frac in self.target_fractions.items():
            current_frac = current_fractions.get(aa, 0.0)
            diff_sum += abs(target_frac - current_frac)

        return diff_sum

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
    def __init__(self, target_value: float, weight: float = 1.0):
        super().__init__(target_value, weight)

    def calculate(self, protein: 'sparrow.Protein') -> float:
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
    def __init__(self, target_value: float, weight: float = 1.0):
        super().__init__(target_value, weight)

    def calculate(self, protein: 'sparrow.Protein') -> float:
        return protein.SHD

class Complexity(ProteinProperty):
    """
    Calculates the Wootton-Federhen complexity of a sequence (also called
    seg complexity, as this the theory used in the classic SEG algorithm.
    """
    def __init__(self, target_value: float, weight: float = 1.0):
        super().__init__(target_value, weight)

    def calculate(self, protein: 'sparrow.Protein') -> float:
        return protein.complexity



class FractionDisorder(ProteinProperty):
    """
    Calculates the amount of disorder in
    the current sequence. 

    Example usage:

    # imports
    import goose
    from sparrow.protein import Protein as pr
    import metapredict as meta

    optimizer = goose.SequenceOptimizer(target_length=200, verbose=True, gap_to_report=100)
    optimizer.add_property(goose.FractionDisorder, weight=1, target_value=1)
    optimizer.set_optimization_params(max_iterations=2000, tolerance=0.01)
    sequence=optimizer.run()
    print(sequence)

    # double check
    disorder=meta.predict_disorder(sequence)
    print((disorder>0.5).sum()/len(disorder))

    """
    def __init__(self, 
                 target_value: float, 
                 weight: float = 1.0,
                 disorder_cutoff = 0.5):
        super().__init__(target_value, weight)
        self.disorder_cutoff = disorder_cutoff

    def calculate(self, protein: 'sparrow.Protein') -> float:
        disorder=meta.predict_disorder(protein.sequence)
        percent_disorder = (disorder>self.disorder_cutoff).sum()/len(disorder)
        return percent_disorder
    


class MatchingResidues(ProteinProperty):
    '''
    Determines the number of residues that match a target sequence in the current sequence.
    '''
    def __init__(self, 
                 target_sequence: str,
                 target_value: float, 
                 weight: float = 1.0):
        super().__init__(target_value, weight)
        self.target_sequence = target_sequence

    def calculate(self, protein: 'sparrow.Protein') -> float:
        return sum([1 for i in range(len(protein.sequence)) if protein.sequence[i] == self.target_sequence[i]])
    
    
class MaxMatchingResidues(ProteinProperty):
    '''
    Determines the number of residues that match a target sequence in the current sequence.
    Penalizes the score if the number of matching residues is greater than the target value.
    '''
    def __init__(self, 
                 target_sequence: str,
                 target_value: float, 
                 weight: float = 1.0):
        super().__init__(target_value, weight)
        self.target_sequence = target_sequence

    def calculate(self, protein: 'sparrow.Protein') -> float:
        matches = sum([1 for i in range(len(protein.sequence)) if protein.sequence[i] == self.target_sequence[i]])
        if matches > self.target_value:
            return matches - self.target_value
        else:
            return self.target_value