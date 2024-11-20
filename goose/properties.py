from abc import ABC, abstractmethod
from typing import Any, Tuple, Dict
import sparrow
import numpy as np
import finches

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

# initialize mpipi frontend
mpipi_frontend=finches.frontend.mpipi_frontend.Mpipi_frontend()


class epsilon_vector_diff(ProteinProperty):
    """
    Calculate the difference in the epsilon attractive and repulsive matrices,
    returns the diff.

    Usage example:
    seq='NGDNFNRTPASSSEMDDGPSRRDHFMKSGFASGRNFGNRDAGECNKRDNTSTMG'
    target='ATQSYGAYPTQPGQGYSQQSSQPYGQQSYSGYSQSTDTSGYGQSSYSSYGQSQNTGYGT'
    original_vectors = mpipi_frontend.epsilon_vectors(seq, target)
    optimizer = sequence_optimization.SequenceOptimizer(target_length=len(seq), verbose=True, gap_to_report=100)
    optimizer.add_property(epsilon_vector_diff, target_value=0.0, weight=1.0, target_interacting_sequence=target, 
                            original_epsilon_vectors=original_vectors, loaded_frontend=mpipi_frontend)
    optimizer.set_optimization_params(max_iterations=10000)
    opt=optimizer.run()


    """
    def __init__(self, 
                 target_value: float, 
                 weight: float = 1.0,
                 target_interacting_sequence: str = None,
                 original_epsilon_vectors: tuple = None,
                 loaded_frontend: finches.frontend.mpipi_frontend.Mpipi_frontend = mpipi_frontend):
        super().__init__(target_value, weight)
        self.target_interacting_sequence = target_interacting_sequence
        self.original_epsilon_vectors = original_epsilon_vectors
        self.mpipi_frontend = mpipi_frontend

    def calculate(self, protein: 'sparrow.Protein') -> float:
        current_vectors = self.mpipi_frontend.epsilon_vectors(protein.sequence, self.target_interacting_sequence)
        attractive_vectors=current_vectors[0]
        repulsive_vectors=current_vectors[1]
        # calculate the difference between the original and current matrix
        attractive_diff = np.abs(self.original_epsilon_vectors[0] - attractive_vectors).sum()
        repulsive_diff = np.abs(self.original_epsilon_vectors[1] - repulsive_vectors).sum()
        return attractive_diff + repulsive_diff


class epsilon_total(ProteinProperty):
    """
    Calculate the difference in total epsilon value between 2 sequences. 

    Usage example:
    seq='MGDEDWEAEINPHMSSYVPIFEKDRYSGENGDNFNRTPASSSEMDDGPSRRDHFMKSGFASGRNFGNRDAGECNKRDNTSTMGGFGVGKSFGNRGFSNSR'
    target='MASNDYTQQATQSYGAYPTQPGQGYSQQSSQPYGQQSYSGYSQSTDTSGYGQSSYSSYGQSQNTGYGTQSTPQGYGSTGGYGSSQSSQSSYGQQSSYPGY'
    original_epsilon = mpipi_frontend.epsilon(seq, target)
    optimizer = sequence_optimization.SequenceOptimizer(target_length=len(seq), verbose=True, gap_to_report=100)
    optimizer.add_property(epsilon_total, target_value=0.0, weight=1.0, target_interacting_sequence=target, 
                            original_epsilon_value=original_epsilon, loaded_frontend=mpipi_frontend)
    optimizer.set_optimization_params(max_iterations=1000, tolerance=1e-3)
    opt=optimizer.run()


    """
    def __init__(self, 
                 target_value: float, 
                 weight: float = 1.0,
                 target_interacting_sequence: str = None,
                 original_epsilon_value: float = None,
                 loaded_frontend: finches.frontend.mpipi_frontend.Mpipi_frontend = mpipi_frontend):
        super().__init__(target_value, weight)
        self.target_interacting_sequence = target_interacting_sequence
        self.original_epsilon_value = original_epsilon_value
        self.mpipi_frontend = mpipi_frontend

    def calculate(self, protein: 'sparrow.Protein') -> float:
        current_epsilon = self.mpipi_frontend.epsilon(protein.sequence, self.target_interacting_sequence)
        # calculate the difference between the original and current matrix
        return np.abs(self.original_epsilon_value-current_epsilon)