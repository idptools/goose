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
                            original_epsilon_vectors=original_vectors, finces_frontend=mpipi_frontend)
    optimizer.set_optimization_params(max_iterations=10000)
    opt=optimizer.run()


    """
    def __init__(self, 
                 target_value: float, 
                 weight: float = 1.0,
                 target_interacting_sequence: str = None,
                 original_epsilon_vectors: tuple = None,
                 model: str = 'mpipi'):
        super().__init__(target_value, weight)
        self.model = model
        self.target_interacting_sequence = target_interacting_sequence
        self.original_epsilon_vectors = original_epsilon_vectors
        self.loaded_model = None

    def load_model(self, model: str = None):
        """
        Load the model to be used for the calculation.

        Parameters:
        model (str): The model to be used for the calculation.
        """
        if model is None:
            model = self.model.lower()
        
        if self.loaded_model==None:
            if model == 'mpipi':
                self.loaded_model = finches.frontend.mpipi_frontend.Mpipi_frontend()
            elif model == 'calvados':
                self.loaded_model = finches.frontend.calvados_frontend.CALVADOS_frontend()
            else:
                raise ValueError(f"Model {model} not supported.")
        return self.loaded_model

    def calculate(self, protein: 'sparrow.Protein') -> float:
        self.loaded_model = self.load_model()
        current_vectors = self.loaded_model.epsilon_vectors(protein.sequence, self.target_interacting_sequence)
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
                            original_epsilon_value=original_epsilon, finces_frontend=mpipi_frontend)
    optimizer.set_optimization_params(max_iterations=1000, tolerance=1e-3)
    opt=optimizer.run()


    """
    def __init__(self, 
                 target_value: float, 
                 weight: float = 1.0,
                 original_sequence: str = None,
                 target_interacting_sequence: str = None,
                 model = 'mpipi'):
        super().__init__(target_value, weight)
        self.target_interacting_sequence : str = target_interacting_sequence
        self.original_sequence : str = original_sequence
        self.original_epsilon = None
        self.loaded_model = None
        self.model = model

    def load_model(self, model: str = None):
        """
        Load the model to be used for the calculation.

        Parameters:
        model (str): The model to be used for the calculation.
        """
        if model is None:
            model = self.model.lower()
        
        if self.loaded_model==None:
            if model == 'mpipi':
                self.loaded_model = finches.frontend.mpipi_frontend.Mpipi_frontend()
            elif model == 'calvados':
                self.loaded_model = finches.frontend.calvados_frontend.CALVADOS_frontend()
            else:
                raise ValueError(f"Model {model} not supported.")
        return self.loaded_model
    
    def get_original_epsilon(self, original_sequence: str = None, target_sequence: str = None):  
        """
        Calculate the original epsilon value of the target sequence.

        Parameters:
        sequence (str): The sequence.
        target_sequence (str): The target sequence.

        Returns:
        float: The original epsilon value of the target sequence.
        """
        if self.original_epsilon is not None:
            return self.original_epsilon
        
        if original_sequence is None:
            original_sequence = self.original_sequence
        if target_sequence is None:
            target_sequence = self.target_interacting_sequence
        self.loaded_model = self.load_model()
        self.original_epsilon = self.loaded_model.epsilon(original_sequence, target_sequence)
        return self.original_epsilon

    def calculate(self, protein: 'sparrow.Protein') -> float:
        self.original_epsilon_value = self.get_original_epsilon()
        self.loaded_model = self.load_model()
        current_epsilon = self.loaded_model.epsilon(protein.sequence, self.target_interacting_sequence)
        # calculate the difference between the original and current matrix
        return np.abs(self.original_epsilon_value-current_epsilon)
    

class chemical_fingerprint(ProteinProperty):
    """
    Uses the chemical foot print from the FINCHES manuscript to generate a sequence with
    a similar chemical fingerprint to the target sequence.

    The chemical fingerprint is calculated by taking the difference between the target
    and the current sequence and summing the differences in the epsilon matrix.
    """
    def __init__(self, 
                 target_value: float, 
                 weight: float = 1.0,
                 target_sequence: str = None,
                 model = 'mpipi'):
        super().__init__(target_value, weight)
        self.target_sequence = target_sequence
        self.loaded_model  = None
        self.target_sequence_fingerprint = None
        self.model = model
        self.chemistries = None

    def set_chemistries(self, model: str=None):
        """
        Set the chemistries to be used for the calculation.

        Parameters:
        model (str): The model to be used for the calculation.
        """
        if self.chemistries is not None:
            return self.chemistries
        if model is None:
            model = self.model.lower()
        if model == 'mpipi':
            self.chemistries = {'c1': 'KRKRKRKRKRKRKRKRKRKR', 'c2': 'RRRRRRRRRRRRRRRRRRRR', 'c3': 'KKKKKKKKKKKKKKKKKKKK', 'c4': 'HRHRHRHRHRHRHRHRHRHR', 'c5': 'HKHKHKHKHKHKHKHKHKHK', 'c6': 'AKAKAKAKAKAKAKAKAKAK', 'c7': 'KQKQKQKQKQKQKQKQKQKQ', 'c8': 'QRQRQRQRQRQRQRQRQRQR', 'c9': 'GRGRGRGRGRGRGRGRGRGR', 'c10': 'MRMRMRMRMRMRMRMRMRMR', 'c11': 'RWRWRWRWRWRWRWRWRWRW', 'c12': 'HHHHHHHHHHHHHHHHHHHH', 'c13': 'HYHYHYHYHYHYHYHYHYHY', 'c14': 'HQHQHQHQHQHQHQHQHQHQ', 'c15': 'HLHLHLHLHLHLHLHLHLHL', 'c16': 'AIAIAIAIAIAIAIAIAIAI', 'c17': 'AGAGAGAGAGAGAGAGAGAG', 'c18': 'ININININININININININ', 'c19': 'GQGQGQGQGQGQGQGQGQGQ', 'c20': 'GSGSGSGSGSGSGSGSGSGS', 'c21': 'DKDKDKDKDKDKDKDKDKDK', 'c22': 'ERERERERERERERERERER', 'c23': 'DHDHDHDHDHDHDHDHDHDH', 'c24': 'QQQQQQQQQQQQQQQQQQQQ', 'c25': 'FGFGFGFGFGFGFGFGFGFG', 'c26': 'QYQYQYQYQYQYQYQYQYQY', 'c27': 'AFAFAFAFAFAFAFAFAFAF', 'c28': 'LWLWLWLWLWLWLWLWLWLW', 'c29': 'FYFYFYFYFYFYFYFYFYFY', 'c30': 'WWWWWWWWWWWWWWWWWWWW', 'c31': 'DWDWDWDWDWDWDWDWDWDW', 'c32': 'EYEYEYEYEYEYEYEYEYEY', 'c33': 'EGEGEGEGEGEGEGEGEGEG', 'c34': 'EPEPEPEPEPEPEPEPEPEP', 'c35': 'AEAEAEAEAEAEAEAEAEAE', 'c36': 'DEDEDEDEDEDEDEDEDEDE'}
        elif model == 'calvados':
            self.chemistries = {'c1': 'DEDEDEDEDEDEDEDEDEDE', 'c2': 'EGEGEGEGEGEGEGEGEGEG', 'c3': 'DSDSDSDSDSDSDSDSDSDS', 'c4': 'EVEVEVEVEVEVEVEVEVEV', 'c5': 'DMDMDMDMDMDMDMDMDMDM', 'c6': 'EWEWEWEWEWEWEWEWEWEW', 'c7': 'AVAVAVAVAVAVAVAVAVAV', 'c8': 'APAPAPAPAPAPAPAPAPAP', 'c9': 'AQAQAQAQAQAQAQAQAQAQ', 'c10': 'GQGQGQGQGQGQGQGQGQGQ', 'c11': 'HQHQHQHQHQHQHQHQHQHQ', 'c12': 'QQQQQQQQQQQQQQQQQQQQ', 'c13': 'DRDRDRDRDRDRDRDRDRDR', 'c14': 'DKDKDKDKDKDKDKDKDKDK', 'c15': 'LMLMLMLMLMLMLMLMLMLM', 'c16': 'IMIMIMIMIMIMIMIMIMIM', 'c17': 'LWLWLWLWLWLWLWLWLWLW', 'c18': 'FMFMFMFMFMFMFMFMFMFM', 'c19': 'AWAWAWAWAWAWAWAWAWAW', 'c20': 'IQIQIQIQIQIQIQIQIQIQ', 'c21': 'LQLQLQLQLQLQLQLQLQLQ', 'c22': 'HLHLHLHLHLHLHLHLHLHL', 'c23': 'AMAMAMAMAMAMAMAMAMAM', 'c24': 'FQFQFQFQFQFQFQFQFQFQ', 'c25': 'HWHWHWHWHWHWHWHWHWHW', 'c26': 'FWFWFWFWFWFWFWFWFWFW', 'c27': 'MRMRMRMRMRMRMRMRMRMR', 'c28': 'RWRWRWRWRWRWRWRWRWRW', 'c29': 'KYKYKYKYKYKYKYKYKYKY', 'c30': 'ARARARARARARARARARAR', 'c31': 'KLKLKLKLKLKLKLKLKLKL', 'c32': 'AKAKAKAKAKAKAKAKAKAK', 'c33': 'KSKSKSKSKSKSKSKSKSKS', 'c34': 'KKKKKKKKKKKKKKKKKKKK', 'c35': 'KRKRKRKRKRKRKRKRKRKR', 'c36': 'RRRRRRRRRRRRRRRRRRRR'}
        else:
            raise ValueError(f"Model {model} not supported.")
        return self.chemistries
        
    
    def load_model(self, model: str = None):
        """
        Load the model to be used for the calculation.

        Parameters:
        model (str): The model to be used for the calculation.
        """
        if model is None:
            model = self.model.lower()
        
        if self.loaded_model==None:
            if model == 'mpipi':
                self.loaded_model = finches.frontend.mpipi_frontend.Mpipi_frontend()
            elif model == 'calvados':
                self.loaded_model = finches.frontend.calvados_frontend.CALVADOS_frontend()
            else:
                raise ValueError(f"Model {model} not supported.")
        return self.loaded_model
    
    def calculate_fingerprint(self, sequence: str = None, chemistries: dict = None):
        """
        Calculate the chemical fingerprint of the target sequence.

        Parameters:
        sequence (str): The sequence.
        chemistries (dict): The chemistries to use for the calculation.

        Returns:
        dict: The chemical fingerprint of the target sequence
        """
        if chemistries is None:    
            chemistries = self.set_chemistries()
        
        # make sure the model is loaded
        self.loaded_model = self.load_model()

        # now for each chemistry, we calculate the epsilon vectors.
        sequence_fingerprint = {}
        for chemistry, chemistry_sequence in chemistries.items():
            cur_vec = self.loaded_model.epsilon_vectors(sequence, chemistry_sequence)
            sequence_fingerprint[chemistry] = {'attractive':cur_vec[0], 'repulsive':cur_vec[1]}
        return sequence_fingerprint



    def get_target_fingerprint(self, target_sequence: str = None, chemistries: dict = None):
        """
        Calculate the chemical fingerprint of the target sequence.

        Parameters:
        target_sequence (str): The target sequence.
        chemistries (dict): The chemistries to use for the calculation.

        Returns:
        dict: The chemical fingerprint of the target sequence
        """
        if self.target_sequence_fingerprint is not None:
            return self.target_sequence_fingerprint
        if target_sequence is None:
            target_sequence = self.target_sequence
        if chemistries is None:
            if self.chemistries is None:
                self.chemistries = self.set_chemistries()
            chemistries = self.chemistries
        self.target_sequence_fingerprint = self.calculate_fingerprint(target_sequence, chemistries)
        return self.target_sequence_fingerprint
    
    
    def calculate(self, protein: 'sparrow.Protein') -> float:
        current_fingerprint = self.calculate_fingerprint(protein.sequence)
        target_fingerprint = self.get_target_fingerprint()
        # calculate the difference between the original and current matrix
        diff = 0
        for key in target_fingerprint.keys():
            diff += np.abs(target_fingerprint[key]['attractive'] - current_fingerprint[key]['attractive']).sum()
            diff += np.abs(target_fingerprint[key]['repulsive'] - current_fingerprint[key]['repulsive']).sum()
        return diff