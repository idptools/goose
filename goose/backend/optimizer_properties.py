from abc import ABC, abstractmethod
from typing import Any, Tuple, Dict, Callable
import inspect
import numpy as np

import sparrow
import metapredict as meta
import finches
from finches.utils.folded_domain_utils import FoldedDomain

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

    Example usage:
    import goose
    from sparrow.protein import Protein as pr
    # initialize optimizer
    optimizer=goose.SequenceOptimizer(target_length=50, 
                                gap_to_report=1000, num_shuffles=1)

    # set optimization parameters
    optimizer.set_optimization_params(max_iterations=50000, 
                                    tolerance=1e-2,
                                    shuffle_interval=5)
    # add IWD as a property to optimize
    optimizer.add_property(goose.IWD, residues=('S'), target_value=3)
    # set best_seq to the optimized sequence
    best_seq=optimizer.run()
    # print the IWD for S for the best_seq and then print it. 
    print(pr(best_seq).compute_iwd('S'))
    print(best_seq)
    """
    def __init__(self, residues: Tuple[str, ...], target_value: float, weight: float = 1.0):
        super().__init__(target_value, weight)
        self.residues = residues

    def calculate(self, protein: 'sparrow.Protein') -> float:
        return protein.compute_iwd(list(self.residues))

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

class MaxFractions(ProteinProperty):
    """
    Compute the difference between the max amino acid fractions and the current amino acid fractions.

    Example usage:
    import goose
    from sparrow.protein import Protein as pr

    # initialize optimizer
    optimizer=goose.SequenceOptimizer(target_length=100, 
                                gap_to_report=1000, num_shuffles=1)

    # set optimization params
    optimizer.set_optimization_params(max_iterations=10000, 
                                    tolerance=1e-2,
                                    shuffle_interval=5)

    # add property
    optimizer.add_property(goose.MaxFractions, {'A':0.1, 'S':0.1, 'G':0.1})
    best_seq=optimizer.run()

    print(pr(best_seq).amino_acid_fractions)
    print(best_seq)    
    """
    def __init__(self, max_fractions: Dict[str, float], weight: float = 1.0):
        """
        Parameters:
        max_fractions (Dict[str, float]): A dictionary where keys are single-letter amino acid codes,
            and values are the max fractions for each amino acid.
        weight (float): The weight of this property in the combined objective function.
        
        """
        super().__init__(0.0, weight)  # Target value is set to 0.0 since we want to minimize the difference
        self.max_fractions = max_fractions

    def calculate(self, protein: 'sparrow.Protein') -> float:
        """
        Calculate the difference between the max amino acid fractions and the current amino acid fractions.

        Parameters:
        protein (sparrow.Protein): The protein instance.

        Returns:
        float: The difference between the max and current amino acid fractions.
        """
        current_fractions = protein.amino_acid_fractions

        diff_sum = 0.0
        for aa, target_frac in self.max_fractions.items():
            current_frac = current_fractions.get(aa, 0.0)
            if current_frac > target_frac:
                diff_sum += abs(target_frac - current_frac)

        return diff_sum

class MinFractions(ProteinProperty):
    """
    Compute the difference between the max amino acid fractions and the current amino acid fractions.

    Example usage:
    import goose
    from sparrow.protein import Protein as pr

    # initialize optimizer
    optimizer=goose.SequenceOptimizer(target_length=100, 
                                gap_to_report=1000, num_shuffles=1)

    # set optimization params
    optimizer.set_optimization_params(max_iterations=10000, 
                                    tolerance=1e-2,
                                    shuffle_interval=5)

    # add property
    optimizer.add_property(goose.MinFractions, {'A':0.1, 'S':0.1, 'G':0.1})
    best_seq=optimizer.run()

    print(pr(best_seq).amino_acid_fractions)
    print(best_seq)    
    """
    def __init__(self, min_fractions: Dict[str, float], weight: float = 1.0):
        """
        Parameters:
        min_fractions (Dict[str, float]): A dictionary where keys are single-letter amino acid codes,
            and values are the min fractions for each amino acid.
        weight (float): The weight of this property in the combined objective function.
        
        """
        super().__init__(0.0, weight)  # Target value is set to 0.0 since we want to minimize the difference
        self.min_fractions = min_fractions

    def calculate(self, protein: 'sparrow.Protein') -> float:
        """
        Calculate the difference between the min amino acid fractions and the current amino acid fractions.

        Parameters:
        protein (sparrow.Protein): The protein instance.

        Returns:
        float: The difference between the min and current amino acid fractions.
        """
        current_fractions = protein.amino_acid_fractions

        diff_sum = 0.0
        for aa, target_frac in self.max_fractions.items():
            current_frac = current_fractions.get(aa, 0.0)
            if current_frac < target_frac:
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
    

class MatchSequenceDisorder(ProteinProperty):
    """
    Matches the disorder of a specified sequence. By default, the 
    specified sequence sets the MINIMUM disorder. 
    You can also try to get the exact disorder.

    Example usage:

    # imports
    import goose
    from sparrow.protein import Protein as pr
    import metapredict as meta

    target_seq='MPHIQKKAQGWCNNPGQQLNPHLQWQWQNQQYIDFDPYNY'
    optimizer = goose.SequenceOptimizer(target_length=200, verbose=True, gap_to_report=100)
    optimizer.add_property(goose.MatchSequenceDisorder, target_sequence=target_seq, weight=1)
    optimizer.set_optimization_params(max_iterations=2000, tolerance=0.01)
    sequence=optimizer.run()
    print(sequence)

    # double check
    disorder=meta.predict_disorder(sequence)
    print((disorder>0.5).sum()/len(disorder))

    """
    def __init__(self, 
                 target_sequence: str, 
                 exact_match: bool = False,
                 target_value : float = 0,
                 weight: float = 1.0):
        super().__init__(target_value, weight)
        self.target_sequence = target_sequence
        self.exact_match = exact_match
        self.target_disorder = None


    def set_initial_disorder(self, target_sequence=None):
        if type(self.target_disorder) != np.array:
            if target_sequence==None:
                target_sequence=self.target_sequence
            self.target_disorder = meta.predict_disorder(target_sequence)
        return self.target_disorder

    def calculate(self, protein: 'sparrow.Protein') -> float:
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
          


class EpsilonVectorBySequence(ProteinProperty):
    """
    Calculate the difference in the epsilon attractive and repulsive matrices,
    returns the diff.

    Usage example:
    # imports
    import numpy as np
    import finches
    import goose
    from sparrow.protein import Protein as pr


    seq='NGDNFNRTPASSSEMDDGPSRRDHFMKSGFASGRNFGNRDAGECNKRDNTSTMG'
    target='ATQSYGAYPTQPGQGYSQQSSQPYGQQSYSGYSQSTDTSGYGQSSYSSYGQSQNTGYGT'
    optimizer = goose.SequenceOptimizer(target_length=len(seq), verbose=True, gap_to_report=100)
    optimizer.add_property(goose.EpsilonVectorBySequence, weight=1.0, 
                            original_sequence=seq,
                            target_interacting_sequence=target, 
                            model='mpipi')
    optimizer.set_optimization_params(max_iterations=5000, tolerance=1)
    best_seq=optimizer.run()
    print(best_seq)

    # double check
    model=finches.frontend.mpipi_frontend.Mpipi_frontend()
    original_vectors=model.epsilon_vectors(seq, target)
    current_vectors = model.epsilon_vectors(best_seq, target)
    attractive_vectors=current_vectors[0]
    repulsive_vectors=current_vectors[1]

    attractive_diff = np.abs(original_vectors[0] - attractive_vectors).sum()
    repulsive_diff = np.abs(original_vectors[1] - repulsive_vectors).sum()

    print(attractive_diff+repulsive_diff)

    """
    def __init__(self, 
                 original_sequence: str,
                 target_interacting_sequence: str,
                 target_value: float = 0, 
                 weight: float = 1.0,
                 model: str = 'mpipi'):
        super().__init__(target_value, weight)
        self.model = model
        self.original_sequence = original_sequence
        self.target_interacting_sequence = target_interacting_sequence
        self.original_epsilon_vectors = None
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

    def get_original_epsilon_vectors(self, original_sequence: str = None, target_sequence: str = None):  
        """
        Calculate the original epsilon value of the target sequence.

        Parameters:
        sequence (str): The sequence.
        target_sequence (str): The target sequence.

        Returns:
        float: The original epsilon value of the target sequence.
        """
        if self.original_epsilon_vectors is not None:
            return self.original_epsilon_vectors
        
        if original_sequence is None:
            original_sequence = self.original_sequence
        if target_sequence is None:
            target_sequence = self.target_interacting_sequence
        self.loaded_model = self.load_model()
        self.original_epsilon_vectors = self.loaded_model.epsilon_vectors(original_sequence, target_sequence)
        return self.original_epsilon_vectors


    def calculate(self, protein: 'sparrow.Protein') -> float:
        # make sure we have the original epsilon vectors
        self.original_epsilon_vectors = self.get_original_epsilon_vectors()
        self.loaded_model = self.load_model()
        current_vectors = self.loaded_model.epsilon_vectors(protein.sequence, self.target_interacting_sequence)
        attractive_vectors=current_vectors[0]
        repulsive_vectors=current_vectors[1]
        # calculate the difference between the original and current matrix
        attractive_diff = np.abs(self.original_epsilon_vectors[0] - attractive_vectors).sum()
        repulsive_diff = np.abs(self.original_epsilon_vectors[1] - repulsive_vectors).sum()
        return attractive_diff + repulsive_diff




class EpsilonByValue(ProteinProperty):
    """
    Make a sequence that interacts with a specific sequence with a specific epsilon value. 

    Usage example:
    import finches
    import goose
    from sparrow.protein import Protein as pr

    seq='MGDEDWEAEINPHMSSYVPIFEKDRYSGENGDNFNRTPASSSEMDDGPSRRDHFMKSGFASGRNFGNRDAGECNKRDNTSTMGGFGVGKSFGNRGFSNSR'
    optimizer = goose.SequenceOptimizer(target_length=len(seq), verbose=True, gap_to_report=100)
    optimizer.add_property(goose.properties.EpsilonByValue, target_value=5, weight=1.0, 
                            target_sequence=seq,
                            model='mpipi')

    optimizer.set_optimization_params(max_iterations=5000, tolerance=1e-2)
    best_seq=optimizer.run()
    print(best_seq)

    # double check
    model=finches.frontend.mpipi_frontend.Mpipi_frontend()
    interaction_value=model.epsilon(best_seq, seq)

    print(interaction_value)


    """
    def __init__(self, 
                 target_value: float, 
                 target_sequence: str,
                 weight: float = 1.0,
                 model = 'mpipi'):
        super().__init__(target_value, weight)
        self.target_sequence : str = target_sequence
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


    def calculate(self, protein: 'sparrow.Protein') -> float:
        self.loaded_model = self.load_model()
        current_epsilon = self.loaded_model.epsilon(protein.sequence, self.target_sequence)
        # calculate the difference between the original and current matrix
        return current_epsilon


class EpsilonBySequence(ProteinProperty):
    """
    Calculate the difference in total epsilon value between 2 sequences. 

    Usage example:

    import finches
    import goose
    from sparrow.protein import Protein as pr

    seq='MGDEDWEAEINPHMSSYVPIFEKDRYSGENGDNFNRTPASSSEMDDGPSRRDHFMKSGFASGRNFGNRDAGECNKRDNTSTMGGFGVGKSFGNRGFSNSR'
    target='MASNDYTQQATQSYGAYPTQPGQGYSQQSSQPYGQQSYSGYSQSTDTSGYGQSSYSSYGQSQNTGYGTQSTPQGYGSTGGYGSSQSSQSSYGQQSSYPGY'

    optimizer = goose.SequenceOptimizer(target_length=len(seq), verbose=True, gap_to_report=100)

    optimizer.add_property(goose.EpsilonBySequence, weight=1.0, 
                            original_sequence=seq,
                            target_interacting_sequence=target,
                            model='mpipi')

    optimizer.set_optimization_params(max_iterations=5000, tolerance=1e-2)
    best_seq=optimizer.run()
    print(best_seq)

    # double check
    model=finches.frontend.mpipi_frontend.Mpipi_frontend()
    orignal_interaction_value = model.epsilon(seq, target)
    designed_seq_interaction_value=model.epsilon(best_seq, target)

    print(orignal_interaction_value, '\n', designed_seq_interaction_value)

    """
    def __init__(self, 
                 original_sequence: str,
                 target_interacting_sequence: str,                 
                 target_value: float = 0, 
                 weight: float = 1.0,
                
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


class SelfEpsilon(ProteinProperty):
    """
    Calculate the self interaction epsilon value of a sequence.

    Example usage:

    import finches
    import goose
    from sparrow.protein import Protein as pr

    optimizer = goose.SequenceOptimizer(target_length=100, verbose=True, gap_to_report=100)
    optimizer.add_property(goose.SelfEpsilon, target_value=5, weight=1.0,
                            model='mpipi')

    optimizer.set_optimization_params(max_iterations=5000, tolerance=1e-2)
    best_seq=optimizer.run()
    print(best_seq)

    # double check
    model=finches.frontend.mpipi_frontend.Mpipi_frontend()
    interaction_value=model.epsilon(best_seq, best_seq)

    print(interaction_value)

    """    
    def __init__(self, 
                 target_value: float, 
                 weight: float = 1.0,
                 model = 'mpipi',
                 preloaded_model=None):
        super().__init__(target_value, weight)
        self.loaded_model = None
        self.model = model

        if preloaded_model != None:
            self.loaded_model = preloaded_model

    def load_model(self, model: str = None):
        """
        Load the model to be used for the calculation.

        Parameters:
        model (str): The model to be used for the calculation.
        """
        if self.loaded_model==None:
            # get self.model
            if model is None:
                model = self.model.lower()
            # make sure is valid, then load
            if model == 'mpipi':
                self.loaded_model = finches.frontend.mpipi_frontend.Mpipi_frontend()
            elif model == 'calvados':
                self.loaded_model = finches.frontend.calvados_frontend.CALVADOS_frontend()
            else:
                raise ValueError(f"Model {model} not supported.")
        return self.loaded_model
    
    
    def calculate(self, protein: 'sparrow.Protein') -> float:
        self.loaded_model = self.load_model()
        current_epsilon = self.loaded_model.epsilon(protein.sequence, protein.sequence)
        return current_epsilon


class FDSurfaceInteractionByValue(ProteinProperty):
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
        



class ChemicalFingerprint(ProteinProperty):
    """
    Uses the chemical foot print from the FINCHES manuscript to generate a sequence with
    a similar chemical fingerprint to the target sequence.

    The chemical fingerprint is calculated by taking the difference between the target
    and the current sequence and summing the differences in the epsilon matrix.

    usage example:
    import goose

    target='MASNDYTQQATQSYGAYPTQPGQGYSQQSSQPYGQQSYSGYSQSTDTSGYGQSSYSSYGQSQNTGYGTQSTPQGYGSTGGYGSSQSSQSSYGQQSSYPGY'
    optimizer = goose.SequenceOptimizer(target_length=len(target), verbose=True, gap_to_report=100)
    optimizer.add_property(goose.ChemicalFingerprint, target_value=0.0, weight=1.0, 
                            target_sequence=target, 
                            model='mpipi')
    optimizer.set_optimization_params(max_iterations=2000, tolerance=100)
    opt=optimizer.run()

    """
    def __init__(self, 
                 target_sequence: str,
                 target_value: float = 0.0, 
                 weight: float = 1.0,
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
    



class Comparison(ProteinProperty):
    """
    This is designed to perform generic comparisons of a protein property to a reference.
    By packaging the property and comparison in a function, it allows users to 
    easily design a whole host of comparisons without the need to write a new
    class for it each time.
    """
    def __init__(self,
                 property_function : Callable[[str],Any],
                 comparison_function : Callable[[Any,Any], float],
                 reference_value : Any,
                 target_value: float, 
                 weight: float = 1.0,
                 check_function_compatibility : bool = True,
                 check_reference_value_compatibility : bool = True):
        
        '''
        
        Parameters
        ----------
        property_function : Callable[[str],Any]
            This is the function that computes some property of a sequence.
            For this function you do not need to specify the return type
            unless you are planning on chekcing the reference values type.
            If this is the case, including the type will prevent unwanted Errors
            from being thrown during the checking process. You MUST specify the 
            input parameter type as a string if you are are going to check function
            compatibility. Whether or not you specify the types, the input to this 
            function MUST be a str in the first parameter and the output should be 
            the same type as the reference value.
        comparison_function : Callable[[Any,Any], float]
            This is the function that compares the output property of interest
            to a reference value. This allows you to shoot for some target value 
            of interest. You do not need to specify the input type of this function
            unless you are checking the reference value for compatibility. You MUST
            specify the return type as float if you are checking the function compatibility.
            Regardless of whether you decide to specify the types, this function
            MUST return a float. It should also take two values of the type returned by
            property_function. The first parameter is reserved for the sequence that will be
            tested by the optimizer. The second parameter will be the reference value.
            During computation time this function gets called like...
            comparison_function(test_sequence_property_result, reference_value).
            The output of this is directly fed into the optimizer.
        reference_value : Any
            This is the property value you will compare the generated sequences against.
            This should be the same type as what is returned by property_function
            and accepted by comparison_function.
        target_value : float
            This is the value returned by the comparison function that you wish to achieve.
        weight : float
            This is a relative weight to give this property with respect to other properties you are also trying
            to optimize for in your design.
        check_function_compatibility : bool
            This determines if the object checks the function(s) compatibility
            with one another in terms of inputs and outputs.
        check_reference_value_compatibility : bool
            This determines if the object will check the reference values
            types against those that are output and accepted by the two user defined
            functions

        Returns
        -------
        Comparison(ProteinProperty)
            This is the object that will be used by the optimizer to find a sequnce
            that suits your needs.
        '''
        #intialize the parent class properties
        super().__init__(target_value, weight)

        #initialize object specific code to None if needed
        self._property_function = None
        self._comparison_function = None
        self._check_function_compatibility = None
        self._check_reference_value_compatibiliity = None

        #set wether to check the functions or not for compatibility
        self._check_function_compatibility = check_function_compatibility
        self._check_reference_value_compatibiliity = check_reference_value_compatibility

        #check that the passed values for the property function and comparison function PASSED are not none
        #This avoid the awkward case where the user initializes with None to get around the compatibility
        #later on in the calculation. AKA they could set the one of the functions later on to None and 
        #this would not trigger a test.
        if property_function is None:
            raise ValueError(f"The property function cannot equal None. Please pass a valid function instead.")
        if comparison_function is None:
            raise ValueError(f"The comparison function cannot be None. Please pass a valid function.")

        #Assign the property function (the check is performed inside the setter)
        self.property_function = property_function
        self.comparison_function = comparison_function

        #set the reference value (THIS MUST HAPPEN AFTER THE FUNCTIONS are assigned - checks fail otherwise)
        self.reference_value = reference_value


    @property
    def check_function_compatibility(self) -> bool:
        '''Returns wether or not to check the functions passed for compatibility'''
        return self._check_function_compatibility
    
    @check_function_compatibility.setter
    def check_function_compatibility(self, val : bool) -> None:
        '''Sets the check_function_compatilibity property and checks that it is a bool.'''
        #ensure that the value being set is a boolean
        if not isinstance(val, bool):
            raise TypeError(f"The value assigned to check_function_compatibility is not a bool ({type(val)}).")
        
        #if it does not raise (is a bool) then set the new value
        self._check_function_compatibility = val
    
    @property 
    def check_reference_value_compatibility(self) -> bool:
        '''Returns whether or not to check the reference value type against the functions'''
        return self._check_reference_value_compatibiliity
    
    @check_reference_value_compatibility.setter
    def check_reference_value_compatibility(self, val : bool) -> None:
        '''Sets the boolean associated with checking the reference values typing. Check that the value is in fact a bool'''
        #check that the value is a bool
        if not isinstance(val, bool):
            raise TypeError(f"The value assigned to check_reference_vale_compatibility is not a bool ({type(val)})")
        
        #set the value if it did not throw a type error
        self._check_reference_value_compatibiliity = val
    
    @property
    def property_function(self) -> Callable[[str], Any]:
        '''Returns the property function'''
        return self._property_function
    
    @property_function.setter
    def property_function(self, func : Callable[[str],Any]) -> None:
        '''Sets the value of the property function and checks for compatibility if requested.'''
        #check if the passed value is callable
        if not callable(func):
            raise ValueError(f"The property function assignment could not be made because '{func.__name__}' was not callable")
        
        #See if the comparison function exists for if we need to check compatibility
        compare_func_exists = self.comparison_function is not None
        #check if we need to see if the functions are compatible
        if self.check_function_compatibility and compare_func_exists:
            self._validate_function_compatibility(func1=func, func2=self.comparison_function)
        
        #if no exceptions are raised assign the value
        self._property_function = func

    @property
    def comparison_function(self) -> Callable[[Any, Any], float]:
        '''Returns the value of the comparison function'''
        return self._comparison_function
    
    @comparison_function.setter
    def comparison_function(self, func : Callable[[Any,Any], float]) -> None:
        '''Sets the value for the comparison function and check for compatibility if requested.'''
        #check if the passed value is callable
        if not callable(func):
            raise ValueError(f"The comparison function assignment could not be made because '{func.__name__}' was not callable")
        
        #See if the property function exists for if we need to check for compatibility
        property_func_exists = self.property_function is not None
        #check if we need to look for function compatibility
        if self.check_function_compatibility and property_func_exists:
            self._validate_function_compatibility(func1=self._property_function,func2=func)

        #if no exceptions are raised assign the value
        self._comparison_function = func

    @property
    def reference_value(self) -> Any:
        '''Returns the reference value that you are comparing against.'''
        return self._reference_value
    
    @reference_value.setter
    def reference_value(self, val : Any) -> None:
        '''Sets the reference value and checks its type if needed'''
        #check the type for the reference if that was specified
        if self.check_reference_value_compatibility:
            self._check_reference_value_type(val)
        
        #set the value if it did not throw an error
        self._reference_value = val
    
    
    def _check_reference_value_type(self, value : Any) -> None:
        '''This function checks if the reference value is the same type as the output of the property function and input of the comparison function
        
        Parameters
        ----------
        value : Any
            This is the value of the reference value to check

        Returns
        -------
        None
            This function only raises an exception if something goes wrong and does NOT
            return a boolean value.
        '''
        #get the output value of the property function and the input type for the comparison function
        out1_type = self._get_function_return_type(self.property_function)
        in1_type = self._get_function_first_parameter_type(self.comparison_function)

        #check if the reference value does not adhere to either type
        if not isinstance(value, out1_type):
            raise TypeError(f"The type of the reference value ({type(self.reference_value)}) does not match the output of the property function ({out1_type}).")
        if not isinstance(value, in1_type):
            raise TypeError(f"The type of the reference value ({type(self.reference_value)}) does not match the output of the property function ({in1_type}).")
        


    def _get_function_details(self, func : Callable) -> Tuple[Dict[str,object],object]:
        '''Check what parameters are to be passed to a function and their classes.
        
        Parameters
        ----------
        func : Callable
            This is the function you wish to learn mre information about

        Returns
        -------
        Tuple[Dict[str,object],Dict[str,object]]
            This tuple contains a dictionary and an object. The dictionary comes first
            It contains all the inputs to the function. The keys for the dictionary 
            are the names of the parameters while the value for the keys are the objects
            that parameter takes as input.
        '''
        # Get the function signature
        signature = inspect.signature(func)
        
        # Extract parameter details
        inputs = {
            name: param.annotation if param.annotation != inspect.Parameter.empty else "Any"
            for name, param in signature.parameters.items()
        }
        
        # Extract return type
        output = signature.return_annotation if signature.return_annotation != inspect.Signature.empty else "Any"
        
        return inputs, output

    def _validate_function_compatibility(self, func1: Callable, func2: Callable) -> None:
        """Validate that the functions are composable and the composition takes a string
        Validates that the output type of `func1` matches the input type of `func2`.
        This function does not return a true or false. It simply raises Exceptions
        if the functions cannot be applied as follows...
        func2(func1(params), other_params) -> float
        This function also validates that the input to func1 is a string and the
        second functions output is a float

        Parameters
        ----------
        func1 : Callable
            This is the innermost function in the composition
        func2: Callable
            This is the outtermost function in the composition

        Returns
        -------
        None
            This function will raise an exception if there is an issue rather
            that return a True/False value.
        """
        # Get the input type on the first function is a string or throw an error
        input1_type = self._get_function_first_parameter_type(func=func1)
        if input1_type is not str:
            raise TypeError(f"The input type to '{func1.__name__}' in the composition must a be a string. It is currently {input1_type}.")
        
        # Check that the output of function 2 is a float
        output2_type = self._get_function_return_type(func=func2)
        if output2_type is not float:
            raise TypeError(f"The output of the function composition must be a float. This means that '{func2.__name__}' must have a return type of float rather than {output2_type}")
        
        # Check that the first functions output matches the second functions input
        output1_type = self._get_function_return_type(func=func1)
        input2_type = self._get_function_first_parameter_type(func=func2)

        # Validate compatibility of the functions outputs and inputs
        if output1_type != input2_type:
            raise TypeError(f"Type mismatch: '{func1.__name__}' returns {output1_type}, which cannot be passed as input to '{func2.__name__}' expecting {input2_type}.")

    def _get_function_first_parameter_type(self, func : Callable) -> object:
        '''Gets the first first parameter type of the function
        
        Parameters
        ----------
        func : Callable
            The function to determine the first parameters type for

        Returns
        -------
        object
            This the type of the first parameter of the function passed
        '''

        # Ensure the function haa annotations for input
        # if "return" not in func_hints:
        #     raise ValueError(f"The function '{func.__name__}' is missing a return type annotation. Please add typing to your function.")

        #obtain the types parameters and types for the output
        func_input_params, a = self._get_function_details(func=func)

        # Get the input parameter types of func2
        if len(func_input_params) == 0:
            raise ValueError(f"The function '{func.__name__}' has no input parameters.")

        #get the first parameter
        func_first_param = next(iter(func_input_params.values()))

        return func_first_param
    

    def _get_function_return_type(self, func : Callable) -> object:
        '''Gets the return type of the function
        
        Parameters
        ----------
        func : Callable
            The function to determine the first parameters type for

        Returns
        -------
        object
            This the type of the return parameter
        '''

        # if not func_hints:
        #     raise ValueError(f"The function '{func.__name__}' is missing parameter type annotations. Please add typing to your function.")
        
        #obtain the function return type
        a, func_return_type = self._get_function_details(func=func)

        return func_return_type


    
    
    def calculate(self, protein: 'sparrow.Protein') -> float:
        '''Calculates the compositional comparison value for a given protein
        
        Parameters
        ----------
        protein : sparrow.Protein
            This is the protein object that the optimizer will pass in order to
            determine if the current stocastic sequence is getting closer to the
            target value for the comparison value
        
        Returns
        -------
        float
            This is the comparison value of interest you are trying to find
        '''
        #since the function property we are interested in must be definable by the 
        #sequence we must first obtain the protein sequence
        protein_sequence = protein.sequence

        #next we will compute the property of interest for the sequence
        protein_property = self.property_function(protein_sequence)

        #after we get the property we will compare against the reference value that was passed at initialization
        comparison_val = self.comparison_function(protein_property, self.reference_value)

        #return the comparison value back to the optimizer
        return comparison_val
    
