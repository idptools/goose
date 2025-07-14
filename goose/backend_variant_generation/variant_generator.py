"""
VariantGenerator class that encapsulates all variant generation functions from variant_sequence_generation.py.
Allows setting common parameters once and reusing them across different methods.
This works as follows:
    1. The variant is attempted to be generated using the specified method.
    2. If the generated sequence does not meet the disorder criteria, it is discarded.
    3. The process is repeated for a specified number of attempts until a valid sequence is found or all attempts are exhausted.
    4. If a valid sequence is found, it is checked for disorder.
    5. If the sequence is disordered, it is returned. If not, the function tries again.
"""
import random
from sparrow.protein import Protein
from goose.backend import parameters
from goose.backend_sequence_generation.sequence_generation import by_properties
from goose.backend_variant_generation.helper_functions import check_variant_disorder_vectorized, hydropathy_range
from goose.backend_variant_generation import variant_sequence_generation as vsg



class VariantGenerator:
    """
    A class to generate protein sequence variants with consistent parameters.
    
    This class encapsulates all variant generation functions and allows you to set
    common parameters once, which will be used across all generation methods.
    """
    
    def __init__(self,
                 num_attempts: int = 200,
                 strict_disorder: bool = False,
                 disorder_cutoff: float = 0.5,
                 metapredict_version: int = 3,
                 hydropathy_tolerance: float = parameters.MAXIMUM_HYDRO_ERROR,
                 kappa_tolerance: float = parameters.MAXIMUM_KAPPA_ERROR):
        """
        Initialize the VariantGenerator with default parameters.
        
        Parameters
        ----------
        num_attempts : int
            Number of attempts to generate a variant.
            Default is 5
        strict_disorder : bool
            If True, all residues must be above the disorder cutoff to be considered disordered.
            Default is False
        disorder_cutoff : float
            Cutoff for disorder. Above this value is considered disordered.
            Default is 0.5
        metapredict_version : int
            Version of MetaPredict to use (1, 2, or 3)
            Default is 3
        hydropathy_tolerance : float
            Acceptable difference between achieved and target hydropathy.
            Default is parameters.MAXIMUM_HYDRO_ERROR
        kappa_tolerance : float
            Acceptable difference between achieved and target kappa.
            Default is parameters.MAXIMUM_KAPPA_ERROR
        """
        self.num_attempts = num_attempts
        self.strict_disorder = strict_disorder
        self.disorder_cutoff = disorder_cutoff
        self.metapredict_version = metapredict_version
        self.hydropathy_tolerance = hydropathy_tolerance
        self.kappa_tolerance = kappa_tolerance
    
    def _check_disorder_and_return(self, input_sequence: str, variant_sequences) -> str:
        """
        Helper method to check disorder and return the first valid sequence.
        
        Parameters
        ----------
        input_sequence : str
            The original input sequence
        variant_sequences : str or list
            The variant sequence(s) to check
            
        Returns
        -------
        str or None
            The first valid disordered sequence, or None if none found
        """
        disordered_sequence = check_variant_disorder_vectorized(
            original_sequence=input_sequence,
            sequences=variant_sequences,
            strict_disorder=self.strict_disorder,
            disorder_cutoff=self.disorder_cutoff,
            metapredict_version=self.metapredict_version,
            return_best_sequence=False
        )
        
        if disordered_sequence is None:
            return None
        
        # Make sure we return a single sequence
        if isinstance(disordered_sequence, list):
            disordered_sequence = disordered_sequence[0]
        
        return disordered_sequence

    
    def shuffle_specific_regions(self,
                                  input_sequence: str,
                                  shuffle_regions: list) -> str:
        """
        Generate a variant sequence by shuffling specified regions of the input sequence.
        
        Parameters
        ----------
        input_sequence : str
            The input sequence to generate a variant from.
        shuffle_regions : list
            List of tuples specifying regions to shuffle. Each tuple should contain
            (start_index, end_index) for the region to shuffle.
            
        Returns
        -------
        str
            A generated variant sequence with specified regions shuffled.
        """
        # validate that shuffle_regions are valid
        if not all(isinstance(region, tuple) and len(region) == 2 for region in shuffle_regions):
            raise ValueError("shuffle_regions must be a list of tuples with (start_index, end_index) pairs.")
        # validate that start_index < end_index for each region
        for region in shuffle_regions:
            if region[0] >= region[1]:
                raise ValueError("Each region's start_index must be less than end_index.")
        # validate that start_index and end_index are within the bounds of the input_sequence
        sequence_length = len(input_sequence)
        for region in shuffle_regions:
            if region[0] < 0 or region[1] > sequence_length:
                raise ValueError("Each region's start_index and end_index must be within the bounds of the input_sequence.")

        for _ in range(self.num_attempts):
            shuffle_vars = []
            for _ in range(10):
                shuffle_vars.append(vsg.shuffle_specific_regions_sequence(input_sequence, shuffle_regions))
            
            disordered_sequence = self._check_disorder_and_return(input_sequence, shuffle_vars)
            if disordered_sequence is not None:
                return disordered_sequence
        
        return None
    
    def shuffle_except_specific_regions(self,
                                        input_sequence: str,
                                        excluded_regions: list) -> str:
        """
        Generate a variant sequence by shuffling all regions except specified ones.
        
        Parameters
        ----------
        input_sequence : str
            The input sequence to generate a variant from.
        excluded_regions : list
            List of tuples specifying regions to exclude from shuffling.
            Each tuple should contain (start_index, end_index) for the region to exclude.
            
        Returns
        -------
        str
            A generated variant sequence with specified regions excluded from shuffling.
        """
        # validate that excluded_regions are valid
        if not all(isinstance(region, tuple) and len(region) == 2 for region in excluded_regions):
            raise ValueError("excluded_regions must be a list of tuples with (start_index, end_index) pairs.")
        # validate that start_index < end_index for each region
        for region in excluded_regions:
            if region[0] >= region[1]:
                raise ValueError("Each region's start_index must be less than end_index.")
        # validate that start_index and end_index are within the bounds of the input_sequence
        sequence_length = len(input_sequence)
        for region in excluded_regions:
            if region[0] < 0 or region[1] > sequence_length:
                raise ValueError("Each region's start_index and end_index must be within the bounds of the input_sequence.")


        for _ in range(self.num_attempts):
            shuffle_vars = []
            for _ in range(10):
                shuffle_vars.append(vsg.shuffle_except_specific_regions_sequence(input_sequence, excluded_regions))
            
            disordered_sequence = self._check_disorder_and_return(input_sequence, shuffle_vars)
            if disordered_sequence is not None:
                return disordered_sequence
        
        return None
    
    def shuffle_specific_residues(self,
                                   input_sequence: str,
                                   target_residues: list) -> str:
        """
        Generate a variant sequence by shuffling residues within a specified target residues.
        
        Parameters
        ----------
        input_sequence : str
            The input sequence to generate a variant from.
        target_residues : list
            List of residues to shuffle within the sequence.
            
        Returns
        -------
        str
            A generated variant sequence with specified residues shuffled.
        """
        # dict of classes that are possible to choose
        classdict={'charged':['D', 'E', 'K', 'R'], 'polar':['Q', 'N', 'S', 'T'], 'aromatic':
        ['F', 'W', 'Y'], 'aliphatic': ['I', 'V', 'L', 'A', 'M'], 'negative':['D', 'E'], 'positive':['K', 'R']}
        
        # if input is a class as a list, pull the first element
        if isinstance(target_residues, list):
            if target_residues[0] in classdict:
                target_residues=target_residues[0] 
        
        # if string is the input, see if its the classdict or amino acids. 
        if isinstance(target_residues, str):
            if target_residues in classdict:
                # if target_aas is a class, get the list of amino acids in that class
                target_residues = classdict[target_residues]
            else:
                # if target_aas is a single amino acid, convert to list
                target_residues = list(target_residues)

        # make sure that there is at least 1 target residue in the input sequence
        if not any(residue in input_sequence for residue in target_residues):
            raise ValueError("No target residues found in the input sequence. Please ensure that at least one of the target residues is present in the input sequence.")

        # get the identity of the target residues in the sequence. 
        target_residues_in_sequence = [residue for residue in target_residues if residue in input_sequence]
        # make sure not only one residue getting shuffled, this does nothing.
        if len(set(target_residues_in_sequence)) < 2:
            raise ValueError("Only one target residue found in the input sequence. This will not change the sequence. Make sure at least two targets in input sequence.")

        for _ in range(self.num_attempts):
            shuff_seqs = []
            for _ in range(10):
                shuff_seqs.append(vsg.shuffle_specific_residues_sequence(input_sequence, target_residues))
            
            disordered_sequence = self._check_disorder_and_return(input_sequence, shuff_seqs)
            if disordered_sequence is not None:
                return disordered_sequence
        
        return None

    def shuffle_except_specific_residues(self,
                                   input_sequence: str,
                                   excluded_residues: list) -> str:
        """
        Generate a variant sequence by shuffling residues not specified in excluded_residues.
        
        Parameters
        ----------
        input_sequence : str
            The input sequence to generate a variant from.
        excluded_residues : list
            List of residues to exclude from shuffling.
            
        Returns
        -------
        str
            A generated variant sequence with specified residues excluded from shuffling.
        """
        # dict of classes that are possible to choose
        classdict={'charged':['D', 'E', 'K', 'R'], 'polar':['Q', 'N', 'S', 'T'], 'aromatic':
        ['F', 'W', 'Y'], 'aliphatic': ['I', 'V', 'L', 'A', 'M'], 'negative':['D', 'E'], 'positive':['K', 'R']}

        # if input is a class as a list, pull the first element
        if isinstance(excluded_residues, list):
            if excluded_residues[0] in classdict:
                excluded_residues=excluded_residues[0] 

        if isinstance(excluded_residues, str):
            if excluded_residues in classdict:
                # if target_aas is a class, get the list of amino acids in that class
                excluded_residues = classdict[excluded_residues]
            else:
                # if target_aas is a single amino acid, convert to list
                excluded_residues = list(excluded_residues) 

        # get all residues to be excluded from the input sequence
        excluded_residues_in_sequence = [residue for residue in excluded_residues if residue in input_sequence]
        # make sure not all residues are excluded, this does nothing.
        if len(set(excluded_residues_in_sequence)) == len(set(input_sequence)):
            raise ValueError("All residues in the input sequence are excluded. This will not change the sequence. Make sure at least one residue is not in excluded_residues.")

        for _ in range(self.num_attempts):
            shuff_seqs = []
            for _ in range(10):
                shuff_seqs.append(vsg.shuffle_except_specific_residues_sequence(
                    input_sequence, 
                    excluded_residues=excluded_residues
                ))
            
            disordered_sequence = self._check_disorder_and_return(input_sequence, shuff_seqs)
            if disordered_sequence is not None:
                return disordered_sequence
        
        return None


    def weighted_shuffle_specific_residues(self,
                                     input_sequence: str,
                                    target_residues: list,
                                    shuffle_weight: float) -> str:
        """
        Generate a variant sequence by shuffling residues with a specified weight.
        
        Parameters
        ----------
        input_sequence : str
            The input sequence to generate a variant from.
        target_residues : list
            List of residues to shuffle within the sequence.
        shuffle_weight : float
            Weight for the shuffling process.
            
        Returns
        -------
        str
            A generated variant sequence with specified residues shuffled.
        """


        # dict of classes that are possible to choose
        classdict={'charged':['D', 'E', 'K', 'R'], 'polar':['Q', 'N', 'S', 'T'], 'aromatic':
        ['F', 'W', 'Y'], 'aliphatic': ['I', 'V', 'L', 'A', 'M'], 'negative':['D', 'E'], 'positive':['K', 'R']}

        # if input is a class as a list, pull the first element
        if isinstance(target_residues, list):
            if target_residues[0] in classdict:
                target_residues=target_residues[0] 

        if isinstance(target_residues, str):
            if target_residues in classdict:
                # if target_aas is a class, get the list of amino acids in that class
                target_residues = classdict[target_residues]
            else:
                # if target_aas is a single amino acid, convert to list
                target_residues = list(target_residues)


        # make sure that there is at least 1 target residue in the input sequence
        if not any(residue in input_sequence for residue in target_residues):
            raise ValueError("No target residues found in the input sequence. Please ensure that at least one of the target residues is present in the input sequence.")

        # get the identity of the target residues in the sequence. 
        target_residues_in_sequence = [residue for residue in target_residues if residue in input_sequence]
        # make sure not only one residue getting shuffled, this does nothing.
        if len(set(target_residues_in_sequence)) < 2:
            raise ValueError("Only one target residue found in the input sequence. This will not change the sequence. Make sure at least two targets in input sequence.")


        for _ in range(self.num_attempts):
            variant_sequence = vsg.weighted_shuffle_specific_residues_sequence(
                input_sequence,
                target_aas=target_residues,
                shuffle_weight=shuffle_weight
            )
            if variant_sequence is None:
                continue
            
            disordered_sequence = self._check_disorder_and_return(input_sequence, variant_sequence)
            if disordered_sequence is not None:
                return disordered_sequence
        
        return None

    def targeted_reposition_specific_residues(self,
                                        input_sequence: str,
                                        target_residues: list) -> str:
        """
        Generate a variant sequence by repositioning specified target residues.
        
        Parameters
        ----------
        input_sequence : str
            The input sequence to generate a variant from.
        target_residues : list
            List of residues to reposition within the sequence.
            
        Returns
        -------
        str
            A generated variant sequence with specified residues repositioned.
        """
        # Define amino acid classes
        amino_acid_classes = {
            'charged': ['D', 'E', 'K', 'R'], 
            'polar': ['Q', 'N', 'S', 'T'], 
            'aromatic': ['F', 'W', 'Y'], 
            'aliphatic': ['I', 'V', 'L', 'A', 'M'], 
            'negative': ['D', 'E'], 
            'positive': ['K', 'R']
        }

        # if input is a class as a list, pull the first element
        if isinstance(target_residues, list):
            if target_residues[0] in amino_acid_classes:
                target_residues = target_residues[0]

        # Process target_residues input
        if isinstance(target_residues, str):
            if target_residues in amino_acid_classes:
                target_residues = amino_acid_classes[target_residues]
            else:
                # Convert single amino acid string to list
                target_residues = list(target_residues.upper())   

        # make sure that there is at least 1 target residue in the input sequence
        if not any(residue in input_sequence for residue in target_residues):
            raise ValueError("No target residues found in the input sequence. Please ensure that at least one of the target residues is present in the input sequence.")
 

        for _ in range(self.num_attempts):
            variant_sequence = vsg.targeted_reposition_specific_residues_sequence(
                input_sequence,
                target_residues=target_residues
            )
            if variant_sequence is None:
                continue
            
            disordered_sequence = self._check_disorder_and_return(input_sequence, variant_sequence)
            if disordered_sequence is not None:
                return disordered_sequence
        
        return None


    def change_residue_asymmetry(self,
                            input_sequence: str,
                            target_residues: list,
                            num_changes: int = 1,
                            increase_or_decrease: str = 'increase') -> str:
        """
        Generate a variant sequence by introducing asymmetry in specified residues.
        
        Parameters
        ----------
        input_sequence : str
            The input sequence to generate a variant from.
        target_residues : list
            List of residues to introduce asymmetry in.
        num_changes : int
            Number of residues to change to make the sequence asymmetric.
            Default is 1
        increase_or_decrease : str
            Whether to increase or decrease the asymmetry.
            Default is 'increase'
            
        Returns
        -------
        str
            A generated variant sequence with specified residues made asymmetric.
        """
        # dict of classes that are possible to choose
        classdict={'charged':['D', 'E', 'K', 'R'], 'polar':['Q', 'N', 'S', 'T'], 'aromatic':
        ['F', 'W', 'Y'], 'aliphatic': ['I', 'V', 'L', 'A', 'M'], 'negative':['D', 'E'], 'positive':['K', 'R']}

        if isinstance(target_residues, list):
            if target_residues[0] in classdict:
                target_residues = classdict[target_residues[0]]

        if isinstance(target_residues, str):
            if target_residues in classdict:
                # if target_aas is a class, get the list of amino acids in that class
                target_residues = classdict[target_residues]
            else:
                # if target_aas is a single amino acid, convert to list
                target_residues = list(target_residues)

        # make sure that there is at least 1 target residue in the input sequence
        if not any(residue in input_sequence for residue in target_residues):
            raise ValueError("No target residues found in the input sequence. Please ensure that at least one of the target residues is present in the input sequence.")

        if isinstance(target_residues, str):
            target_residues = list(target_residues)
        
        for _ in range(self.num_attempts):
            variant_sequence = vsg.change_residue_asymmetry_sequence(
                input_sequence,
                target_residues=target_residues,
                num_changes=num_changes,
                increase_or_decrease=increase_or_decrease
            )
            if variant_sequence is None:
                continue
            
            disordered_sequence = self._check_disorder_and_return(input_sequence, variant_sequence)
            if disordered_sequence is not None:
                return disordered_sequence
        
        return None
    
    
    def constant_properties(self,
                       input_sequence: str,
                       exclude_residues: list = None) -> str:
        """
        Generate a new variant where the sequence is constrained by hydropathy, 
        NCPR, FCR, and kappa values.
        
        Parameters
        ----------
        input_sequence : str
            The input sequence to generate a variant from.
        exclude_residues : list
            List of residues to exclude from the variant generation.
            If None, no residues are excluded.
            
        Returns
        -------
        str
            A generated variant sequence that meets the disorder criteria.
        """
        protein = Protein(input_sequence)
        original_hydropathy = protein.hydrophobicity
        original_kappa = protein.kappa
        original_FCR = protein.FCR
        original_NCPR = protein.NCPR
        length = len(input_sequence)
        
        for _ in range(self.num_attempts):
            seq = by_properties(
                length, 
                fcr=original_FCR, 
                ncpr=original_NCPR, 
                hydropathy=original_hydropathy, 
                kappa=original_kappa, 
                exclude_residues=exclude_residues,
                num_attempts=5000,
                hydropathy_tolerance=self.hydropathy_tolerance,
                kappa_tolerance=self.kappa_tolerance,
                metapredict_version=self.metapredict_version,
                return_all_sequences=True,
                check_sequence_disorder=False,
                batch_size=200
            )
            if seq is None:
                continue
            
            disordered_sequence = self._check_disorder_and_return(input_sequence, seq)
            if disordered_sequence is not None:
                return disordered_sequence
        
        return None
    

    def constant_residues_and_properties(self,
                                   input_sequence: str,
                                   constant_residues: list) -> str:
        """
        Generate a variant sequence by keeping specified residues constant.
        Hydropathy, kappa, ncpr and fcr are also held constant.
        
        Parameters
        ----------
        input_sequence : str
            The input sequence to generate a variant from.
        constant_residues : list
            List of residues to keep constant in the variant sequence.
            
        Returns
        -------
        str
            A generated variant sequence that meets the disorder criteria.
        """
        # make sure that there is at least 1 target residue in the input sequence
        if not any(residue in input_sequence for residue in constant_residues):
            raise ValueError("No constant residues found in the input sequence. Please ensure that at least one of the constant residues is present in the input sequence.")


        for _ in range(self.num_attempts):
            variant_sequence = vsg.constant_residues_and_properties_sequence(
                input_sequence,
                constant_residues=constant_residues,
                hydropathy_tolerance=self.hydropathy_tolerance
            )
            if variant_sequence is None:
                continue
            
            disordered_sequence = self._check_disorder_and_return(input_sequence, variant_sequence)
            if disordered_sequence is not None:
                return disordered_sequence
        
        return None
    

    def constant_properties_and_class(self, input_sequence: str) -> str:
        """
        Generate a new variant sequence by constant class mutation.
        
        Parameters
        ----------
        input_sequence : str
            The input sequence to generate a variant from.
            
        Returns
        -------
        str
            A generated variant sequence that meets the disorder criteria.
        """
        
        for _ in range(self.num_attempts):
            variant_sequence = vsg.constant_properties_and_class_sequence(
                input_sequence,
                kappa_tolerance=self.kappa_tolerance,
                hydropathy_tolerance=self.hydropathy_tolerance
            )
            if variant_sequence is None:
                continue
            
            disordered_sequence = self._check_disorder_and_return(input_sequence, variant_sequence)
            if disordered_sequence is not None:
                return disordered_sequence
        
        return None

    def constant_properties_and_class_by_order(self, input_sequence: str) -> str:
        """
        Generate a constant class variant of the input sequence.
        
        Parameters
        ----------
        input_sequence : str
            The input sequence to generate a variant from.
            
        Returns
        -------
        str
            A generated variant sequence that meets the disorder criteria.
        """
        for _ in range(self.num_attempts):
            variant_sequence = vsg.constant_properties_and_class_by_order_sequence(input_sequence)
            if variant_sequence is None:
                continue
            
            disordered_sequence = self._check_disorder_and_return(input_sequence, variant_sequence)
            if disordered_sequence is not None:
                return disordered_sequence
        
        return None


    def change_hydropathy_constant_class(self,
                                   input_sequence: str,
                                   target_hydropathy: float) -> str:
        """
        Generate a hydropathy class variant of the input sequence.
        
        Parameters
        ----------
        input_sequence : str
            The input sequence to generate a variant from.
        target_hydropathy : float
            Target mean hydropathy value to achieve.
            
        Returns
        -------
        str
            A generated variant sequence that meets the disorder criteria.
        """
        # determine the possible hydropathy range
        prot_object = Protein(input_sequence)
        min_hydropathy, max_hydropathy = hydropathy_range(
            len(input_sequence),
            prot_object.FCR,
            prot_object.NCPR
        )

        # check if target hydropathy is within the range
        if target_hydropathy < min_hydropathy or target_hydropathy > max_hydropathy:
            raise ValueError(f"Target hydropathy {target_hydropathy} is out of bounds based on required FCR and NCPR values for input sequence."
                             f"Valid range is [{min_hydropathy}, {max_hydropathy}].")

        for _ in range(self.num_attempts):
            variant_sequence = vsg.change_hydropathy_constant_class_sequence(
                input_sequence,
                target_hydropathy=target_hydropathy,
                hydropathy_tolerance=self.hydropathy_tolerance
            )
            if variant_sequence is None:
                continue
            
            disordered_sequence = self._check_disorder_and_return(input_sequence, variant_sequence)
            if disordered_sequence is not None:
                return disordered_sequence
        
        return None
    
    def change_fcr_minimize_class_changes(self,
                            input_sequence: str,
                            target_FCR: float) -> str:
        """
        Generate a FCR class variant of the input sequence.
        
        Parameters
        ----------
        input_sequence : str
            The input sequence to generate a variant from.
        target_FCR : float
            Target FCR value to achieve.
            
        Returns
        -------
        str
            A generated variant sequence that meets the disorder criteria.
        """
        for _ in range(self.num_attempts):
            variant_sequence = vsg.change_fcr_minimize_class_changes_sequence(
                input_sequence,
                target_FCR=target_FCR,
                hydropathy_tolerance=self.hydropathy_tolerance,
                kappa_tolerance=self.kappa_tolerance
            )
            if variant_sequence is None:
                continue
            
            disordered_sequence = self._check_disorder_and_return(input_sequence, variant_sequence)
            if disordered_sequence is not None:
                return disordered_sequence
        
        return None
    
    def change_ncpr_constant_class(self,
                             input_sequence: str,
                             target_NCPR: float) -> str:
        """
        Generate a NCPR class variant of the input sequence.
        
        Parameters
        ----------
        input_sequence : str
            The input sequence to generate a variant from.
        target_NCPR : float
            Target NCPR value to achieve.
            
        Returns
        -------
        str
            A generated variant sequence that meets the disorder criteria.
        """
        # make sure NCPR is not greater than FCR
        protein = Protein(input_sequence)
        if target_NCPR > protein.FCR:
            raise ValueError("Target NCPR cannot be greater than FCR. Please provide a valid target NCPR value.")
        
        for _ in range(self.num_attempts):
            variant_sequence = vsg.change_ncpr_constant_class_sequence(
                input_sequence,
                target_NCPR=target_NCPR,
                hydropathy_tolerance=self.hydropathy_tolerance,
                kappa_tolerance=self.kappa_tolerance
            )
            if variant_sequence is None:
                continue
            
            disordered_sequence = self._check_disorder_and_return(input_sequence, variant_sequence)
            if disordered_sequence is not None:
                return disordered_sequence
        
        return None
    
    def change_kappa(self,
                        input_sequence: str,
                        target_kappa: float) -> str:
        """
        Generate a kappa class variant of the input sequence.
        
        Parameters
        ----------
        input_sequence : str
            The input sequence to generate a variant from.
        target_kappa : float
            Target kappa value to achieve.
            
        Returns
        -------
        str
            A generated variant sequence that meets the kappa criteria.
        """
        # make sure input sequence has FCR > 0.
        protein = Protein(input_sequence)
        if protein.FCR <= 0:
            raise ValueError("Input sequence must have charged residues (FCR > 0) to make kappa variant.")
        # make sure FCR != NCPR
        if protein.FCR == protein.NCPR:
            raise ValueError("Input sequence must have at least one negatively charged and one positively charged residue to make kappa variant.")
        # make a copy of the input sequence so we can go back to it if we need to
        original_sequence = input_sequence
        # iterate.
        for cur_iter in range(self.num_attempts):
            variant_sequence = vsg.change_kappa_sequence(
                input_sequence,
                target_kappa=target_kappa,
                kappa_tolerance=self.kappa_tolerance
            )
            if variant_sequence is None:
                continue
            
            disordered_sequence = self._check_disorder_and_return(input_sequence, variant_sequence)
            if disordered_sequence is not None:
                return disordered_sequence
            else:
                # some times sequences get hung on not being disordered enough and the function gets kappa
                # before it changes the sequence significantly. This can get around that. 
                # only shuffle if cur_iter % 2 == 0
                if cur_iter % 2 == 0:
                    input_sequence = list(input_sequence)
                    random.shuffle(input_sequence)
                    input_sequence = ''.join(input_sequence)
                else:
                    # go back to the original sequence
                    input_sequence = original_sequence
        
        return None

    
    def change_properties_minimze_differences(self,
                           input_sequence: str,
                           target_hydropathy: float = None,
                           target_kappa: float = None,
                           target_FCR: float = None,
                           target_NCPR: float = None) -> str:
        """
        Generate a minimal variant of the input sequence.
        
        Parameters
        ----------
        input_sequence : str
            The input sequence to generate a variant from.
        target_hydropathy : float
            Target mean hydropathy value to achieve.
            If None, hydropathy optimization is not performed.
        target_kappa : float
            Target kappa value to achieve.
            If None, kappa optimization is not performed.
        target_FCR : float
            Target FCR value to achieve.
            If None, FCR optimization is not performed.
        target_NCPR : float
            Target NCPR value to achieve.
            If None, NCPR optimization is not performed.
            
        Returns
        -------
        str
            A generated variant sequence that meets the disorder criteria.
        """
        for _ in range(self.num_attempts):
            variant_sequence = vsg.change_properties_minimize_differences_sequence(
                input_sequence,
                target_hydropathy=target_hydropathy,
                target_kappa=target_kappa,
                target_FCR=target_FCR,
                target_NCPR=target_NCPR,
                hydropathy_tolerance=self.hydropathy_tolerance,
                kappa_tolerance=self.kappa_tolerance
            )
            if variant_sequence is None:
                continue
            
            disordered_sequence = self._check_disorder_and_return(input_sequence, variant_sequence)
            if disordered_sequence is not None:
                return disordered_sequence
        
        return None


    def change_any_properties(self,
                                  input_sequence: str,
                                  target_FCR: float=None,
                                  target_NCPR: float=None,
                                  target_kappa: float=None,
                                  target_hydropathy: float=None) -> str:
        """
        Generate a variant sequence by constraining FCR, NCPR, kappa, and hydropathy.
        
        Parameters
        ----------
        input_sequence : str
            The input sequence to generate a variant from.
        target_FCR : float
            Target FCR value to achieve.
        target_NCPR : float
            Target NCPR value to achieve.
        target_kappa : float
            Target kappa value to achieve.
        target_hydropathy : float
            Target hydropathy value to achieve.
            
        Returns
        -------
        str
            A generated variant sequence that meets the disorder criteria.
        """
        for _ in range(self.num_attempts):
            variant_sequence = vsg.change_any_properties_sequence(
                input_sequence,
                target_FCR=target_FCR,
                target_NCPR=target_NCPR,
                target_kappa=target_kappa,
                target_hydropathy=target_hydropathy,
                hydropathy_tolerance=self.hydropathy_tolerance,
                kappa_tolerance=self.kappa_tolerance
            )
            if variant_sequence is None:
                continue
            
            disordered_sequence = self._check_disorder_and_return(input_sequence, variant_sequence)
            if disordered_sequence is not None:
                return disordered_sequence
        
        return None


    def change_dimensions(self,
                          input_sequence: str,
                          increase_or_decrease: str,
                          rg_or_re: str,
                          num_dim_attempts: int=5,
                          allowed_error = parameters.MAXIMUM_RG_RE_ERROR,
                          reduce_pos_charged: bool = False,
                          exclude_aas=None
                          ) -> str:
        """
        Generate a variant sequence by constraining radius of gyration (Rg) and end-to-end distance (Re).

        Parameters
        ----------
        sequence : str
            The input amino acid sequence to modify
        increase_or_decrease : {'increase', 'decrease'}
            Whether to increase or decrease the Rg
        rg_or_re : {'rg', 're'}
            Whether to optimize for radius of gyration (Rg) or radius of elongation (Re)
        allowed_error : float, optional
            If specified, the maximum allowed error for the generated Rg or Re compared to the original value
        reduce_pos_charged : bool, default=False
            Whether to reduce positively charged residues (K, R) in the sequence.
        exclude_aas : list, optional
            A list of amino acids to exclude from the optimization process.

        Returns
        -------
        str or list of str
        The generated variant(s) of the input sequence
        """ 
        if isinstance(rg_or_re, str):
            rg_or_re = rg_or_re.lower()
        else:
            raise ValueError("rg_or_re must be a string, either 'rg' or 're'.")
        
        if rg_or_re not in ['rg', 're']:
            raise ValueError("rg_or_re must be either 'rg' or 're'.")

        # iterate
        for _ in range(self.num_attempts):
            variant_sequence = vsg.change_dimensions_sequence(
                input_sequence,
                increase_or_decrease=increase_or_decrease,
                rg_or_re=rg_or_re,
                num_attempts=num_dim_attempts,
                allowed_error=allowed_error,
                reduce_pos_charged=reduce_pos_charged,
                exclude_aas=exclude_aas
            )
            if variant_sequence is None:
                continue
            
            disordered_sequence = self._check_disorder_and_return(input_sequence, variant_sequence)
            if disordered_sequence is not None:
                return disordered_sequence
        
        return None