"""
VariantGenerator class that encapsulates all variant generation functions.
Allows setting common parameters once and reusing them across different methods.
"""
from sparrow.protein import Protein
from goose.backend_sequence_generation.sequence_generation_vectorized import generate_seq_by_props
from goose.backend_variant_generation.helper_functions import check_variant_disorder_vectorized
from goose.backend_variant_generation import variant_generation_functions as vgf
from goose.backend import parameters


class VariantGenerator:
    """
    A class to generate protein sequence variants with consistent parameters.
    
    This class encapsulates all variant generation functions and allows you to set
    common parameters once, which will be used across all generation methods.
    """
    
    def __init__(self,
                 num_attempts: int = 5,
                 strict_disorder: bool = False,
                 disorder_cutoff: float = 0.5,
                 metapredict_version: int = 3,
                 hydropathy_tolerance: float = parameters.HYDRO_ERROR,
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
            Default is parameters.HYDRO_ERROR
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
    
    def gen_constant_class_variant(self, input_sequence: str) -> str:
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
            variant_sequence = vgf.gen_constant_class_variant(input_sequence)
            if variant_sequence is None:
                continue
            
            disordered_sequence = self._check_disorder_and_return(input_sequence, variant_sequence)
            if disordered_sequence is not None:
                return disordered_sequence
        
        return None
    
    def gen_minimal_variant(self,
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
            variant_sequence = vgf.generate_minimal_variant(
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
    
    def gen_region_shuffle_variant(self,
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
        for _ in range(self.num_attempts):
            shuffle_vars = []
            for _ in range(10):
                shuffle_vars.append(vgf.generate_region_shuffle_variant(input_sequence, shuffle_regions))
            
            disordered_sequence = self._check_disorder_and_return(input_sequence, shuffle_vars)
            if disordered_sequence is not None:
                return disordered_sequence
        
        return None
    
    def gen_excluded_shuffle_variant(self,
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
        for _ in range(self.num_attempts):
            shuff_seqs = []
            for _ in range(10):
                shuff_seqs.append(vgf.generate_excluded_shuffle_variant(
                    input_sequence, 
                    excluded_residues=excluded_residues
                ))
            
            disordered_sequence = self._check_disorder_and_return(input_sequence, shuff_seqs)
            if disordered_sequence is not None:
                return disordered_sequence
        
        return None
    
    def gen_targeted_shuffle_variant(self,
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
        for _ in range(self.num_attempts):
            shuff_seqs = []
            for _ in range(10):
                shuff_seqs.append(vgf.generate_targeted_shuffle_variant(input_sequence, target_residues))
            
            disordered_sequence = self._check_disorder_and_return(input_sequence, shuff_seqs)
            if disordered_sequence is not None:
                return disordered_sequence
        
        return None
    
    def gen_new_var_constant_class(self, input_sequence: str) -> str:
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
        # Use a higher number of attempts for this method (as in original)
        attempts = max(self.num_attempts, 50)
        
        for _ in range(attempts):
            variant_sequence = vgf.generate_new_seq_constant_class_variant(
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
    
    def gen_new_variant(self,
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
            seq = generate_seq_by_props(
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
                check_disorder=False,
                batch_size=200
            )
            if seq is None:
                continue
            
            disordered_sequence = self._check_disorder_and_return(input_sequence, seq)
            if disordered_sequence is not None:
                return disordered_sequence
        
        return None
    
    def gen_constant_residue_variant(self,
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
        # Use a higher number of attempts for this method (as in original)
        attempts = max(self.num_attempts, 50)
        
        for _ in range(attempts):
            variant_sequence = vgf.generate_constant_residue_variant(
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
    
    def gen_asymmetry_variant(self,
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
        # Use a higher number of attempts for this method (as in original)
        attempts = max(self.num_attempts, 50)
        
        for _ in range(attempts):
            variant_sequence = vgf.generate_asymmetry_variant(
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
    
    def gen_hydropathy_class_variant(self,
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
        for _ in range(self.num_attempts):
            variant_sequence = vgf.generate_hydro_class_variant(
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
    
    def gen_fcr_class_variant(self,
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
            variant_sequence = vgf.generate_fcr_class_variant(
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
    
    def gen_ncpr_class_variant(self,
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
        for _ in range(self.num_attempts):
            variant_sequence = vgf.generate_ncpr_class_variant(
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
    
    def gen_kappa_variant(self,
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
        for _ in range(self.num_attempts):
            variant_sequence = vgf.generate_kappa_variant(
                input_sequence,
                target_kappa=target_kappa,
                kappa_tolerance=self.kappa_tolerance
            )
            if variant_sequence is None:
                continue
            
            disordered_sequence = self._check_disorder_and_return(input_sequence, variant_sequence)
            if disordered_sequence is not None:
                return disordered_sequence
        
        return None
    
    def gen_all_props_class_variant(self,
                                  input_sequence: str,
                                  target_FCR: float,
                                  target_NCPR: float,
                                  target_kappa: float,
                                  target_hydropathy: float) -> str:
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
            variant_sequence = vgf.generate_all_props_class_var(
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

    def gen_weighted_shuffle_variant(self,
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
        for _ in range(self.num_attempts):
            variant_sequence = vgf.gen_weighted_shuffle_variant(
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

    def gen_targeted_reposition_variant(self,
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
        for _ in range(self.num_attempts):
            variant_sequence = vgf.generate_targeted_reposition_variant(
                input_sequence,
                target_residues=target_residues
            )
            if variant_sequence is None:
                continue
            
            disordered_sequence = self._check_disorder_and_return(input_sequence, variant_sequence)
            if disordered_sequence is not None:
                return disordered_sequence
        
        return None



    def gen_rg_re_variant(self,
                          input_sequence: str,
                          increase_or_decrease: str,
                          rg_or_re: str,
                          return_all: bool = False,
                          return_all_interval: float = 0.2,
                          include_original: bool = False,
                          num_dim_attempts: int=5,
                          allowed_error = None,
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
        return_all : bool, default=False
            Whether to return all generated variants or just the best one
        return_all_interval : float, default=0.2
            Interval at which to return all generated variants
        include_original : bool, default=False
            Whether to include the original sequence in the output
        num_attempts : int, default=5
            Number of attempts to generate a valid variant
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
        for _ in range(self.num_attempts):
            variant_sequence = vgf.gen_dimensions_variant(
                input_sequence,
                increase_or_decrease=increase_or_decrease,
                rg_or_re=rg_or_re,
                return_all=return_all,
                return_all_interval=return_all_interval,
                include_original=include_original,
                num_dim_attempts=num_dim_attempts,
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