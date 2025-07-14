"""
This was originally build before SPARROW. Now it is
asically a wrapper over SPARROW. 
"""

import random

from sparrow import Protein as SparrowProtein
from goose.goose_exceptions import GooseInputError


class Protein:
    """
    Class that holds the properties of an amino acid sequence. This class makes
    use of functionality encoded in sparrow
    """
    def __init__(self, seq):
        """
        Constructor takes an input amino acid sequence and sets up a 

        """
        
        # make the sequence all uppercase
        self.seq = seq.upper()
        
        # Check if all amino acids in the sequence are standard
        standard_amino_acids = ['A', 'C', 'D', 'E', 'F', 'G', 'H', 'I', 'K', 'L', 'M', 'N', 'P', 'Q', 'R', 'S', 'T', 'V', 'W', 'Y']
        if not set(seq).issubset(standard_amino_acids):
            raise GooseInputError("Invalid amino acid detected. Make sure sequence only has canonical amino acids.")

        # build a new sparrow Protein object which can then be called as needed
        self._sparrow_protein = SparrowProtein(seq)


    # ......................................................................
    #
    @property
    def sequence(self):
        return self.seq

    
    # ......................................................................
    #    
    @property
    def length(self):
        return len(self.seq)

    
    # ......................................................................
    #    
    @property
    def fractions(self):
        """
        Returns a dictionary with the fractions o
        """

        # get the raw fractions
        fraction_amino_acids = self._sparrow_protein.amino_acid_fractions

        # Use dictionary comprehension for faster rounding
        return {amino_acid: round(frac, 5) for amino_acid, frac in fraction_amino_acids.items()}

    
    # ......................................................................
    #
    @property
    def FCR(self):
        return self._sparrow_protein.FCR

    
    # ......................................................................
    #
    @property
    def NCPR(self):
        return self._sparrow_protein.NCPR


    # ......................................................................
    #
    @property
    def sigma(self):
        # can't be calculated in SPARROW any more.
        if self.FCR ==  0:
            return 0
        else:
            return round(((self.NCPR**2) / self.FCR), 6)

    

    # ......................................................................
    #    
    @property
    def SCD(self):
        return self._sparrow_protein.SCD

    @property
    def kappa(self):
        return self._sparrow_protein.kappa

    @property
    def hydro(self):
        return self._sparrow_protein.hydrophobicity

    @property
    def hydropathy(self):
        return self._sparrow_protein.hydrophobicity

    @property
    def percent_polar(self):
        return self._sparrow_protein.compute_residue_fractions(['Q','N','S','T'])

    @property
    def percent_aliphatic(self):
        return self._sparrow_protein.compute_residue_fractions(['I','L','V','A','M'])

    @property
    def percent_aromatic(self):
        return self._sparrow_protein.compute_residue_fractions(['Y','W','F'])

    @property
    def properties(self):
        properties_dict = {}
        
        properties_dict['length']     = self.length
        properties_dict['FCR']        = round(self.FCR,5)
        properties_dict['NCPR']       = round(self.NCPR,5)
        properties_dict['hydropathy'] = round(self.hydropathy,5)
        properties_dict['sigma']      = round(self.sigma,5)
        properties_dict['SCD']        = round(self.SCD,5)
        properties_dict['kappa']      = round(self.kappa,5)
        
        return properties_dict

    def calc_basic_properties(self):
        properties_dict = {}

        properties_dict['length']     = self.length
        properties_dict['FCR']        = round(self.FCR,5)
        properties_dict['NCPR']       = round(self.NCPR,5)
        properties_dict['hydropathy'] = round(self.hydropathy,5)

        return properties_dict

    def auto_name(self):
        standard_amino_acids = ['A', 'C', 'D', 'E', 'F', 'G', 'H', 'I', 'K', 'L', 'M', 'N', 'P', 'Q', 'R', 'S', 'T', 'V', 'W', 'Y']
        # generates an autoname based on the sequence properties
        return '>' + ''.join(f"{random.choice(list(standard_amino_acids))}{random.randint(0, 9)}" for _ in range(5))


