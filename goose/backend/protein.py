"""
code that uses the AminoAcid class to calculate protein parameter
values based off of the amino acids
"""

import random
import math

from goose.backend.amino_acids import AminoAcid
from sparrow import Protein as SP


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
        
        #check sequence
        for i in seq:
            if i in AminoAcid.standard_amino_acids:
                continue
            else:
                raise Exception ("Invalid amino acid detected. Make sure sequence only has canonical amino acids.")

        # build a new sparrow Protein object which can then be called as needed
        self._sparrow_protein = SP(seq)


    # ......................................................................
    #
    @property
    def sequence(self):
        return self._sparrow_protein.sequence

    
    # ......................................................................
    #    
    @property
    def length(self):
        return len(self._sparrow_protein)

    
    # ......................................................................
    #    
    @property
    def fractions(self):
        """
        Returns a dictionary with the fractions o
        """

        # get the raw fractions
        fraction_amino_acids = self._sparrow_protein.amino_acid_fractions

        for amino_acid in fraction_amino_acids:
            fraction_amino_acids[amino_acid] = round(fraction_amino_acids[amino_acid],5)

        return fraction_amino_acids

    
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
        # generates an autoname based on the sequence properties
        random_name = '>'
        amino_acids = ['A', 'C', 'D', 'E', 'F', 'G', 'H', 'I', 'K', 'L', 'M', 'N', 'P', 'Q', 'R', 'S', 'T', 'V', 'W', 'Y']
        for i in range(0, 5):
            random_name += amino_acids[random.randint(0, len(amino_acids)-1)]
            random_name += str(random.randint(0, 9))
        return random_name



            
