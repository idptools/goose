"""
code that holds the properties of amino acids 
"""

class AminoAcid:
    """
    Class that holds properties for individual amino acids.
    """
    
    def __init__(self, amino_acid):
        #make the sequence all uppercase
        self.amino_acid = amino_acid.upper()
        #check value is an amino acid
        if self.amino_acid not in AminoAcid.standard_amino_acids:
            raise Exception ("Invalid amino acid detected. Make sure value is a canonical amino acid.")

        self.AA = self.amino_acid
        self.hydropathy = AminoAcid.hydro(amino_acid)
        self.charge = AminoAcid.charge_value(amino_acid)
        self.aromatic = AminoAcid.aromatic_check(amino_acid)
        self.polar = AminoAcid.polar_check(amino_acid)
        self.aliphatic = AminoAcid.aliphatic_check(amino_acid)   
        self.AA_class = AminoAcid.return_AA_class(amino_acid)     
    
    #standard amino acids
    standard_amino_acids = ['A', 'C', 'D', 'E', 'F', 'G', 'H', 'I', 'K', 'L', 'M', 'N', 'P', 'Q', 'R', 'S', 'T', 'V', 'W', 'Y']

    def hydro(amino_acid):
        #function that returns the hydropathy of an amino acid


        # KYTE-DOOLITTLE SCALES
        # References
        # A simple method for displaying the hydropathic character of a protein. 
        # Kyte J, Doolittle RF. J Mol Biol. 1982 May 5;157(1):105-32.
        # Why are "natively unfolded" proteins unstructured under physiological conditions?
        # Valdimir N. Uversky, Joel R. Gillespie, and Anthony L. Frink
        # Protines: Structure, function, and genetics 41:415-427 (2000)
        # Main hydrophobicity scale

        AA_hydro = {"A": 6.3,
        "R": 0.0,
        "N": 1.0,
        "D": 1.0,
        "C": 7.0,
        "Q": 1.0,
        "E": 1.0,
        "G": 4.1,
        "H": 1.3,
        "I": 9.0,
        "L": 8.3,
        "K": 0.6,
        "M": 6.4,
        "F": 7.3,
        "P": 2.9,
        "S": 3.7,
        "T": 3.8,
        "W": 3.6,
        "Y": 3.2,
        "V": 8.7
        }
        return AA_hydro[amino_acid]

    def charge_value(amino_acid):
        #function that returns the charge of an amino acid
        if amino_acid == "D" or amino_acid == "E":
            charge = -1
        elif amino_acid == "K" or amino_acid == "R":
            charge = 1
        else:
            charge = 0
        return charge

    def aromatic_check(amino_acid):
        #function to see if an amino acid is aromatic
        aromatics = ["F", "W", "Y"]
        if amino_acid in aromatics:
            return True
        else:
            return False
            
    def polar_check(amino_acid):
        #function to see if an amino acid is a polar amino acid
        polar_aas = ["Q", "N", "S", "T"]
        if amino_acid in polar_aas:
            return True
        else:
            return False

    def aliphatic_check(amino_acid):
        # function to see if an amino acid is amphipathic
        aliphatic_aas = ["I", "V", "L", "A", "M"]
        if amino_acid in aliphatic_aas:
            return True
        else:
            return False

    def return_AA_class(amino_acid):
        # returns which class the amino acid is in
        aromatics = ['F', 'W', 'Y']
        polar = ['Q', 'N', 'S', 'T']
        hydrophobics = ['I', 'V', 'L', 'A', 'M']
        positive = ['K', 'R']
        negative = ['D', 'E']

        if amino_acid in aromatics:
            return 'aromatic'
        elif amino_acid in polar:
            return 'polar'
        elif amino_acid in hydrophobics:
            return 'hydrophobic'
        elif amino_acid in positive:
            return 'positive'
        elif amino_acid in negative:
            return 'negative'
        else:
            return amino_acid



