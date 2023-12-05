'''
backend code for generating sequence variants where it is necessary to keep the 
classes of amino acids in the sequence the same
'''


import math
import random
from random import randint

from goose.backend.amino_acids import AminoAcid


from metapredict import meta

# standard amino acids
amino_acids = ['A', 'C', 'D', 'E', 'F', 'G', 'H', 'I', 'K', 'L', 'M', 'N', 'P', 'Q', 'R', 'S', 'T', 'V', 'W', 'Y']


#function that returns a random amino acid from a specified list.
def random_amino_acid(seq_list):

    """

    Function that returns a random amino acid from
    a list of amino acids.
    

    Parameters
    -------------

    seq_list : List
        List of amino acids from which to select a random residue.


    Returns
    ---------
    
    String
        Returns a single amino acid from the selected list as a string.

    """

    return seq_list[randint(0, len(seq_list)-1)]


def generate_random_seq(length, seq_list = []):
    """

    Function that returns a random amino acid sequence
    of a specified length. Mainly for testing.
    

    Parameters
    -------------

    length : Int
        Integer value representing the length of the desired sequence


    Returns
    ---------
    
    String
        Returns an amino acid sequence of desired length

    """
    # make empty string to hold amino acids
    sequence = ""

    # if no seq list provided, provide standard amino acids
    if seq_list == []:
        seq_list = amino_acids

    # add residues to the sequence
    for i in range(0, length):
        sequence += random.choice(amino_acids)

    # return the sequence 
    return sequence



# function to figure out the max changes in 
# hydropathy if the desire is to keep classes
# of amino acids similar in generated variants


def hydro_range_constant_classes(sequence):
    """

    Function to calculate the minimum and maximum possible hydropathy values
    that a sequence can have if a variant is to be made that keeps the
    number of amino acids in each class the same. The amino acid classes
    are as follows
        Aromatics : F W Y
        Polar : Q N S T
        Hydrophobic : I V L A M
        Positive charge : K R
        Negative Charge : D E
        Special cases (not in any class)
            Glycine (G)
            Proline (P)
            Histidine (H)    
            Cystine (C)

    Parameters
    -------------
    sequence : String
        Amino acid sequence as a string


    Returns
    ---------
    List
        Returns a list where the first value is the lowest possible hydropathy
        for the variant and the second value is the highest possible hydropathy

    """    

    # set min and max values
    min_possible = 0
    max_possible = 0

    # min possible hydro values from each class
    AA_hydro_min = {"A": 6.3,
    "R": 0.0,
    "N": 1.0,
    "D": 1.0,
    "C": 7.0,
    "Q": 1.0,
    "E": 1.0,
    "G": 4.1,
    "H": 1.3,
    "I": 6.3,
    "L": 6.3,
    "K": 0.0,
    "M": 6.3,
    "F": 3.2,
    "P": 2.9,
    "S": 1.0,
    "T": 1.0,
    "W": 3.2,
    "Y": 3.2,
    "V": 6.3
    }

    # max possible hydro values for each class
    AA_hydro_max = {"A": 9.0,
    "R": 0.6,
    "N": 3.8,
    "D": 1.0,
    "C": 7.0,
    "Q": 3.8,
    "E": 1.0,
    "G": 4.1,
    "H": 1.3,
    "I": 9.0,
    "L": 9.0,
    "K": 0.6,
    "M": 9.0,
    "F": 7.3,
    "P": 2.9,
    "S": 3.8,
    "T": 3.8,
    "W": 7.3,
    "Y": 7.3,
    "V": 9.0
    }

    # for each of the amino acids, figure out the value to add for each class
    for aa in sequence:
        min_possible += AA_hydro_min[aa]
        max_possible += AA_hydro_max[aa]

    # return the max and min values dividued by sequence length in a list 
    return[round(min_possible/len(sequence), 6), round(max_possible/len(sequence), 6)]





