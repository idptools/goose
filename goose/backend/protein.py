"""
code that uses the AminoAcid class to calculate protein parameter
values based off of the amino acids
"""

import random
import math

from goose.backend.amino_acids import AminoAcid

class Protein:
    """
    Class that holds the properties of an amino acid sequence.
    """
    def __init__(self, seq):
        #make the sequence all uppercase
        self.seq = seq.upper()
        #check sequence
        for i in seq:
            if i in AminoAcid.standard_amino_acids:
                continue
            else:
                raise Exception ("Invalid amino acid detected. Make sure sequence only has canonical amino acids.")
                
        self.sequence = seq
        self.length = len(seq)
        self.fractions = Protein.calc_frac(seq)
        self.FCR = Protein.calc_FCR(seq)
        self.NCPR = Protein.calc_NCPR(seq)
        self.sigma = Protein.calc_sigma(seq)
        self.delta = Protein.calc_delta(seq)
        self.scd = Protein.calc_SCD(seq)
        self.hydro = Protein.calc_mean_hydro(seq)
        self.hydropathy = Protein.calc_mean_hydro(seq)
        self.properties = Protein.all_properties(seq)
        self.percent_polar = Protein.calc_percent_polar(seq)
        self.percent_aliphatic = Protein.calc_percent_aliphatic(seq)
        self.basic_properties = Protein.basic_properties(seq)

    #function that returns the fraction of each amino acid in a sequence    
    def calc_frac(seq):
        fraction_amino_acids = {
        'A': 0, 'C': 0, 'D': 0, 'E': 0, 'F': 0, 'G': 0, 'H': 0, 
        'I': 0, 'K': 0, 'L': 0,'M': 0,'N': 0, 'P': 0,
        'Q': 0,'R': 0,'S': 0,'T': 0,'V': 0,'W': 0,'Y': 0}
        for amino_acid in seq:
            fraction_amino_acids[amino_acid] += 1
        for amino_acid in fraction_amino_acids:
            fraction_amino_acids[amino_acid] = round((fraction_amino_acids[amino_acid] / len(seq)),5)
        return fraction_amino_acids

    def calc_FCR(seq):
        charged_res = 0
        for i in seq:
            if AminoAcid(i).charge != 0:
                charged_res += 1
        return round(charged_res / len(seq), 6)

    #function that returns the fraction of positive values
    def calc_NCPR(seq):
        charge = 0
        for i in seq:
            charge += AminoAcid(i).charge
        return round(charge / len(seq), 6)

    #function to calculate sigma of a sequence
    def calc_sigma(seq):
        FCR_seq = Protein.calc_FCR(seq)
        NCPR_seq = Protein.calc_NCPR(seq)
        if FCR_seq ==  0:
            Sigma = 0
        else:
            Sigma = round(((NCPR_seq**2) / FCR_seq), 6)
        return Sigma

    def delta_pre(seq, bloblen = 5):
        sigma_seq = Protein.calc_sigma(seq)
        nblobs = len(seq) - bloblen + 1
        final_delta = 0
        total_blob_sigma = 0
        if len(seq) < 5:
            print("can't calculate for length less than 6")
        else:
            arb_val = 0
            for i in range(0, nblobs):
                cur_blob = seq[arb_val: (arb_val + bloblen)]
                cur_sigma = Protein.calc_sigma(cur_blob)
                add_sigma = ((sigma_seq - cur_sigma) ** 2) / nblobs
                total_blob_sigma += add_sigma
                arb_val += 1
        return total_blob_sigma

    #function to calculate delta with bloblen 5.5
    def calc_delta(seq):
        return round(((Protein.delta_pre(seq, bloblen=5) + Protein.delta_pre(seq, bloblen=6)) /2), 6)

    # function to calculate sequence charge decoration
    def calc_SCD(sequence):
        total = 0
        for m in range(1, len(sequence)):
            m_val = m+1
            for n in range(0, m-1):
                n_val = n+1

                cur_m_res = sequence[m]
                cur_n_res = sequence[n]
                
                if cur_m_res == "D" or cur_m_res == "E":
                    cur_m_charge = -1
                elif cur_m_res == "K" or cur_m_res == "R":
                    cur_m_charge = 1
                else:
                    cur_m_charge = 0

                if cur_n_res == "D" or cur_n_res == "E":
                    cur_n_charge = -1
                elif cur_n_res == "K" or cur_n_res == "R":
                    cur_n_charge = 1
                else:
                    cur_n_charge = 0

                charge_val = cur_m_charge * cur_n_charge
                final_val = charge_val * (math.sqrt((m_val)-(n_val)))

                total += final_val

        return round(total * (1/len(sequence)), 5)


    #function to calculate the mean hydropathy of a sequence
    def calc_mean_hydro(seq):
        total_hydro = 0
        for i in seq:
            total_hydro += AminoAcid(i).hydropathy
        return round(total_hydro / len(seq), 6)

    def calc_percent_aromatic(seq):
        #function to determine percent of amino acids that are aromatic
        N=0
        for i in seq:
            if AminoAcid(i).aromatic_check:
                N+=1
            else:
                continue
        if N > 0:
            percent_aromatic = N/len(seq)
        else:
            percent_aromatic = 0
        return percent_aromatic

    def calc_percent_polar(seq):
        #function to determine percent of amino acids that are polar
        N=0
        for i in seq:
            if AminoAcid(i).polar_check:
                N+=1
            else:
                continue
        if N > 0:
            percent_polar = N/len(seq)
        else:
            percent_polar = 0
        return percent_polar

    def calc_percent_aliphatic(seq):
        N=0
        for i in seq:
            if AminoAcid(i).aliphatic_check:
                N+=1
            else:
                continue
        if N > 0:
            percent_aliphatic = N/len(seq)
        else:
            percent_aliphatic = 0
        return percent_aliphatic

    def auto_name():
        # generates an autoname based on the sequence properties
        random_name = '>'
        amino_acids = ['A', 'C', 'D', 'E', 'F', 'G', 'H', 'I', 'K', 'L', 'M', 'N', 'P', 'Q', 'R', 'S', 'T', 'V', 'W', 'Y']
        for i in range(0, 5):
            random_name += amino_acids[random.randint(0, len(amino_acids)-1)]
            random_name += str(random.randint(0, 9))
        return random_name

    def calc_all_properties(seq, return_seq = False):
        properties_dict = {}
        # properties_dict['Name'] = Protein.auto_name()
        properties_dict['length'] = len(seq)
        properties_dict['FCR'] = Protein.calc_FCR(seq)
        properties_dict['NCPR'] = Protein.calc_NCPR(seq)
        properties_dict['hydropathy'] = Protein.calc_mean_hydro(seq)
        properties_dict['sigma'] = Protein.calc_sigma(seq)
        # properties_dict['delta'] = Protein.calc_delta(seq)
        properties_dict['SCD'] = Protein.calc_SCD(seq)
        properties_dict['fractions'] = Protein.calc_frac(seq)
        properties_dict['sequence'] = seq
        return properties_dict

    def calc_basic_properties(seq):
        properties_dict = {}
        properties_dict['length'] = len(seq)
        properties_dict['FCR'] = Protein.calc_FCR(seq)
        properties_dict['NCPR'] = Protein.calc_NCPR(seq)
        properties_dict['hydropathy'] = Protein.calc_mean_hydro(seq)
        properties_dict['SCD'] = Protein.calc_SCD(seq)
            
