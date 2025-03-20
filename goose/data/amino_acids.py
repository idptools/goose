'''
.py file for amino acids. Contains a dictionary of amino acids and various
properties of amino acids.

Using this instead of the .pkl file because .pkl files are falling out of favor
due to their being a security risk.
'''

amino_acid_dat={
                'L' : {'FCR' : 0.0, 'NCPR' : 0.0, 'hydrophobicity' : 8.3, 'aliphatic' :1,'aromatic' : 0, 
                    'polar' : 0, 'negative' : 0, 'positive' : 0, 'pro' : 0, 'his' : 0,'gly' : 0, 'cys' : 0, 'special' : 0},
                'M' : {'FCR' : 0.0, 'NCPR' : 0.0, 'hydrophobicity' : 6.4, 'aliphatic' :1,'aromatic' : 0, 
                    'polar' : 0, 'negative' : 0, 'positive' : 0, 'pro' : 0, 'his' : 0,'gly' : 0, 'cys' : 0, 'special' : 0},
                'A' : {'FCR' : 0.0, 'NCPR' : 0.0, 'hydrophobicity' : 6.3, 'aliphatic' :1,'aromatic' : 0, 
                    'polar' : 0, 'negative' : 0, 'positive' : 0, 'pro' : 0, 'his' : 0,'gly' : 0, 'cys' : 0, 'special' : 0},
                'V' : {'FCR' : 0.0, 'NCPR' : 0.0, 'hydrophobicity' : 8.7, 'aliphatic' :1,'aromatic' : 0, 
                    'polar' : 0, 'negative' : 0, 'positive' : 0, 'pro' : 0, 'his' : 0,'gly' : 0, 'cys' : 0, 'special' : 0},
                'I' : {'FCR' : 0.0, 'NCPR' : 0.0, 'hydrophobicity' : 9.0, 'aliphatic' :1,'aromatic' : 0, 
                    'polar' : 0, 'negative' : 0, 'positive' : 0, 'pro' : 0, 'his' : 0,'gly' : 0, 'cys' : 0, 'special' : 0},
                'Q' : {'FCR' : 0.0, 'NCPR' : 0.0, 'hydrophobicity' : 1.0, 'aliphatic' : 0,'aromatic' : 0, 
                    'polar' : 1, 'negative' : 0, 'positive' : 0, 'pro' : 0, 'his' : 0,'gly' : 0, 'cys' : 0, 'special' : 0},
                'N' : {'FCR' : 0.0, 'NCPR' : 0.0, 'hydrophobicity' : 1.0, 'aliphatic' : 0,'aromatic' : 0, 
                    'polar' : 1, 'negative' : 0, 'positive' : 0, 'pro' : 0, 'his' : 0,'gly' : 0, 'cys' : 0, 'special' : 0},
                'S' : {'FCR' : 0.0, 'NCPR' : 0.0, 'hydrophobicity' : 3.7, 'aliphatic' : 0,'aromatic' : 0, 
                    'polar' : 1, 'negative' : 0, 'positive' : 0, 'pro' : 0, 'his' : 0,'gly' : 0, 'cys' : 0, 'special' : 0},
                'T' : {'FCR' : 0.0, 'NCPR' : 0.0, 'hydrophobicity' : 3.8, 'aliphatic' : 0,'aromatic' : 0, 
                    'polar' : 1, 'negative' : 0, 'positive' : 0, 'pro' : 0, 'his' : 0,'gly' : 0, 'cys' : 0, 'special' : 0},
                'D' : {'FCR' : 1.0, 'NCPR' : -1.0, 'hydrophobicity' : 1.0, 'aliphatic' : 0,'aromatic' : 0, 
                    'polar' : 0, 'negative' : 1, 'positive' : 0, 'pro' : 0, 'his' : 0,'gly' : 0, 'cys' : 0, 'special' : 0},
                'E' : {'FCR' : 1.0, 'NCPR' : -1.0, 'hydrophobicity' : 1.0, 'aliphatic' : 0,'aromatic' : 0, 
                    'polar' : 0, 'negative' : 1, 'positive' : 0, 'pro' : 0, 'his' : 0,'gly' : 0, 'cys' : 0, 'special' : 0},
                'K' : {'FCR' : 1.0, 'NCPR' : 1.0, 'hydrophobicity' : 0.6, 'aliphatic' : 0,'aromatic' : 0, 
                    'polar' : 0, 'negative' : 0, 'positive' : 1, 'pro' : 0, 'his' : 0,'gly' : 0, 'cys' : 0, 'special' : 0},
                'R' : {'FCR' : 1.0, 'NCPR' : 1.0, 'hydrophobicity' : 0.0, 'aliphatic' : 0,'aromatic' : 0, 
                    'polar' : 0, 'negative' : 0, 'positive' : 1, 'pro' : 0, 'his' : 0,'gly' : 0, 'cys' : 0, 'special' : 0},
                'F' : {'FCR' : 0.0, 'NCPR' : 0.0, 'hydrophobicity' : 7.3, 'aliphatic' : 0,'aromatic' : 1,
                    'polar' : 0, 'negative' : 0, 'positive' : 0, 'pro' : 0, 'his' : 0,'gly' : 0, 'cys' : 0, 'special' : 0},
                'W' : {'FCR' : 0.0, 'NCPR' : 0.0, 'hydrophobicity' : 3.6, 'aliphatic' : 0,'aromatic' : 1, 
                    'polar' : 0, 'negative' : 0, 'positive' : 0, 'pro' : 0, 'his' : 0,'gly' : 0, 'cys' : 0, 'special' : 0},
                'Y' : {'FCR' : 0.0, 'NCPR' : 0.0, 'hydrophobicity' : 3.2, 'aliphatic' : 0,'aromatic' : 1, 
                    'polar' : 0, 'negative' : 0, 'positive' : 0, 'pro' : 0, 'his' : 0,'gly' : 0, 'cys' : 0, 'special' : 0},
                'P' : {'FCR' : 0.0, 'NCPR' : 0.0, 'hydrophobicity' : 2.9, 'aliphatic' : 0,'aromatic' : 0, 
                    'polar' : 0, 'negative' : 0, 'positive' : 0, 'pro' : 1, 'his' : 0,'gly' : 0, 'cys' : 0, 'special' : 1},
                'H' : {'FCR' : 0.0, 'NCPR' : 0.0, 'hydrophobicity' : 1.3, 'aliphatic' : 0,'aromatic' : 0, 
                    'polar' : 0, 'negative' : 0, 'positive' : 0, 'pro' : 0, 'his' : 1,'gly' : 0, 'cys' : 0, 'special' :0},
                'G' : {'FCR' : 0.0, 'NCPR' : 0.0, 'hydrophobicity' : 4.1, 'aliphatic' : 0,'aromatic' : 0, 
                    'polar' : 0, 'negative' : 0, 'positive' : 0, 'pro' : 0, 'his' : 0,'gly' : 1, 'cys' : 0, 'special' : 1},
                'C' : {'FCR' : 0.0, 'NCPR' : 0.0, 'hydrophobicity' : 7.0, 'aliphatic' : 0,'aromatic' : 0, 
                    'polar' : 0, 'negative' : 0, 'positive' : 0, 'pro' : 0, 'his' : 0,'gly' : 0, 'cys' : 1, 'special' : 0}
    }
    
    