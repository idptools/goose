'''
Where the amino acid classes are defined. This will make it easy
to change the classes in the future if needed.
'''

aa_classes = {
    'aliphatic': ['A', 'I', 'L', 'M', 'V'],
    'polar': ['S', 'T', 'N', 'Q'],
    'positive': ['K', 'R'],
    'negative': ['D', 'E'],
    'aromatic': ['F', 'W', 'Y'],
    'glycine': ['G'],
    'proline': ['P'],
    'cysteine': ['C'],
    'histidine': ['H'],
}

'''
The amino acids in the class of the given amino acid
'''
aa_classes_by_aa = {
    'A': ['I', 'L', 'M', 'V'],
    'I': ['A', 'L', 'M', 'V'],
    'L': ['A', 'I', 'M', 'V'],
    'M': ['A', 'I', 'L', 'V'],
    'V': ['A', 'I', 'L', 'M'],
    'S': ['T', 'N', 'Q'],
    'T': ['S', 'N', 'Q'],
    'N': ['S', 'T', 'Q'],
    'Q': ['S', 'T', 'N'],
    'K': ['R'],
    'R': ['K'],
    'D': ['E'],
    'E': ['D'],
    'F': ['W', 'Y'],
    'W': ['F', 'Y'],
    'Y': ['F', 'W'],
    'G': ['G'],
    'P': ['P'],
    'C': ['C'],
    'H': ['H']
}

aas_to_class_nums = {
    'A': 0,
    'I': 0,
    'L': 0,
    'M': 0,
    'V': 0,
    'S': 1,
    'T': 1,
    'N': 1,
    'Q': 1,
    'K': 2,
    'R': 2,
    'D': 3,
    'E': 3,
    'F': 4,
    'W': 4,
    'Y': 4,
    'G': 5,
    'P': 6,
    'C': 7,
    'H': 8
}

# Define amino acid classes (using indices)
aa_class_indices = {
    # Hydrophobic: A, I, L, M, V
    0: [7, 9, 10, 17],   # A -> A, I, L, M, V
    7: [0, 9, 10, 17],   # I -> A, I, L, M, V  
    9: [0, 7, 10, 17],   # L -> A, I, L, M, V
    10: [0, 7, 9, 17],  # M -> A, I, L, M, V
    17: [0, 7, 9, 10],  # V -> A, I, L, M, V
    
    # Aromatic: F, W, Y
    4: [18, 19],         # F -> F, W, Y
    18: [4, 19],        # W -> F, W, Y
    19: [4, 18],        # Y -> F, W, Y
    
    # Polar: S, T, N, Q
    15: [16, 11, 13],   # S -> S, T, N, Q
    16: [15, 11, 13],   # T -> S, T, N, Q
    11: [15, 16, 13],   # N -> S, T, N, Q
    13: [15, 16, 11],   # Q -> S, T, N, Q
    
    # Positive: D, E
    2: [3],              # D -> D, E
    3: [3],              # E -> D, E
    
    # Negative: K, R
    8: [14],             # K -> K, R
    14: [8],            # R -> K, R
    
    # Immutable residues
    5: [5],                 # G -> G (glycine)
    1: [1],                 # C -> C (cysteine)
    12: [12],               # P -> P (proline)
    6: [6],                 # H -> H (histidine)
}


max_class_hydro = {
    'A': 9,
    'I': 9,
    'L': 9,
    'M': 9,
    'V': 9,
    'S': 3.8,
    'T': 3.8,
    'N': 3.8,
    'Q': 3.8,
    'K': 0.6,
    'R': 0.6,
    'D': 1.0,
    'E': 1.0,
    'F': 7.3,
    'W': 7.3,
    'Y': 7.3,
    'G': 4.1,
    'P': 2.9,
    'C': 7.0,
    'H': 1.3
    }

min_class_hydro = {
    'A': 6.3,
    'I': 6.3,
    'L': 6.3,
    'M': 6.3,
    'V': 6.3,
    'S': 1,
    'T': 1,
    'N': 1,
    'Q': 1,
    'K': 0,
    'R': 0,
    'D': 1.0,
    'E': 1.0,
    'F': 3.2,
    'W': 3.2,
    'Y': 3.2,
    'G': 4.1,
    'P': 2.9,
    'C': 7.0,
    'H': 1.3
    }
