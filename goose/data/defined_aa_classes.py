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