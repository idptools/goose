'''
code to verify that variants maintain expected properties. 
'''
from goose import goose_exceptions
from sparrow.protein import Protein
from goose.data.defined_aa_classes import aa_classes, aas_to_class_nums
from goose.backend import parameters
from goose.backend_variant_generation.helper_functions import check_variant_disorder_vectorized
from goose.backend_property_optimization.optimize_dimensions import predict_re, predict_rg

def verify_same_number_by_class(original_seq, variant_seq):
    """
    Verify that the number of amino acids in each class remains the same
    between the original sequence and the variant sequence.
    """
    for cur_class in aa_classes:
        original_count = sum([1 for a in original_seq if a in aa_classes[cur_class]])
        variant_count =  sum([1 for a in variant_seq if a in aa_classes[cur_class]])
        if original_count != variant_count:
            return False
    return True

def verify_same_order_by_class(original_seq, variant_seq):
    """
    verify that the order of amino acids by class is the same between 2 sequences
    """
    original_class_order = [aas_to_class_nums[a] for a in original_seq]
    variant_class_order = [aas_to_class_nums[a] for a in variant_seq]
    return original_class_order == variant_class_order

def verify_constant_properties(original_seq, variant_seq,
                               kappa_tolerance=parameters.MAXIMUM_KAPPA_ERROR,
                               hydropathy_tolerance=parameters.MAXIMUM_HYDRO_ERROR,
                               FCR_tolerance = 0.001):
    """
    Verify that the variant sequence maintains the same properties as the original sequence.
    This includes checking the number of amino acids in each class and the order of classes.
    """
    original_protein = Protein(original_seq)
    variant_protein = Protein(variant_seq)

    # Verify kappa and hydropathy values
    if abs(original_protein.kappa - variant_protein.kappa) > kappa_tolerance:
        return False
    if abs(original_protein.hydrophobicity - variant_protein.hydrophobicity) > hydropathy_tolerance:
        return False
    # verify hydrophobicity
    if abs(original_protein.hydrophobicity - variant_protein.hydrophobicity) > hydropathy_tolerance:
        return False
    # verify FCR
    if abs(original_protein.FCR - variant_protein.FCR) > FCR_tolerance: 
        return False
    # verify NCPR. NCPR may vary by 1 due to prioritizing FCR.
    if abs(original_protein.NCPR - variant_protein.NCPR) > (1/len(original_seq))+1e-6:
        return False
    return True

def verify_same_residue_count(original_seq, variant_seq):
    """
    Verify that the number of each amino acid in the original sequence
    is the same as in the variant sequence.
    """
    original_count = {a: original_seq.count(a) for a in set(original_seq)}
    variant_count = {a: variant_seq.count(a) for a in set(variant_seq)}
    
    for aa in original_count:
        if original_count[aa] != variant_count.get(aa, 0):
            return False
    return True

def verify_constant_hydropathy(original_seq, variant_seq, hydropathy_tolerance=parameters.MAXIMUM_HYDRO_ERROR):
    """
    Verify that the hydropathy of the original sequence is the same as the variant sequence.
    """
    original_protein = Protein(original_seq)
    variant_protein = Protein(variant_seq)

    return abs(round(original_protein.hydrophobicity,8) - round(variant_protein.hydrophobicity,8)) < hydropathy_tolerance

def verify_constant_kappa(original_seq, variant_seq, kappa_tolerance=parameters.MAXIMUM_KAPPA_ERROR):
    """
    Verify that the kappa value of the original sequence is the same as the variant sequence.
    """
    original_protein = Protein(original_seq)
    variant_protein = Protein(variant_seq)

    return abs(round(original_protein.kappa,8) - round(variant_protein.kappa,8)) < kappa_tolerance

def verify_constant_FCR(original_seq, variant_seq, FCR_tolerance=0.001):
    """
    Verify that the FCR value of the original sequence is the same as the variant sequence.
    """
    original_protein = Protein(original_seq)
    variant_protein = Protein(variant_seq)

    return abs(round(original_protein.FCR,8) - round(variant_protein.FCR,8)) < FCR_tolerance

def verify_constant_NCPR(original_seq, variant_seq):
    """
    Verify that the NCPR value of the original sequence is the same as the variant sequence.
    NCPR may vary by 1 due to prioritizing FCR.
    """
    original_protein = Protein(original_seq)
    variant_protein = Protein(variant_seq)

    return abs(round(original_protein.NCPR,8) - round(variant_protein.NCPR,8)) <= (1/len(original_seq)) + 1e-6

def verify_target_hydropathy(variant_seq, hydropathy_target, hydropathy_tolerance=parameters.MAXIMUM_HYDRO_ERROR):
    """
    Verify that the hydropathy of the variant sequence is within the specified tolerance of the target hydropathy.
    """
    variant_protein = Protein(variant_seq)
    return abs(round(variant_protein.hydrophobicity,8) - round(hydropathy_target,8)) <= hydropathy_tolerance

def verify_target_kappa(variant_seq, kappa_target, kappa_tolerance=parameters.MAXIMUM_KAPPA_ERROR):
    """
    Verify that the kappa value of the variant sequence is within the specified tolerance of the target kappa.
    """
    variant_protein = Protein(variant_seq)
    return abs(round(variant_protein.kappa,8) - round(kappa_target,8)) <= kappa_tolerance

def verify_target_FCR(variant_seq, FCR_target):
    """
    Verify that the FCR value of the variant sequence is within the specified tolerance of the target FCR.
    """
    FCR_tolerance = 1/len(variant_seq) + 0.0001 # Tolerance based on sequence length
    variant_protein = Protein(variant_seq)
    return abs(round(variant_protein.FCR,8) - round(FCR_target,8)) <= FCR_tolerance

def verify_target_NCPR(variant_seq, NCPR_target):
    """
    Verify that the NCPR value of the variant sequence is within the specified tolerance of the target NCPR.
    NCPR may vary by 1 due to prioritizing FCR.
    """
    variant_protein = Protein(variant_seq)
    return abs(round(variant_protein.NCPR,8) - round(NCPR_target,8)) <= (1/len(variant_seq)) + 0.0001

def verify_constant_region(original_seq, variant_seq, region):
    """
    Verify that the constant region of the original sequence is the same as the variant sequence.
    The constant region is defined by the region parameter, which is a tuple of start and end indices.
    """
    start, end = region
    return original_seq[start:end] == variant_seq[start:end]

def verify_region_changed(original_seq, variant_seq, region):
    """
    Verify that the specified region of the original sequence is different from the variant sequence.
    The region is defined by the region parameter, which is a tuple of start and end indices.
    """
    start, end = region
    # if there is only 1 amino acid, even if the region was shuffled it would not change. Thus,
    # we return True for that to avoid false negatives.
    if len(set(original_seq[start:end]))==1 and len(set(variant_seq[start:end])) == 1:
        return True
    return original_seq[start:end] != variant_seq[start:end]

def verify_constant_residue_positions(original_seq, variant_seq, residues):
    """
    Verify that the specified residue positions in the original sequence are the same as the variant sequence.
    The residue positions are defined by the residues parameter, which is a list of indices.
    """
    for res in residues:
        # get the indices for the residue in original_sequence
        original_indices = [i for i, a in enumerate(original_seq) if a == res]
        # get the indices for the residue in variant_sequence
        variant_indices = [i for i, a in enumerate(variant_seq) if a == res]
        
        # if the indices are not the same, return False
        if original_indices != variant_indices:
            return False
    # if all residues are in the same positions, return True
    return True

def verify_disorder(original_seq, variant_seq, strict_disorder=False):
    """
    Verify that the disorder vector of the original sequence is the same as the variant sequence.
    This is done using a vectorized approach for efficiency.
    """
    if check_variant_disorder_vectorized(original_seq, variant_seq, strict_disorder=strict_disorder) == None:
        return False
    return True

def verify_changed_iwd(original_seq, variant_seq, increase_or_decrease, residues):
    """
    Verify that the IWD (Intrinsic Disorder) of the variant sequence has changed as expected.
    The increase_or_decrease parameter should be 'increase' or 'decrease'.
    """
    original_protein = Protein(original_seq)
    variant_protein = Protein(variant_seq)

    classdict={'charged':['D', 'E', 'K', 'R'], 'polar':['Q', 'N', 'S', 'T'], 'aromatic':
        ['F', 'W', 'Y'], 'aliphatic': ['I', 'V', 'L', 'A', 'M'], 'negative':['D', 'E'], 'positive':['K', 'R']}
    if isinstance(residues, list):
        if residues[0] in classdict:
            residues = classdict[residues[0]]
    if isinstance(residues, str):
        if residues in classdict:
            # if residues is a class, get the list of amino acids in that class
            residues = classdict[residues]
        else:
            residues=list(residues)

    if increase_or_decrease == 'increase':
        return variant_protein.compute_iwd(residues) > original_protein.compute_iwd(residues)
    elif increase_or_decrease == 'decrease':
        return variant_protein.compute_iwd(residues) < original_protein.compute_iwd(residues)
    else:
        raise goose_exceptions.GooseException("increase_or_decrease must be 'increase' or 'decrease'.")
    
def verify_dimensions(original_seq, variant_seq, rg_or_re, increase_or_decrease):
    """
    Verify that the dimensions (Rg or Re) of the variant sequence have changed as expected.
    The rg_or_re parameter should be 'Rg' or 'Re', and increase_or_decrease should be 'increase' or 'decrease'.
    """
    if rg_or_re == 'rg':
        predict_function = predict_rg
    elif rg_or_re == 're':
        predict_function = predict_re
    else:
        raise goose_exceptions.GooseException("rg_or_re must be 'rg' or 're'.")
    
    original_dimension = predict_function(original_seq)
    variant_dimension = predict_function(variant_seq)

    if increase_or_decrease == 'increase':
        return variant_dimension > original_dimension
    elif increase_or_decrease == 'decrease':
        return variant_dimension < original_dimension
    else:
        raise goose_exceptions.GooseException("increase_or_decrease must be 'increase' or 'decrease'.")

def verify_same_length(original_seq, variant_seq):
    """
    Verify that the original sequence and the variant sequence have the same length.
    """
    return len(original_seq) == len(variant_seq)