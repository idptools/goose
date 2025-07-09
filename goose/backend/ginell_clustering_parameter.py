"""
This method was developed by Garrett Ginell, a form member of the
Holehouse Lab. 

"""

import numpy as np

###############################
## Bionary Masking Functions ##
###############################

# ------------------------------------------------------------------ 
#
def mask_sequence(sequence, target_residues):
    """
    Function converts sequence to a binary mask (list where each
    position is either 1 or 0) based on the residues passed in
    target_residues.
    Parameters
    --------------
    sequence : str
        Input amino acid sequence
    target_residues : list 
        A list of residues which will be used to to mask the sequence
        into 1s and zeros
    
    Returns
    ---------------
    Mask : list
        List where every position is either a 0 or a 1, depending on if 
        the original residue was in the target_residue list or not.
    """

    # Convert sequence to numpy array for vectorized operations
    seq_array = np.array(list(sequence))
    target_set = set(target_residues)  # Use set for O(1) lookup
    
    # Vectorized operation using numpy's isin function
    binary_mask = np.isin(seq_array, list(target_set)).astype(int)
    
    return binary_mask.tolist()


##------------------------------------------------------------------ 
#
def filter_out_islands(mask, island_threshold=0):
    """
    Takes a binary mask and filters out 1s based on parameter input removing 
    1s in the mask that are above to island threshold 
    Parameters
    --------------
    mask : list 
        mask of 1s and 0s 
    island_threshold : int 
        define maximum allowed distances betweens hits in a sequence
                      _ 
        IE. if mask = 100000000000000000000101000100001000000000000000101000000000000
        and island_threshold=20 
                      _
           new mask = 000000000000000000000101000100001000000000000000101000000000000
    Returns
    ---------------
    Mask : list
        returns new vector mask 
    """

    # find hits in mask - use vectorized numpy operation
    mask_array = np.array(mask)
    hits = np.where(mask_array == 1)[0]

    # if only one hit just return the original mask as there are no pairs to compare
    if len(hits) < 2:
        return mask

    # Vectorized distance calculation using broadcasting
    # Create a matrix of distances between all pairs of hits
    distance_matrix = np.abs(hits[:, np.newaxis] - hits[np.newaxis, :])
    
    # Set diagonal to infinity to ignore self-distances
    np.fill_diagonal(distance_matrix, np.inf)
    
    # Find minimum distance for each hit (excluding self)
    min_distances = np.min(distance_matrix, axis=1)

    # Create new mask and remove hits with distances above threshold
    new_mask = mask.copy()
    hits_to_remove = hits[min_distances > island_threshold]
    
    for hit_idx in hits_to_remove:
        new_mask[hit_idx] = 0

    return new_mask


##------------------------------------------------------------------ 
#
def extract_fragments(mask, max_separation=1):
    """
    Converts a binary mask of 1s and 0s to fragments based on the max_separation
    (largest gap between two fragments).
    For example, if max_separation = 1 then:
    
       In:  [0, 0, 1, 1, 1, 0, 1, 1, 0, 0, 0, 1, 0, 0, 1, 1, 0, 0, 1, 0, 1, 1, 0, 1]
        
       Out: ['111011', '1', '11', '11', '101101']
    
    Parameters
    --------------
    mask : list 
       Binary mask of 0s and 1s
    max_separation : int 
        Define maximum number of 0s betweens 1s before a new fragment is identified 
    Returns
    ---------------
    list
        Returns list of strings where each element is a fragment extracted from the
        input mask
    """

    mask = [str(s) for s in mask]
    spliter= ''.join(['0'] * (max_separation+1))
    frag_list = ''.join(mask).split(spliter)
    
    return [f.strip('0') for f in frag_list if f ]


###################################
##  Non-Binary Masking Functions ##
###################################

##------------------------------------------------------------------ 
#
def MASK_n_closest_nearest_neighbors(mask, max_separation=1, max_distance=4):
    """
    Takes in mask and converts this residues to into none binary mask 
    based on the relitive posision for the hit residues to each other 
    Parameters
    --------------
    mask : list 
        Mask of 1s and 0s 
    max_separation : int 
        Define maximum number of 0s betweens 1s that is allowed when counting the
        nearest neighbors to the target residue (default = 1) 
    max_distance : int 
        define maximum linear distance that is summed when counting the
        nearest neighbors to the target residue ie the window size (default = 4) 
    Returns
    ---------------
    mask : list
        returns new vector mask where target residues are assigned the sum of there
        nearest neibors
    """
    w_half = max_distance    
    
    # extract fragments per aliphatic residue (IE split mask by mask seperator, iterate by aliphatic and get sum within fragment 
    # for centered window around aliphatic residue
    frags = extract_fragments(mask,max_separation=max_separation)
    per_ali_frags=[] 
    
    for frag in frags:
        l = len(frag)
        # handle residues where nearest neighbor fragment is less than the window size
        if l <= w_half:
            for r in frag:
                if r == '1':
                    per_ali_frags.append(frag)
        else: 
            # parse frags larger than window size
            l_mask = frag
            out_mask=[]
            for i,r in enumerate(l_mask):
                if r == '1':
                    if w_half >= i and i+w_half <= l:
                        per_ali_frags.append(l_mask[:i+w_half+1])
                    elif  i+w_half > l:
                        per_ali_frags.append(l_mask[i-w_half:])
                    elif w_half <= i <= l-w_half: 
                        per_ali_frags.append(l_mask[i-w_half:i+w_half+1])
                    else:
                        raise Exception('Parsing ERROR')      
    
    # get sum of neighbors in fragments
    nearest_neighbor_sums_mask = []
    ali_count=0
    for i in mask:
        if i == 1:
            nearest_neighbor_sums_mask.append(sum([int(r) for r in per_ali_frags[ali_count]]))
            ali_count+=1
        else:
            nearest_neighbor_sums_mask.append(i)
    
    return nearest_neighbor_sums_mask



# ------------------------------------------------------------------
def find_all_indices(input_list, boolean_function):
    """
    Function that returns the list indices where values in the input_list
    are True as assessed by the boolean function.
    Parameters
    -------------
    input_list : list
        List of elements, where each element will be evaluated by the 
        boolean_function()
    boolean_function : function
        Function which takes each element from input_list and must return
        either True or False. 
    Returns
    --------
    list
        Returns a list of indices where the input_list elements 
        where True
    """

    # Convert to numpy array for vectorized operations
    input_array = np.array(input_list)
    
    # Apply boolean function vectorized using np.vectorize for complex functions
    # or direct operations for simple comparisons
    try:
        # Try direct vectorized operation first (faster for simple functions)
        if hasattr(boolean_function, '__name__') and boolean_function.__name__ == '<lambda>':
            # For lambda functions, use np.vectorize
            vectorized_func = np.vectorize(boolean_function)
            mask = vectorized_func(input_array)
        else:
            # For simple functions, try direct application
            mask = boolean_function(input_array)
    except (TypeError, ValueError):
        # Fallback to element-wise application for complex functions
        vectorized_func = np.vectorize(boolean_function)
        mask = vectorized_func(input_array)
    
    # Return indices where condition is True
    return np.where(mask)[0].tolist()



# ------------------------------------------------------------------
#
def calculate_average_inverse_distance_from_sequence(sequence, target_residues):
    """
    Function that takes an amino acid sequence and a set of target residues and
    computes the inverse-weighted distance (IWD) between residues in the 
    target_residue group.
    Parameters
    ----------------
    sequence : str
        Valid amino acid sequence
    target_residues : list
        List of amino acids to be used as the residues of interest with 
        respect to clustering. Practically, the IWD metric is calculated
        by applying a mask to the sequence where residues found inside the
        target_residues list are set to 1 and all other residues set to 0
    Returns
    ----------------
    float
        Returns the IWD value, where is a number between 0 and some 
        positive value.
    """

    # build binary mask
    binary_mask = mask_sequence(sequence, target_residues)

    # compute IWD 
    return __compute_IWD_from_binary_mask(binary_mask)



def __compute_IWD_from_binary_mask(binary_mask):
    """
    Internal function that actually computes the inverse weighted distance
    from a binary mask. This function should not be called directly, but
    intsead should be called by functions inside this module. This is the
    single place where the monovariate IWD should be calculated .
    Parameters
    ----------------
    binary_mask : list
        A list where each element is either 0 or 1. This mask is then used
        to calculate the inverse weighted distance.
    Returns
    ----------------
    float
        Returns the IWD value, where is a number between 0 and some 
        positive value.
    """

    # Convert to numpy array and find hit indices using vectorized operation
    mask_array = np.array(binary_mask)
    hit_indices = np.where(mask_array == 1)[0]
    
    # Return 0 if no hits found
    if len(hit_indices) == 0:
        return 0
    
    # Return 0 if only one hit (no distances to calculate)
    if len(hit_indices) == 1:
        return 0
            
    # Vectorized distance calculation using broadcasting
    # Create distance matrix: hit_indices[:, None] - hit_indices[None, :]
    distance_matrix = np.abs(hit_indices[:, np.newaxis] - hit_indices[np.newaxis, :])
    
    # Convert to float to handle infinity values
    distance_matrix = distance_matrix.astype(float)
    
    # Set diagonal to infinity to avoid division by zero (self-distances)
    np.fill_diagonal(distance_matrix, np.inf)
    
    # Calculate inverse distances, avoiding division by zero
    with np.errstate(divide='ignore', invalid='ignore'):
        inverse_distance_matrix = 1.0 / distance_matrix
        inverse_distance_matrix[distance_matrix == 0] = 0  # Handle any remaining zeros
        inverse_distance_matrix[np.isinf(inverse_distance_matrix)] = 0  # Handle infinities
    
    # Sum inverse distances for each hit (excluding self)
    inverse_distance_sums = np.sum(inverse_distance_matrix, axis=1)
    
    # Return average inverse distance
    return np.mean(inverse_distance_sums)
        

# ------------------------------------------------------------------
#
def calculate_average_inverse_distance_from_array(input_array, binerize_function=lambda a: a == 1):
    """
    Function that returns the average inverse weighted distance (IWD)
    calculated from an input array (e.g. a list, string, any other 
    iterable data type that can work with the syntax:
        for i in input_array:
            # do something with i
    The input list or string should report on per-position information,
    while the binerize_function takes that list or string and binerizes
    it, converting it into a list of 0s or 1s.
    For example, the input_array could be the linear disorder score, and
    the binerize_function a function that defines a threshold above which
    disorder is converted to 1 and below disorder converted to 0.
    This function will check that after the binerize_function has been 
    applied, a list of the same length consisting of only 0s and 1s is 
    generated.
    NOTE: For working directly with amino acid sequences where you have 
          specific residues you want to cluster we recommend using the 
          calculate_average_inverse_distance_from_sequence()
    Parameters
    -------------
    input_array : iterable
        Input data which is iterated over and each element evaluated by
        the boolean_function. 
    binerize_function : function
        Function which takes each element from input_array and performs
        a logical operation to cast that element to be either 0 or 1.
       This can be passed as a lambda function, e.g.::
            'lambda a: a < 0' : 
                where positions less than a value of zero in the mask are returned 
            'lambda a: a if a in ['Y','F','W'] else None' : 
                where postions of specific residues are returned 
        or simply any named function - e.g::
            def bool_test(i):
                if i > 0.5:
                    return 1
                return 0
            calculate_average_inverse_distance_from_array(input_array, bool_test)
    Returns
    --------
    float
        Returns the inverse weighted distance (IWD) clustering value, as computed
        based on the binerized array generated by the binerize_function
    """

    # Convert input to numpy array for vectorized operations
    input_array = np.array(input_array)
    
    try:
        # Try vectorized operation first (much faster)
        vectorized_func = np.vectorize(binerize_function)
        binary_mask = vectorized_func(input_array)
    except Exception as e:
        # Fallback to element-wise application if vectorization fails
        print(f'Vectorization failed, falling back to element-wise processing. Error: {e}')
        binary_mask = []
        for idx, i in enumerate(input_array):
            try:
                binary_mask.append(binerize_function(i))
            except Exception as inner_e:
                print(f'Error when applying binerize_function to element {idx} in input_array. See exception below...')
                raise(inner_e)

    # Convert to numpy array if it isn't already
    binary_mask = np.array(binary_mask)

    # Check that this is a true binary array using vectorized operations
    unique_values = np.unique(binary_mask)
    invalid_values = unique_values[~np.isin(unique_values, [0, 1])]
    
    if len(invalid_values) > 0:
        raise Exception(f'After parsing the input array with the binerize function, one or more elements generated was {invalid_values[0]}. The binerize_function must ONLY generate 1s or 0s (ints)')

    # Convert to list for compatibility with existing function
    binary_mask = binary_mask.tolist()

    # finally, compute IWD
    return __compute_IWD_from_binary_mask(binary_mask)

