"""
This method was developed by Garrett Ginell, a member of the
Holehouse Lab. 

Code written by Holehouse Lab member Garrett Ginell. This
code allows for calculating the clustering of any amino acid, 
which will be used in GOOSE to calculate approximate residue
clustering. Thanks Garrett!
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

    value_track = []

    #iterate sequence and build track 
    for i in sequence:
        if str(i) in target_residues:
            value_track.append(1)
        else:
            value_track.append(0)      
            
    #check if track matches len of sequence, if not throw error 
    if len(value_track) != len(sequence):
        raise Exception('Masking ERROR - mask legnth does not match sequence length')
    
    return value_track


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

    # find hits in mask
    hits = np.array(find_all_indicies(mask, 1))

    # if only one hit just return the original mask as there are no pairs to compare
    if len(hits) < 2:
        return mask

    # generate list of list where each sublist are tuples of all pairs of hits except self-pairing
    pair_sets = [list(map(lambda p1: (p, p1), hits[np.arange(len(hits))!=i])) for i,p in enumerate(hits)]

    # get the minimum abs value of subtracted pairs for each hit 
    min_distances = [min(abs(pair[0]-pair[1]) for pair in pair_set) for pair_set in pair_sets]

    # define new mask
    new_mask = mask.copy() 

    # generate new mask by copying old mask and removing hits with minimum pair distance above threshold
    for i, distance in enumerate(min_distances):
        if distance > island_threshold:
            new_mask[hits[i]] = 0

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

    indices = []

    # cycle over each element in input_list
    for idx, value in enumerate(input_list):

        # if an element is True (as assessed by the boolean_function)
        # then add that index to the indices list
        if boolean_function(value):
            indices.append(idx)

    return indices



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

    # get the indices where the binary mask is 1
    hit_indices = []
    all_hits = {}
    
    for idx, value in enumerate(binary_mask):
        if value == 1:
            hit_indices.append(idx)
            all_hits[idx] = 0

    # convert to an array so we can reference multiple positions with 
    # slice notation
    hit_indices = np.array(hit_indices)
            
    # now cycle over each position that was '1' and calculate the index
    # position between 
    for index, position in enumerate(hit_indices):
        resi_distances = 0
        
        # for every OTHER hit index  
        for other_position in hit_indices[np.arange(len(hit_indices)) != index]:
            resi_distances = resi_distances + 1/(np.abs(other_position - position))

        all_hits[position] = resi_distances

    if len(all_hits) > 0:
        return sum(all_hits.values())/len(all_hits)
    else:
        return 0
        

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

    binary_mask = []
    try:
        for idx, i in enumerate(input_array):
            binary_mask.append(binerize_function(i))

    except Exception as e:
        print(f'Error when applying binerize_function to element {idx} in input_array. See exception below...')
        raise(e)


    # check that this is a true binary array
    for i in list(set(binary_mask)):
        if i not in [0,1]:
            raise Exception(f'After parsing the input array with the binerize function, one or more elements generated was {i}. The binerize_function must ONLY generate 1s or 0s (ints)')

    # finally, compute IWD
    return __compute_IWD_from_binary_mask(binary_mask)



# ------------------------------------------------------------------
#
def calculate_average_bivariate_inverse_distance(mask, delimiter_function1, delimiter_function2):
    """
    Function that returns the average inverse distance BETWEEN TWO specific sets of residues in
    a sequence mask. Each delimiter function determines a set of residues, the average 
    bivariate inverse distance is calculated between the two sets of resitues. This metric can 
    be used to quantify the relitive clustering between two sets of residues in an amino acid 
    sequence. This bivariate clustering can also be thought of as 'how interdespersed are the two
    sets of residues relitive to each other?'. 
    NOTE: housetools.sequence_tools.sequence_masking is written specifcally to assit
          with creating sequence masks from Amino Acid sequences. 
    Parameters
    -------------
    mask : list
        Mask of sequence as list, which is pasted to both delimiter_functions to  
        idenify which positions will be evaluated in the average inverse distance
        calculation. This mask and each delimiter function is passed to find_all_indices.  
        
    boolean_function1 : function
        Function which takes each element from input mask and must return either True or False,
        to indicate whether index should be included in the first index group. Defines how to extract 
        group1 possitions for clustering from the mask. 
    boolean_function2 : function
        Function which takes each element from input mask and must return either True or False,
        to indicate whether index should be included in the second index group. Defines how to extract 
        group2 possitions for clustering from the mask.
        Example functions could be: 
            'lambda a: a < 0' :  
                where positions less than a value of zero in the mask are returned 
            'lambda a: a if a in ['Y','F','W'] else None' : 
                where postions of specific residues are returned 
    Returns
    --------
    float
        float of average bivariate IWD value calculated between the two sets of indecies
        determined by delimiter_function1 and delimiter_function2 
        NOTE: Intution of this parameter is nuanced in that the higher the average bivariate IWD the more 
        interdisperced the two groups of residues are relitive to each other.   
    """
    
    # dictionary of empty values for each index point for both catagories
    all_1_hits = {i:0 for i in find_all_indices(mask, delimiter_function1)} 
    hits1 = np.array([i for i in all_1_hits.keys()])
    
    all_2_hits = {i:0 for i in find_all_indices(mask, delimiter_function2)}
    hits2 = np.array([i for i in all_2_hits.keys()])
    
    
    # make sure both catigories are populated with hits 
    if len(hits1) == 0 or len(hits2) == 0:
        return 0
    
    else:
    
        # iterate through index
        for i, p in enumerate(hits1):
            resi_distances= 0
            
            # iterate through pairs for that index 
            # but do not compare to self as self could be in both catagories
            for p1 in [i for i in hits2 if i != p]:
                resi_distances += 1 / np.abs(p1-p)

            # save to all_hits
            all_1_hits[p] = resi_distances
            
        # return bivariate average IWD of pairs 
        return sum(all_1_hits.values())/(len(hits1)+len(hits2))


# ------------------------------------------------------------------
#
def calculate_average_inverse_distance_charge(mask, sequence, charge=['-','+']):
    """
    Function that returns the charge weighted average inverse distance of either the positive or
    negative residues in the sequence. For more information on the average_inverse_distance see 
    the calculate_average_inverse_distance function. The only difference is here the per residue sum 
    of the inverse distances is weighted by the absolute value of the charge at that posision. IE 
    residues with higher charge hold more more weight when evaluating how clustered that specific 
    charge residue is. This metric can be used to quantify the local clustering of either positive 
    or negitive charge residues in an amino acid sequence.
    Parameters
    -------------
    mask : list
        Linear net charge per residue (NCPR) as list calculated accross the sequence. For 
        more information on how to calculate the NCPR see SPARROW
    sequence : list
        Amino acid sequence spit into list, where element 0 is the first residue in sequence and 
        element N-1 is the Nth residue in the sequence
    
    charge : string  options=['-','+']
        Pass '-' to quantify the clustering of negitive residues.
        Pass '+' to quantify the clustering of positive residues.
    Returns
    --------
    float 
        Returns average charge weighted IWD value for the sequence based on the passed,
        NCPR mask, sequence, and which prefered charge to calculate the charge-weighted average IWD    
    """
    
    # dictionary to select delimiter function based on passed charge prefference 
    delimiter_function = {'-':lambda a: a in ['D','E'], '+':lambda a: a in ['R','K']}[charge]

    # dictionary to test map prefered charge to sign of NCPR     
    charge_tester = {'-':-1.0,'+':1.0}
    
    # dictionary of empty values for each index point 
    # (extracts all hit residues based on charge definition)
    all_hits = {i:0 for i in find_all_indices(sequence, delimiter_function)} 
    hits = np.array([i for i in all_hits.keys()])
    
    # iterate through index (here p is an index in the sequence)
    for i, p in enumerate(hits):
        resi_distances= 0
        # iterate through pairs for that index 
        for p1 in hits[np.arange(len(hits))!=i]:
            resi_distances += 1 / np.abs(p1-p)

        # multiply the residue charge by the resi_distances
        # if NCPR of residue in the mask is opposite of the charge being evaluated set charge to 0 
        if np.sign(mask[p]) == charge_tester[charge]:
            all_hits[p] = abs(mask[p])*resi_distances
        else:
            all_hits[p] = 0*resi_distances
        
    if len(hits) > 0:
        return sum(all_hits.values())/len(hits)
    else:
        return 0
    

# ------------------------------------------------------------------
#
def calculate_average_bivariate_inverse_distance_charge(mask, sequence):
    """
    Function returns the charge-weight average inverse distance BETWEEN positive and negative
    residues in an amino acid sequence. This function is similar to the 
    calculate_average_bivariate_inverse_distance in its calculation except it is specific to charge
    residues and is weigthed by the difference in NCPR between charge pairs. 
    The per residue pair bivariate inverse distance is weighted by the difference in charge 
    between the two oppositly charged residues of intrest: 
    (charge2 - charge1) / (distance2 - distance1)
    This metric can be used to quantify the relitive clustering between charge residues or 
    the interdespersion of charge residues in a sequence relitive to each other. 
    Parameters
    -------------
    mask : list
        Linear net charge per residue (NCPR) as list calculated accross the sequence. For 
        more information on how to calculate the NCPR see SPARROW
    sequence : list
        Amino acid sequence spit into list, where element 0 is the first residue in sequence and 
        element N-1 is the Nth residue in the sequence
  
    Returns
    --------
    float 
        Returns average bivariate charge-weighted IWD value for the sequence based on the passed,
        NCPR mask and sequence. NOTE the intuition of bivariate clustering the parameter is nuanced 
        in that the higher the value the more interspersed the the positive and negative residues are. 
    """
    
    # dictionary to select delimiter function based on passed charge prefference 
    delimiter_function = {'-':lambda a: a in ['D','E'], '+':lambda a: a in ['R','K']}

    #  dictionary to test map prefered charge to sign of NCPR     
    charge_tester = {-1.0:['D','E'], 1.0:['R','K'], 0:[None]}
    
    # adjust charge mask such that charged residues with a charge opposite to that of the 
    # their defined charge are automatically set to zero. This is important to do so that during 
    # the charge difference calculation for weighting the difference of between charges is always
    # calculated a the difference between a - & + charge 
    l_mask = [c if sequence[i] in charge_tester[np.sign(c)] else 0 for i,c in enumerate(mask)]
    
    # dictionary of empty values for each index point for both catagories
    all_neg_hits = {i:0 for i in find_all_indices(sequence, delimiter_function['-'])} 
    neg_hits = np.array([i for i in all_neg_hits.keys()])
    
    all_pos_hits = {i:0 for i in find_all_indices(sequence, delimiter_function['+'])}
    pos_hits = np.array([i for i in all_pos_hits.keys()])
    
    
    # make sure both catigories are populated with hits 
    if len(neg_hits) == 0 or len(neg_hits) == 0:
        return 0
    
    else:
    
        # iterate through index
        for i, p in enumerate(neg_hits):
            resi_distances= 0
            # iterate through pairs for that index 
            for p1 in pos_hits:
                # sum over inverse distance
                # resi_distances += 1 / np.abs(p1-p)
                
                # sum over charge difference over linear distance
                resi_distances += np.abs(l_mask[p1]-l_mask[p]) / np.abs(p1-p)

            # save to all_hits
            all_neg_hits[p] = resi_distances
        
        
        # return bivariate average IWD of pairs 
        return sum(all_neg_hits.values())/(len(neg_hits)+len(pos_hits))


