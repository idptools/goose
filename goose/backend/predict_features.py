'''
Integrating the sparrow predictors into GOOSE. 
This is for running predictions and then formatting them. 
This will allow easier updating of networks as they become better over time.
'''


from sparrow import Protein as pr


def get_helical_regions(sequence):
    '''
    function that returns amino acid coordinates 
    of helical regions

    Parameters
    ----------
    sequence : string
        the amino acid sequence as a string
    Returns
    -------
    helical_regions : list
        list of lists of helical regions
    '''
    helicity_scores=pr(sequence).predictor.dssp_helicity()
    helical_regions = []
    subregion=[]
    for aa in range(0, len(helicity_scores)):
        if helicity_scores[aa] == 1:
            subregion.append(aa)
        else:
            if len(subregion) > 0:
                helical_regions.append([subregion[0]+1, subregion[-1]+1])
                subregion = []
    return helical_regions



def predict_mitochondrial_localization(sequence, return_raw=False):
    '''
    function for predicting mitochondrial localization sequences.

    Parameters
    ----------
    sequence : string
        the amino acid sequence as a string

    return_raw : bool
        whether to return raws scores
        default is false.

    Returns
    -------
    mito_loc : dict
        dictionary of mitochondrial localization sequences and their locations
    '''
    seq=pr(sequence)
    mitochondrial_targeting = seq.predictor.mitochondrial_targeting_sequence()
    if return_raw==True:
        return mitochondrial_targeting
    mito_loc = {}
    if mitochondrial_targeting.count(1) > 10:
        cur_consecutive = 0
        best = 0
        mito_seq = ''
        best_mito = ''
        for i in range(0, len(sequence)):
            cur_mito = mitochondrial_targeting[i]
            cur_aa = sequence[i]
            if cur_mito == 1:
                cur_consecutive += 1
                mito_seq += cur_aa
            else:
                if cur_consecutive > best:
                    best = cur_consecutive
                    best_mito = mito_seq
                cur_consecutive = 0
                mito_seq = ''
        if best_mito != '':
            sec_loc = sequence.index(best_mito)+1
            end_loc = sec_loc+len(best_mito)
            mito_loc = {best_mito: [sec_loc, end_loc]}
    return mito_loc

def predict_NES(sequence, return_raw=False):
    '''
    function for predicting NES

    Parameters
    ----------
    sequence : string
        the amino acid sequence as a string

    return_raw : bool
        whether to return raws scores
        default is false.

    Returns
    -------
    nes_loc : dict
        dictionary of NES and their locations in the sequence
    '''
    seq=pr(sequence)
    nes_loc={}
    num_consec = 0
    potential_nes = ''
    nes_seqs = []
    nes_targeting = seq.predictor.nuclear_export_signal()
    if return_raw==True:
        return nes_targeting
    for prob in range(0, len(nes_targeting)):
        cur_prob = nes_targeting[prob]
        cur_aa = sequence[prob]
        if cur_prob > 0.5:
            num_consec += 1
            potential_nes += cur_aa
        else:
            if num_consec > 4:
                nes_seqs.append(potential_nes)
            potential_nes = ''
            num_consec = 0
    for nes in nes_seqs:
        cur_start = sequence.index(nes)+1
        cur_end = cur_start + len(nes)
        nes_loc[nes] = [cur_start, cur_end]
    return nes_loc

def predict_NLS(sequence, return_raw=False):
    '''
    function for predicting NLS

    Parameters
    ----------
    sequence : string
        the amino acid sequence as a string

    return_raw : bool
        whether to return raws scores
        default is false.

    Returns
    -------
    nes_loc : dict
        dictionary of NLS and their locations in the sequence
    '''
    seq=pr(sequence)    
    num_consec = 0
    nls_loc={}
    potential_nls = ''
    nls_seqs = []
    nls_targeting=seq.predictor.nuclear_import_signal()
    if return_raw==True:
        return nls_targeting
    for prob in range(0, len(nls_targeting)):
        cur_prob = nls_targeting[prob]
        cur_aa = sequence[prob]
        if cur_prob > 0.5:
            num_consec += 1
            potential_nls += cur_aa
        else:
            if num_consec > 4:
                nls_seqs.append(potential_nls)
            potential_nls = ''
            num_consec = 0
    for nls in nls_seqs:
        cur_start = sequence.index(nls)+1
        cur_end = cur_start + len(nls)
        nls_loc[nls] = [cur_start, cur_end]
    return nls_loc

def predict_phosphosites(sequence, return_sites=True):
    '''
    function for predicting phosphosites

    Parameters
    ----------
    sequence : string
        the amino acid sequence as a string

    return_sites : bool
        whether to return the sites or return raw scores
        default is to return sites

    Returns
    -------
    dict
        dictionary of S, T, and Y phosphosites and their locations in the sequence
    '''
    seq=pr(sequence)
    final_returned_values={}
    S=seq.predictor.serine_phosphorylation(return_sites_only=return_sites)
    T=seq.predictor.threonine_phosphorylation(return_sites_only=return_sites)
    Y=seq.predictor.tyrosine_phosphorylation(return_sites_only=return_sites)
    return {'S':S,'T':T,'Y':Y}

def predict_TAD(sequence, return_raw=False):
    '''
    function for predicting TADs

    Parameters
    ----------
    sequence : string
        the amino acid sequence as a string

    return_raw : bool
        whether to return raws scores
        default is false.
    Returns
    -------
    tad_locs : dict
        dictionary of TADs and their locations in the sequence
        
    '''
    seq=pr(sequence)
    num_consec = 0
    potential_tad = ''
    tad_seqs = []
    tad_locs = seq.predictor.transactivation_domains()
    if return_raw==True:
        return tad_locs
    for prob in range(0, len(tad_locs)):
        cur_prob = tad_locs[prob]
        cur_aa = sequence[prob]
        if cur_prob > 0.5:
            num_consec += 1
            potential_tad += cur_aa
        else:
            if num_consec > 5:
                tad_seqs.append(potential_tad)
            potential_tad = ''
            num_consec = 0
    tad_locs = {}
    if tad_seqs != []:
        for tad in tad_seqs:
            start_seq = sequence.index(tad)+1
            end_seq = start_seq + len(tad)
            tad_locs[tad] = [start_seq, end_seq]
    return tad_locs

def predict_polymer_props(sequence):
    '''
    function for predicting polymer properties

    Parameters
    ----------
    sequence : string
        the amino acid sequence as a string

    Returns
    -------
    dict
        dictionary of polymer properties, specifically the Rg and Re
    '''
    seq=pr(sequence)
    Rg=round(seq.predictor.radius_of_gyration(),4)
    Re=round(seq.predictor.end_to_end_distance(),4)
    return {'Rg':Rg,'Re':Re}

