'''
User facing functionality for analysis of generated proteins.

Limited for now, planning on expanding on in the future
'''

from goose.backend.protein import Protein
from goose.backend.predict_features import predict_mitochondrial_localization, predict_NES, predict_NLS, predict_phosphosites, predict_TAD, predict_polymer_props, get_helical_regions

def properties(sequence, fractions=True):
    '''
    analyzes various sequence properties and returns them to the user

    parameters
    -----------
    sequence : str
        the amino acid sequence as a string
            
    fractions : bool
        whether or not to include the specific fractions 
        of all amino acids for the sequence

    returns
    -------
    props_dict : dict
            returns a dictionary of the properties. 
    '''

    tmp = Protein(sequence)
    props_dict = tmp.calc_basic_properties()
    props_dict['kappa'] = round(tmp.kappa,4)
    
    # if fractions are wanted, add to dict
    if fractions == True:
        props_dict['fractions'] = tmp.fractions

    # return the dict
    return props_dict

def polymer_properties(sequence): 
    '''
    returns the polymer properties of the sequence
    specifically predicted Rg and Re
    '''
    return predict_polymer_props(sequence)

def phosphosites(sequence, raw_vals=False):
    '''
    Function that predicts possible phosphosites
    along the sequence.

    Parameters
    --------------
    sequence : str
        Amino acid sequence

    raw_vals : bool
        Flag which, if set to true, means the dictionary
        returns values 

    '''
    if raw_vals==False:
        return predict_phosphosites(sequence)
    else:
        return predict_phosphosites(sequence, return_sites=True)


def cellular_localization(sequence):
    '''
    returns whether the sequence may have any 
    cellular localization signals

    parameters
    ---------
    sequence : string
            amino acid sequence as a string

    returns 
    -------
    all_locs : dict
            a dictionary of all predicted mitochondrial, NLS, 
            and NES sequences. The dict has three keys
            1. mitochondria
            2. NES
            3. NLS
            within each of those key pair values is another dictionary that has the 
            sequence of the predicted feature as well as the sequence coordinates.
    '''

    all_locs = {}
    mito_loc = predict_mitochondrial_localization(sequence)
    nes_loc = predict_NES(sequence)
    nls_loc = predict_NLS(sequence)

    if mito_loc != {}:
        all_locs['mitochondria'] = mito_loc
    else:
        all_locs['mitochondria'] = 'No mitochondrial targeting sequences predicted.'

    if nes_loc != {}:
        all_locs['NES'] = nes_loc
    else:
        all_locs['NES'] = 'No NES sequences predicted.'

    if nls_loc != {}:
        all_locs['NLS'] = nls_loc
    else:
        all_locs['NLS'] = 'No NLS targeting sequences predicted.'

    return all_locs


def transcriptional_activation(sequence):
    '''
    function to return any potential transcriptional activation domains

    parameters
    ---------
    sequence : string
            amino acid sequence as a string

    returns 
    -------
    all_locs : dict
            a dictionary of the tad sequences and their coordinates in the seq.
    '''
    TADs = predict_TAD(sequence)
    if TADs=={}:
        return {'Predicted TADs':'No transcriptional activation domains predicted.'}
    else:
        return TADs


def everything(sequence, split_predictions = False, just_predictions=False):
    '''
    for when you want all the info in one sweet go.

    parameters
    ---------
    sequence : string
            amino acid sequence as a string

    split_predictions : bool
        whether to split the predictions further. Instead of phosphorylation
        as a single key in the dict, S phosphorylation will be independent
        as will others.

    just_predictions : bool
        whether to return just the predictions

    returns 
    -------
    all_info : dict
            a dictionary of all of the info.

    '''

    if just_predictions == False:
        all_info = properties(sequence)
    else:
        all_info={}

    # now get the rest
    if split_predictions == False:
        all_info['helical regions'] = get_helical_regions(sequence)
        all_info['predicted phosphosites'] = phosphosites(sequence)
        all_info['predicted cellular localization'] = cellular_localization(sequence)
        all_info['predicted transcriptional activation'] = transcriptional_activation(sequence)
        all_info['predicted polymer properties'] = predict_polymer_props(sequence)
    else:
        all_phosphosites = phosphosites(sequence)
        all_info['S phosphorylation'] = all_phosphosites['S']
        all_info['T phosphorylation'] = all_phosphosites['T']
        all_info['Y phosphorylation'] = all_phosphosites['Y']
        all_localization = cellular_localization(sequence)
        pol_props=(sequence)
        all_info['Re']=round(pol_props['Re'],4)
        all_info['Rg']=round(pol_props['Rg'],4)
        all_info['NLS'] = all_localization['NLS']
        all_info['NES'] = all_localization['NES']
        all_info['mitochondrial'] = all_localization['mitochondria']
        all_TADs = transcriptional_activation(sequence)
        all_info['helical regions'] = get_helical_regions(sequence)
        if all_TADs.keys() != 'No TAD sequences':
            list_TADs=[]
            start_num=1
            for TAD in all_TADs.keys():
                cur_seq = TAD
                cur_loc = all_TADs[TAD]
                list_TADs.append(f'{cur_seq} : {cur_loc}')
                start_num+=1
            all_info['predicted transcriptional activation'] = list_TADs
        else:
            all_info['predicted transcriptional activation'] = 'No predicted TADs.'

    if just_predictions==False:
        all_info['fractions'] = Protein(sequence).fractions
        all_info['fraction aromatic'] = round(Protein(sequence).percent_aromatic,5)
        all_info['fraction polar'] = round(Protein(sequence).percent_polar,5)
        all_info['fraction aliphatic']=round(Protein(sequence).percent_aliphatic,5)
        all_info['sequence'] = sequence

    return all_info



def prediction_diffs(sequence_1, sequence_2):
    '''
    Use to show differences in predicted features
    including phosphosites, localization, and TADs.

    parameters
    ----------
    sequence_1 : string
        The amino acid sequence as a string for the first sequence

    sequence_2 : string
        The amino acid sequence as a string for the second sequence

    returns
    -------
    dict_diffs : dict
        returns a dict of differences between sequence 1 and sequnce 2.
    '''
    # get all predictions
    all_predictions_seq1 = everything(sequence_1, split_predictions=True, just_predictions=True)
    all_predictions_seq2 = everything(sequence_2, split_predictions=True, just_predictions=True)

    # make a dict to hold differences
    dict_diffs = {}

    # iterate through predictions
    for prediction in list(all_predictions_seq1.keys()):
        cur_prediction_1 = all_predictions_seq1[prediction]
        cur_prediction_2 = all_predictions_seq2[prediction]
        if prediction == 'predicted transcriptional activation':
            TAD_diffs = []
            for preds in cur_prediction_1:
                if preds not in cur_prediction_2:
                    TAD_diffs.append(f'Sequence 1 predicted TAD - {preds} not in sequence 2')
            for preds in cur_prediction_2:
                if preds not in cur_prediction_1:
                    TAD_diffs.append(f'Sequence 2 predicted TAD - {preds} not in sequence 1')
            if TAD_diffs == []:
                dict_diffs[prediction] = 'No differences.'     
            else:
                dict_diffs[prediction] = TAD_diffs     
        else:
            if cur_prediction_1 != cur_prediction_2:
                dict_diffs[prediction] = f'sequence 1: {cur_prediction_1}, sequence 2: {cur_prediction_2}'
            else:
                dict_diffs[prediction] = f'No differences.'
    return dict_diffs


