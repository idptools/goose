import os
import numpy as np
from goose.backend.predictors import py_predictor_v2



def predict_phospho_S(sequence):

    # get path to network
    PATH = os.path.dirname(os.path.realpath(__file__))
    
    # selcet the chosen network, kept as separate line of code in 
    used_predictor = 'S_phospho_network.pt'
    
    # set location of chosen network
    predictor_path = f'{PATH}/networks/{used_predictor}'
    
    # set predictor value using py_predictor_V2
    my_predictor = py_predictor_v2.Predictor(predictor_path, 
                                            dtype="residues")    
    # get values of prediction
    values = my_predictor.predict(sequence)
    
    # make list of phosphorylation probability
    s_phospho_probability = []

    # get just the NLS probabilites, nothing else.
    for val in values:
        s_phospho_probability.append(round(val[1],5))

    # return the NLS probabilities as a list that has 
    # the probability per amino acid across the seq
    return s_phospho_probability




def predict_phospho_T(sequence):

    # get path to network
    PATH = os.path.dirname(os.path.realpath(__file__))
    
    # selcet the chosen network, kept as separate line of code in 
    used_predictor = 'T_phospho_network.pt'
    
    # set location of chosen network
    predictor_path = f'{PATH}/networks/{used_predictor}'
    
    # set predictor value using py_predictor_V2
    my_predictor = py_predictor_v2.Predictor(predictor_path, 
                                            dtype="residues")    
    # get values of prediction
    values = my_predictor.predict(sequence)
    
    # make list of phosphorylation probability
    t_phospho_probability = []

    # get just the NLS probabilites, nothing else.
    for val in values:
        t_phospho_probability.append(round(val[1],5))

    # return the NLS probabilities as a list that has 
    # the probability per amino acid across the seq
    return t_phospho_probability



def predict_phospho_Y(sequence):

    # get path to network
    PATH = os.path.dirname(os.path.realpath(__file__))
    
    # selcet the chosen network, kept as separate line of code in 
    used_predictor = 'Y_phospho_network.pt'
    
    # set location of chosen network
    predictor_path = f'{PATH}/networks/{used_predictor}'
    
    # set predictor value using py_predictor_V2
    my_predictor = py_predictor_v2.Predictor(predictor_path, 
                                            dtype="residues")    
    # get values of prediction
    values = my_predictor.predict(sequence)
    
    # make list of phosphorylation probability
    y_phospho_probability = []

    # get just the NLS probabilites, nothing else.
    for val in values:
        y_phospho_probability.append(round(val[1],5))

    # return the NLS probabilities as a list that has 
    # the probability per amino acid across the seq
    return y_phospho_probability


# predict all 3 sites

def predict_phosphorylation(sequence):
    '''
    returns a dict where the amino acid letter for each
    site is the key to the phosphorylation probability values
    '''
    # phophorlation dictionary
    phospho_dict = {}
    phospho_dict['Y'] = predict_phospho_Y(sequence)
    phospho_dict['S'] = predict_phospho_S(sequence)
    phospho_dict['T'] = predict_phospho_T(sequence)

    return phospho_dict

