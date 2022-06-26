import os
import numpy as np
from goose.backend.predictors import py_predictor_v2



def predict_nls_seq(sequence):

    # get path to network
    PATH = os.path.dirname(os.path.realpath(__file__))
    
    # selcet the chosen network, kept as separate line of code in 
    used_predictor = 'NLS_predictor_200e_nl1_v1.pt'
    
    # set location of chosen network
    predictor_path = f'{PATH}/networks/{used_predictor}'
    
    # set predictor value using py_predictor_V2
    my_predictor = py_predictor_v2.Predictor(predictor_path, 
                                            dtype="residues")    
    # get values of prediction
    values = my_predictor.predict(sequence)
    
    # make list of NLS probability
    NLS_probability = []

    # get just the NLS probabilites, nothing else.
    for val in values:
        NLS_probability.append(round(val[1],5))

    # return the NLS probabilities as a list that has 
    # the probability per amino acid across the seq
    return NLS_probability

