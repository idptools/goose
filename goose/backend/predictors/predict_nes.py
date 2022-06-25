#nes_predictor.py
# for predicting NES. WARNING - this was not trained on a ton of data,
# I was on a plane... so I used what I had and trained it on my laptop.
# better predictor coming soon!

import os
import numpy as np
from goose.backend.predictors import py_predictor_v2


def predict_nes_seq(sequence):

    # get path to network
    PATH = os.path.dirname(os.path.realpath(__file__))
    
    # selcet the chosen network, kept as separate line of code in 
    used_predictor = 'basic_nes_predictor_v1_100e_nl2_hs20.pt'
    
    # set location of chosen network
    predictor_path = f'{PATH}/networks/{used_predictor}'
    
    # set predictor value using py_predictor_V2
    my_predictor = py_predictor_v2.Predictor(predictor_path, 
                                            dtype="residues")    
    # get values of prediction
    values = my_predictor.predict(sequence)
    
    # make list of NES probability
    NES_probability = []

    # get just the NES probabilites, nothing else.
    for val in values:
        NES_probability.append(round(val[1],4))

    # return the NES probabilities as a list that has 
    # the probability per amino acid across the seq
    return NES_probability



