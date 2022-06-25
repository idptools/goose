import os
from goose.backend.predictors import py_predictor_v2


def predict_tad_seq(sequence):

    # get path to network
    PATH = os.path.dirname(os.path.realpath(__file__))
    
    # selcet the chosen network, kept as separate line of code in 
    used_predictor = 'AD_prediction_network.pt'
    
    # set location of chosen network
    predictor_path = f'{PATH}/networks/{used_predictor}'
    
    # set predictor value using py_predictor_V2
    my_predictor = py_predictor_v2.Predictor(predictor_path, 
                                            dtype="residues")    
    # get values of prediction
    values = my_predictor.predict(sequence)
    
    # make list of phosphorylation probability
    tad_prob = []

    # get just the NLS probabilites, nothing else.
    for val in values:
        tad_prob.append(round(val[1],5))

    # return the NLS probabilities as a list that has 
    # the probability per amino acid across the seq
    return tad_prob


