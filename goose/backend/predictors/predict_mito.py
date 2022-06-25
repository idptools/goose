'''
need to snag the sparrow network first
'''


import os
import numpy as np
from goose.backend.predictors import py_predictor_v2
from goose.backend.predictors import encode_sequence
from goose.backend.predictors import brnn_architecture
import torch



def predict_mitochondrial_targeting(seq):
    """
    Prediction function. seq should be a valid amino acid sequence.

    NOTE that this assumes mitochondrial targeting sequences (MTSs) are 
    N-terminal, so truncates anything over 168 residues. This threshold
    was empyrically determined based on the set of annottated MTSs.

    Parameters
    ------------
    seq : str
        Valid amino acid sequence

    Returns
    ----------
    np.ndarray
        Returns a 1D np.ndarray the length of the sequence where each position
        is the transient helicity at that position.

    """

    # convert sequence to uppercase
    seq = seq.upper()

    # truncate all but 168 - if shorter than this just gets everything
    if len(seq)>=168:
        sub_seq = seq[0:168]
    else:
        sub_seq = seq

    # get path to network
    PATH = os.path.dirname(os.path.realpath(__file__))

    # selcet the chosen network, kept as separate line of code in 
    used_predictor = 'mitochondrial_targeting_predictor_network_v1.pt'

    # set location of chosen network
    predictor_path = f'{PATH}/networks/{used_predictor}'

    saved_weights = predictor_path

    loaded_model = torch.load(saved_weights, map_location=torch.device('cpu'))

    # Dynamically read in correct hyperparameters:
    num_layers = 0
    while True:
        s = f'lstm.weight_ih_l{num_layers}'
        try:
            temp = loaded_model[s]
            num_layers += 1
        except KeyError:
            break
                 
    ##  determine the number of classes; note you may need to change the key names here no leading
    # module. in ther
    number_of_classes = np.shape(loaded_model['fc.bias'])[0]
    input_size = 20 # (hardcoded at 20 for 20 amino acids)

    hidden_vector_size = int(np.shape(loaded_model['lstm.weight_ih_l0'])[0] / 4)

    # Instantiate network weights into object
    network = brnn_architecture.BRNN_MtM(input_size, hidden_vector_size, num_layers, number_of_classes, 'cpu')                                          
    network.load_state_dict(loaded_model)


    # Convert to one-hot sequence vector
    seq_vector = encode_sequence.one_hot(sub_seq)
    seq_vector = seq_vector.view(1, len(seq_vector), -1)  # formatting

    # Forward pass  -this is specific for classication
    prediction = network(seq_vector.float()).detach().numpy()

    int_vals = []
    for row in prediction[0]:
        int_vals.append(np.argmax(row))

    prediction = int_vals

    # append empty 0s for remainder of sequence
    extra = [0]*(len(seq)-len(sub_seq))

    prediction.extend(extra)
    # return prediction + extra zeros
    return prediction

