'''
Functions that are in beta testing. Most will probably
never be implemented.
'''

import goose.backend.predict_disorder
from goose.backend.predict_disorder import calculate_disorder_for_AA
from goose.backend.amino_acid_lists import disordered_list
from goose.backend.constant_class_functionality import generate_random_seq, random_amino_acid



# function to try to predict each value for a sequence in real time using a seed sequence
# and adding the amino acid to the sequence as it goes.

def generate_disordered_seq(length, cutoff = 0.9):
    '''
    predicts disorder at every step for all amino acids
    just out of curiousity
    '''

    # make a seed sequence
    seed_seq = generate_random_seq(50, disordered_list)
    
    # make the sequence
    final_sequence = ""

    # add to final sequence
    for i in range(0, length):

        # calculate disorder scores for all amino acids
        if len(final_sequence) < 20:
            disorder_scores = calculate_disorder_for_AA(seed_seq)
        else:
            disorder_scores = calculate_disorder_for_AA(final_sequence)
        
        # possible amino acid list
        possible_amino_acids = []
        # iterate through dict and get corresponding aas
        for i in disorder_scores.keys():
            if disorder_scores[i] >= cutoff:
                possible_amino_acids.append(i)

        # if the list is empty grab the top 3
        if possible_amino_acids == []:
            possible_scores = sorted(disorder_scores.values(), reverse=True)[0:2]
            for i in disorder_scores.keys():
                if disorder_scores[i] in possible_scores:
                    possible_amino_acids.append(i)
                
        # add amino acids to the final sequence and seed sequence
        chosen_amino_acid = random_amino_acid(possible_amino_acids)
        final_sequence += chosen_amino_acid
        seed_seq += chosen_amino_acid

    return final_sequence

generated_seqs = []
