'''
Test for GOOSE functionality in calculating the specific 
properties.
'''
import pytest
import sys
import os

from sparrow import Protein as pr

from goose.backend.protein import Protein 
<<<<<<< HEAD
from goose.tests.local_data import test_prop_dicts
=======
#from goose.tests.local_data import test_prop_dicts
from ..tests.local_data import test_prop_dicts
>>>>>>> origin/sparrow_predictions
from goose.backend.predictors.predict_mito import predict_mitochondrial_targeting
from goose.backend.predictors.predict_nes import predict_nes_seq
from goose.backend.predictors.predict_nls import predict_nls_seq
from goose.backend.predictors.predict_phosphosite import predict_phospho_S
from goose.backend.predictors.predict_phosphosite import predict_phospho_T
from goose.backend.predictors.predict_phosphosite import predict_phospho_Y
from goose.backend.predictors.predict_tad import predict_tad_seq
from goose.analyze import phosphosites
from goose.analyze import cellular_localization
from goose.analyze import transcriptional_activation
from goose.analyze import everything

# test the calculation of basic sequence properties
def test_goose_property_calculations():
    # iterate through precalculated sequence properties
    for seq_and_props in test_prop_dicts:
        # grab current sequence
        current_sequence = seq_and_props['sequence']
<<<<<<< HEAD
        # test FCR
        assert Protein.calc_FCR(current_sequence) == seq_and_props['FCR']
        # test NCPR
        assert Protein.calc_NCPR(current_sequence) == seq_and_props['NCPR']
        # test mean hydropathy
        assert Protein.calc_mean_hydro(current_sequence) == seq_and_props['hydro']        
        # test kappa
        assert pr(current_sequence).kappa == seq_and_props['kappa']
=======
        current_sequence = Protein(current_sequence)
        # test FCR
        assert current_sequence.FCR == seq_and_props['FCR']
        # test NCPR
        assert current_sequence.NCPR == seq_and_props['NCPR']
        # test mean hydropathy
        assert current_sequence.hydropathy == seq_and_props['hydro']        
        # test kappa
        assert current_sequence.kappa == seq_and_props['kappa']
>>>>>>> origin/sparrow_predictions


def within_spec(list_vals1, list_vals2, acceptable_diff = 0.000001):
    '''
    function to make sure that the difference between
    the values in 2 lists is less than a specified
    value. Using becasue np.float32 keeps persisting
    despite my best efforts. FFS.

    parameters
    ----------
    list_vals1 : list
        A list of vals to check that should be equivalent to
        list_vals2

    list_vals2 : list
        A list of vals to check that should be equivalent to
        list_vals1

    acceptable_diff : Float
        The acceptable difference between each value in the list
    '''

    for val in range(0, len(list_vals1)):
        if abs(list_vals1[val]-list_vals2[val]) > acceptable_diff:
            return False
    return True

# test the predictors as far as the raw values returned for each
def test_goose_raw_predictions():
    # iterate through precalculated sequence properties
    for seq_and_props in test_prop_dicts:
        # grab current sequence
        current_sequence = seq_and_props['sequence']
        # test raw mitochondrial targeting scores
        assert predict_mitochondrial_targeting(current_sequence) == seq_and_props['mitochondrial_targeting_raw']
        # test raw NES scores
        assert within_spec(predict_nes_seq(current_sequence), seq_and_props['NES_raw'])==True
        # test raw NLS scores
        assert within_spec(predict_nls_seq(current_sequence), seq_and_props['NLS_raw'])==True
        # test raw S phoshphorylation scores
        assert within_spec(predict_phospho_S(current_sequence), seq_and_props['phosphorylation_S_raw'])==True
        # test raw T phosphorylation scores
        assert within_spec(predict_phospho_T(current_sequence), seq_and_props['phosphorylation_T_raw'])==True
        # test raw Y phosphorylation scores
        assert within_spec(predict_phospho_Y(current_sequence), seq_and_props['phosphorylation_Y_raw'])==True
        # test raw TAD scores
        assert within_spec(predict_tad_seq(current_sequence), seq_and_props['TAD_raw'])==True

# test the analyze module everything() function
def test_goose_analyze_predictions():
    # iterate through precalculated sequence properties
    for seq_and_props in test_prop_dicts:
        # grab current sequence
        current_sequence = seq_and_props['sequence']
<<<<<<< HEAD
        # make sure all predictions in the analyze.everything function come back correct
        assert everything(current_sequence) == seq_and_props['complete_analysis']
=======
        # get current predictions
        current_predictions = everything(current_sequence)
        # make sure all predictions in the analyze.everything function come back correct
        known_predictions = seq_and_props['complete_analysis']
        # tests
        testing_factors = list(current_predictions.keys())
        # need to round because sparrow implementation...
        round_me = ['FCR', 'NCPR', 'hydropathy', 'sigma', 'delta']
        for factor in testing_factors:
            if factor in round_me:
                assert round(current_predictions[factor], 6) == round(known_predictions[factor], 6)

>>>>>>> origin/sparrow_predictions






