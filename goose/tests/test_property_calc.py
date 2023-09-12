'''
Test for GOOSE functionality in calculating the specific 
properties.
'''
import pytest
import sys
import os

from sparrow import Protein as pr

from goose.backend.protein import Protein
from ..tests.local_data import test_prop_dicts


def test_goose_property_calculations():
    # iterate through precalculated sequence properties
    for seq_and_props in test_prop_dicts:
        # grab current sequence
        current_sequence = seq_and_props['sequence']
        current_sequence = Protein(current_sequence)
        # test FCR
        assert current_sequence.FCR == seq_and_props['FCR']
        # test NCPR
        assert current_sequence.NCPR == seq_and_props['NCPR']
        # test mean hydropathy
        assert current_sequence.hydropathy == seq_and_props['hydropathy']
        # test kappa
        assert current_sequence.kappa == seq_and_props['kappa']


def within_spec(list_vals1, list_vals2, acceptable_diff=0.000001, round_to=4):
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
        if abs(round(list_vals1[val], round_to)-round(list_vals2[val], round_to)) > acceptable_diff:
            return False
    return True



