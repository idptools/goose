'''
Test for GOOSE functionality to generate sequences.
'''

import pytest
import sys
import os
import metapredict as meta
from sparrow import Protein as pr

from goose import create
from goose import analyze
from goose.backend.protein import Protein
from goose.backend import parameters
from goose.backend.sequence_generation import check_disorder
from goose.goose_exceptions import GooseError, GooseInputError


def test_sequence_by_properties_generation():
    '''
    Testing the generation of sequences by different parameters
    '''
    # simple test of just specifying length
    seq_len = create.sequence(100)
    assert len(seq_len) == 100

    # specfying length and FCR
    seq_len_FCR = create.sequence(90, FCR=0.2)
    assert len(seq_len_FCR) == 90
    assert round(Protein(seq_len_FCR).FCR, 4) == 0.2

    # specifying length and NCPR
    seq_len_NCPR = create.sequence(50, NCPR=-0.24)
    assert len(seq_len_NCPR) == 50
    assert round(Protein(seq_len_NCPR).NCPR, 4) == -0.24

    # specifying length and hydropathy
    seq_len_hydropathy = create.sequence(150, hydropathy=3.1)
    assert (len(seq_len_hydropathy)) == 150
    assert abs(Protein(seq_len_hydropathy).hydropathy-3.1) < parameters.HYDRO_ERROR

    # specifying length, ncpr, and FCR
    seq_len_ncpr_fcr = create.sequence(90, FCR=0.4, NCPR=-0.2)
    assert len(seq_len_ncpr_fcr) == 90
    assert round(Protein(seq_len_ncpr_fcr).FCR, 4) == 0.4
    assert round(Protein(seq_len_ncpr_fcr).NCPR, 4) == -0.2

    # specifying length, hydropathy, and FCR
    seq_len_fcr_hydro = create.sequence(40, FCR=0.2, hydropathy=3.5)
    assert len(seq_len_fcr_hydro) == 40
    assert round(Protein(seq_len_fcr_hydro).FCR, 4) == 0.2
    assert abs(Protein(seq_len_fcr_hydro).hydropathy-3.5) < parameters.HYDRO_ERROR

    # specifying length, ncpr, and hydropathy
    seq_len_ncpr_hydro = create.sequence(240, NCPR=-0.2, hydropathy=3.9)
    assert len(seq_len_ncpr_hydro) == 240
    assert round(Protein(seq_len_ncpr_hydro).NCPR, 4) == -0.2
    assert abs(Protein(seq_len_ncpr_hydro).hydropathy-3.9) < parameters.HYDRO_ERROR

    # specifying length, fcr, ncpr, and hydropathy
    seq_len_ncpr_fcr_hydro = create.sequence(80, FCR=0.3, NCPR=-0.2, hydropathy=2.9)
    assert len(seq_len_ncpr_fcr_hydro) == 80
    assert round(Protein(seq_len_ncpr_fcr_hydro).FCR, 4) == 0.3
    assert round(Protein(seq_len_ncpr_fcr_hydro).NCPR, 4) == -0.2
    assert abs(Protein(seq_len_ncpr_fcr_hydro).hydropathy-2.9) < parameters.HYDRO_ERROR

    # specifying length, fcr, ncpr, kappa, and hydropathy
    seq_len_ncpr_fcr_hydro_kappa = create.sequence(140, FCR=0.3, NCPR=-0.2, hydropathy=2.95, kappa=0.15)
    assert len(seq_len_ncpr_fcr_hydro_kappa) == 140
    assert round(Protein(seq_len_ncpr_fcr_hydro_kappa).FCR, 4) == 0.3
    assert round(Protein(seq_len_ncpr_fcr_hydro_kappa).NCPR, 4) == -0.2
    assert abs(Protein(seq_len_ncpr_fcr_hydro_kappa).hydropathy-2.95) < parameters.HYDRO_ERROR
    assert abs(Protein(seq_len_ncpr_fcr_hydro_kappa).kappa - 0.15) < parameters.MAXIMUM_KAPPA_ERROR

    # Finally, verify that the disorder for all sequencs made are within spec.
    all_generated_seqs = [seq_len, seq_len_FCR, seq_len_NCPR, seq_len_hydropathy, seq_len_ncpr_fcr,
                          seq_len_fcr_hydro, seq_len_ncpr_hydro, seq_len_ncpr_fcr_hydro, seq_len_ncpr_fcr_hydro_kappa]
    for sequence in all_generated_seqs:
        assert check_disorder(sequence) == True


def test_sequence_by_fractions_generation():
    '''
    Testing the generation of sequences by different fractions
    '''
    # specify one amino acid above 0
    seq_one_fraction = create.seq_fractions(100, K=0.3)
    assert Protein(seq_one_fraction).fractions['K'] == 0.3
    assert (len(seq_one_fraction)) == 100

    # specify one amino acid at fraction=0
    seq_one_fraction_zero = create.seq_fractions(79, Q=0)
    assert Protein(seq_one_fraction_zero).fractions['Q'] == 0
    assert (len(seq_one_fraction_zero)) == 79

    # specify one amino acid above 0 and one below 0
    seq_one_fraction_one_not = create.seq_fractions(100, K=0.3, Q=0)
    fracs = Protein(seq_one_fraction_one_not).fractions
    assert fracs['K'] == 0.3
    assert fracs['Q'] == 0
    assert (len(seq_one_fraction_one_not)) == 100

    # specify multiple amino acids as zero
    seq_multiple_zero_fractions = create.seq_fractions(193, K=0.0, R=0.0, Q=0.0, A=0.0)
    fracs = Protein(seq_multiple_zero_fractions).fractions
    assert fracs['K'] == 0
    assert fracs['Q'] == 0
    assert fracs['R'] == 0
    assert fracs['A'] == 0
    assert (len(seq_multiple_zero_fractions)) == 193

    # specify multiple non-zero amino acid values
    seq_multiple_specified = create.seq_fractions(200, R=0.05, K=0.1, D=0.01, E=0.01, Q=0.2, A=0.1, N=0.15)
    fracs = Protein(seq_multiple_specified).fractions
    assert fracs['R'] == 0.05
    assert fracs['K'] == 0.1
    assert fracs['D'] == 0.01
    assert fracs['E'] == 0.01
    assert fracs['Q'] == 0.2
    assert fracs['A'] == 0.1
    assert fracs['N'] == 0.15
    assert (len(seq_multiple_specified)) == 200

    # specify fraction values up to one
    seq_all_fractions = create.seq_fractions(100, A=0.02, C=0.01, D=0.03, E=0.01, F=0.01, G=0.06, H=0.1,
                                             I=0.01, K=0.05, L=0.06, M=0.01, N=0.13, P=0.1, Q=0.1, R=0.06, S=0.1, T=0.1, V=0.02, W=0.01, Y=0.01)
    assert len(seq_all_fractions) == 100
    fracs = Protein(seq_all_fractions).fractions
    assert fracs['A'] == 0.02
    assert fracs['C'] == 0.01
    assert fracs['D'] == 0.03
    assert fracs['E'] == 0.01
    assert fracs['F'] == 0.01
    assert fracs['G'] == 0.06
    assert fracs['H'] == 0.1
    assert fracs['I'] == 0.01
    assert fracs['K'] == 0.05
    assert fracs['L'] == 0.06
    assert fracs['M'] == 0.01
    assert fracs['N'] == 0.13
    assert fracs['P'] == 0.1
    assert fracs['Q'] == 0.1
    assert fracs['R'] == 0.06
    assert fracs['S'] == 0.1
    assert fracs['T'] == 0.1
    assert fracs['V'] == 0.02
    assert fracs['W'] == 0.01
    assert fracs['Y'] == 0.01

    # Finally, verify that the disorder for all sequencs made are within spec.
    all_generated_seqs = [seq_one_fraction, seq_one_fraction_zero, seq_one_fraction_one_not,
                          seq_multiple_zero_fractions, seq_multiple_specified, seq_all_fractions]
    for sequence in all_generated_seqs:
        assert check_disorder(sequence) == True

    # test max_aa_fractions over-ride
    new_seq = create.seq_fractions(100, Q=0.9, max_aa_fractions={'Q': 0.95})
    assert Protein(new_seq).fractions['Q'] == 0.9


def test_sequence_by_properties_generation_with_bad_values():

    # check that pasing a Q does not work...
    with pytest.raises(GooseInputError):
        seq_len = create.sequence(100, Q=0.2)


def test_sequence_by_fraction_generation_with_bad_values():
    # check that passing an overloaded number of fractions doesn't work
    with pytest.raises(GooseInputError):
        seq_len = create.seq_fractions(100, Q=1.2)

    # this should fail
    with pytest.raises(GooseInputError):
        seq_len = create.seq_fractions(100, W=0.9)
