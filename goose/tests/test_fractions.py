"""
A few tests for the seq_by_fractions function. 
"""

import pytest
from goose import create
from goose import goose_exceptions
from goose.backend import parameters
from sparrow.protein import Protein
import re
from goose.data.defined_aa_classes import aa_classes


def test_seq_by_fractions():
    """Test seq_by_fractions function."""
    for i in range(500):
        seq = create.seq_by_fractions(100, G=0, W=0, P=0.1, E=0.1)
        assert len(seq) == 100
        assert seq.count('G') == 0
        assert seq.count('W') == 0
        assert seq.count('P') == 10
        assert seq.count('E') == 10

if __name__ == '__main__':
    pytest.main([__file__])