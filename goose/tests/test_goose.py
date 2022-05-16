"""
Unit and regression test for the goose package.
"""

# Import package, test suite, and other packages as needed
import sys

import pytest

import goose


def test_goose_imported():
    """Sample test, will always pass so long as import statement worked."""
    assert "goose" in sys.modules
