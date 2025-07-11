
"""
Comprehensive pytest tests for create.sequence variant functions.

This module contains thorough tests for the main sequence variant functions.
Tests should cover all parameters, edge cases, and error conditions...
"""

import pytest
from goose import create
from goose import goose_exceptions
from goose.backend import parameters
from sparrow.protein import Protein
from goose.data.defined_aa_classes import aa_classes
from goose.backend_variant_generation.verify_variants import (
    verify_same_number_by_class,
    verify_same_order_by_class,
    verify_constant_properties,
    verify_same_residue_count,
    verify_constant_hydropathy,
    verify_constant_kappa,
    verify_constant_FCR,
    verify_constant_NCPR,
    verify_target_hydropathy,
    verify_target_kappa,
    verify_target_FCR,
    verify_target_NCPR,
    verify_constant_region,
    verify_region_changed,
    verify_constant_residue_positions,
    verify_disorder)

'''
For each of the functions, we need to verify specific things, going to list that here:



'''



if __name__ == '__main__':
    pytest.main([__file__])
