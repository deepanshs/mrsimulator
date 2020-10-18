# -*- coding: utf-8 -*-
import numpy as np
from mrsimulator import Site
from mrsimulator import SpinSystem
from mrsimulator.method.utils import query_permutations

__author__ = "Maxwell Venetos"
__email__ = "mvenetos@berkeley.edu"


def basic_transition_query_tests(iso):
    # Single site tests
    test_1 = query_permutations(
        query={"P": {"channel-1": [[-1]]}},
        isotope=iso[0].get_isotopes(),
        channel=["1H"],
    )
    assert test_1 == [np.array([-1.0])]
    test_2 = query_permutations(
        query={"P": {"channel-1": [[-1], [1]]}},
        isotope=iso[0].get_isotopes(),
        channel=["1H"],
    )
    assert test_2 == [np.array([-1.0]), np.array([1.0])]

    # Multi sites same channel tests
    test_3 = query_permutations(
        query={"P": {"channel-1": [[-1, 1]]}},
        isotope=iso[1].get_isotopes(),
        channel=["1H"],
    )
    test_3_check = [
        np.array([-1.0, 1.0, 0.0]),
        np.array([-1.0, 0.0, 1.0]),
        np.array([1.0, -1.0, 0.0]),
        np.array([1.0, 0.0, -1.0]),
        np.array([0.0, -1.0, 1.0]),
        np.array([0.0, 1.0, -1.0]),
    ]
    assert np.array_equal(test_3, test_3_check)

    # Multi sites same channel tests
    test_4 = query_permutations(
        query={"P": {"channel-1": [[-1]]}},
        isotope=iso[2].get_isotopes(),
        channel=["17O"],
    )
    assert np.array_equal(test_4, [np.array([0.0, -1.0, 0.0])])


def test_transition_query():
    H1 = Site(isotope="1H", shielding_symmetric=dict(zeta=50, eta=0))
    Si29 = Site(isotope="29Si", shielding_symmetric=dict(zeta=50, eta=0))
    O17 = Site(isotope="17O", quadrupolar=dict(Cq=50, eta=0))

    iso = [
        SpinSystem(sites=[H1]),
        SpinSystem(sites=[H1, H1, H1]),
        SpinSystem(sites=[Si29, O17, Si29]),
    ]
    basic_transition_query_tests(iso)
