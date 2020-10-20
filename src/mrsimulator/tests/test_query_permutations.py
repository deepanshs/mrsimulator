# -*- coding: utf-8 -*-
import numpy as np
import pytest
from mrsimulator import Site
from mrsimulator import SpinSystem
from mrsimulator.method.utils import query_permutations

__author__ = "Maxwell Venetos"
__email__ = "mvenetos@berkeley.edu"

H1 = Site(isotope="1H", shielding_symmetric=dict(zeta=50, eta=0))
Si29 = Site(isotope="29Si", shielding_symmetric=dict(zeta=50, eta=0))
O17 = Site(isotope="17O", quadrupolar=dict(Cq=50, eta=0))


def basic_transition_query_tests(iso):
    # Single site tests
    test_1 = query_permutations(
        query={"P": {"channel-1": [[-1]]}},
        isotope=iso[0].get_isotopes(),
        channel=["1H"],
    )
    assert np.array_equal(test_1, [[-1.0]])

    test_2 = query_permutations(
        query={"P": {"channel-1": [[-1], [1]]}},
        isotope=iso[0].get_isotopes(),
        channel=["1H"],
    )
    assert np.array_equal(test_2, [[-1.0], [1.0]])

    # Multi sites same channel tests
    test_3 = query_permutations(
        query={"P": {"channel-1": [[-1, 1]]}},
        isotope=iso[1].get_isotopes(),
        channel=["1H"],
    )
    test_3_check = [
        [-1.0, 1.0, 0.0],
        [-1.0, 0.0, 1.0],
        [1.0, -1.0, 0.0],
        [1.0, 0.0, -1.0],
        [0.0, -1.0, 1.0],
        [0.0, 1.0, -1.0],
    ]
    for item in test_3_check:
        assert item in test_3

    # Multi sites same channel tests
    test_4 = query_permutations(
        query={"P": {"channel-1": [[-1, -1]]}},
        isotope=iso[1].get_isotopes(),
        channel=["1H"],
    )
    test_4_check = [
        [-1.0, -1.0, 0.0],
        [-1.0, 0.0, -1.0],
        [0.0, -1.0, -1.0],
    ]
    for item in test_4_check:
        assert item in test_4

    # Multi sites same channel tests
    test_5 = query_permutations(
        query={"P": {"channel-1": [[-1]]}},
        isotope=iso[2].get_isotopes(),
        channel=["17O"],
    )
    assert np.array_equal(test_5, [[0.0, -1.0, 0.0]])

    test_6 = query_permutations(
        query={"P": {"channel-1": [[-1]], "channel-2": [[2]]}},
        isotope=iso[2].get_isotopes(),
        channel=["29Si", "17O"],
    )
    test_6_check = [[-1.0, 2.0, 0.0], [0.0, 2.0, -1.0]]
    for item in test_6_check:
        assert item in test_6


def test_transition_query():
    iso = [
        SpinSystem(sites=[H1]),
        SpinSystem(sites=[H1, H1, H1]),
        SpinSystem(sites=[Si29, O17, Si29]),
    ]
    basic_transition_query_tests(iso)


def test_two_site():
    # sys = SpinSystem(sites=[H1])
    with pytest.raises(ValueError, match=".*The length of the transition query*"):
        query_permutations(
            query={"P": {"channel-1": [[-1, 1]]}},
            isotope=["1H"],
            channel=["1H"],
        )
