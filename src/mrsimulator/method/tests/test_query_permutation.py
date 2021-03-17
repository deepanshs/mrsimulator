# -*- coding: utf-8 -*-
import numpy as np
from mrsimulator.method.transition_query import TransitionQuery

__author__ = "Maxwell Venetos"
__email__ = "mvenetos@berkeley.edu"

# H1 = Site(isotope="1H", shielding_symmetric=dict(zeta=50, eta=0))
# Si29 = Site(isotope="29Si", shielding_symmetric=dict(zeta=50, eta=0))
# O17 = Site(isotope="17O", quadrupolar=dict(Cq=50, eta=0))


def check_equal(query, isotopes, channels, res):
    test = TransitionQuery(**query).permutation(isotopes, channels)
    assert np.allclose(test["P"], res[0])
    assert np.allclose(test["D"], res[1])


def tests_00():
    # Single site tests
    # P = -1 on A system, channel A
    query = {"ch1": {"P": [-1]}}
    res = [[-1]], []
    check_equal(query, ["A"], ["A"], res)

    # P = -1 on A A system, channel A
    res = [[-1, 0], [0, -1]], []
    check_equal(query, ["A", "A"], ["A"], res)

    # P = -1 on A A B system, channel A
    res = [[-1, 0, 0], [0, -1, 0]], []
    check_equal(query, ["A", "A", "B"], ["A"], res)

    # P = -1 on A A B system, channel A, B
    check_equal(query, ["A", "A", "B"], ["A", "B"], res)

    # P = -1 on A A B system, channel A, B, C
    check_equal(query, ["A", "A", "B"], ["A", "B", "C"], res)

    # P = -1 on A B B A system, channel A
    res = [[-1, 0, 0, 0], [0, 0, 0, -1]], []
    check_equal(query, ["A", "B", "B", "A"], ["A"], res)

    # P = -1 on A B B A system, channel A, B
    check_equal(query, ["A", "B", "B", "A"], ["A", "B"], res)

    # P = -1 -1 on A B B A system, channel A
    query = {"ch1": {"P": [-1, -1]}}
    res = [[-1, 0, 0, -1]], []
    check_equal(query, ["A", "B", "B", "A"], ["A"], res)

    # P = -1 -1 on A B B A system, channel A, B
    check_equal(query, ["A", "B", "B", "A"], ["A", "B"], res)

    # P = -1 D = 1 on A B B A system, channel A
    query = {"ch1": {"P": [-1], "D": [1]}}
    res = [[-1, 0, 0, 0], [0, 0, -1, 0]], [[1, 0, 0, 0], [0, 0, 1, 0]]
    check_equal(query, ["A", "B", "A", "B"], ["A"], res)

    # P = -1 D=1 on A B B A system, channel A, B
    check_equal(query, ["A", "B", "A", "B"], ["A", "B"], res)


def tests_01():
    # two channels
    # A B B A system, P = -1 channel A| P = -1 channel B
    query = {"ch1": {"P": [-1]}, "ch2": {"P": [-1]}}
    res = [[-1, -1, 0, 0], [-1, 0, -1, 0], [0, -1, 0, -1], [0, 0, -1, -1]], []
    check_equal(query, ["A", "B", "B", "A"], ["A", "B"], res)

    # A B B A system, P = -1 channel A| P = -1 channel B, but active channels = A.
    res = [[-1, 0, 0, 0], [0, 0, 0, -1]], []
    check_equal(query, ["A", "B", "B", "A"], ["A"], res)

    # A B B A system, P = -1 channel A| P = -1, D = 2 channel B.
    query = {"ch1": {"P": [-1], "D": [0]}, "ch2": {"P": [-1], "D": [2]}}
    res = (
        [[-1, -1, 0, 0], [-1, 0, -1, 0], [0, -1, 0, -1], [0, 0, -1, -1]],
        [[0, 0, 2, 0], [0, 2, 0, 0]],
    )
    check_equal(query, ["A", "B", "B", "A"], ["A", "B"], res)

    # A B B A system, P = -1 D = 1 channel A| P = -1, 1 channel B.
    query = {"ch1": {"P": [-1], "D": [1]}, "ch2": {"P": [-1, 1], "D": [0]}}
    res = (
        [[-1, -1, 1, 0], [-1, 1, -1, 0], [0, -1, 1, -1], [0, 1, -1, -1]],
        [[1, 0, 0, 0], [0, 0, 0, 1]],
    )
    check_equal(query, ["A", "B", "B", "A"], ["A", "B"], res)

    # query = {"ch1": {"P": [-1]}}, {"ch1": None, "ch2": {"P": [-1]}}
    # res = ([[-1, 0, 0, 0], [0, 0, 0, -1], [0, -1, 0, 0], [0, 0, -1, 0]], [])
    # check_equal(query, ["A", "B", "B", "A"], ["A", "B"], res)


#     test_2 = query_permutations(
#         query={"P": {"channel-1": [[-1], [1]]}},
#         isotope=iso[0].get_isotopes(symbol=True),
#         channel=["1H"],
#     )
#     check_equal(test_2, [[-1], [1]])
#     # Multi sites same channel tests
#     test_3 = query_permutations(
#         query={"P": {"channel-1": [[-1, 1]]}},
#         isotope=iso[1].get_isotopes(symbol=True),
#         channel=["1H"],
#     )
#     test_3_check = [
#         [-1, 0, 1],
#         [0, -1, 1],
#         [1, 0, -1],
#         [-1, 1, 0],
#         [1, -1, 0],
#         [0, 1, -1],
#     ]
#     check_equal(test_3, test_3_check)
#     # Multi sites same channel tests
#     test_4 = query_permutations(
#         query={"P": {"channel-1": [[-1, -1]]}},
#         isotope=iso[1].get_isotopes(symbol=True),
#         channel=["1H"],
#     )
#     test_4_check = [
#         [0, -1, -1],
#         [-1, -1, 0],
#         [-1, 0, -1],
#     ]
#     check_equal(test_4, test_4_check)
#     # Multi sites same channel tests
#     test_5 = query_permutations(
#         query={"P": {"channel-1": [[-1]]}},
#         isotope=iso[2].get_isotopes(symbol=True),
#         channel=["17O"],
#     )
#     check_equal(test_5, [[0, -1, 0]])
#     test_6 = query_permutations(
#         query={"P": {"channel-1": [[-1]], "channel-2": [[2]]}},
#         isotope=iso[2].get_isotopes(symbol=True),
#         channel=["29Si", "17O"],
#     )
#     test_6_check = [[-1, 2, 0], [0, 2, -1]]
#     check_equal(test_6, test_6_check)
#     # test by swapping the channels and channel query
#     test_7 = query_permutations(
#         query={"P": {"channel-1": [[-1, -1], [1]], "channel-2": [[2]]}},
#         isotope=iso[2].get_isotopes(symbol=True),
#         channel=["29Si", "17O"],
#     )
#     test_7_check = [[-1, 2, -1], [1, 2, 0], [0, 2, 1]]
#     check_equal(test_7, test_7_check)
#     test_7 = query_permutations(
#         query={"P": {"channel-1": [[2]], "channel-2": [[-1, -1], [1]]}},
#         isotope=iso[2].get_isotopes(symbol=True),
#         channel=["17O", "29Si"],
#     )
#     check_equal(test_7, test_7_check)
#     test_7 = query_permutations(
#         query={"P": {"channel-2": [[-1, -1], [1]], "channel-1": [[2]]}},
#         isotope=iso[2].get_isotopes(symbol=True),
#         channel=["17O", "29Si"],
#     )
#     check_equal(test_7, test_7_check)
#     test_7 = query_permutations(
#         query={"P": {"channel-2": [[2]], "channel-1": [[-1, -1], [1]]}},
#         isotope=iso[2].get_isotopes(symbol=True),
#         channel=["29Si", "17O"],
#     )
#     check_equal(test_7, test_7_check)
#     test_8 = query_permutations(
#         query={"P": {"channel-1": [[-1, -1], [1]]}},
#         isotope=iso[2].get_isotopes(symbol=True),
#         channel=["29Si", "27Al"],
#     )
#     test_8_check = []
#     check_equal(test_8, test_8_check)
#     test_9 = query_permutations(
#         query={"P": {"channel-1": [[-1, -1], [1]], "channel-2": [[2]]}},
#         isotope=iso[2].get_isotopes(symbol=True),
#         channel=["29Si", "27Al"],
#     )
#     test_9_check = []
#     check_equal(test_9, test_9_check)


# def test_transition_query():
#     iso = [
#         SpinSystem(sites=[H1]),
#         SpinSystem(sites=[H1, H1, H1]),
#         SpinSystem(sites=[Si29, O17, Si29]),
#     ]
#     basic_transition_query_tests(iso)


# def test_two_site():
#     # sys = SpinSystem(sites=[H1])
#     with pytest.raises(ValueError, match=".*The length of the transition query*"):
#         query_permutations(
#             query={"P": {"channel-1": [-1, 1]}},
#             isotope=["1H"],
#             channel=["1H"],
#         )
