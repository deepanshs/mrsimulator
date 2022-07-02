# -*- coding: utf-8 -*-
import numpy as np
from mrsimulator.method.query import TransitionQuery

__author__ = "Maxwell Venetos"
__email__ = "mvenetos@berkeley.edu"

NAN = np.nan


def check_equal(query, isotopes, channels, desired):
    test = TransitionQuery(**query).combination(isotopes, channels)

    for element in ["P", "D"]:
        assert len(test[element]) == len(desired[element])
        for test_item in test[element]:
            assert np.any(
                [
                    np.allclose(test_item, desired_item, equal_nan=True)
                    for desired_item in desired[element]
                ]
            ), f"Error in element {element}"


def test_00():
    # Single site tests
    # P = -1 on A system, channel A
    query = {"ch1": {"P": [-1]}}
    desired = {"P": [[-1]], "D": []}
    check_equal(query, ["A"], ["A"], desired)

    # P = -1 on A A system, channel A
    desired = {"P": [[-1, 0], [0, -1]], "D": []}
    check_equal(query, ["A", "A"], ["A"], desired)

    # P = -1 on A A B system, channel A
    desired = {"P": [[-1, 0, 0], [0, -1, 0]], "D": []}
    check_equal(query, ["A", "A", "B"], ["A"], desired)

    # P = -1 on A A B system, channel A, B
    check_equal(query, ["A", "A", "B"], ["A", "B"], desired)

    # P = -1 on A A B system, channel A, B, C
    check_equal(query, ["A", "A", "B"], ["A", "B", "C"], desired)

    # P = -1 on A B B A system, channel A
    desired = {"P": [[-1, 0, 0, 0], [0, 0, 0, -1]], "D": []}
    check_equal(query, ["A", "B", "B", "A"], ["A"], desired)

    # P = -1 on A B B A system, channel A, B
    check_equal(query, ["A", "B", "B", "A"], ["A", "B"], desired)

    # P = -1 -1 on A B B A system, channel A
    query = {"ch1": {"P": [-1, -1]}}
    desired = {"P": [[-1, 0, 0, -1]], "D": []}
    check_equal(query, ["A", "B", "B", "A"], ["A"], desired)

    # P = -1 -1 on A B B A system, channel A, B
    check_equal(query, ["A", "B", "B", "A"], ["A", "B"], desired)

    # P = -1 -1 on A B B A system, channel A
    desired = {"P": [], "D": []}
    check_equal(query, ["A"], ["A"], desired)

    # P = -1 D = 1 on A B B A system, channel A
    query = {"ch1": {"P": [-1], "D": [1]}}
    desired = {
        "P": [[-1, 0, 0, 0], [0, 0, -1, 0]],
        "D": [[1, NAN, 00, NAN], [00, NAN, 1, NAN]],
    }
    check_equal(query, ["A", "B", "A", "B"], ["A"], desired)

    # P = -1 D=1 on A B B A system, channel A, B
    check_equal(query, ["A", "B", "A", "B"], ["A", "B"], desired)


def test_01():
    query = {"ch1": {"P": [3], "D": [0]}}
    desired = {"P": [[3, 0]], "D": [[0, NAN]]}
    check_equal(query, ["A", "B"], ["A"], desired)

    query = {"ch1": {"P": [-1], "D": [0]}}
    desired = {"P": [[-1, 0]], "D": [[0, NAN]]}
    check_equal(query, ["A", "B"], ["A"], desired)


def test_02():
    # two channels
    # A B B A system, P = -1 channel A| P = -1 channel B
    query = {"ch1": {"P": [-1]}, "ch2": {"P": [-1]}}
    desired = {
        "P": [[-1, -1, 0, 0], [-1, 0, -1, 0], [0, -1, 0, -1], [0, 0, -1, -1]],
        "D": [],
    }
    check_equal(query, ["A", "B", "B", "A"], ["A", "B"], desired)

    # A B B A system, P = -1 channel A| P = -1 channel B, but active channels = A.
    desired = {"P": [[-1, 0, 0, 0], [0, 0, 0, -1]], "D": []}
    check_equal(query, ["A", "B", "B", "A"], ["A"], desired)

    # A B B A system, P = -1 channel A| P = -1, D = 2 channel B.
    query = {"ch1": {"P": [-1], "D": [0]}, "ch2": {"P": [-1], "D": [2]}}
    desired = {
        "P": [[-1, -1, 0, 0], [-1, 0, -1, 0], [0, -1, 0, -1], [0, 0, -1, -1]],
        "D": [[0, 0, 2, 0], [0, 2, 0, 0]],
    }
    check_equal(query, ["A", "B", "B", "A"], ["A", "B"], desired)

    # A B B A system, P = -1 D = 1 channel A| P = -1, 1 channel B.
    query = {"ch1": {"P": [-1], "D": [1]}, "ch2": {"P": [-1, 1], "D": [0]}}
    desired = {
        "P": [[-1, -1, 1, 0], [-1, 1, -1, 0], [0, -1, 1, -1], [0, 1, -1, -1]],
        "D": [[1, 0, 0, 0], [0, 0, 0, 1]],
    }
    check_equal(query, ["A", "B", "B", "A"], ["A", "B"], desired)

    query = {}
    desired = {"P": [[0, 0, 0, 0]], "D": []}
    check_equal(query, ["A", "B", "B", "A"], ["A", "B"], desired)

    query = {"ch1": {"P": [-1]}, "ch2": {"D": [0]}}
    desired = {"P": [[-1, 0]], "D": [[NAN, 0]]}
    check_equal(query, ["A", "B"], ["A", "B"], desired)

    query = {"ch1": {"P": [-1]}, "ch2": {"D": [1]}}
    desired = {
        "P": [[-1, 0, 0, 0], [0, 0, 0, -1]],
        "D": [[NAN, 1, 0, NAN], [NAN, 0, 1, NAN]],
    }
    check_equal(query, ["A", "B", "B", "A"], ["A", "B"], desired)

    query = {"ch1": {"P": [-1], "D": [1]}, "ch2": {"P": [-1, 1], "D": [0]}}
    desired = {"P": [], "D": [[1, 0, NAN, 0], [0, 0, NAN, 1]]}
    check_equal(query, ["A", "B", "C", "A"], ["A", "B", "C"], desired)

    query = {"ch1": {"P": [-1, 1], "D": [-1, 1]}, "ch2": {"P": [-1]}}
    desired = {
        "P": [
            [-1, 0, -1, 1],
            [0, -1, -1, 1],
            [1, 0, -1, -1],
            [-1, 1, -1, 0],
            [1, -1, -1, 0],
            [0, 1, -1, -1],
        ],
        "D": [
            [-1, 0, NAN, 1],
            [0, -1, NAN, 1],
            [1, 0, NAN, -1],
            [-1, 1, NAN, 0],
            [1, -1, NAN, 0],
            [0, 1, NAN, -1],
        ],
    }
    check_equal(query, ["A", "A", "B", "A"], ["A", "B"], desired)
