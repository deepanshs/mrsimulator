import numpy as np
import pytest
from mrsimulator.utils.extra import _check_dict_value_lengths
from mrsimulator.utils.extra import broadcast_dict_values
from mrsimulator.utils.extra import zip_dict


__author__ = ["Matthew D. Giammar"]
__email__ = ["giammar.7@buckeyemail.osu.edu"]


def test_broadcast_dict_values():
    d = {
        "foo": [1, 2, 3, 4],
        "bar": [5, 6, 7, 8],
        "baz": 9,
    }

    test_d = broadcast_dict_values(d)
    assert np.array_equal(test_d["foo"], np.asarray([1, 2, 3, 4]))
    assert np.array_equal(test_d["bar"], np.asarray([5, 6, 7, 8]))
    assert np.array_equal(test_d["baz"], np.asarray([9, 9, 9, 9]))

    e = ".*Length of lists in passed dict did not match `n_sites`.*"
    with pytest.raises(ValueError, match=e):
        broadcast_dict_values(d, 10)


def test_zip_dict():
    d = {
        "foo": [1, 2, 3, 4],
        "bar": [5, 6, 7, 8],
        "baz": 9,
    }

    test_lst = zip_dict(d)
    assert test_lst == [
        {"foo": 1, "bar": 5, "baz": 9},
        {"foo": 2, "bar": 6, "baz": 9},
        {"foo": 3, "bar": 7, "baz": 9},
        {"foo": 4, "bar": 8, "baz": 9},
    ]


def test_check_dict_value_lengths():
    d = {
        "foo": 1,
        "bar": 5,
        "baz": 9,
    }
    assert _check_dict_value_lengths(d) == 1

    d = {
        "foo": [1, 2, 3, 4],
        "bar": [5, 6, 7, 8],
        "baz": 9,
    }
    assert _check_dict_value_lengths(d) == 4

    d["foo"] = [1] * 10
    e = "Values with multiple lengths passed"
    with pytest.raises(ValueError, match=e):
        _check_dict_value_lengths(d)
