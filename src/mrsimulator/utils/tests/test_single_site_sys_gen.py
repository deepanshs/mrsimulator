# -*- coding: utf-8 -*-
import numpy as np
import pytest
from mrsimulator.utils.collection import _check_lengths
from mrsimulator.utils.collection import _extend_dict_values
from mrsimulator.utils.collection import _extend_to_nparray
from mrsimulator.utils.collection import _fix_item
from mrsimulator.utils.collection import _zip_dict
from mrsimulator.utils.collection import single_site_system_generator

# from mrsimulator.utils.collection import generate_site_list


__author__ = ["Deepansh Srivastava", "Matthew D. Giammar"]
__email__ = ["srivastava.89@osu.edu", "giammar.7@buckeyemail.osu.edu"]


# def test_get_default_lists():
#     a = 2.0
#     assert np.all(_get_default_lists(a, 5) == np.asarray([2.0] * 5))

#     a = "32G"
#     assert np.all(_get_default_lists(a, 5) == np.asarray(["32G"] * 5))


# def test_get_length():
#     a = 2.0
#     assert _get_length(a) == 0

#     a = [1, 2, 3, 4, 5]
#     assert _get_length(a) == 5

#     a = np.arange(43)
#     assert _get_length(a) == 43


def test_fix_item():
    assert _fix_item("str") == "str"

    item = [1, 2, 3]

    assert np.array_equal(_fix_item(item), np.array([1, 2, 3]))

    item = np.array(item)

    # 3d array
    item = np.array(
        [
            [
                [-12, -11, -10, -9],
                [-8, -7, -6, -5],
                [-4, -3, -2, -1],
            ],
            [
                [0, 1, 2, 3],
                [4, 5, 6, 7],
                [8, 9, 10, 11],
            ],
        ]
    )
    check = np.arange(-12, 12, 1)

    assert np.array_equal(_fix_item(item), check)


def test_zip_dict():
    dictionary = {
        "key1": [0, 1, 2, 3, 4],
        "key2": [5, 6, 7, 8, 9],
        "key3": [10, 11, 12, 13, 14],
        "key4": [15, 16, 17, 18, 19],
    }
    zipped = _zip_dict(dictionary)
    for i, row in enumerate(zipped):
        assert row == {"key1": i, "key2": 5 + i, "key3": 10 + i, "key4": 15 + i}

    dictionary["key2"] = [None] * 5
    zipped = _zip_dict(dictionary)
    for i, row in enumerate(zipped):
        assert row == {"key1": i, "key2": None, "key3": 10 + i, "key4": 15 + i}

    dictionary = {
        "key1": [0, None, 2, 3, 4],
        "key2": [5, None, 7, 8, 9],
        "key3": [10, None, 12, 13, 14],
        "key4": [15, None, 17, 18, 19],
    }
    zipped = _zip_dict(dictionary)
    for i, row in enumerate(zipped):
        if i == 1:
            assert row is None
        else:
            assert row == {"key1": i, "key2": 5 + i, "key3": 10 + i, "key4": 15 + i}


def test_extend_to_nparray():
    item = ["foo", "bar", "baz", "spam", "ham", "eggs"]
    extended = _extend_to_nparray(item, 10)

    assert isinstance(extended, np.ndarray)
    assert extended.ndim == 1
    assert extended.size == 6

    item = np.array(item)
    extended = _extend_to_nparray(item, 10)

    assert isinstance(extended, np.ndarray)
    assert extended.size == 6

    extended = _extend_to_nparray("foo", 10)

    assert isinstance(extended, np.ndarray)
    assert extended.size == 10

    extended_none = _extend_to_nparray(None, 15)

    assert isinstance(extended_none, np.ndarray)
    assert extended_none.size == 15


def test_extend_dict_values():
    def numpy_dict_equality(d1, d2):
        for key1, val1 in d1.items():
            if key1 not in d2:
                return False
            val2 = d2[key1]
            if not np.all(val1 == val2):
                return False
        return True

    _dict = {"key1": 1, "key2": 2, "key3": 3}

    # Single length dict remains unchanged with n_sites of 1
    assert numpy_dict_equality(_extend_dict_values(_dict, 1)[0], _dict)
    assert numpy_dict_equality(_extend_dict_values(_dict, 5)[0], _dict)

    _dict = {"key1": [1], "key2": [2], "key3": [3]}
    check_dict = {"key1": 1, "key2": 2, "key3": 3}

    # Single length dict remains unchanged with n_sites of 1
    assert numpy_dict_equality(_extend_dict_values(_dict, 1)[0], check_dict)
    assert numpy_dict_equality(_extend_dict_values(_dict, 5)[0], check_dict)

    _dict = {"key1": [1] * 5, "key2": 2, "key3": 3}
    check_list = [{"key1": 1, "key2": 2, "key3": 3}] * 5

    assert _extend_dict_values(_dict, 1) == (check_list, 5)
    assert _extend_dict_values(_dict, 5) == (check_list, 5)

    error = ".*A list in a dictionary was misshapen.*"
    with pytest.raises(ValueError, match=error):
        _extend_dict_values(_dict, 4)

    _dict = {"key1": [1] * 5, "key2": [2] * 4, "key3": 3}

    error = ".*An array or list was either too short or too long.*"
    with pytest.raises(ValueError, match=error):
        _extend_dict_values(_dict, 1)


def test_check_lengths():
    lists = [[0]] * 5

    assert _check_lengths(lists) == 1

    lists.append([2] * 7)

    assert _check_lengths(lists) == 7

    lists = [[i for i in range(10)] for _ in range(5)]

    assert _check_lengths(lists) == 10

    lists = [np.zeros(10) for _ in range(5)]

    assert _check_lengths(lists) == 10

    lists = [[0], [0], [0], [3, 3, 3], [4, 4, 4, 4]]

    error = ".*An array or list was either too short or too long.*"
    with pytest.raises(ValueError, match=error):
        _check_lengths(lists)

    lists = [np.zeros(10) for _ in range(5)] + [np.ones(9)]

    with pytest.raises(ValueError, match=error):
        _check_lengths(lists)


def test_unbalanced_lists():
    isotopes = ["13C", "71Ga", "15N"]
    shifts = np.arange(10)

    # error = ".*Each entry can either be a single item or a list of items.*"
    error = ".*An array or list was either too short or too long.*"
    with pytest.raises(ValueError, match=error):
        single_site_system_generator(isotope=isotopes, isotropic_chemical_shift=shifts)


def test_shielding_01():
    isotopes = ["13C", "71Ga", "15N", "14N", "27Al", "29Si", "1H", "17O", "33S", "31P"]

    sys = single_site_system_generator(isotope=isotopes)
    for i in range(10):
        assert sys[i].sites[0].isotope.symbol == isotopes[i]
        assert sys[i].sites[0].isotropic_chemical_shift == 0
        assert sys[i].sites[0].shielding_symmetric is None
        assert sys[i].sites[0].quadrupolar is None


def test_shielding_02():
    isotropics = np.arange(20)
    sys = single_site_system_generator(
        isotope="71Ga", isotropic_chemical_shift=isotropics
    )

    for i in range(20):
        assert sys[i].sites[0].isotope.symbol == "71Ga"
        assert sys[i].sites[0].isotropic_chemical_shift == isotropics[i]
        assert sys[i].sites[0].shielding_symmetric is None
        assert sys[i].sites[0].quadrupolar is None


def test_shielding_03():
    zeta_dist = np.arange(10)
    eta_dist = np.ones(10) * 0.5
    sys = single_site_system_generator(
        isotope="13C", shielding_symmetric={"zeta": zeta_dist, "eta": eta_dist}
    )

    for i in range(10):
        assert sys[i].sites[0].isotope.symbol == "13C"
        assert sys[i].sites[0].isotropic_chemical_shift == 0
        assert sys[i].sites[0].shielding_symmetric.zeta == zeta_dist[i]
        assert sys[i].sites[0].shielding_symmetric.eta == eta_dist[i]
        assert sys[i].sites[0].shielding_symmetric.alpha is None
        assert sys[i].sites[0].shielding_symmetric.beta is None
        assert sys[i].sites[0].shielding_symmetric.gamma is None
        assert sys[i].sites[0].quadrupolar is None


def test_shielding_04():
    zeta_dist = np.arange(10)
    eta_dist = np.random.rand(10)
    beta_dist = np.random.rand(10) * 3.1415
    sys = single_site_system_generator(
        isotope="13C",
        shielding_symmetric={"zeta": zeta_dist, "eta": eta_dist, "beta": beta_dist},
    )

    for i in range(10):
        assert sys[i].sites[0].isotope.symbol == "13C"
        assert sys[i].sites[0].isotropic_chemical_shift == 0
        assert sys[i].sites[0].shielding_symmetric.zeta == zeta_dist[i]
        assert sys[i].sites[0].shielding_symmetric.eta == eta_dist[i]
        assert sys[i].sites[0].shielding_symmetric.alpha is None
        assert sys[i].sites[0].shielding_symmetric.beta == beta_dist[i]
        assert sys[i].sites[0].shielding_symmetric.gamma is None
        assert sys[i].sites[0].quadrupolar is None


def test_shielding_05():
    iso_dist = np.random.normal(0, 10, 10)
    zeta_dist = np.arange(10)
    eta_dist = np.random.rand(10)
    gamma_dist = np.random.rand(10) * 3.1415
    sys = single_site_system_generator(
        isotope="13C",
        isotropic_chemical_shift=iso_dist,
        shielding_symmetric={"zeta": zeta_dist, "eta": eta_dist, "gamma": gamma_dist},
    )

    for i in range(10):
        assert sys[i].sites[0].isotope.symbol == "13C"
        assert sys[i].sites[0].isotropic_chemical_shift == iso_dist[i]
        assert sys[i].sites[0].shielding_symmetric.zeta == zeta_dist[i]
        assert sys[i].sites[0].shielding_symmetric.eta == eta_dist[i]
        assert sys[i].sites[0].shielding_symmetric.alpha is None
        assert sys[i].sites[0].shielding_symmetric.beta is None
        assert sys[i].sites[0].shielding_symmetric.gamma == gamma_dist[i]
        assert sys[i].sites[0].quadrupolar is None


def assertion_quad(sys, isotope, iso_dist, Cq_dist, eta_dist, abundances):
    for i in range(Cq_dist.size):
        assert sys[i].sites[0].isotope.symbol == isotope
        assert sys[i].sites[0].isotropic_chemical_shift == iso_dist[i]
        assert sys[i].sites[0].shielding_symmetric is None
        assert sys[i].sites[0].quadrupolar.Cq == Cq_dist[i]
        assert sys[i].sites[0].quadrupolar.eta == eta_dist[i]
        assert sys[i].sites[0].quadrupolar.alpha is None
        assert sys[i].sites[0].quadrupolar.beta is None
        assert sys[i].sites[0].quadrupolar.gamma is None
        assert sys[i].abundance == abundances


def test_quad_01():
    Cq_dist = np.arange(10)
    eta_dist = np.ones(10) * 0.5
    sys = single_site_system_generator(
        isotope="27Al", quadrupolar={"Cq": Cq_dist, "eta": eta_dist}
    )

    iso_dict = np.zeros(10)
    assertion_quad(sys, "27Al", iso_dict, Cq_dist, eta_dist, 0.1)


def test_quad_02():
    iso_dist = np.random.normal(0, 10, 10)
    Cq_dist = np.arange(10)
    eta_dist = np.random.rand(10)
    gamma_dist = np.random.rand(10) * 3.1415
    sys = single_site_system_generator(
        isotope="17O",
        isotropic_chemical_shift=iso_dist,
        quadrupolar={"Cq": Cq_dist, "eta": eta_dist, "gamma": gamma_dist},
    )

    for i in range(10):
        assert sys[i].sites[0].isotope.symbol == "17O"
        assert sys[i].sites[0].isotropic_chemical_shift == iso_dist[i]
        assert sys[i].sites[0].shielding_symmetric is None
        assert sys[i].sites[0].quadrupolar.Cq == Cq_dist[i]
        assert sys[i].sites[0].quadrupolar.eta == eta_dist[i]
        assert sys[i].sites[0].quadrupolar.alpha is None
        assert sys[i].sites[0].quadrupolar.beta is None
        assert sys[i].sites[0].quadrupolar.gamma == gamma_dist[i]


def test_quad_shield():
    iso_dist = np.random.normal(0, 10, 10)
    zeta_dist = np.random.rand(10)
    eta_dist_s = np.random.rand(10)
    Cq_dist = np.random.rand(10)
    eta_dist_q = np.random.rand(10)
    gamma_dist = np.random.rand(10) * 3.1415
    sys = single_site_system_generator(
        isotope="17O",
        isotropic_chemical_shift=iso_dist,
        shielding_symmetric={"zeta": zeta_dist, "eta": eta_dist_s},
        quadrupolar={"Cq": Cq_dist, "eta": eta_dist_q, "gamma": gamma_dist},
    )

    for i in range(10):
        assert sys[i].sites[0].isotope.symbol == "17O"
        assert sys[i].sites[0].isotropic_chemical_shift == iso_dist[i]

        assert sys[i].sites[0].shielding_symmetric.zeta == zeta_dist[i]
        assert sys[i].sites[0].shielding_symmetric.eta == eta_dist_s[i]
        assert sys[i].sites[0].shielding_symmetric.alpha is None
        assert sys[i].sites[0].shielding_symmetric.beta is None
        assert sys[i].sites[0].shielding_symmetric.gamma is None

        assert sys[i].sites[0].quadrupolar.Cq == Cq_dist[i]
        assert sys[i].sites[0].quadrupolar.eta == eta_dist_q[i]
        assert sys[i].sites[0].quadrupolar.alpha is None
        assert sys[i].sites[0].quadrupolar.beta is None
        assert sys[i].sites[0].quadrupolar.gamma == gamma_dist[i]


def test_abundance_01():
    Cq_dist = np.arange(10)
    eta_dist = np.ones(10) * 0.5
    abundances = 0.6
    sys = single_site_system_generator(
        isotope="27Al",
        quadrupolar={"Cq": Cq_dist, "eta": eta_dist},
        abundance=abundances,
    )

    iso_dict = np.zeros(10)
    assertion_quad(sys, "27Al", iso_dict, Cq_dist, eta_dist, 0.6)
    # for i in range(10):
    #     assert sys[i].sites[0].isotope.symbol == "27Al"
    #     assert sys[i].sites[0].isotropic_chemical_shift == 0
    #     assert sys[i].sites[0].shielding_symmetric is None
    #     assert sys[i].sites[0].quadrupolar.Cq == Cq_dist[i]
    #     assert sys[i].sites[0].quadrupolar.eta == eta_dist[i]
    #     assert sys[i].sites[0].quadrupolar.alpha is None
    #     assert sys[i].sites[0].quadrupolar.beta is None
    #     assert sys[i].sites[0].quadrupolar.gamma is None
    #     assert sys[i].abundance == 0.6


def test_abundance_02():
    Cq_dist = np.arange(10)
    eta_dist = np.ones(10) * 0.5
    abundances = np.zeros(10)
    sys = single_site_system_generator(
        isotope="27Al",
        quadrupolar={"Cq": Cq_dist, "eta": eta_dist},
        abundance=abundances,
    )
    assert sys == []


def test_abundance_03():
    Cq_dist = np.arange(10)
    eta_dist = np.ones(10) * 0.5
    gamma = np.random.rand(10)
    abundances = np.zeros(10)

    indexes = [2, 5, 8]
    for i in indexes:
        abundances[i] = 1

    sys = single_site_system_generator(
        isotope="27Al",
        shielding_symmetric={"gamma": gamma},
        quadrupolar={"Cq": Cq_dist, "eta": eta_dist},
        abundance=abundances,
    )
    assert len(sys) == 3

    for i, j in enumerate(indexes):
        assert sys[i].sites[0].isotope.symbol == "27Al"
        assert sys[i].sites[0].isotropic_chemical_shift == 0
        assert sys[i].sites[0].shielding_symmetric.zeta is None
        assert sys[i].sites[0].shielding_symmetric.eta is None
        assert sys[i].sites[0].shielding_symmetric.alpha is None
        assert sys[i].sites[0].shielding_symmetric.beta is None
        assert sys[i].sites[0].shielding_symmetric.gamma == gamma[j]
        assert sys[i].sites[0].quadrupolar.Cq == Cq_dist[j]
        assert sys[i].sites[0].quadrupolar.eta == eta_dist[j]
        assert sys[i].sites[0].quadrupolar.alpha is None
        assert sys[i].sites[0].quadrupolar.beta is None
        assert sys[i].sites[0].quadrupolar.gamma is None
        assert sys[i].abundance == abundances[j]
