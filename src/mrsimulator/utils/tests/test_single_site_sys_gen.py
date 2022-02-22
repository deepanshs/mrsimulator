# -*- coding: utf-8 -*-
import numpy as np
import pytest
from mrsimulator.utils.collection import _check_lengths_of_args
from mrsimulator.utils.collection import _flatten_item
from mrsimulator.utils.collection import single_site_system_generator
from mrsimulator.utils.collection import site_generator


__author__ = ["Deepansh Srivastava", "Matthew D. Giammar"]
__email__ = ["srivastava.89@osu.edu", "giammar.7@buckeyemail.osu.edu"]


def test_flatten_item():
    # Check scalar value
    item = 1

    assert _flatten_item(item) == 1

    # Check 1d list
    lst = [0, 1, 2, 3, 4]
    flat_lst = _flatten_item(lst)

    assert isinstance(flat_lst, np.ndarray)
    assert np.array_equal(flat_lst, np.array([0, 1, 2, 3, 4]))

    # Check list flattenning
    lst = [
        [0, 1, 2, 3],
        [4, 5, 6, 7],
        [8, 9, 10, 11],
    ]
    flat_lst = _flatten_item(lst)

    assert isinstance(flat_lst, np.ndarray)
    assert np.array_equal(flat_lst, np.arange(12))

    # Check passing dicts
    d = {"scalar": 0, "flat": np.arange(50), "3d": np.arange(50).reshape((5, 5, 2))}
    flat_d = _flatten_item(d)

    assert flat_d["scalar"] == 0
    assert np.array_equal(flat_d["flat"], np.arange(50))
    assert np.array_equal(flat_d["3d"], np.arange(50))


def test_unbalanced_lists():
    isotopes = ["13C", "71Ga", "15N"]
    shifts = np.arange(10)

    error = ".*Not all arrays/lists passed were of the same length.*"
    with pytest.raises(ValueError, match=error):
        single_site_system_generator(isotope=isotopes, isotropic_chemical_shift=shifts)

    shielding = {"zeta": np.arange(10)}

    error = ".*Not all arrays/lists passed were of the same length.*"
    with pytest.raises(ValueError, match=error):
        single_site_system_generator(isotope=isotopes, shielding_symmetric=shielding)

    # Calling function directly
    assert _check_lengths_of_args(isotopes, [1, 2, 3], "foo") == 3
    error = ".*Not all arrays/lists passed were of the same length.*"
    with pytest.raises(ValueError, match=error):
        _check_lengths_of_args(isotopes, shielding)


def test_site_generator():
    iso_dist = np.random.normal(0, 10, 10)
    zeta_dist = np.arange(10)
    eta_dist = np.random.rand(10)
    gamma_dist = np.random.rand(10) * 3.1415
    sites = site_generator(
        isotope="13C",
        isotropic_chemical_shift=iso_dist,
        shielding_symmetric={"zeta": zeta_dist, "eta": eta_dist, "gamma": gamma_dist},
    )

    assert isinstance(sites, list)

    for i in range(10):
        assert sites[i].isotope.symbol == "13C"
        assert sites[i].isotropic_chemical_shift == iso_dist[i]
        assert sites[i].shielding_symmetric.zeta == zeta_dist[i]
        assert sites[i].shielding_symmetric.eta == eta_dist[i]
        assert sites[i].shielding_symmetric.alpha is None
        assert sites[i].shielding_symmetric.beta is None
        assert sites[i].shielding_symmetric.gamma == gamma_dist[i]
        assert sites[i].quadrupolar is None


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


def test_shielding_06():
    iso_dist = np.random.normal(0, 10, 10)
    zeta_dist = np.arange(10)
    alpha_dist = np.random.rand(10)
    beta_dist = np.random.rand(10) * 3.1415
    sys = single_site_system_generator(
        isotope="13C",
        isotropic_chemical_shift=iso_dist,
        shielding_antisymmetric={
            "zeta": zeta_dist,
            "alpha": alpha_dist,
            "beta": beta_dist,
        },
    )

    for i in range(10):
        assert sys[i].sites[0].isotope.symbol == "13C"
        assert sys[i].sites[0].isotropic_chemical_shift == iso_dist[i]
        assert sys[i].sites[0].shielding_antisymmetric.zeta == zeta_dist[i]
        assert sys[i].sites[0].shielding_antisymmetric.alpha == alpha_dist[i]
        assert sys[i].sites[0].shielding_antisymmetric.beta == beta_dist[i]
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
    abundances = [0] * 10

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
