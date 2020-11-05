# -*- coding: utf-8 -*-
import numpy as np
import pytest
from mrsimulator.utils.collection import _get_default_lists
from mrsimulator.utils.collection import _get_length
from mrsimulator.utils.collection import single_site_system_generator


def test_get_default_lists():
    a = 2.0
    assert np.all(_get_default_lists(a, 5) == np.asarray([2.0] * 5))

    a = "32G"
    assert np.all(_get_default_lists(a, 5) == np.asarray(["32G"] * 5))


def test_get_length():
    a = 2.0
    assert _get_length(a) == 0

    a = [1, 2, 3, 4, 5]
    assert _get_length(a) == 5

    a = np.arange(43)
    assert _get_length(a) == 43


def test_unbalanced_lists():
    isotopes = ["13C", "71Ga", "15N"]
    shifts = np.arange(10)

    error = ".*Each entry can either be a single item or a list of items.*"
    with pytest.raises(ValueError, match=error):
        single_site_system_generator(
            isotopes=isotopes, isotropic_chemical_shifts=shifts
        )


def test_shielding_01():
    isotopes = ["13C", "71Ga", "15N", "14N", "27Al", "29Si", "1H", "17O", "33S", "31P"]

    sys = single_site_system_generator(isotopes=isotopes)
    for i in range(10):
        assert sys[i].sites[0].isotope.symbol == isotopes[i]
        assert sys[i].sites[0].isotropic_chemical_shift == 0
        assert sys[i].sites[0].shielding_symmetric is None
        assert sys[i].sites[0].quadrupolar is None


def test_shielding_02():
    isotropics = np.arange(20)
    sys = single_site_system_generator(
        isotopes="71Ga", isotropic_chemical_shifts=isotropics
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
        isotopes="13C", shielding_symmetric={"zeta": zeta_dist, "eta": eta_dist}
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
        isotopes="13C",
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
        isotopes="13C",
        isotropic_chemical_shifts=iso_dist,
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


def test_quad_01():
    Cq_dist = np.arange(10)
    eta_dist = np.ones(10) * 0.5
    sys = single_site_system_generator(
        isotopes="27Al", quadrupolar={"Cq": Cq_dist, "eta": eta_dist}
    )

    for i in range(10):
        assert sys[i].sites[0].isotope.symbol == "27Al"
        assert sys[i].sites[0].isotropic_chemical_shift == 0
        assert sys[i].sites[0].shielding_symmetric is None
        assert sys[i].sites[0].quadrupolar.Cq == Cq_dist[i]
        assert sys[i].sites[0].quadrupolar.eta == eta_dist[i]
        assert sys[i].sites[0].quadrupolar.alpha is None
        assert sys[i].sites[0].quadrupolar.beta is None
        assert sys[i].sites[0].quadrupolar.gamma is None


def test_quad_02():
    iso_dist = np.random.normal(0, 10, 10)
    Cq_dist = np.arange(10)
    eta_dist = np.random.rand(10)
    gamma_dist = np.random.rand(10) * 3.1415
    sys = single_site_system_generator(
        isotopes="17O",
        isotropic_chemical_shifts=iso_dist,
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
        isotopes="17O",
        isotropic_chemical_shifts=iso_dist,
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
    abundance = 0.6
    sys = single_site_system_generator(
        isotopes="27Al",
        quadrupolar={"Cq": Cq_dist, "eta": eta_dist},
        abundance=abundance,
    )

    for i in range(10):
        assert sys[i].sites[0].isotope.symbol == "27Al"
        assert sys[i].sites[0].isotropic_chemical_shift == 0
        assert sys[i].sites[0].shielding_symmetric is None
        assert sys[i].sites[0].quadrupolar.Cq == Cq_dist[i]
        assert sys[i].sites[0].quadrupolar.eta == eta_dist[i]
        assert sys[i].sites[0].quadrupolar.alpha is None
        assert sys[i].sites[0].quadrupolar.beta is None
        assert sys[i].sites[0].quadrupolar.gamma is None
        assert sys[i].abundance == 0.6


def test_abundance_02():
    Cq_dist = np.arange(10)
    eta_dist = np.ones(10) * 0.5
    abundance = np.zeros(10)
    sys = single_site_system_generator(
        isotopes="27Al",
        quadrupolar={"Cq": Cq_dist, "eta": eta_dist},
        abundance=abundance,
    )
    assert sys == []


def test_abundance_03():
    Cq_dist = np.arange(10)
    eta_dist = np.ones(10) * 0.5
    gamma = np.random.rand(10)
    abundance = np.zeros(10)

    print(gamma)
    indexes = [2, 5, 8]
    for i in indexes:
        abundance[i] = 1

    sys = single_site_system_generator(
        isotopes="27Al",
        shielding_symmetric={"gamma": gamma},
        quadrupolar={"Cq": Cq_dist, "eta": eta_dist},
        abundance=abundance,
    )
    assert len(sys) == 3

    for i, j in enumerate(indexes):
        print(sys[i].sites[0].shielding_symmetric.gamma, gamma[j])
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
        assert sys[i].abundance == abundance[j]
