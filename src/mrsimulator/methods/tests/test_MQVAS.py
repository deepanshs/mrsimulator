# -*- coding: utf-8 -*-
import numpy as np
import pytest
from mrsimulator.method.query import TransitionQuery
from mrsimulator.methods import FiveQ_VAS
from mrsimulator.methods import SevenQ_VAS
from mrsimulator.methods import ThreeQ_VAS

__author__ = "Deepansh J. Srivastava"
__email__ = "srivastava.89@osu.edu"

methods = [ThreeQ_VAS, FiveQ_VAS, SevenQ_VAS]
names = ["ThreeQ_VAS", "FiveQ_VAS", "SevenQ_VAS"]


def sample_test_output(n):
    return {
        "magnetic_flux_density": "9.4 T",
        "rotor_angle": "0.9553166181245 rad",
        "rotor_frequency": "1000000000000.0 Hz",
        "spectral_dimensions": [
            {
                "count": 1024,
                "spectral_width": "25000.0 Hz",
                "events": [{"transition_query": [{"ch1": {"P": [n], "D": [0]}}]}],
            },
            {
                "count": 1024,
                "spectral_width": "25000.0 Hz",
                "events": [{"transition_query": [{"ch1": {"P": [-1], "D": [0]}}]}],
            },
        ],
    }


def test_MQ_VAS_rotor_freq():
    e = "`rotor_frequency=1e12 Hz` is fixed for all NamedMethod"
    isotopes = ["87Rb", "27Al", "51V"]
    for iso, method in zip(isotopes, methods):
        with pytest.raises(ValueError, match=f".*{e}.*"):
            method(channels=[iso], rotor_frequency=10, spectral_dimensions=[{}, {}])


def test_MQ_VAS_affine():
    sites = ["87Rb", "27Al", "51V"]
    spins = [1.5, 2.5, 3.5]
    k_MQ_MAS = {
        3: {1.5: 21 / 27, 2.5: 114 / 72, 3.5: 303 / 135, 4.5: 546 / 216},
        5: {2.5: 150 / 72, 3.5: 165 / 135, 4.5: 570 / 216},
        7: {3.5: 483 / 135, 4.5: 84 / 216},
        9: {4.5: 1116 / 216},
    }
    for j, method in enumerate(methods):
        for i, isotope in zip(spins[j:], sites[j:]):
            meth = method(channels=[isotope])
            k = k_MQ_MAS[3 + 2 * j][i]
            assert meth.spectral_dimensions[0].events[0].fraction == 1
            assert meth.spectral_dimensions[1].events[0].fraction == 1
            assert np.allclose(meth.affine_matrix, [1 / (k + 1), k / (k + 1), 0, 1])


def test_3Q_VAS_general():
    """3Q-VAS method test"""
    mth = ThreeQ_VAS(channels=["87Rb"], spectral_dimensions=[{}, {}])
    assert mth.name == "ThreeQ_VAS"
    assert mth.description == "Simulate a 3Q variable-angle spinning spectrum."
    assert mth.spectral_dimensions[0].events[0].transition_query == [
        TransitionQuery(ch1={"P": [-3], "D": [0]})
    ]
    assert mth.spectral_dimensions[1].events[0].transition_query == [
        TransitionQuery(ch1={"P": [-1], "D": [0]})
    ]
    assert ThreeQ_VAS.parse_dict_with_units(mth.json()) == mth

    assert np.allclose(mth.affine_matrix, [0.5625, 0.4375, 0.0, 1.0])

    serialize = mth.json()
    _ = serialize.pop("affine_matrix")
    assert serialize == {
        "channels": ["87Rb"],
        "description": "Simulate a 3Q variable-angle spinning spectrum.",
        "name": "ThreeQ_VAS",
        **sample_test_output(-3),
    }


def test_5Q_VAS_general():
    """5Q-VAS method test"""
    mth = FiveQ_VAS(channels=["17O"], spectral_dimensions=[{}, {}])

    assert mth.name == "FiveQ_VAS"
    assert mth.description == "Simulate a 5Q variable-angle spinning spectrum."
    assert mth.spectral_dimensions[0].events[0].transition_query == [
        TransitionQuery(ch1={"P": [-5], "D": [0]})
    ]
    assert mth.spectral_dimensions[1].events[0].transition_query == [
        TransitionQuery(ch1={"P": [-1], "D": [0]})
    ]
    assert FiveQ_VAS.parse_dict_with_units(mth.json()) == mth

    assert np.allclose(
        mth.affine_matrix, [0.3243243243243243, 0.6756756756756757, 0.0, 1.0]
    )

    serialize = mth.json()
    _ = serialize.pop("affine_matrix")

    assert serialize == {
        "channels": ["17O"],
        "description": "Simulate a 5Q variable-angle spinning spectrum.",
        "name": "FiveQ_VAS",
        **sample_test_output(-5),
    }


def test_7Q_VAS_general():
    """7Q-VAS method test"""
    mth = SevenQ_VAS(channels=["51V"], spectral_dimensions=[{}, {}])

    assert mth.name == "SevenQ_VAS"
    assert mth.description == "Simulate a 7Q variable-angle spinning spectrum."
    assert mth.spectral_dimensions[0].events[0].transition_query == [
        TransitionQuery(ch1={"P": [-7], "D": [0]})
    ]
    assert mth.spectral_dimensions[1].events[0].transition_query == [
        TransitionQuery(ch1={"P": [-1], "D": [0]})
    ]
    assert SevenQ_VAS.parse_dict_with_units(mth.json()) == mth

    assert np.allclose(mth.affine_matrix, [0.2184466, 0.7815534, 0.0, 1.0])

    serialize = mth.json()
    _ = serialize.pop("affine_matrix")

    assert serialize == {
        "channels": ["51V"],
        "description": "Simulate a 7Q variable-angle spinning spectrum.",
        "name": "SevenQ_VAS",
        **sample_test_output(-7),
    }
