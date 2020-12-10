# -*- coding: utf-8 -*-
import numpy as np
import pytest
from mrsimulator.method import Method
from mrsimulator.method.transition_query import TransitionQuery
from mrsimulator.methods import FiveQ_VAS
from mrsimulator.methods import SevenQ_VAS
from mrsimulator.methods import ThreeQ_VAS

__author__ = "Deepansh J. Srivastava"
__email__ = "srivastava.89@osu.edu"

methods = [ThreeQ_VAS, FiveQ_VAS, SevenQ_VAS]
names = ["ThreeQ_VAS", "FiveQ_VAS", "SevenQ_VAS"]


def test_MQ_VAS_rotor_freq():
    def error(name):
        return f"`rotor_frequency` attribute cannot be modified for {name} method."

    for name, method in zip(names, methods):
        e = error(name)
        with pytest.raises(AttributeError, match=f".*{e}.*"):
            method(rotor_frequency=10, spectral_dimensions=[{}, {}])


def test_MQ_VAS_setting_transition_query():
    def error(name):
        return f"`transition_query` attribute cannot be modified for {name} method."

    for name, method in zip(names, methods):
        e = error(name)
        with pytest.raises(AttributeError, match=f".*{e}.*"):
            method(
                spectral_dimensions=[{"events": [{"transition_query": {"P": [-1]}}]}],
            )


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
    assert mth.spectral_dimensions[0].events[0].transition_query == TransitionQuery(
        P={"channel-1": [[-3]]}, D={"channel-1": [[0]]}
    )
    assert mth.spectral_dimensions[1].events[0].transition_query == TransitionQuery(
        P={"channel-1": [[-1]]}, D={"channel-1": [[0]]}
    )
    assert Method.parse_dict_with_units(mth.json()) == mth


def test_5Q_VAS_general():
    """5Q-VAS method test"""
    mth = FiveQ_VAS(channels=["17O"], spectral_dimensions=[{}, {}])

    assert mth.name == "FiveQ_VAS"
    assert mth.description == "Simulate a 5Q variable-angle spinning spectrum."
    assert mth.spectral_dimensions[0].events[0].transition_query == TransitionQuery(
        P={"channel-1": [[-5]]}, D={"channel-1": [[0]]}
    )
    assert mth.spectral_dimensions[1].events[0].transition_query == TransitionQuery(
        P={"channel-1": [[-1]]}, D={"channel-1": [[0]]}
    )
    assert Method.parse_dict_with_units(mth.json()) == mth


def test_7Q_VAS_general():
    """7Q-VAS method test"""
    mth = SevenQ_VAS(channels=["51V"], spectral_dimensions=[{}, {}])

    assert mth.name == "SevenQ_VAS"
    assert mth.description == "Simulate a 7Q variable-angle spinning spectrum."
    assert mth.spectral_dimensions[0].events[0].transition_query == TransitionQuery(
        P={"channel-1": [[-7]]}, D={"channel-1": [[0]]}
    )
    assert mth.spectral_dimensions[1].events[0].transition_query == TransitionQuery(
        P={"channel-1": [[-1]]}, D={"channel-1": [[0]]}
    )
    assert Method.parse_dict_with_units(mth.json()) == mth
