# -*- coding: utf-8 -*-
import numpy as np
import pytest
from mrsimulator.method.transition_query import TransitionQuery
from mrsimulator.methods import FiveQ_VAS
from mrsimulator.methods import ThreeQ_VAS


def test_VAS_rotor_freq():
    error = " `rotor_frequency` cannot be modified for ThreeQ_VAS class."
    with pytest.raises(AttributeError, match=f".*{error}.*"):
        ThreeQ_VAS(channels=["87Rb"], rotor_frequency=10, spectral_dimensions=[{}, {}])

    error = " `rotor_frequency` cannot be modified for FiveQ_VAS class."
    with pytest.raises(AttributeError, match=f".*{error}.*"):
        FiveQ_VAS(channels=["87Rb"], rotor_frequency=10, spectral_dimensions=[{}, {}])


def test_VAS_setting_transition_query():
    error = "`transition_query` attribute cannot be modified for ThreeQ_VAS class."
    with pytest.raises(AttributeError, match=f".*{error}.*"):
        ThreeQ_VAS(
            channels=["87Rb"],
            spectral_dimensions=[{"events": [{"transition_query": {"P": [-1]}}]}, {}],
        )

    error = "`transition_query` attribute cannot be modified for FiveQ_VAS class."
    with pytest.raises(AttributeError, match=f".*{error}.*"):
        FiveQ_VAS(
            channels=["87Rb"],
            spectral_dimensions=[{"events": [{"transition_query": {"P": [-1]}}]}, {}],
        )


def test_3Q_VAS_fractions():
    sites = ["87Rb", "27Al"]
    spins = [1.5, 2.5]
    k_MQ_MAS = {
        3: {1.5: 21 / 27, 2.5: 114 / 72, 3.5: 303 / 135, 4.5: 546 / 216},
        5: {2.5: 150 / 72, 3.5: 165 / 135, 4.5: 570 / 216},
        7: {3.5: 483 / 135, 4.5: 84 / 216},
        9: {4.5: 1116 / 216},
    }
    for i, isotope in zip(spins, sites):
        meth = ThreeQ_VAS(channels=[isotope])
        k = k_MQ_MAS[3][i]
        assert meth.spectral_dimensions[0].events[0].fraction == 1
        assert meth.spectral_dimensions[1].events[0].fraction == 1
        assert np.allclose(meth.affine_matrix, [1 / (k + 1), k / (k + 1), 0, 1])


def test_3Q_VAS_general():
    """MQMAS method declaration"""
    mth = ThreeQ_VAS(
        channels=["87Rb"],
        magnetic_flux_density=9.4,  # in T
        # rotor_angle=54.735 * np.pi / 180,
        spectral_dimensions=[
            {
                "count": 1024,
                "spectral_width": 5e4,  # in Hz
                "reference_offset": 0,  # in Hz
            },
            {
                "count": 1024,
                "spectral_width": 5e4,  # in Hz
                "reference_offset": 0,  # in Hz
            },
        ],
    )
    assert mth.name == "ThreeQ_VAS"

    assert mth.description == "Simulate a 3Q variable-angle spinning spectrum."

    assert mth.spectral_dimensions[0].events[0].transition_query == TransitionQuery(
        P={"channel-1": [[-3]]}, D={"channel-1": [[0]]},
    )
    assert mth.spectral_dimensions[1].events[0].transition_query == TransitionQuery(
        P={"channel-1": [[-1]]}, D={"channel-1": [[0]]},
    )


def test_5Q_VAS_general():
    """Five quantum variable-angle spinning method"""
    mth = FiveQ_VAS(
        channels=["17O"],
        magnetic_flux_density=9.4,  # in T
        # rotor_angle=54.735 * np.pi / 180,
        spectral_dimensions=[
            {
                "count": 1024,
                "spectral_width": 5e4,  # in Hz
                "reference_offset": 0,  # in Hz
            },
            {
                "count": 1024,
                "spectral_width": 5e4,  # in Hz
                "reference_offset": 0,  # in Hz
            },
        ],
    )
    assert mth.name == "FiveQ_VAS"

    assert mth.description == "Simulate a 5Q variable-angle spinning spectrum."

    assert mth.spectral_dimensions[0].events[0].transition_query == TransitionQuery(
        P={"channel-1": [[-5]]}, D={"channel-1": [[0]]},
    )
    assert mth.spectral_dimensions[1].events[0].transition_query == TransitionQuery(
        P={"channel-1": [[-1]]}, D={"channel-1": [[0]]},
    )
