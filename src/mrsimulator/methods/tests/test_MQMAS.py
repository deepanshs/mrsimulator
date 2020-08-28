# -*- coding: utf-8 -*-
import pytest
from mrsimulator.method.transition_query import TransitionQuery
from mrsimulator.methods import ThreeQ_MAS


def test_3QMAS_rotor_freq():
    error = " `rotor_frequency` cannot be modified for ThreeQ_MAS method."
    with pytest.raises(AttributeError, match=f".*{error}.*"):
        ThreeQ_MAS(channels=["87Rb"], rotor_frequency=10, spectral_dimensions=[{}, {}])


def test_3QMAS_rotor_amgle():
    error = "rotor_angle` is fixed to the magic-angle and cannot be modified"
    with pytest.raises(AttributeError, match=f".*{error}.*"):
        ThreeQ_MAS(channels=["87Rb"], rotor_angle=0, spectral_dimensions=[{}, {}])


def test_3QMAS_fractions():
    sites = ["87Rb", "27Al"]
    spins = [1.5, 2.5]
    k_MQ_MAS = {
        3: {1.5: 21 / 27, 2.5: 114 / 72, 3.5: 303 / 135, 4.5: 546 / 216},
        5: {2.5: 150 / 72, 3.5: 165 / 135, 4.5: 570 / 216},
        7: {3.5: 483 / 135, 4.5: 84 / 216},
        9: {4.5: 1116 / 216},
    }
    for i, isotope in zip(spins, sites):
        meth = ThreeQ_MAS(channels=[isotope])
        k = k_MQ_MAS[3][i]
        assert meth.spectral_dimensions[0].events[0].fraction == 1 / (1 + k)
        assert meth.spectral_dimensions[0].events[1].fraction == k / (1 + k)
        assert meth.spectral_dimensions[1].events[0].fraction == 1


def test_3QMAS_general():
    """MQMAS method declaration"""
    mth = ThreeQ_MAS(
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
    assert mth.name == "ThreeQ_MAS"

    assert mth.description == "Simulate a 3Q magic-angle spinning spectrum."

    assert mth.spectral_dimensions[0].events[0].transition_query == TransitionQuery(
        P={"channel-1": [[-3]]}, D={"channel-1": [[0]]},
    )

    assert mth.spectral_dimensions[0].events[1].transition_query == TransitionQuery(
        P={"channel-1": [[-1]]}, D={"channel-1": [[0]]},
    )

    assert mth.spectral_dimensions[1].events[0].transition_query == TransitionQuery(
        P={"channel-1": [[-1]]}, D={"channel-1": [[0]]},
    )
