# -*- coding: utf-8 -*-
import numpy as np
import pytest
from mrsimulator.method.transition_query import TransitionQuery
from mrsimulator.methods import SSB2D

__author__ = "Deepansh J. Srivastava"
__email__ = "srivastava.89@osu.edu"


def test_SSB_rotor_freq():
    e = "`rotor_frequency` cannot be zero for SSB2D method."
    with pytest.raises(ValueError, match=f".*{e}.*"):
        SSB2D(spectral_dimensions=[{}, {}])

    with pytest.raises(ValueError, match=f".*{e}.*"):
        SSB2D(rotor_frequency=0, spectral_dimensions=[{}, {}])


def test_SSB_setting_events():
    e = "`events` attribute cannot be modified for SSB2D method."
    with pytest.raises(AttributeError, match=f".*{e}.*"):
        SSB2D(spectral_dimensions=[{"events": [{}]}])


def test_SSB_affine():
    meth = SSB2D(channels=["13C"], rotor_frequency=1200)
    np.allclose(meth.affine_matrix, [1, -1, 0, 0])


def test_SSB_general():
    """Inner satellite-transition variable-angle spinning method"""
    mth = SSB2D(
        channels=["87Rb"],
        magnetic_flux_density=9.4,  # in T
        rotor_frequency=1200,
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
    assert mth.name == "SSB2D"

    assert mth.description == "Simulate a 2D sideband separation method."

    # test transition query
    tq = TransitionQuery(P={"channel-1": [[-1]]}, D={"channel-1": [[0]]})
    assert mth.spectral_dimensions[0].events[0].transition_query == tq
    assert mth.spectral_dimensions[1].events[0].transition_query == tq

    # test rotor_frequency
    assert mth.spectral_dimensions[0].events[0].rotor_frequency == 1200
    assert mth.spectral_dimensions[1].events[0].rotor_frequency == 1e9
