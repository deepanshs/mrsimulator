# -*- coding: utf-8 -*-
from os import path

import numpy as np
import pytest
from monty.serialization import loadfn
from mrsimulator.methods import BlochDecaySpectrum
from mrsimulator.methods import Method2D
from mrsimulator.methods import ThreeQ_VAS

__author__ = "Deepansh J. Srivastava"
__email__ = "srivastava.89@osu.edu"

MODULE_DIR = path.dirname(path.abspath(__file__))
TESTDATA = loadfn(path.join(MODULE_DIR, "test_data.json"))


def test_more_spectral_dimensions():
    error = "The method allows 1 spectral dimension"
    with pytest.raises(ValueError, match=f".*{error}.*"):
        BlochDecaySpectrum(spectral_dimensions=[{}, {}])


def test_01():
    error = "method requires exactly 1 channel"
    with pytest.raises(ValueError, match=f".*{error}.*"):
        BlochDecaySpectrum(channels=["1H", "29Si"])


def test_02():
    e = "`rotor_frequency` attribute cannot be modified for Method2D method."
    with pytest.raises(AttributeError, match=f".*{e}.*"):
        Method2D(channels=["87Rb"], rotor_frequency=10, spectral_dimensions=[{}, {}])


def test_03():
    """generic method declaration"""
    mth = Method2D(
        channels=["87Rb"],
        magnetic_flux_density=9.4,  # in T
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

    assert TESTDATA["generic"] == mth.to_dict_with_units()


def test_04():
    """SAS method declaration"""
    mth = Method2D(
        channels=["87Rb"],
        magnetic_flux_density=9.4,  # in T
        spectral_dimensions=[
            {
                "count": 512,
                "spectral_width": 5e4,  # in Hz
                "reference_offset": 0,  # in Hz
                "events": [
                    {
                        "rotor_angle": 70.12 * np.pi / 180,
                        "transition_query": {"P": [-1], "D": [0]},
                    },
                ],
            },
            {
                "count": 128,
                "spectral_width": 5e4,  # in Hz
                "reference_offset": 0,  # in Hz
                "events": [
                    {
                        "rotor_angle": 54.735 * np.pi / 180,
                        "transition_query": {"P": [-1], "D": [0]},
                    },
                ],
            },
        ],
    )

    assert TESTDATA["SAS"] == mth.to_dict_with_units()


def test_05():
    """Satellite to central correlation method declaration"""
    mth = Method2D(
        channels=["87Rb"],
        magnetic_flux_density=9.4,  # in T
        spectral_dimensions=[
            {
                "count": 512,
                "spectral_width": 5e6,  # in Hz
                "reference_offset": 0,  # in Hz
                "events": [
                    {
                        "rotor_angle": 0 * np.pi / 180,
                        "transition_query": {"P": [-1], "D": [2, -2]},
                    },
                ],
            },
            {
                "count": 128,
                "spectral_width": 5e4,  # in Hz
                "reference_offset": 0,  # in Hz
                "events": [
                    {
                        "rotor_angle": 54.735 * np.pi / 180,
                        "transition_query": {"P": [-1], "D": [0]},
                    },
                ],
            },
        ],
    )

    assert TESTDATA["STMAS"] == mth.to_dict_with_units()


def test_3QMAS():
    """3Q MAS correlation method declaration"""
    mth = ThreeQ_VAS(
        channels=["87Rb"],
        magnetic_flux_density=9.4,  # in T
        spectral_dimensions=[
            {
                "count": 512,
                "spectral_width": 5e6,  # in Hz
                "reference_offset": 0,  # in Hz
            },
            {
                "count": 128,
                "spectral_width": 5e4,  # in Hz
                "reference_offset": 0,  # in Hz
            },
        ],
    )

    assert np.allclose(mth.affine_matrix, [0.5625, 0.4375, 0, 1])


def test_06():
    """Test with order"""
    test_03()
    test_05()
    test_04()


def test_methods():
    das = Method2D(
        channels=["87Rb"],
        magnetic_flux_density=4.2,  # in T
        spectral_dimensions=[
            {
                "count": 256,
                "spectral_width": 2e4,  # in Hz
                "reference_offset": -5e3,  # in Hz
                "label": "70.12 dimension",
                "events": [
                    {
                        "rotor_angle": 70.12 * 3.14159 / 180,  # in rads
                        "transition_query": {"P": [-1], "D": [0]},
                    },
                    {
                        "rotor_angle": 30.12 * 3.14159 / 180,  # in rads
                        "transition_query": {"P": [-1], "D": [0]},
                    },
                ],
            },
            # The last spectral dimension block is the direct-dimension
            {
                "count": 256,
                "spectral_width": 3e4,  # in Hz
                "reference_offset": -7e3,  # in Hz
                "label": "MAS dimension",
                "events": [
                    {
                        "rotor_angle": 54.735 * 3.14159 / 180,  # in rads
                        "transition_query": {"P": [-1], "D": [0]},
                    }
                ],
            },
        ],
    )

    assert das.affine_matrix is None
    assert das.to_dict_with_units() == TESTDATA["DAS"]
