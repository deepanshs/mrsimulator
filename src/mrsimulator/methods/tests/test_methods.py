# -*- coding: utf-8 -*-
from os import path

import numpy as np
import pytest
from monty.serialization import loadfn
from mrsimulator.method.frequency_contrib import freq_default
from mrsimulator.methods import BlochDecaySpectrum
from mrsimulator.methods import Method1D
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


def test_method_1D():
    error = "Expecting a 1x1 affine matrix."
    with pytest.raises(ValueError, match=f".*{error}.*"):
        Method1D(spectral_dimensions=[{}], affine_matrix=[1, 2, 3, 4])

    error = r"The method allows 1 spectral dimension\(s\), 2 given."
    with pytest.raises(ValueError, match=f".*{error}.*"):
        Method1D(spectral_dimensions=[{}, {}])

    # parse dict with units test
    dict_1d = {
        "channels": ["87Rb"],
        "magnetic_flux_density": "7 T",  # in T
        "rotor_angle": "54.735 deg",
        "rotor_frequency": "1e9 Hz",
        "spectral_dimensions": [
            {
                "count": 1024,
                "spectral_width": "10 kHz",  # in Hz
                "reference_offset": "-4 kHz",  # in Hz
                "events": [
                    {"fraction": 27 / 17, "freq_contrib": ["Quad2_0"]},
                    {"fraction": 1, "freq_contrib": ["Quad2_4"]},
                ],
            }
        ],
    }
    method1a = Method1D.parse_dict_with_units(dict_1d)

    method1b = Method1D(
        channels=["87Rb"],
        magnetic_flux_density=7,  # in T
        rotor_angle=54.735 * np.pi / 180,
        rotor_frequency=1e9,
        spectral_dimensions=[
            {
                "count": 1024,
                "spectral_width": 1e4,  # in Hz
                "reference_offset": -4e3,  # in Hz
                "events": [
                    {"fraction": 27 / 17, "freq_contrib": ["Quad2_0"]},
                    {"fraction": 1, "freq_contrib": ["Quad2_4"]},
                ],
            }
        ],
    )

    assert method1a == method1b


def test_method_2D():

    # parse dict with units test
    dict_1d = {
        "channels": ["87Rb"],
        "magnetic_flux_density": "7 T",  # in T
        "rotor_angle": "54.735 deg",
        "spectral_dimensions": [
            {
                "count": 1024,
                "spectral_width": "10 kHz",  # in Hz
                "events": [
                    {"fraction": 27 / 17, "freq_contrib": ["Quad2_0"]},
                    {"fraction": 1, "freq_contrib": ["Quad2_4"]},
                ],
            },
            {
                "count": 1024,
                "spectral_width": "10 kHz",  # in Hz
                "events": [
                    {"fraction": 27 / 17, "freq_contrib": ["Quad2_0"]},
                    {"fraction": 1, "freq_contrib": ["Quad2_4"]},
                ],
            },
        ],
    }
    method1a = Method2D.parse_dict_with_units(dict_1d)

    method1b = Method2D(
        channels=["87Rb"],
        magnetic_flux_density=7,  # in T
        rotor_angle=54.735 * np.pi / 180,
        spectral_dimensions=[
            {
                "count": 1024,
                "spectral_width": 1e4,  # in Hz
                "events": [
                    {"fraction": 27 / 17, "freq_contrib": ["Quad2_0"]},
                    {"fraction": 1, "freq_contrib": ["Quad2_4"]},
                ],
            },
            {
                "count": 1024,
                "spectral_width": 1e4,  # in Hz
                "events": [
                    {"fraction": 27 / 17, "freq_contrib": ["Quad2_0"]},
                    {"fraction": 1, "freq_contrib": ["Quad2_4"]},
                ],
            },
        ],
    )

    assert method1a == method1b


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


def test_BlochDecaySpectrum():
    # test-1
    m1 = BlochDecaySpectrum()

    event_dictionary_ = {
        "fraction": 1.0,
        "freq_contrib": freq_default,
        "magnetic_flux_density": "9.4 T",
        "rotor_frequency": "0.0 Hz",
        "rotor_angle": "0.9553166 rad",
        "transition_query": {"P": {"channel-1": [[-1.0]]}},
        "user_variables": ["magnetic_flux_density", "rotor_frequency", "rotor_angle"],
    }
    dimension_dictionary_ = {
        "count": 1024,
        "spectral_width": "25000.0 Hz",
        "reference_offset": "0.0 Hz",
        "events": [event_dictionary_],
    }

    should_be = {
        "name": "BlochDecaySpectrum",
        "channels": ["1H"],
        "spectral_dimensions": [dimension_dictionary_],
    }
    dict_ = m1.to_dict_with_units()
    dict_.pop("description")
    assert dict_ == should_be

    # test-2
    m2_dict = {
        "channels": ["29Si"],
        "magnetic_flux_density": "11.7 T",
        "rotor_angle": "90 deg",
        "spectral_dimensions": [{}],
    }
    m2 = BlochDecaySpectrum.parse_dict_with_units(m2_dict)

    angle = 90 * np.pi / 180
    event_dictionary_ = {
        "fraction": 1.0,
        "freq_contrib": freq_default,
        "magnetic_flux_density": "11.7 T",
        "rotor_frequency": "0.0 Hz",
        "rotor_angle": f"{angle} rad",
        "transition_query": {"P": {"channel-1": [[-1.0]]}},
        "user_variables": ["magnetic_flux_density", "rotor_frequency", "rotor_angle"],
    }
    dimension_dictionary_ = {
        "count": 1024,
        "spectral_width": "25000.0 Hz",
        "reference_offset": "0.0 Hz",
        "events": [event_dictionary_],
    }

    should_be = {
        "name": "BlochDecaySpectrum",
        "channels": ["29Si"],
        "spectral_dimensions": [dimension_dictionary_],
    }

    dict_ = m2.to_dict_with_units()
    dict_.pop("description")
    assert dict_ == should_be


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
