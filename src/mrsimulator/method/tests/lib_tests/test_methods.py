# -*- coding: utf-8 -*-
from os import path

import numpy as np
import pytest
from monty.serialization import loadfn
from mrsimulator.method import Method
from mrsimulator.method.lib import BlochDecaySpectrum
from mrsimulator.method.lib import ThreeQ_VAS

__author__ = "Deepansh J. Srivastava"
__email__ = "srivastava.89@osu.edu"

MODULE_DIR = path.dirname(path.abspath(__file__))
TESTDATA = loadfn(path.join(MODULE_DIR, "test_data.json"))

MAS_DIM = {
    "count": 128,
    "spectral_width": 5e4,  # in Hz
    "reference_offset": 0,  # in Hz
    "events": [
        {
            "rotor_angle": 54.735 * np.pi / 180,
            "transition_query": [{"ch1": {"P": [-1], "D": [0]}}],
        },
    ],
}


sample_method_dict = dict(
    magnetic_flux_density=9.4,
    spectral_dimensions=[
        {"count": 1024, "spectral_width": 5e4},
        {"count": 1024, "spectral_width": 5e4},
    ],
)

sq_tq = {"transition_query": [{"ch1": {"P": [-1]}}]}
sample_test_output = {
    "count": 1024,
    "events": [
        {"fraction": 27 / 17, "freq_contrib": ["Quad2_0"], **sq_tq},
        {"fraction": 1, "freq_contrib": ["Quad2_4"], **sq_tq},
    ],
}


def test_more_spectral_dimensions():
    error = "Method requires exactly 1 spectral dimensions, given 2."
    with pytest.raises(ValueError, match=f".*{error}.*"):
        BlochDecaySpectrum(spectral_dimensions=[{}, {}])


def test_03():
    """generic method declaration"""
    mth = Method(
        name="generic",
        channels=["87Rb"],
        magnetic_flux_density=9.4,  # in T
        rotor_frequency=1000000000000.0,
        spectral_dimensions=[
            {"count": 1024, "spectral_width": 5e4, "events": [sq_tq]},
            {"count": 1024, "spectral_width": 5e4, "events": [sq_tq]},
        ],
    )

    assert mth.json() == TESTDATA["generic"]
    assert Method.parse_dict_with_units(mth.json()) == mth


def test_1D_method():
    # Method1D replaced with generic Method object
    error = "Expecting a 1x1 affine matrix."
    with pytest.raises(ValueError, match=f".*{error}.*"):
        Method(spectral_dimensions=[{}], affine_matrix=[1, 2, 3, 4])


def test_2D_method():
    # Method2D replaced with generic Method object
    error = "The first element of the affine matrix cannot be zero."
    with pytest.raises(ValueError, match=f".*{error}.*"):
        Method(spectral_dimensions=[{}, {}], affine_matrix=[0, 1, 2, 3])


def test_04():
    """SAS method declaration"""
    mth = Method(
        name="SAS",
        channels=["87Rb"],
        magnetic_flux_density=9.4,  # in T
        rotor_angle=70.12 * np.pi / 180,
        rotor_frequency=1000000000000,
        spectral_dimensions=[
            {
                "count": 512,
                "spectral_width": 5e4,  # in Hz
                "events": [
                    {"transition_query": [{"ch1": {"P": [-1], "D": [0]}}]},
                ],
            },
            MAS_DIM.copy(),
        ],
    )

    assert mth.json() == TESTDATA["SAS"]
    assert Method.parse_dict_with_units(mth.json()) == mth


def test_BlochDecaySpectrum():
    # test-1
    m1 = BlochDecaySpectrum(channels=["1H"])

    dimension_dictionary_ = {
        "count": 1024,
        "spectral_width": "25000.0 Hz",
        "events": [{"transition_query": [{"ch1": {"P": [-1]}}]}],
    }

    should_be = {
        "name": "BlochDecaySpectrum",
        "channels": ["1H"],
        "magnetic_flux_density": "9.4 T",
        "rotor_frequency": "0.0 Hz",
        "rotor_angle": "0.9553166181245 rad",
        "spectral_dimensions": [dimension_dictionary_],
    }
    dict_ = m1.json()
    assert Method.parse_dict_with_units(dict_) == m1

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
    dimension_dictionary_ = {
        "count": 1024,
        "spectral_width": "25000.0 Hz",
        "events": [{"transition_query": [{"ch1": {"P": [-1]}}]}],
    }

    should_be = {
        "name": "BlochDecaySpectrum",
        "channels": ["29Si"],
        "magnetic_flux_density": "11.7 T",
        "rotor_frequency": "0.0 Hz",
        "rotor_angle": f"{angle} rad",
        "spectral_dimensions": [dimension_dictionary_],
    }

    dict_ = m2.json()
    assert Method.parse_dict_with_units(dict_) == m2

    dict_.pop("description")
    assert dict_ == should_be


def test_05():
    """Satellite to central correlation method declaration"""
    sp0 = dict(
        count=512,
        spectral_width=5e6,
        events=[
            {
                "rotor_angle": 0 * np.pi / 180,
                "transition_query": [
                    {"ch1": {"P": [-1], "D": [2]}},
                    {"ch1": {"P": [-1], "D": [-2]}},
                ],
            }
        ],
    )
    mth = Method(
        name="STMAS",
        channels=["87Rb"],
        magnetic_flux_density=9.4,  # in T
        rotor_frequency=1000000000000,
        spectral_dimensions=[sp0, MAS_DIM.copy()],
    )

    assert mth.json() == TESTDATA["STMAS"]
    assert Method.parse_dict_with_units(TESTDATA["STMAS"]) == mth


def test_3QMAS():
    """3Q MAS correlation method declaration"""
    mth = ThreeQ_VAS(
        channels=["87Rb"],
        magnetic_flux_density=9.4,  # in T
        spectral_dimensions=[
            {"count": 512, "spectral_width": 5e6},
            {"count": 128, "spectral_width": 5e4},
        ],
    )

    assert np.allclose(mth.affine_matrix, [0.5625, 0.4375, 0, 1])
    assert ThreeQ_VAS.parse_dict_with_units(mth.json()) == mth


def test_06():
    """Test with order"""
    test_03()
    test_05()
    test_04()


def test_methods():
    das = Method(
        name="DAS",
        channels=["87Rb"],
        magnetic_flux_density=4.2,  # in T
        rotor_angle=54.735 * 3.14159 / 180,  # in rads
        rotor_frequency=1000000000000,
        spectral_dimensions=[
            {
                "count": 256,
                "spectral_width": 2e4,  # in Hz
                "reference_offset": -5e3,  # in Hz
                "label": "70.12 dimension",
                "events": [
                    {
                        "fraction": 0.5,
                        "rotor_angle": 70.12 * 3.14159 / 180,  # in rads
                        "transition_query": [{"ch1": {"P": [-1], "D": [0]}}],
                    },
                    {
                        "fraction": 0.5,
                        "rotor_angle": 30.12 * 3.14159 / 180,  # in rads
                        "transition_query": [{"ch1": {"P": [-1], "D": [0]}}],
                    },
                ],
            },
            # The last spectral dimension block is the direct-dimension
            {
                "count": 256,
                "spectral_width": 3e4,  # in Hz
                "reference_offset": -7e3,  # in Hz
                "label": "MAS dimension",
                "events": [{"transition_query": [{"ch1": {"P": [-1], "D": [0]}}]}],
            },
        ],
    )

    assert das.affine_matrix is None
    assert das.json() == TESTDATA["DAS"]
    assert Method.parse_dict_with_units(das.json()) == das
