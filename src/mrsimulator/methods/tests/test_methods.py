# -*- coding: utf-8 -*-
import numpy as np
import pytest
from mrsimulator.methods import BlochDecaySpectrum
from mrsimulator.methods import Method2D
from mrsimulator.methods import ThreeQ_VAS


def test_more_spectral_dimensions():
    error = "The method allows 1 spectral dimension"
    with pytest.raises(ValueError, match=f".*{error}.*"):
        BlochDecaySpectrum(spectral_dimensions=[{}, {}])


def test_01():
    error = "method requires exactly 1 channel"
    with pytest.raises(ValueError, match=f".*{error}.*"):
        BlochDecaySpectrum(channels=["1H", "29Si"])


def test_02():
    error = " `rotor_frequency` cannot be modified for Method2D class."
    with pytest.raises(AttributeError, match=f".*{error}.*"):
        Method2D(channels=["87Rb"], rotor_frequency=10, spectral_dimensions=[{}, {}])


def test_03():
    """MQMAS method declaration"""
    mth = Method2D(
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

    serialize = {
        "name": "Method2D",
        "description": "Simulate a generic two-dimensional correlation spectrum",
        "spectral_dimensions": [
            {
                "count": 1024,
                "spectral_width": "50000.0 Hz",
                "reference_offset": "0.0 Hz",
                "events": [
                    {
                        "fraction": 1.0,
                        "magnetic_flux_density": "9.4 T",
                        "rotor_frequency": "1000000000000.0 Hz",
                        "rotor_angle": "0.9553166 rad",
                        "transition_query": {
                            "P": {"channel-1": [[-1]]},
                            "D": {"channel-1": [[0]]},
                        },
                        "user_variables": ["magnetic_flux_density", "rotor_angle"],
                    }
                ],
            },
            {
                "count": 1024,
                "spectral_width": "50000.0 Hz",
                "reference_offset": "0.0 Hz",
                "events": [
                    {
                        "fraction": 1.0,
                        "magnetic_flux_density": "9.4 T",
                        "rotor_frequency": "1000000000000.0 Hz",
                        "rotor_angle": "0.9553166 rad",
                        "transition_query": {
                            "P": {"channel-1": [[-1]]},
                            "D": {"channel-1": [[0]]},
                        },
                        "user_variables": ["magnetic_flux_density", "rotor_angle"],
                    }
                ],
            },
        ],
        "channels": ["87Rb"],
    }

    assert serialize == mth.to_dict_with_units()


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

    serialize = {
        "name": "Method2D",
        "description": "Simulate a generic two-dimensional correlation spectrum",
        "spectral_dimensions": [
            {
                "count": 512,
                "spectral_width": "50000.0 Hz",
                "reference_offset": "0.0 Hz",
                "events": [
                    {
                        "fraction": 1.0,
                        "magnetic_flux_density": "9.4 T",
                        "rotor_frequency": "1000000000000.0 Hz",
                        "rotor_angle": "1.223824871498424 rad",
                        "transition_query": {
                            "P": {"channel-1": [[-1]]},
                            "D": {"channel-1": [[0]]},
                        },
                        "user_variables": ["magnetic_flux_density", "rotor_angle"],
                    }
                ],
            },
            {
                "count": 128,
                "spectral_width": "50000.0 Hz",
                "reference_offset": "0.0 Hz",
                "events": [
                    {
                        "fraction": 1.0,
                        "magnetic_flux_density": "9.4 T",
                        "rotor_frequency": "1000000000000.0 Hz",
                        "rotor_angle": "0.9553059660790962 rad",
                        "transition_query": {
                            "P": {"channel-1": [[-1]]},
                            "D": {"channel-1": [[0]]},
                        },
                        "user_variables": ["magnetic_flux_density", "rotor_angle"],
                    }
                ],
            },
        ],
        "channels": ["87Rb"],
    }

    assert serialize == mth.to_dict_with_units()


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

    serialize = {
        "name": "Method2D",
        "description": "Simulate a generic two-dimensional correlation spectrum",
        "spectral_dimensions": [
            {
                "count": 512,
                "spectral_width": "5000000.0 Hz",
                "reference_offset": "0.0 Hz",
                "events": [
                    {
                        "fraction": 1.0,
                        "magnetic_flux_density": "9.4 T",
                        "rotor_frequency": "1000000000000.0 Hz",
                        "rotor_angle": "0.0 rad",
                        "transition_query": {
                            "P": {"channel-1": [[-1]]},
                            "D": {"channel-1": [[2], [-2]]},
                        },
                        "user_variables": ["magnetic_flux_density", "rotor_angle"],
                    }
                ],
            },
            {
                "count": 128,
                "spectral_width": "50000.0 Hz",
                "reference_offset": "0.0 Hz",
                "events": [
                    {
                        "fraction": 1.0,
                        "magnetic_flux_density": "9.4 T",
                        "rotor_frequency": "1000000000000.0 Hz",
                        "rotor_angle": "0.9553059660790962 rad",
                        "transition_query": {
                            "P": {"channel-1": [[-1]]},
                            "D": {"channel-1": [[0]]},
                        },
                        "user_variables": ["magnetic_flux_density", "rotor_angle"],
                    }
                ],
            },
        ],
        "channels": ["87Rb"],
    }

    assert serialize == mth.to_dict_with_units()


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
    res = {
        "name": "Method2D",
        "description": "Simulate a generic two-dimensional correlation spectrum",
        "spectral_dimensions": [
            {
                "count": 256,
                "spectral_width": "20000.0 Hz",
                "reference_offset": "-5000.0 Hz",
                "label": "70.12 dimension",
                "events": [
                    {
                        "fraction": 1.0,
                        "magnetic_flux_density": "4.2 T",
                        "rotor_frequency": "1000000000000.0 Hz",
                        "rotor_angle": "1.223788777777778 rad",
                        "transition_query": {
                            "P": {"channel-1": [[-1]]},
                            "D": {"channel-1": [[0]]},
                        },
                        "user_variables": ["magnetic_flux_density", "rotor_angle"],
                    },
                    {
                        "fraction": 1.0,
                        "magnetic_flux_density": "4.2 T",
                        "rotor_frequency": "1000000000000.0 Hz",
                        "rotor_angle": "0.5256776666666667 rad",
                        "transition_query": {
                            "P": {"channel-1": [[-1]]},
                            "D": {"channel-1": [[0]]},
                        },
                        "user_variables": ["magnetic_flux_density", "rotor_angle"],
                    },
                ],
            },
            {
                "count": 256,
                "spectral_width": "30000.0 Hz",
                "reference_offset": "-7000.0 Hz",
                "label": "MAS dimension",
                "events": [
                    {
                        "fraction": 1.0,
                        "magnetic_flux_density": "4.2 T",
                        "rotor_frequency": "1000000000000.0 Hz",
                        "rotor_angle": "0.9552777916666667 rad",
                        "transition_query": {
                            "P": {"channel-1": [[-1]]},
                            "D": {"channel-1": [[0]]},
                        },
                        "user_variables": ["magnetic_flux_density", "rotor_angle"],
                    }
                ],
            },
        ],
        "channels": ["87Rb"],
    }

    method = Method2D(
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
                        "rotor_angle": 70.12 * 3.1415 / 180,  # in rads
                        "transition_query": {"P": [-1], "D": [0]},
                    },
                    {
                        "rotor_angle": 30.12 * 3.1415 / 180,  # in rads
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
                        "rotor_angle": 54.735 * 3.1415 / 180,  # in rads
                        "transition_query": {"P": [-1], "D": [0]},
                    }
                ],
            },
        ],
    )

    assert method.to_dict_with_units() == res
