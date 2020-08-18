# -*- coding: utf-8 -*-
import numpy as np
import pytest
from mrsimulator.methods import BlochDecaySpectrum
from mrsimulator.methods import MQVAS


def test_01():
    error = "method requires exactly 1 channel"
    with pytest.raises(ValueError, match=f".*{error}.*"):
        BlochDecaySpectrum(channels=["1H", "29Si"])


def test_02():
    error = " {'rotor_frequency'} connot be modified for MQVAS method."
    with pytest.raises(AttributeError, match=f".*{error}.*"):
        MQVAS(channels=["87Rb"], rotor_frequency=10, spectral_dimensions=[{}, {}])


def test_03():
    """MQMAS method declaration"""
    mth = MQVAS(
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
        "name": "MQVAS",
        "description": "Multi-quantum variable-angle spinning",
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
                            "P": {"channel-1": [[-3]]},
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
    mth = MQVAS(
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
        "name": "MQVAS",
        "description": "Multi-quantum variable-angle spinning",
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
    mth = MQVAS(
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
        "name": "MQVAS",
        "description": "Multi-quantum variable-angle spinning",
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


def test_06():
    """Test with order"""
    test_03()
    test_05()
    test_04()
