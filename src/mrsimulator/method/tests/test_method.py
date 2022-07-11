"""Test for the base Dimension class."""
from copy import deepcopy

import csdmpy as cp
import numpy as np
import pytest
from mrsimulator.method import Method
from mrsimulator.method import SpectralDimension
from mrsimulator.method import SpectralEvent
from mrsimulator.method.frequency_contrib import freq_default
from mrsimulator.spin_system.isotope import Isotope
from mrsimulator.utils.error import MissingSpectralDimensionError
from pydantic import ValidationError

__author__ = "Deepansh J. Srivastava"
__email__ = "srivastava.89@osu.edu"

event_dictionary = {
    "fraction": 0.5,
    "freq_contrib": freq_default,
    "magnetic_flux_density": "9.6 T",
    "rotor_frequency": "1 kHz",
    "rotor_angle": "54.735 deg",
}
dimension_dictionary = {
    "count": 1024,
    "spectral_width": "100 Hz",
    "reference_offset": "0 GHz",
    "events": [event_dictionary, event_dictionary],
}


def basic_method_tests(the_method):
    assert the_method != "r"
    assert the_method.name == "test-1-d"
    the_method.name = "test worked"
    assert the_method.name == "test worked"

    assert the_method.description == "Test-1"
    the_method.description = "test worked again"
    assert the_method.description == "test worked again"

    # spectral width test
    assert the_method.channels == [Isotope(symbol="29Si")]
    the_method.channels = ["1H", "17O"]
    assert the_method.channels == [Isotope(symbol="1H"), Isotope(symbol="17O")]

    with pytest.raises(ValidationError, match=".*value is not a valid list.*"):
        the_method.channels = "6Li"

    dimension = SpectralDimension.parse_dict_with_units(dimension_dictionary)
    assert the_method.spectral_dimensions[0] == dimension

    the_method2 = deepcopy(the_method)
    the_method2.affine_matrix = [1]
    assert the_method2 != the_method

    # json()

    evt = [{"fraction": 0.5, "transition_query": [{"ch1": {"P": [0]}}]}] * 2
    serialize = {
        "name": "test worked",
        "description": "test worked again",
        "channels": ["1H", "17O"],
        "magnetic_flux_density": "9.6 T",
        "rotor_frequency": "1000.0 Hz",
        "rotor_angle": "0.9553059660790962 rad",
        "spectral_dimensions": [
            {
                "count": 1024,
                "spectral_width": "100.0 Hz",
                "events": evt,
            }
        ],
    }
    assert the_method.json() == serialize

    # json(units=False)
    assert the_method.json(units=False) == {
        "name": "test worked",
        "description": "test worked again",
        "channels": ["1H", "17O"],
        "magnetic_flux_density": 9.6,
        "rotor_frequency": 1000.0,
        "rotor_angle": 0.9553059660790962,
        "spectral_dimensions": [
            {
                "count": 1024,
                "spectral_width": 100.0,
                "events": evt,
            }
        ],
    }


def test_method():
    # test-1 single dimension method

    # parse dict with units
    method_dictionary = {
        "name": "test-1-d",
        "description": "Test-1",
        "channels": ["29Si"],
        "magnetic_flux_density": "9.6 T",
        "rotor_frequency": "1 kHz",
        "rotor_angle": "54.735 deg",
        "spectral_dimensions": [dimension_dictionary],
    }
    the_method = Method.parse_dict_with_units(method_dictionary)
    basic_method_tests(the_method)

    # test-2 two dimensional two events method

    # parse dict with units
    method_dictionary = {
        "name": "test-1-d",
        "description": "Test-1",
        "channels": ["29Si"],
        "magnetic_flux_density": "9.6 T",
        "rotor_frequency": "1 kHz",
        "rotor_angle": "54.735 deg",
        "spectral_dimensions": [dimension_dictionary, dimension_dictionary],
    }
    the_method = Method.parse_dict_with_units(method_dictionary)

    # test experiment assignment
    assert the_method.experiment is None

    with pytest.raises(ValidationError, match="Unable to read the dataset."):
        the_method.experiment = "test"

    data = np.random.rand(100).reshape(10, 10)
    csdm_dataset = cp.as_csdm(data)

    csdm_dataset.x[0] *= cp.ScalarQuantity("Hz")
    csdm_dataset.x[1] *= cp.ScalarQuantity("Hz")

    the_method.experiment = csdm_dataset
    the_method.simulation = csdm_dataset
    assert isinstance(the_method.experiment, cp.CSDM)
    assert isinstance(the_method.simulation, cp.CSDM)

    csdm_dict = csdm_dataset.to_dict()
    the_method.experiment = csdm_dict
    assert isinstance(the_method.experiment, cp.CSDM)
    assert the_method.experiment == csdm_dataset

    the_method.simulation = csdm_dict
    assert isinstance(the_method.simulation, cp.CSDM)
    assert the_method.simulation == csdm_dataset

    # json()
    event_dictionary_ = {"fraction": 0.5, "transition_query": [{"ch1": {"P": [0]}}]}
    dimension_dictionary_ = {
        "count": 1024,
        "spectral_width": "100.0 Hz",
        "events": [event_dictionary_, event_dictionary_],
    }
    method_dictionary_ = {
        "name": "test-1-d",
        "description": "Test-1",
        "channels": ["29Si"],
        "magnetic_flux_density": "9.6 T",
        "rotor_frequency": "1000.0 Hz",
        "rotor_angle": "0.9553059660790962 rad",
        "spectral_dimensions": [dimension_dictionary_, dimension_dictionary_],
        "simulation": csdm_dataset.to_dict(),
        "experiment": csdm_dataset.to_dict(),
    }
    serialize = the_method.json()
    serialize["simulation"]["csdm"].pop("timestamp")
    assert serialize == method_dictionary_

    # json(units=False)
    event_dictionary_ = {"fraction": 0.5, "transition_query": [{"ch1": {"P": [0]}}]}
    dimension_dictionary_ = {
        "count": 1024,
        "spectral_width": 100.0,
        "events": [event_dictionary_, event_dictionary_],
    }
    method_dictionary_ = {
        "name": "test-1-d",
        "description": "Test-1",
        "channels": ["29Si"],
        "magnetic_flux_density": 9.6,
        "rotor_frequency": 1000.0,
        "rotor_angle": 0.9553059660790962,
        "spectral_dimensions": [dimension_dictionary_, dimension_dictionary_],
        "simulation": csdm_dataset.to_dict(),
        "experiment": csdm_dataset.to_dict(),
    }
    serialize = the_method.json(units=False)
    serialize["simulation"]["csdm"].pop("timestamp")
    assert serialize == method_dictionary_


def test_rotor_frequency():
    """Ensures only 1 non-zero finite spinning speed in method"""
    # Good method, should not throw error
    Method(
        channels=["1H"],
        spectral_dimensions=[
            SpectralDimension(
                events=[
                    SpectralEvent(fraction=0.5, rotor_frequency=123),
                    SpectralEvent(fraction=0.5, rotor_frequency=0),
                ]
            ),
            SpectralDimension(events=[SpectralEvent(rotor_frequency=np.inf)]),
        ],
    )

    # Bad method, should throw error for multiple finite speeds
    for cls in [Method]:
        with pytest.raises(NotImplementedError):
            cls(
                channels=["1H"],
                spectral_dimensions=[
                    SpectralDimension(
                        events=[
                            SpectralEvent(fraction=0.5, rotor_frequency=123),
                            SpectralEvent(fraction=0.5, rotor_frequency=456),
                        ]
                    )
                ],
            )

    with pytest.raises(NotImplementedError):
        Method(
            channels=["1H"],
            spectral_dimensions=[
                SpectralDimension(
                    events=[
                        SpectralEvent(fraction=0.5, rotor_frequency=123),
                        SpectralEvent(fraction=0.5, rotor_frequency=np.inf),
                    ]
                ),
                SpectralDimension(
                    events=[
                        SpectralEvent(fraction=0.5, rotor_frequency=0),
                        SpectralEvent(fraction=0.5, rotor_frequency=456),
                    ]
                ),
            ],
        )

    with pytest.raises(NotImplementedError):
        Method(
            channels=["27Al"],
            rotor_frequency=np.inf,
            spectral_dimensions=[
                SpectralDimension(
                    events=[
                        SpectralEvent(fraction=0.5, rotor_frequency=0.1),
                        SpectralEvent(fraction=0.5, rotor_frequency=456),
                    ]
                )
            ],
        )


def test_empty_spectral_dimensions():
    e = ".*Method requires at least one SpectralDimension, none found.*"
    with pytest.raises(MissingSpectralDimensionError, match=e):
        Method(channels=["1H"])
