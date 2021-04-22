# -*- coding: utf-8 -*-
"""Test for the base Dimension class."""
from copy import deepcopy

import csdmpy as cp
import numpy as np
import pytest
from mrsimulator.method import Method
from mrsimulator.method.frequency_contrib import freq_default
from mrsimulator.method.spectral_dimension import SpectralDimension
from mrsimulator.spin_system.isotope import Isotope
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
    "events": [event_dictionary],
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

    # to_dict_with_unit()
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
                "reference_offset": "0.0 Hz",
                "events": [
                    {
                        "fraction": 0.5,
                        "transition_query": {"P": {"channel-1": [[-1.0]]}},
                    }
                ],
            }
        ],
    }
    assert the_method.json() == serialize

    # reduced_dict()
    assert the_method.reduced_dict() == {
        "name": "test worked",
        "description": "test worked again",
        "channels": ["1H", "17O"],
        "spectral_dimensions": [dimension.reduced_dict()],
    }


def test_method():
    # test-1 single dimension method

    # parse dict with units
    method_dictionary = {
        "name": "test-1-d",
        "description": "Test-1",
        "channels": ["29Si"],
        "spectral_dimensions": [dimension_dictionary],
    }
    the_method = Method.parse_dict_with_units(method_dictionary)
    basic_method_tests(the_method)

    # test-2 two dimensional two events method

    # parse dict with units
    dimension_dictionary["events"].append(event_dictionary)

    method_dictionary = {
        "name": "test-1-d",
        "description": "Test-1",
        "channels": ["29Si"],
        "spectral_dimensions": [dimension_dictionary, dimension_dictionary],
    }
    the_method = Method.parse_dict_with_units(method_dictionary)

    # test experiment assignment
    assert the_method.experiment is None

    with pytest.raises(ValidationError, match="Unable to read the data."):
        the_method.experiment = "test"

    data = np.random.rand(100).reshape(10, 10)
    csdm_data = cp.as_csdm(data)

    csdm_data.x[0] *= cp.ScalarQuantity("Hz")
    csdm_data.x[1] *= cp.ScalarQuantity("Hz")

    the_method.experiment = csdm_data
    the_method.simulation = csdm_data
    assert isinstance(the_method.experiment, cp.CSDM)
    assert isinstance(the_method.simulation, cp.CSDM)

    csdm_dict = csdm_data.to_dict()
    the_method.experiment = csdm_dict
    assert isinstance(the_method.experiment, cp.CSDM)
    assert the_method.experiment == csdm_data

    the_method.simulation = csdm_dict
    assert isinstance(the_method.simulation, cp.CSDM)
    assert the_method.simulation == csdm_data

    # to_dict_with_unit()
    event_dictionary_ = {
        "fraction": 0.5,
        "transition_query": {"P": {"channel-1": [[-1.0]]}},
    }
    dimension_dictionary_ = {
        "count": 1024,
        "spectral_width": "100.0 Hz",
        "reference_offset": "0.0 Hz",
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
        "simulation": csdm_data.to_dict(),
        "experiment": csdm_data.to_dict(),
    }
    serialize = the_method.json()
    serialize["simulation"]["csdm"].pop("timestamp")
    assert serialize == method_dictionary_

    # reduced_dict()
    event_dictionary_ = {
        "fraction": 0.5,
        "freq_contrib": freq_default,
        "magnetic_flux_density": 9.6,
        "rotor_frequency": 1000.0,
        "rotor_angle": 0.9553059660790962,
        "transition_query": {"P": {"channel-1": [[-1.0]]}},
    }
    dimension_dictionary_ = {
        "count": 1024,
        "spectral_width": 100.0,
        "reference_offset": 0.0,
        "events": [event_dictionary_, event_dictionary_],
    }
    method_dictionary_ = {
        "name": "test-1-d",
        "description": "Test-1",
        "channels": ["29Si"],
        "spectral_dimensions": [dimension_dictionary_, dimension_dictionary_],
        "simulation": csdm_data.to_dict(),
        "experiment": csdm_data.to_dict(),
    }
    serialize = the_method.reduced_dict()
    serialize["simulation"]["csdm"].pop("timestamp")
    assert serialize == method_dictionary_

    # update_spectral_dimension_attributes_from_experiment
    the_method.update_spectral_dimension_attributes_from_experiment()
    for i in range(2):
        assert the_method.spectral_dimensions[i].count == 10
        assert the_method.spectral_dimensions[i].spectral_width == 10
        assert the_method.spectral_dimensions[i].reference_offset == 0
        assert the_method.spectral_dimensions[i].origin_offset == 0
