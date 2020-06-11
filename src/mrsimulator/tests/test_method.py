# -*- coding: utf-8 -*-
"""Test for the base Dimension class."""
import csdmpy as cp
import numpy as np
import pytest
from mrsimulator.isotope import Isotope
from mrsimulator.method import Method
from mrsimulator.method.event import Event
from mrsimulator.method.spectral_dimension import SpectralDimension
from pydantic import ValidationError


def basic_event_tests(the_event):
    # fraction
    assert the_event.fraction == 0.5
    the_event.fraction = 1.2
    assert the_event.fraction == 1.2
    with pytest.raises(ValidationError, match="value is not a valid float"):
        the_event.fraction = "test"

    # magnetic flux density
    assert the_event.magnetic_flux_density == 9.6
    the_event.magnetic_flux_density = 11.7
    assert the_event.magnetic_flux_density == 11.7
    with pytest.raises(ValidationError, match="value is not a valid float"):
        the_event.magnetic_flux_density = "test"
    # ensure the default value is T
    assert the_event.property_units["magnetic_flux_density"] == "T"

    # rotor frequency
    assert the_event.rotor_frequency == 1000
    the_event.rotor_frequency = 2.5e4
    assert the_event.rotor_frequency == 25000
    with pytest.raises(ValidationError, match="value is not a valid float"):
        the_event.rotor_frequency = "test"
    # ensure the default value is Hz
    assert the_event.property_units["rotor_frequency"] == "Hz"

    # rotor angle
    assert the_event.rotor_angle == 54.735 * np.pi / 180
    the_event.rotor_angle = 90 * np.pi / 180
    assert the_event.rotor_angle == 90 * np.pi / 180
    with pytest.raises(ValidationError, match="value is not a valid float"):
        the_event.rotor_angle = "test"
    # ensure the default value is radians
    assert the_event.property_units["rotor_angle"] == "rad"

    # to dict with units
    angle = 90 * np.pi / 180
    should_be = dict(
        fraction=1.2,
        magnetic_flux_density="11.7 T",
        rotor_frequency="25000.0 Hz",
        rotor_angle=f"{angle} rad",
        transition_query={"P": {"channel-1": [[-1.0]]}},
    )
    assert should_be == the_event.to_dict_with_units()

    # reduced_dict()
    assert the_event.reduced_dict() == {
        "fraction": 1.2,
        "magnetic_flux_density": 11.7,
        "rotor_frequency": 25000,
        "rotor_angle": angle,
        "transition_query": {"P": {"channel-1": [[-1.0]]}},
    }


def test_events():
    # parse dict with units
    event_dictionary = {
        "fraction": 0.5,
        "magnetic_flux_density": "9.6 T",
        "rotor_frequency": "1 kHz",
        "rotor_angle": "54.735 deg",
    }
    the_event = Event.parse_dict_with_units(event_dictionary)
    basic_event_tests(the_event)

    # direct initialization
    magic_angle_in_rad = 54.735 * np.pi / 180
    the_event = Event(
        fraction=0.5,
        magnetic_flux_density=9.6,
        rotor_frequency=1000,
        rotor_angle=magic_angle_in_rad,
    )
    basic_event_tests(the_event)


def basic_spectral_dimension_tests(the_dimension):
    assert the_dimension.count == 1024
    error = "ensure this value is greater than 0"
    with pytest.raises(ValidationError, match=".*{0}.*".format(error)):
        the_dimension.count = -1024
    with pytest.raises(ValidationError, match="value is not a valid integer"):
        the_dimension.count = "test"

    # spectral width test
    assert the_dimension.spectral_width == 100
    with pytest.raises(ValidationError, match=".*{0}.*".format(error)):
        the_dimension.spectral_width = 0
    # ensure the default value is Hz
    assert the_dimension.property_units["spectral_width"] == "Hz"

    # reference offset test
    assert the_dimension.reference_offset == 0
    the_dimension.reference_offset = -120
    assert the_dimension.reference_offset == -120
    # ensure the default value is Hz
    assert the_dimension.property_units["spectral_width"] == "Hz"
    with pytest.raises(ValidationError, match=".*value is not a valid float.*"):
        the_dimension.reference_offset = "-120 Hz"

    # label
    assert the_dimension.label == ""
    the_dimension.label = "This is a spectral dimensions"
    assert the_dimension.label == "This is a spectral dimensions"
    the_dimension.label = 45.0
    assert the_dimension.label == "45.0"

    # coordinate in Hz and ppm
    the_dimension.count = 32
    the_dimension.reference_offset = 0
    the_dimension.spectral_width = 32

    coordinates = np.arange(32) - 16
    assert np.allclose(the_dimension.coordinates_Hz, coordinates)

    with pytest.warns(
        UserWarning, match=".*The coordinates along the dimension without an origin.*"
    ):
        the_dimension.coordinates_ppm

    the_dimension.reference_offset = -4
    the_dimension.origin_offset = 5e6
    assert np.allclose(the_dimension.coordinates_Hz, coordinates - 4)
    assert np.allclose(the_dimension.coordinates_ppm, (coordinates - 4) / (5 - 4e-6))

    the_dimension.count = 31
    the_dimension.reference_offset = 0
    the_dimension.spectral_width = 31
    coordinates = np.arange(31) - 15
    assert np.allclose(the_dimension.coordinates_Hz, coordinates)
    assert np.allclose(the_dimension.coordinates_ppm, coordinates / 5)

    # CSDM dimension
    csdm_dimension = the_dimension.to_csdm_dimension()
    assert np.allclose(csdm_dimension.coordinates.to("Hz").value, coordinates)
    csdm_dimension.to("ppm", "nmr_frequency_ratio")
    assert np.allclose(csdm_dimension.coordinates.value, the_dimension.coordinates_ppm)

    # to dict with units
    should_be = dict(
        count=31,
        spectral_width="31.0 Hz",
        reference_offset="0.0 Hz",
        origin_offset="5000000.0 Hz",
        label="45.0",
        events=[
            {
                "fraction": 0.5,
                "magnetic_flux_density": "9.6 T",
                "rotor_angle": "0.9553059660790962 rad",
                "rotor_frequency": "1000.0 Hz",
                "transition_query": {"P": {"channel-1": [[-1.0]]}},
            }
        ],
    )
    assert should_be == the_dimension.to_dict_with_units()

    # reduced_dict()
    assert the_dimension.reduced_dict() == dict(
        count=31,
        spectral_width=31,
        reference_offset=0.0,
        origin_offset=5000000.0,
        label="45.0",
        events=[
            {
                "fraction": 0.5,
                "magnetic_flux_density": 9.6,
                "rotor_angle": 0.9553059660790962,
                "rotor_frequency": 1000,
                "transition_query": {"P": {"channel-1": [[-1.0]]}},
            }
        ],
    )


def test_spectral_dimension():
    # test-1 single event spectral dimensions

    # parse dict with units
    event_dictionary = {
        "fraction": 0.5,
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
    the_dimension = SpectralDimension.parse_dict_with_units(dimension_dictionary)
    basic_spectral_dimension_tests(the_dimension)

    # direct initialization
    the_dimension = SpectralDimension(
        count=1024,
        spectral_width=100,
        events=[Event.parse_dict_with_units(event_dictionary)],
    )
    basic_spectral_dimension_tests(the_dimension)

    # test-2: two event spectral dimension
    # parse dict with units
    event_dictionary = {
        "fraction": 0.5,
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
    the_dimension = SpectralDimension.parse_dict_with_units(dimension_dictionary)

    # direct initialization
    the_dimension2 = SpectralDimension(
        count=1024,
        spectral_width=100,
        events=[Event.parse_dict_with_units(event_dictionary) for _ in range(2)],
    )

    # to dict with units
    should_be = dict(
        count=1024,
        spectral_width="100.0 Hz",
        reference_offset="0.0 Hz",
        events=[
            {
                "fraction": 0.5,
                "magnetic_flux_density": "9.6 T",
                "rotor_angle": "0.9553059660790962 rad",
                "rotor_frequency": "1000.0 Hz",
                "transition_query": {"P": {"channel-1": [[-1.0]]}},
            },
            {
                "fraction": 0.5,
                "magnetic_flux_density": "9.6 T",
                "rotor_angle": "0.9553059660790962 rad",
                "rotor_frequency": "1000.0 Hz",
                "transition_query": {"P": {"channel-1": [[-1.0]]}},
            },
        ],
    )
    assert should_be == the_dimension.to_dict_with_units()
    assert should_be == the_dimension2.to_dict_with_units()

    # reduced_dict()
    reduced_dict = dict(
        count=1024,
        spectral_width=100,
        reference_offset=0.0,
        label="",
        events=[
            {
                "fraction": 0.5,
                "magnetic_flux_density": 9.6,
                "rotor_angle": 0.9553059660790962,
                "rotor_frequency": 1000,
                "transition_query": {"P": {"channel-1": [[-1.0]]}},
            },
            {
                "fraction": 0.5,
                "magnetic_flux_density": 9.6,
                "rotor_angle": 0.9553059660790962,
                "rotor_frequency": 1000,
                "transition_query": {"P": {"channel-1": [[-1.0]]}},
            },
        ],
    )
    assert the_dimension.reduced_dict() == reduced_dict
    assert the_dimension2.reduced_dict() == reduced_dict


def basic_method_tests(the_method):
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

    event_dictionary = {
        "fraction": 0.5,
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
    dimension = SpectralDimension.parse_dict_with_units(dimension_dictionary)
    assert the_method.spectral_dimensions[0] == dimension

    # to_dict_with_unit()
    serialize = {
        "name": "test worked",
        "description": "test worked again",
        "channels": ["1H", "17O"],
        "spectral_dimensions": [dimension.to_dict_with_units()],
    }
    assert the_method.to_dict_with_units() == serialize

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
    event_dictionary = {
        "fraction": 0.5,
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
    event_dictionary = {
        "fraction": 0.5,
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

    csdm_data.dimensions[0] *= cp.ScalarQuantity("Hz")
    csdm_data.dimensions[1] *= cp.ScalarQuantity("Hz")
    the_method.experiment = csdm_data

    assert isinstance(the_method.experiment, cp.CSDM)

    csdm_dict = csdm_data.to_dict()
    the_method.experiment = csdm_dict
    assert isinstance(the_method.experiment, cp.CSDM)
    assert the_method.experiment == csdm_data

    # to_dict_with_unit()
    event_dictionary_ = {
        "fraction": 0.5,
        "magnetic_flux_density": "9.6 T",
        "rotor_frequency": "1000.0 Hz",
        "rotor_angle": "0.9553059660790962 rad",
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
        "spectral_dimensions": [dimension_dictionary_, dimension_dictionary_],
        "experiment": csdm_data.to_dict(),
    }
    assert the_method.to_dict_with_units() == method_dictionary_

    # reduced_dict()
    event_dictionary_ = {
        "fraction": 0.5,
        "magnetic_flux_density": 9.6,
        "rotor_frequency": 1000.0,
        "rotor_angle": 0.9553059660790962,
        "transition_query": {"P": {"channel-1": [[-1.0]]}},
    }
    dimension_dictionary_ = {
        "count": 1024,
        "spectral_width": 100.0,
        "reference_offset": 0.0,
        "label": "",
        "events": [event_dictionary_, event_dictionary_],
    }
    method_dictionary_ = {
        "name": "test-1-d",
        "description": "Test-1",
        "channels": ["29Si"],
        "spectral_dimensions": [dimension_dictionary_, dimension_dictionary_],
        "experiment": csdm_data.to_dict(),
    }
    assert the_method.reduced_dict() == method_dictionary_

    # update_spectral_dimension_attributes_from_experiment
    the_method.update_spectral_dimension_attributes_from_experiment()
    for i in range(2):
        assert the_method.spectral_dimensions[i].count == 10
        assert the_method.spectral_dimensions[i].spectral_width == 10
        assert the_method.spectral_dimensions[i].reference_offset == 0
        assert the_method.spectral_dimensions[i].origin_offset == 0
