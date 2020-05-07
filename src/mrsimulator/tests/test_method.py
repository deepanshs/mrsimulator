# -*- coding: utf-8 -*-
"""Test for the base Dimension class."""
import numpy as np
import pytest
from mrsimulator.method import Event
from mrsimulator.method import SpectralDimension
from pydantic import ValidationError


def basic_spectral_dimension_tests(the_method):
    assert the_method.count == 1024
    error = "ensure this value is greater than 0"
    with pytest.raises(ValidationError, match=".*{0}.*".format(error)):
        the_method.count = -1024
    with pytest.raises(ValidationError, match="value is not a valid integer"):
        the_method.count = "test"

    # spectral width test
    assert the_method.spectral_width == 100
    with pytest.raises(ValidationError, match=".*{0}.*".format(error)):
        the_method.spectral_width = 0
    # ensure the default value is Hz
    assert the_method.property_units["spectral_width"] == "Hz"

    # reference offset test
    assert the_method.reference_offset == 0
    the_method.reference_offset = -120
    assert the_method.reference_offset == -120
    # ensure the default value is Hz
    assert the_method.property_units["spectral_width"] == "Hz"
    with pytest.raises(ValidationError, match=".*value is not a valid float.*"):
        the_method.reference_offset = "-120 Hz"

    # label
    assert the_method.label == ""
    the_method.label = "This is a spectral dimensions"
    assert the_method.label == "This is a spectral dimensions"
    the_method.label = 45.0
    assert the_method.label == "45.0"


def test_spectral_dimension():
    # parse dict with units
    event_dictionary = {
        "fraction": 0.5,
        "magnetic_flux_density": "9.6 T",
        "rotor_frequency": "1 kHz",
        "rotor_angle": "54.735 deg",
    }
    dictionary = {
        "count": 1024,
        "spectral_width": "100 Hz",
        "reference_offset": "0 GHz",
        "events": [event_dictionary],
    }
    the_method = SpectralDimension.parse_dict_with_units(dictionary)
    basic_spectral_dimension_tests(the_method)

    # direct initialization
    the_method = SpectralDimension(
        count=1024,
        spectral_width=100,
        events=[Event.parse_dict_with_units(event_dictionary)],
    )
    basic_spectral_dimension_tests(the_method)

    # coordinate in Hz
    the_method.count = 32
    the_method.reference_offset = 0
    the_method.spectral_width = 32

    coordinates = np.arange(32) - 16
    assert np.allclose(the_method.coordinates_Hz, coordinates)

    the_method.reference_offset = -4
    assert np.allclose(the_method.coordinates_Hz, coordinates - 4)

    the_method.count = 31
    the_method.reference_offset = 0
    the_method.spectral_width = 31
    coordinates = np.arange(31) - 15
    assert np.allclose(the_method.coordinates_Hz, coordinates)

    # CSDM dimension
    csdm_dimension = the_method.to_csdm_dimension()
    assert np.allclose(csdm_dimension.coordinates.to("Hz").value, coordinates)

    # to dict with units
    should_be = dict(
        count=31,
        spectral_width="31.0 Hz",
        reference_offset="0.0 Hz",
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
    assert should_be == the_method.to_dict_with_units()


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


def test_events():
    # parse dict with units
    dictionary = {
        "fraction": 0.5,
        "magnetic_flux_density": "9.6 T",
        "rotor_frequency": "1 kHz",
        "rotor_angle": "54.735 deg",
    }
    the_event = Event.parse_dict_with_units(dictionary)
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
