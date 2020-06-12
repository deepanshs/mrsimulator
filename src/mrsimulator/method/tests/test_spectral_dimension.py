# -*- coding: utf-8 -*-
import numpy as np
import pytest
from mrsimulator.method.event import Event
from mrsimulator.method.spectral_dimension import SpectralDimension
from pydantic import ValidationError


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
