# -*- coding: utf-8 -*-
import numpy as np
import pytest
from mrsimulator.method.event import Event
from mrsimulator.method.frequency_contrib import freq_default
from pydantic import ValidationError

__author__ = "Deepansh J. Srivastava"
__email__ = "srivastava.89@osu.edu"


def test_freq_contrib():
    event = Event(freq_contrib=["Quad2_4", "Quad2_0"])
    assert event.json()["freq_contrib"] == ["Quad2_4", "Quad2_0"]
    assert event.dict()["freq_contrib"] == ["Quad2_4", "Quad2_0"]
    assert event.reduced_dict()["freq_contrib"] == ["Quad2_4", "Quad2_0"]
    assert event.get_value_int() == [0, 0, 0, 1, 0, 1]

    event = Event()
    assert event.dict()["freq_contrib"] == freq_default
    assert event.reduced_dict()["freq_contrib"] == freq_default
    assert event.get_value_int() == [1, 1, 1, 1, 1, 1]


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
    assert should_be == the_event.json()

    # reduced_dict()
    assert the_event.reduced_dict() == {
        "fraction": 1.2,
        "freq_contrib": freq_default,
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
