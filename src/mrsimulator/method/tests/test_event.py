# -*- coding: utf-8 -*-
import numpy as np
import pytest
from mrsimulator.method.event import BaseEvent
from mrsimulator.method.event import ConstantDurationEvent
from mrsimulator.method.event import MixingEvent
from mrsimulator.method.event import SpectralEvent
from mrsimulator.method.frequency_contrib import freq_default
from mrsimulator.method.query import MixingEnum
from pydantic import ValidationError

__author__ = "Deepansh J. Srivastava"
__email__ = "srivastava.89@osu.edu"


def test_freq_contrib():
    event = BaseEvent(freq_contrib=["Quad2_4", "Quad2_0"])
    assert event.json()["freq_contrib"] == ["Quad2_4", "Quad2_0"]
    assert event.dict()["freq_contrib"] == ["Quad2_4", "Quad2_0"]
    assert event.json(units=False)["freq_contrib"] == ["Quad2_4", "Quad2_0"]
    assert np.all(event._freq_contrib_flags() == [0, 0, 0, 1, 0, 1])

    event = BaseEvent(freq_contrib=["Shielding1_2"])
    assert np.all(event._freq_contrib_flags() == [0, 1, 0, 0, 0, 0])

    event = BaseEvent()
    assert event.dict()["freq_contrib"] == freq_default
    assert np.all(event._freq_contrib_flags() == [1, 1, 1, 1, 1, 1])


def basic_spectral_and_constant_time_event_tests(the_event, type_="spectral"):
    # fraction
    if type_ == "spectral":
        assert the_event.fraction == 0.5
        the_event.fraction = 1.2
        assert the_event.fraction == 1.2
        with pytest.raises(ValidationError, match="value is not a valid float"):
            the_event.fraction = "test"

    # duration
    if type_ == "constant_duration":
        assert the_event.duration == 0.5
        the_event.duration = 1.2
        assert the_event.duration == 1.2
        with pytest.raises(ValidationError, match="value is not a valid float"):
            the_event.duration = "test"

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

    should_be_units = dict(
        magnetic_flux_density="11.7 T",
        rotor_frequency="25000.0 Hz",
        rotor_angle=f"{angle} rad",
        transition_query=[{"ch1": {"P": [-1]}}],
    )
    should_be = {
        "magnetic_flux_density": 11.7,
        "rotor_frequency": 25000.0,
        "rotor_angle": angle,
    }

    if type_ == "spectral":
        should_be_units = dict(fraction=1.2, **should_be_units)
        assert the_event.json() == should_be_units
        assert the_event.json(units=False) == {
            "fraction": 1.2,
            **should_be,
            "transition_query": [{"ch1": {"P": [-1]}}],
        }

    if type_ == "constant_duration":
        should_be_units = dict(duration="1.2 µs", **should_be_units)
        assert the_event.json() == should_be_units
        assert the_event.json(units=False) == {
            "duration": 1.2,
            **should_be,
            "transition_query": [{"ch1": {"P": [-1.0]}}],
        }


def test_spectral_and_constant_time_events():
    # parse dict with units
    base_event_dictionary = {
        "magnetic_flux_density": "9.6 T",
        "rotor_frequency": "1 kHz",
        "rotor_angle": "54.735 deg",
    }
    evt_dict = {"fraction": 0.5, **base_event_dictionary}
    the_event = SpectralEvent.parse_dict_with_units(evt_dict)
    basic_spectral_and_constant_time_event_tests(the_event, type_="spectral")

    evt_dict = {"duration": "0.5 µs", **base_event_dictionary}
    the_event = ConstantDurationEvent.parse_dict_with_units(evt_dict)
    basic_spectral_and_constant_time_event_tests(the_event, type_="constant_duration")

    # direct initialization
    magic_angle_in_rad = 54.735 * np.pi / 180
    base_event_dict = dict(
        magnetic_flux_density=9.6, rotor_frequency=1000, rotor_angle=magic_angle_in_rad
    )
    the_event = SpectralEvent(fraction=0.5, **base_event_dict)
    basic_spectral_and_constant_time_event_tests(the_event, type_="spectral")

    the_event = ConstantDurationEvent(duration=0.5, **base_event_dict)
    basic_spectral_and_constant_time_event_tests(the_event, type_="constant_time")


def basic_mixing_event_tests(the_event):
    mix = the_event.query.ch1

    # tip angle
    assert mix.angle == np.pi / 2
    mix.angle = 3.2123
    assert mix.angle == 3.2123
    with pytest.raises(ValidationError, match="value is not a valid float"):
        mix.angle = "test"
    # ensure the default value is rad
    assert mix.property_units["angle"] == "rad"

    # phase
    assert mix.phase == 0.0
    mix.phase = 1.745
    assert mix.phase == 1.745
    with pytest.raises(ValidationError, match="value is not a valid float"):
        mix.phase = "test"
    # ensure the default value is rad
    assert mix.property_units["phase"] == "rad"

    # json()
    should_be_units = dict(ch1=dict(angle="3.2123 rad", phase="1.745 rad"))
    should_be = dict(ch1=dict(angle=3.2123, phase=1.745))

    should_be_units = dict(query=should_be_units)
    assert the_event.json() == should_be_units
    assert the_event.json(units=False) == {"query": should_be}


def test_Mixing_event():
    mix_event_dict = {"query": {"ch1": {"angle": "90 degree", "phase": "0 rad"}}}
    the_event = MixingEvent.parse_dict_with_units(mix_event_dict)
    basic_mixing_event_tests(the_event)


def test_total_and_no_mixing():
    no_mix = MixingEvent(query=MixingEnum.NoMixing)
    assert no_mix.query.value.ch1.angle == 0
    assert no_mix.query.value.ch1.phase == 0
    assert no_mix.query.value.ch2.angle == 0
    assert no_mix.query.value.ch2.phase == 0
    assert no_mix.query.value.ch3.angle == 0
    assert no_mix.query.value.ch3.phase == 0

    total_mix = MixingEvent(query=MixingEnum.TotalMixing)
    assert total_mix.query.value == "TotalMixing"


def check_equal(query, isotopes, channels, res):
    test = SpectralEvent(transition_query=query).permutation(isotopes, channels)
    for i, item in enumerate(res):
        for item2 in item[0]:
            assert item2 in test[i]["P"]

        for item2 in item[1]:
            assert item2 in test[i]["D"]


def test_BaseEvent_permutation():
    # P = -1 D = -1 on A B B A system, channel A, B
    # P = +1 D = -1 on A B B A system, channel A, B
    query = [
        {"ch1": {"P": [-1], "D": [1]}},
        {"ch1": {"P": [1], "D": [-1]}},
    ]
    res = [
        [
            [-1, 0, 0, 0],
            [0, 0, -1, 0],
        ],
        [
            [1, 0, 0, 0],
            [0, 0, 1, 0],
        ],
    ], [[[1, 0, 0, 0], [0, 0, 1, 0]], [[-1, 0, 0, 0], [0, 0, -1, 0]]]
    check_equal(query, ["A", "B", "A", "B"], ["A", "B"], res)
