# -*- coding: utf-8 -*-
"""Tests for the base Parseable pattern"""
from typing import ClassVar

import pytest
from mrsimulator.parseable import enforce_units
from mrsimulator.parseable import Parseable


class ParseableTestClass(Parseable):
    """
    Dummy test class for Parseable pattern
    """

    foo: float = 0
    bar: int = 0
    property_unit_types: ClassVar = {
        "foo": "angle",
        "bar": ["dimensionless", "frequency"],
    }
    property_default_units: ClassVar = {"foo": "rad", "bar": ["pct", "Hz"]}


# Test Enforce Units
def test_good_units():
    enforce_units("300 Hz", "frequency", "Hz")


def test_bad_units():
    with pytest.raises(Exception):
        enforce_units("300 Hz", "angle", "rad")


# Test Parseable pattern


def test_blank_init():
    pr = ParseableTestClass()
    pr


def test_parse_json():

    good_json = {"foo": "300 rad", "bar": "300 ppm"}
    pr = ParseableTestClass.parse_dict_with_units(good_json)
    pr

    good_json2 = {"foo": "300 rad", "bar": "300 Hz"}
    pr = ParseableTestClass.parse_dict_with_units(good_json2)
    pr

    bad_json = {"foo": "300 Hz", "bar": "300 ppm"}

    with pytest.raises(Exception) as err:
        ParseableTestClass.parse_dict_with_units(bad_json)
    assert (
        str(err.value) == "Error enforcing units for foo: 300 Hz\n"
        "A angle value is required but got a frequency instead"
    )
