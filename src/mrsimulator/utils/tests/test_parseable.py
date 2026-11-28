"""Tests for the base Parseable pattern"""
from typing import ClassVar

import pytest
from mrsimulator.utils.parseable import enforce_units
from mrsimulator.utils.parseable import Parseable

default = {"name": None, "description": None, "label": None}


class ParseableTestClass(Parseable):
    """Dummy test class for Parseable pattern"""

    foo: float = 0
    bar: float = 0
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
    assert pr.foo == 0
    assert pr.bar == 0


def test_parse_json():

    good_json = {"foo": "300 rad", "bar": "300 ppm"}
    pr = ParseableTestClass.parse_dict_with_units(good_json)
    assert pr.dict() == {
        "foo": 300.0,
        "bar": 0.03,
        **default,
        "property_units": {"foo": "rad", "bar": "pct"},
    }

    good_json2 = {"foo": "300 rad", "bar": "300 kHz"}
    pr = ParseableTestClass.parse_dict_with_units(good_json2)
    assert pr.dict() == {
        "foo": 300.0,
        "bar": 300000.0,
        **default,
        "property_units": {"foo": "rad", "bar": "Hz"},
    }

    bad_json = {"foo": "300 Hz", "bar": "300 ppm"}

    with pytest.raises(Exception) as err:
        ParseableTestClass.parse_dict_with_units(bad_json)
    assert (
        str(err.value) == "Error enforcing units for foo: 300 Hz\n"
        "A angle value is required but got a frequency instead"
    )

    bad_json = {"foo": "300 rad", "bar": "300 Pa"}

    with pytest.raises(Exception) as err:
        ParseableTestClass.parse_dict_with_units(bad_json)
    assert str(err.value) == "Could not enforce any units on bar"


def test_reduced_dict():
    obj = {"foo": "300 rad", "bar": "300 ppm"}
    pr = ParseableTestClass.parse_dict_with_units(obj)
    assert pr.reduced_dict() == {"foo": 300, "bar": 0.03}


def test_parse_dict_with_units_does_not_modify_the_input():
    """parse_dict_with_units leaves the given dictionary unchanged."""
    from copy import deepcopy

    from mrsimulator.method.query import Rotation
    from mrsimulator.spin_system.tensors import AntisymmetricTensor
    from mrsimulator.spin_system.tensors import SymmetricTensor

    for cls, py_dict in [
        (SymmetricTensor, {"zeta": "10 ppm", "eta": 0.1, "alpha": "30 deg"}),
        (AntisymmetricTensor, {"zeta": "10 ppm", "alpha": "30 deg"}),
        (Rotation, {"angle": "90 deg", "phase": "10 deg"}),
    ]:
        original = deepcopy(py_dict)
        cls.parse_dict_with_units(py_dict)
        assert py_dict == original, cls.__name__
