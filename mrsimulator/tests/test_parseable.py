# coding: utf-8
"""
Tests for the base Parseable pattern
"""
import os
import unittest
from typing import ClassVar

from mrsimulator.parseable import Parseable, enforce_units


class ParseableTestClass(Parseable):
    """
	Dummy test class for Parseable pattern
	"""

    foo: float = 0
    bar: str = ""
    property_unit_types: ClassVar = {"foo": "angle"}
    property_default_units: ClassVar = {"foo": "rad"}


class TestEnforceUnits(unittest.TestCase):
    def test_good_units(self):
        enforce_units("300 Hz", "frequency", "Hz")

    def test_bad_units(self):
        with self.assertRaises(Exception):
            enforce_units("300 Hz", "angle", "rad")


class TestParseable(unittest.TestCase):
    def setUp(self):
        pass

    def test_blank_init(self):
        pr = ParseableTestClass()

    def test_parse_json(self):

        good_json = {"foo": "300 rad", "bar": "stuff"}
        pr = ParseableTestClass.parse_json_with_units(good_json)

        bad_json = {"foo": "300 Hz", "bar": "stuff"}

        with self.assertRaises(Exception) as err:
            pr = ParseableTestClass.parse_json_with_units(bad_json)
        self.assertEqual(
            str(err.exception),
            "A angle value is required but got a frequency instead",
        )


if __name__ == "__main__":
    unittest.main()
