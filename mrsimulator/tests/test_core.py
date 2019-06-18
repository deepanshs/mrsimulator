# coding: utf-8
"""
Tests for the base Parseable pattern
"""
import os
import unittest
from typing import ClassVar

import mrsimulator.unit
from mrsimulator import Site, Isotopomer, Spectrum


class TestSite(unittest.TestCase):
    def test_direct_init(self):
        Site(
            nucleus="29Si",
            isotropic_chemical_shift=10,
            shift_anisotropy=10,
            shift_asymmetry=0.1,
            alpha=0.5,
            beta=0.5,
            gamma=0.5,
        )

    def test_parse_json(self):
        good_json = {
            "isotope_symbol": "1H",
            "isotropic_chemical_shift": "0 ppm",
            "shielding_symmetric": {
                "anisotropy": "13.89 ppm",
                "asymmetry": 0.25,
            },
        }

        good_json_2 = {
            "isotope_symbol": "1H",
            "isotropic_chemical_shift": "0 ppm",
            "anisotropy": "13.89 ppm",
            "asymmetry": 0.25,
        }

        bad_json = {
            "isotope_symbol": "1H",
            "isotropic_chemical_shift": "0 ppm",
            "anisotropy": "13.89 Hz",
            "asymmetry": 0.25,
        }

        Site.parse_json_with_units(good_json)
        Site.parse_json_with_units(good_json_2)

        with self.assertRaises(Exception):
            Site.parse_json_with_units(bad_json)


class TestIsotopomer(unittest.TestCase):
    def test_direct_init(self):
        Isotopomer(sites=[], abundance=10)
        test_site = Site(
            nucleus="29Si",
            isotropic_chemical_shift=10,
            shift_anisotropy=10,
            shift_asymmetry=0.1,
            alpha=0.5,
            beta=0.5,
            gamma=0.5,
        )

        Isotopomer(sites=[test_site], abundance=10)
        Isotopomer(sites=[test_site, test_site], abundance=10)

    def test_parse__isojson(self):
        good_json = {"sites": [], "abundance": "10"}

        good_json2 = {
            "sites": [
                {
                    "isotope_symbol": "1H",
                    "isotropic_chemical_shift": "0 ppm",
                    "anisotropy": "13.89 ppm",
                    "asymmetry": 0.25,
                }
            ],
            "abundance": "10",
        }

        bad_json = {"sites": [], "abundance": "10 Hz"}

        Isotopomer.parse_json_with_units(good_json)
        Isotopomer.parse_json_with_units(good_json2)

        with self.assertRaises(Exception):
        	Isotopomer.parse_json_with_units(bad_json)


if __name__ == "__main__":
    unittest.main()
