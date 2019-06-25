# coding: utf-8
"""
Tests for the base Parseable pattern
"""
import os.path
import unittest
from typing import ClassVar

from monty.serialization import loadfn
from mrsimulator import Site, Isotopomer, Spectrum
from mrsimulator.tests import TEST_FOLDER


class TestSite(unittest.TestCase):
    def test_direct_init(self):
        Site(nucleus="29Si", isotropic_chemical_shift=10)

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
        }

        bad_json = {
            "isotope_symbol": "1H",
            "isotropic_chemical_shift": "0 rad",
        }

        Site.parse_json_with_units(good_json)
        Site.parse_json_with_units(good_json_2)

        with self.assertRaises(Exception):
            Site.parse_json_with_units(bad_json)


class TestIsotopomer(unittest.TestCase):
    def test_direct_init(self):
        Isotopomer(sites=[], abundance=10)
        test_site = Site(nucleus="29Si", isotropic_chemical_shift=10)

        Isotopomer(sites=[test_site], abundance=10)
        Isotopomer(sites=[test_site, test_site], abundance=10)

    def test_parse_json(self):
        good_json = {"sites": [], "abundance": "10"}

        good_json2 = {
            "sites": [{
                "isotope_symbol": "1H",
                "isotropic_chemical_shift": "0 ppm"
            }],
            "abundance": "10",
        }

        bad_json = {"sites": [], "abundance": "10 Hz"}

        Isotopomer.parse_json_with_units(good_json)
        Isotopomer.parse_json_with_units(good_json2)

        with self.assertRaises(Exception):
            Isotopomer.parse_json_with_units(bad_json)


class TestSpectrum(unittest.TestCase):
    def test_direct_init(self):
        Spectrum()
        Spectrum(
            number_of_points=1024,
            spectral_width=100,
            reference_offset=0,
            magnetic_flux_density=9.4,
            rotor_frequency=0,
            rotor_angle=0.9553,  # 54.935 degrees in radians
            rotor_phase=0,
            nucleus="1H",
            spin=1,
            natural_abundance=0.04683,
            gyromagnetic_ratio=-8.465,
        )

    def test_parse_json(self):
        good_json = {
            "number_of_points": 1024,
            "spectral_width": "100 Hz",
            "reference_offset": "0 Hz",
            "magnetic_flux_density": "9.4 T",
            "rotor_frequency": "0 Hz",
            "rotor_angle": "0.9553 rad",  # 54.935 degrees in radians
            "rotor_phase": "0 rad",
            "nucleus": "1H",
        }

        spec = Spectrum.parse_json_with_units(good_json)
        self.assertIn("spin", spec.dict())
        self.assertEqual(spec.spin, 1)


class TestJSONData(unittest.TestCase):
    def setUp(self):

        self.mas_data = loadfn(os.path.join(TEST_FOLDER, "mas.json"))
        self.static_data = loadfn(os.path.join(TEST_FOLDER, "static.json"))

    def test_parsing(self):
        Spectrum.parse_json_with_units(self.mas_data["spectrum"])
        Spectrum.parse_json_with_units(self.static_data["spectrum"])

        [Isotopomer.parse_json_with_units(isotopomer) for isotopomer in self.mas_data["isotopomers"]]

        [Isotopomer.parse_json_with_units(isotopomer) for isotopomer in self.static_data["isotopomers"]]


if __name__ == "__main__":
    unittest.main()
