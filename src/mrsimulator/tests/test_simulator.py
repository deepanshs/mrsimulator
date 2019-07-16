# coding: utf-8
"""
Tests for the base Parseable pattern
"""
import os.path
import unittest
from typing import ClassVar

from monty.serialization import loadfn
from mrsimulator import Site, Isotopomer, Spectrum, Simulator
from mrsimulator.tests import TEST_FOLDER


class TestSite(unittest.TestCase):
    def setUp(self):
        spectrum = Spectrum(
            number_of_points=1024,
            spectral_width=100,
            reference_offset=0,
            magnetic_flux_density=9.4,
            rotor_frequency=0,
            rotor_angle=0.9553,  # 54.935 degrees in radians
            rotor_phase=0,
            isotope="1H",
            spin=1,
            natural_abundance=0.04683,
            gyromagnetic_ratio=-8.465,
        )
        isotopomers = [
            Isotopomer(
                sites=[
                    Site(
                        isotope="29Si",
                        isotropic_chemical_shift=10,
                    )
                ],
                abundance=10)
        ]
        self.simulator = Simulator(isotopomers,spectrum)

    def test_allowed_isotopes(self):
        self.assertEqual(set(Simulator.allowed_isotopes()),
            {'19F', '31P', '129Xe', '1H', '57Fe', '13C', '15N', '29Si'})

    def test_all_isotopes(self):
        self.assertEqual(set(self.simulator.all_isotopes),
            {'29Si'})

    def test_valid_isotope_list(self):
        self.assertEqual(set(self.simulator.valid_isotope_list),
            {'29Si'})

    def test_one_d_spectrum(self):
        self.simulator.one_d_spectrum

    def test_run(self):
        pass

if __name__ == "__main__":
    unittest.main()
