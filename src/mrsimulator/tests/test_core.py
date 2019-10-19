# -*- coding: utf-8 -*-
"""
Tests for the base Parseable pattern
"""
# import os.path
import pytest
import numpy as np

# from typing import ClassVar

from mrsimulator import Site, Isotopomer, Dimension

# from mrsimulator.examples import mas_data, static_data


# Test Site
def test_direct_init_site():
    the_site = Site(isotope="29Si", isotropic_chemical_shift=10)
    assert the_site.isotope == "29Si"
    assert the_site.isotropic_chemical_shift == 10
    # assert the_site.property_units["isotropic_chemical_shift"] == "Hz"
    assert the_site.shielding_antisymmetric is None
    assert the_site.quadrupolar is None
    assert the_site.shielding_symmetric is None


def test_parse_json_site():
    good_json = {
        "isotope": "1H",
        "isotropic_chemical_shift": "0 ppm",
        "shielding_symmetric": {"zeta": "13.89 ppm", "eta": 0.25},
    }

    the_site = Site.parse_dict_with_units(good_json)
    assert the_site.isotope == "1H"
    assert the_site.isotropic_chemical_shift == 0
    assert the_site.property_units["isotropic_chemical_shift"] == "ppm"
    assert the_site.shielding_antisymmetric is None
    assert the_site.quadrupolar is None
    assert the_site.shielding_symmetric.zeta == 13.89
    assert the_site.shielding_symmetric.property_units["zeta"] == "ppm"
    assert the_site.shielding_symmetric.eta == 0.25
    assert the_site.shielding_symmetric.alpha is None
    assert the_site.shielding_symmetric.beta is None
    assert the_site.shielding_symmetric.gamma is None

    good_json_2 = {"isotope": "14N", "isotropic_chemical_shift": "-10 ppm"}

    the_site = Site.parse_dict_with_units(good_json_2)
    assert the_site.isotope == "14N"
    assert the_site.isotropic_chemical_shift == -10
    assert the_site.property_units["isotropic_chemical_shift"] == "ppm"
    assert the_site.shielding_antisymmetric is None
    assert the_site.quadrupolar is None
    assert the_site.shielding_symmetric is None

    result = {
        "isotope": "14N",
        "isotropic_chemical_shift": -10.0,
        "property_units": {"isotropic_chemical_shift": "ppm"},
        "quadrupolar": None,
        "shielding_symmetric": None,
        "shielding_antisymmetric": None,
    }
    assert the_site.dict() == result

    bad_json = {"isotope": "1H", "isotropic_chemical_shift": "0 rad"}

    with pytest.raises(Exception):
        Site.parse_dict_with_units(bad_json)


def test_direct_init_isotopomer():
    the_isotopomer = Isotopomer(sites=[], abundance=10)

    assert the_isotopomer.sites == []
    assert the_isotopomer.abundance == 10.0

    test_site = Site(isotope="29Si", isotropic_chemical_shift=10)

    assert test_site.isotope == "29Si"
    assert test_site.isotropic_chemical_shift == 10.0
    # assert test_site.property_units["isotropic_chemical_shift"] == "Hz"

    the_isotopomer = Isotopomer(sites=[test_site], abundance=10)
    assert isinstance(the_isotopomer.sites[0], Site)
    assert the_isotopomer.abundance == 10.0

    the_isotopomer = Isotopomer(sites=[test_site, test_site], abundance=10)
    assert isinstance(the_isotopomer.sites[0], Site)
    assert isinstance(the_isotopomer.sites[1], Site)
    assert id(the_isotopomer.sites[0] != the_isotopomer.sites[1])
    assert the_isotopomer.abundance == 10.0

    # This should raise an error, but it is not.
    the_isotopomer.sites[0] = "Trash"
    assert the_isotopomer.sites[0] == "Trash"


def test_parse_json_isotopomer():
    good_json = {"sites": [], "abundance": "10"}

    good_json2 = {
        "sites": [{"isotope": "1H", "isotropic_chemical_shift": "0 ppm"}],
        "abundance": "10",
    }

    bad_json = {"sites": [], "abundance": "10 Hz"}

    Isotopomer.parse_dict_with_units(good_json)
    Isotopomer.parse_dict_with_units(good_json2)

    with pytest.raises(Exception):
        Isotopomer.parse_dict_with_units(bad_json)


def test_direct_init_spectrum():
    the_spectrum = Dimension(number_of_points=1024, spectral_width=100)
    assert the_spectrum.number_of_points == 1024
    assert the_spectrum.spectral_width == 100
    # assert the_spectrum.property_units["spectral_width"] == 'Hz'

    Dimension(
        number_of_points=1024,
        spectral_width=100,
        reference_offset=0,
        magnetic_flux_density=9.4,
        rotor_frequency=0,
        rotor_angle=0.9553,  # 54.935 degrees in radians
        rotor_phase=0,
        isotope="1H",
    )


def test_parse_json_spectrum():
    good_json = {
        "number_of_points": 1024,
        "spectral_width": "100 Hz",
        "reference_offset": "0 Hz",
        "magnetic_flux_density": "9.4 T",
        "rotor_frequency": "0 Hz",
        "rotor_angle": "54.935 degree",  # 54.935 degrees in radians
        "rotor_phase": "0 rad",
        "isotope": "1H",
    }

    spec = Dimension.parse_dict_with_units(good_json)
    assert spec.spin == 1
    assert spec.isotope == "1H"
    assert np.allclose(spec.rotor_angle, 0.95879662)


# def test_parsing(mas_data, static_data):
#     mas = Dimension.parse_dict_with_units(mas_data["spectrum"])
#     static = Dimension.parse_dict_with_units(static_data["spectrum"])

#     [
#         Isotopomer.parse_dict_with_units(isotopomer)
#         for isotopomer in mas_data["isotopomers"]
#     ]

#     [
#         Isotopomer.parse_dict_with_units(isotopomer)
#         for isotopomer in static_data["isotopomers"]
#     ]

#     assert static.rotor_frequency == 0
#     assert mas.rotor_frequency == 1000
