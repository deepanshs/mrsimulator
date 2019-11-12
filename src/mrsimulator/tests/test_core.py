# -*- coding: utf-8 -*-
"""
Tests for the base Parseable pattern
"""
# import os.path
import numpy as np
import pytest
from mrsimulator import Dimension
from mrsimulator import Isotopomer
from mrsimulator import Site
from pydantic import ValidationError


# Test Site ===========================================================================
def test_direct_init_site1():
    # test 1 --------------------------------------------------------------------------
    the_site = Site(isotope="29Si", isotropic_chemical_shift=10)
    assert the_site.isotope == "29Si"
    assert the_site.isotropic_chemical_shift == 10
    assert the_site.property_units["isotropic_chemical_shift"] == "ppm"

    assert the_site.shielding_antisymmetric is None
    assert the_site.quadrupolar is None
    assert the_site.shielding_symmetric is None

    assert the_site.atomic_number == 14
    assert the_site.gyromagnetic_ratio == -8.465499
    assert the_site.natural_abundance == 4.683
    assert the_site.quadrupole_moment == 0.0
    assert the_site.spin == 0.5

    # test 2 --------------------------------------------------------------------------
    the_site = Site(
        isotope="29Si",
        isotropic_chemical_shift=10,
        shielding_symmetric={"zeta": 12.1, "eta": 0.1},
    )
    assert the_site.isotope == "29Si"
    assert the_site.isotropic_chemical_shift == 10
    assert the_site.property_units["isotropic_chemical_shift"] == "ppm"

    assert the_site.shielding_antisymmetric is None

    assert the_site.quadrupolar is None

    assert the_site.shielding_symmetric.zeta == 12.1
    assert the_site.shielding_symmetric.property_units["zeta"] == "ppm"
    assert the_site.shielding_symmetric.eta == 0.1

    # site properties
    assert the_site.atomic_number == 14
    assert the_site.gyromagnetic_ratio == -8.465499
    assert the_site.natural_abundance == 4.683
    assert the_site.quadrupole_moment == 0.0
    assert the_site.spin == 0.5

    # test 3 --------------------------------------------------------------------------
    error = "ensure this value is less than or equal to 1"
    with pytest.raises(ValidationError, match=".*{0}.*".format(error)):
        the_site = Site(
            isotope="29Si",
            isotropic_chemical_shift=10,
            shielding_symmetric={"zeta": 12.1, "eta": 1.5},
        )

    error = r"Site\nisotope\n  field required"
    with pytest.raises(ValidationError, match=".*{0}.*".format(error)):
        Site()


def test_parse_json_site():
    # test 1 --------------------------------------------------------------------------
    good_json_site = {
        "isotope": "1H",
        "isotropic_chemical_shift": "0 ppm",
        "shielding_symmetric": {"zeta": "13.89 ppm", "eta": 0.25},
    }

    # testing method parse_dict_with_units()
    the_site = Site.parse_dict_with_units(good_json_site)
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

    # site properties
    assert the_site.atomic_number == 1
    assert the_site.gyromagnetic_ratio == 42.57748
    assert the_site.natural_abundance == 99.985
    assert the_site.quadrupole_moment == 0.0
    assert the_site.spin == 0.5

    # test 2 --------------------------------------------------------------------------
    good_json_2 = {
        "isotope": "14N",
        "isotropic_chemical_shift": "-10 ppm",
        "quadrupolar": {"Cq": "5.12 MHz", "eta": 0.5},
    }

    the_site = Site.parse_dict_with_units(good_json_2)
    assert the_site.isotope == "14N"
    assert the_site.isotropic_chemical_shift == -10
    assert the_site.property_units["isotropic_chemical_shift"] == "ppm"

    assert the_site.shielding_antisymmetric is None
    assert the_site.shielding_symmetric is None

    assert the_site.quadrupolar.Cq == 5120000.0
    assert the_site.quadrupolar.eta == 0.5

    # site properties
    assert the_site.atomic_number == 7
    assert the_site.gyromagnetic_ratio == 3.077706
    assert the_site.natural_abundance == 99.634
    assert the_site.quadrupole_moment == 0.0193
    assert the_site.spin == 1

    # test 3 bad input ----------------------------------------------------------------
    bad_json = {"isotope": "1H", "isotropic_chemical_shift": "0 rad"}

    with pytest.raises(Exception):
        Site.parse_dict_with_units(bad_json)

    error = "Error enforcing units for isotropic_chemical_shift: 10 MHz"
    with pytest.raises(Exception, match=".*{0}.*".format(error)):
        the_site = Site.parse_dict_with_units(
            {
                "isotope": "29Si",
                "isotropic_chemical_shift": "10 MHz",
                "shielding_symmetric": {"zeta": "12.1 ppm", "eta": 0.5},
            }
        )


def test_site_object_methods():
    good_json_2 = {"isotope": "14N", "isotropic_chemical_shift": "-10 ppm"}
    the_site = Site.parse_dict_with_units(good_json_2)

    # testing method dict()
    result = {
        "isotope": "14N",
        "isotropic_chemical_shift": -10.0,
        "property_units": {"isotropic_chemical_shift": "ppm"},
        "quadrupolar": None,
        "shielding_symmetric": None,
        "shielding_antisymmetric": None,
    }
    assert the_site.dict() == result, "Failed Site.dict()"

    # testing method to_freq_dict()
    result = {
        "isotope": "14N",
        "isotropic_chemical_shift": -1 * 3.077706 * 9.4 * -10.0,  # -gamma * B0 * iso
        "quadrupolar": None,
        "shielding_symmetric": None,
        "shielding_antisymmetric": None,
    }
    assert the_site.to_freq_dict(B0=9.4) == result, "Failed Site.to_freq_dict()"

    # testing method to_dict_with_units()
    result = {"isotope": "14N", "isotropic_chemical_shift": "-10.0 ppm"}
    assert the_site.to_dict_with_units() == result, "Failed Site.to_dict_with_units()"

    result = {
        "isotope": "27Al",
        "isotropic_chemical_shift": "10.0 ppm",
        "shielding_symmetric": {"zeta": "12.1 ppm", "eta": 0.1, "alpha": "2.1 rad"},
        "shielding_antisymmetric": {
            "zeta": "-1.1 ppm",
            "alpha": "0.1 rad",
            "beta": "2.5 rad",
        },
        "quadrupolar": {"Cq": "10000000.0 Hz", "eta": 0.6},
    }
    the_site = Site(
        isotope="27Al",
        isotropic_chemical_shift=10,
        shielding_symmetric={"zeta": 12.1, "eta": 0.1, "alpha": 2.1},
        shielding_antisymmetric={"zeta": -1.1, "alpha": 0.1, "beta": 2.5},
        quadrupolar={"Cq": 10e6, "eta": 0.6},
    )
    assert the_site.to_dict_with_units() == result, "Failed Site.to_dict_with_units()"

    larmor_freq = -1 * 11.10309 * 9.4  # -gamma * B0
    result = {
        "isotope": "27Al",
        "isotropic_chemical_shift": 10.0 * larmor_freq,  # larmor_freq * iso
        "shielding_symmetric": {
            "zeta": 12.1 * larmor_freq,
            "eta": 0.1,
            "alpha": 2.1,
            "beta": None,
            "gamma": None,
        },
        "shielding_antisymmetric": {
            "zeta": -1.1 * larmor_freq,
            "alpha": 0.1,
            "beta": 2.5,
        },
        "quadrupolar": {
            "Cq": 10000000.0,
            "zeta": None,
            "eta": 0.6,
            "alpha": None,
            "beta": None,
            "gamma": None,
        },
    }
    assert the_site.to_freq_dict(9.4) == result, "Failed Site.to_dict_with_units()"


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
    assert id(the_isotopomer.sites[0]) != id(the_isotopomer.sites[1])
    assert the_isotopomer.abundance == 10.0

    the_isotopomer = Isotopomer(
        name="Just a test",
        description="The same",
        sites=[
            {"isotope": "1H", "isotropic_chemical_shift": 0},
            {
                "isotope": "17O",
                "isotropic_chemical_shift": -10,
                "quadrupolar": {"Cq": 5.1e6, "eta": 0.5},
            },
        ],
        abundance=4.23,
    )
    assert the_isotopomer.name == "Just a test"
    assert the_isotopomer.description == "The same"
    assert the_isotopomer.sites[0].isotope == "1H"
    assert the_isotopomer.sites[0].isotropic_chemical_shift == 0
    assert the_isotopomer.sites[1].isotope == "17O"
    assert the_isotopomer.sites[1].isotropic_chemical_shift == -10
    assert the_isotopomer.sites[1].quadrupolar.Cq == 5.1e6
    assert the_isotopomer.sites[1].quadrupolar.eta == 0.5
    assert the_isotopomer.abundance == 4.23


def test_parse_json_isotopomer():
    good_json = {"sites": [], "abundance": "10%"}

    good_json2 = {
        "sites": [{"isotope": "1H", "isotropic_chemical_shift": "0 ppm"}],
        "abundance": "10%",
    }

    iso1 = Isotopomer.parse_dict_with_units(good_json)
    assert len(iso1.sites) == 0
    assert iso1.abundance == 10

    iso2 = Isotopomer.parse_dict_with_units(good_json2)
    assert len(iso2.sites) == 1
    assert iso2.sites[0].isotope == "1H"
    assert iso2.sites[0].isotropic_chemical_shift == 0
    assert iso2.abundance == 10

    bad_json = {"sites": [], "abundance": "10 Hz"}
    with pytest.raises(Exception):
        Isotopomer.parse_dict_with_units(bad_json)


def test_isotopomer_methods():
    good_json2 = {
        "sites": [{"isotope": "1H", "isotropic_chemical_shift": "2 ppm"}],
        "abundance": "10%",
    }

    # to_freq_dict()
    iso1 = Isotopomer.parse_dict_with_units(good_json2).to_freq_dict(9.4)
    result = {
        "name": "",
        "description": "",
        "sites": [
            {
                "isotope": "1H",
                "isotropic_chemical_shift": -1
                * 42.57748
                * 9.4
                * 2,  # -gamma * B0 * iso
                "quadrupolar": None,
                "shielding_antisymmetric": None,
                "shielding_symmetric": None,
            }
        ],
        "abundance": 10,
    }
    assert iso1 == result

    # to_dict_with_units()
    iso1 = Isotopomer.parse_dict_with_units(good_json2).to_dict_with_units()
    result = {
        "name": "",
        "description": "",
        "sites": [{"isotope": "1H", "isotropic_chemical_shift": "2.0 ppm"}],
        "abundance": "10.0%",
    }
    assert iso1 == result


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
    assert spec.gyromagnetic_ratio == 42.57748
    assert spec.natural_abundance == 99.985
    assert spec.quadrupole_moment == 0.0
    assert spec.atomic_number == 1
    assert np.allclose(spec.rotor_angle, 0.95879662)


def test_parse_json_spectrum2():
    good_json = {
        "number_of_points": 1024,
        "spectral_width": "100 Hz",
        "reference_offset": "0 Hz",
        "magnetic_flux_density": "9.4 T",
        "rotor_frequency": "0 Hz",
        "rotor_angle": "54.935 degree",  # 54.935 degrees in radians
        "rotor_phase": "0 rad",
        "isotope": "29Si",
    }

    spec = Dimension.parse_dict_with_units(good_json)
    assert spec.spin == 1
    assert spec.isotope == "29Si"
    assert spec.gyromagnetic_ratio == -8.465499
    assert spec.natural_abundance == 4.683
    assert spec.quadrupole_moment == 0.0
    assert spec.atomic_number == 14
    assert np.allclose(spec.rotor_angle, 0.95879662)
