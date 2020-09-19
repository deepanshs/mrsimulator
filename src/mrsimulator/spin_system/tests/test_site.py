# -*- coding: utf-8 -*-
"""Test for the base Site class."""
import pytest
from mrsimulator import Site
from mrsimulator.spin_system.isotope import Isotope
from pydantic import ValidationError


# Test Site ===========================================================================
def test_direct_init_site1():
    # test 1 --------------------------------------------------------------------------
    the_site = Site(isotope="29Si", isotropic_chemical_shift=10)
    assert the_site.isotope == Isotope(symbol="29Si")
    assert the_site.isotropic_chemical_shift == 10
    assert the_site.property_units["isotropic_chemical_shift"] == "ppm"

    assert the_site.shielding_antisymmetric is None
    assert the_site.quadrupolar is None
    assert the_site.shielding_symmetric is None

    # test 2 --------------------------------------------------------------------------
    the_site = Site(
        isotope="29Si",
        isotropic_chemical_shift=10,
        shielding_symmetric={"zeta": 12.1, "eta": 0.1},
    )
    assert the_site.isotope == Isotope(symbol="29Si")
    assert the_site.isotropic_chemical_shift == 10
    assert the_site.property_units["isotropic_chemical_shift"] == "ppm"

    assert the_site.shielding_antisymmetric is None

    assert the_site.quadrupolar is None

    assert the_site.shielding_symmetric.zeta == 12.1
    assert the_site.shielding_symmetric.property_units["zeta"] == "ppm"
    assert the_site.shielding_symmetric.eta == 0.1

    # test 3 --------------------------------------------------------------------------
    error = "ensure this value is less than or equal to 1"
    with pytest.raises(ValidationError, match=f".*{error}.*"):
        Site(
            isotope="29Si",
            isotropic_chemical_shift=10,
            shielding_symmetric={"zeta": 12.1, "eta": 1.5},
        )

    assert Site().isotope == Isotope(symbol="1H")

    error = ["with spin quantum number", "does not allow quadrupolar tensor"]
    with pytest.raises(ValidationError, match=".*{0}.*{1}.*".format(*error)):
        Site(quadrupolar={"Cq": 5.1e6})

    with pytest.raises(ValidationError, match=".*{0}.*{1}.*".format(*error)):
        Site.parse_dict_with_units(dict(quadrupolar={"Cq": "5.1 MHz"}))

    ax = Site.parse_dict_with_units({"isotope": "29Si"})
    assert ax.to_dict_with_units() == {
        "isotope": "29Si",
        "isotropic_chemical_shift": "0 ppm",
    }


def test_parse_json_site():
    # test 1 --------------------------------------------------------------------------
    good_json_site = {
        "isotope": "1H",
        "isotropic_chemical_shift": "0 ppm",
        "shielding_symmetric": {"zeta": "13.89 ppm", "eta": 0.25},
    }

    # testing method parse_dict_with_units()
    the_site = Site.parse_dict_with_units(good_json_site)
    assert the_site.isotope == Isotope(symbol="1H")
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

    # test 2 --------------------------------------------------------------------------
    good_json_2 = {
        "isotope": "14N",
        "isotropic_chemical_shift": "-10 ppm",
        "quadrupolar": {"Cq": "5.12 MHz", "eta": 0.5},
    }

    the_site = Site.parse_dict_with_units(good_json_2)
    assert the_site.isotope == Isotope(symbol="14N")
    assert the_site.isotropic_chemical_shift == -10
    assert the_site.property_units["isotropic_chemical_shift"] == "ppm"

    assert the_site.shielding_antisymmetric is None
    assert the_site.shielding_symmetric is None

    assert the_site.quadrupolar.Cq == 5120000.0
    assert the_site.quadrupolar.eta == 0.5

    # test 3 bad input ----------------------------------------------------------------
    bad_json = {"isotope": "1H", "isotropic_chemical_shift": "0 rad"}

    with pytest.raises(Exception):
        Site.parse_dict_with_units(bad_json)

    error = "Error enforcing units for isotropic_chemical_shift: 10 MHz"
    with pytest.raises(Exception, match=f".*{error}.*"):
        Site.parse_dict_with_units(
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
        "isotope": {"symbol": "14N"},
        "isotropic_chemical_shift": -10.0,
        "property_units": {"isotropic_chemical_shift": "ppm"},
        "name": None,
        "label": None,
        "description": None,
        "quadrupolar": None,
        "shielding_symmetric": None,
        "shielding_antisymmetric": None,
    }
    assert the_site.dict() == result, "Failed Site.dict()"

    # testing method to_freq_dict()
    result = {
        "isotope": "14N",
        "isotropic_chemical_shift": -1 * 3.077706 * 9.4 * -10.0,  # -gamma * B0 * iso
        "name": None,
        "label": None,
        "description": None,
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
        "name": None,
        "label": None,
        "description": None,
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


def test_site_quad_set_to_None():
    a = Site(isotope="27Al", quadrupolar=None)
    assert a.isotope.symbol == "27Al"
    assert a.quadrupolar is None


def test_equality():
    a = Site(isotope="1H")
    b = Site(isotope="1H")
    assert a == b

    c = Site(isotope="1H", isotropic_chemical_shift=16)
    assert a != c
