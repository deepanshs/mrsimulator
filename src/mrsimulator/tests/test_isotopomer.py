# -*- coding: utf-8 -*-
"""Test for the base Isotopomers class."""
from random import randint

import pytest
from mrsimulator import Isotopomer
from mrsimulator import Site
from mrsimulator.isotopomer import allowed_isotopes
from pydantic import ValidationError


def test_direct_init_isotopomer():
    the_isotopomer = Isotopomer(sites=[], abundance=10)

    assert the_isotopomer.sites == []
    assert the_isotopomer.abundance == 10.0

    test_site = Site(isotope="29Si", isotropic_chemical_shift=10)

    assert test_site.isotope == "29Si"
    assert test_site.isotropic_chemical_shift == 10.0
    assert test_site.property_units["isotropic_chemical_shift"] == "ppm"

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
                "isotropic_chemical_shift": -2 * 42.57748 * 9.4,  # -gamma * B0 * iso
                "name": None,
                "quadrupolar": None,
                "shielding_antisymmetric": None,
                "shielding_symmetric": None,
            }
        ],
        "abundance": 10,
        "transitions": None,
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


def get_isotopomer_list():
    isotopes = ["19F", "31P", "2H", "6Li", "14N", "27Al", "25Mg", "45Sc", "87Sr"]
    sites = []
    for isotope in isotopes:
        for _ in range(randint(1, 3)):
            sites.append(Site(isotope=isotope))
    return Isotopomer(sites=sites)


def test_get_isotopes():
    isotopes = {"19F", "31P", "2H", "6Li", "14N", "27Al", "25Mg", "45Sc", "87Sr"}
    isotopomer = get_isotopomer_list()
    assert isotopomer.get_isotopes() == isotopes
    assert isotopomer.get_isotopes(I=0.5) == {"19F", "31P"}
    assert isotopomer.get_isotopes(I=1) == {"2H", "6Li", "14N"}
    assert isotopomer.get_isotopes(I=1.5) == set()
    assert isotopomer.get_isotopes(I=2.5) == {"27Al", "25Mg"}
    assert isotopomer.get_isotopes(I=3.5) == {"45Sc"}
    assert isotopomer.get_isotopes(I=4.5) == {"87Sr"}


def test_allowed_isotopes():
    assert {"19F", "31P", "129Xe", "1H", "57Fe", "13C", "15N", "29Si"}.issubset(
        set(allowed_isotopes(I=0.5))
    )
    assert {"2H", "6Li", "14N"}.issubset(set(allowed_isotopes(I=1)))
    assert {"7Li", "9Be", "11B", "21Ne", "23Na", "33S", "37Cl", "41K"}.issubset(
        set(allowed_isotopes(I=1.5))
    )
    assert {"17O", "25Mg", "27Al", "47Ti", "55Mn", "67Zn"}.issubset(
        set(allowed_isotopes(I=2.5))
    )
    assert {"43Ca", "45Sc", "49Ti", "51V", "59Co", "123Sb", "133Cs"}.issubset(
        set(allowed_isotopes(I=3.5))
    )
    assert {"73Ge", "83Kr", "87Sr", "93Nb", "113In"}.issubset(
        set(allowed_isotopes(I=4.5))
    )


def test_bad_assignments():
    error = "value is not a valid list"
    with pytest.raises(ValidationError, match=".*{0}.*".format(error)):
        Isotopomer(sites=Site())


def test_equality():
    a = Isotopomer(sites=[Site(isotope="1H")])
    b = Isotopomer(sites=[Site(isotope="1H")])
    assert a == b

    c = Isotopomer(sites=[Site(isotope="1H", isotropic_chemical_shift=16)])
    assert a != c
