# -*- coding: utf-8 -*-
"""Test for the base SpinSystem class."""
import numpy as np
import pytest
from mrsimulator import Coupling
from mrsimulator import Site
from mrsimulator import SpinSystem
from mrsimulator.spin_system import allowed_isotopes
from mrsimulator.spin_system.isotope import Isotope
from pydantic import ValidationError

__author__ = "Deepansh Srivastava"
__email__ = "srivastava.89@osu.edu"


def test_failure():
    site = Site(isotope="29Si")
    e = "All entries must be of type `Site`."
    with pytest.raises(ValueError, match=f".*{e}.*"):
        SpinSystem(sites=np.asarray(["a", site]))

    e = "All entries must be of type `Coupling`."
    with pytest.raises(ValueError, match=f".*{e}.*"):
        SpinSystem(sites=[site, site], couplings=np.asarray(["a"]))


def test_direct_init_spin_system():
    # test-1 # empty spin system
    the_spin_system = SpinSystem(sites=[], abundance=10, transition_pathways=None)

    assert the_spin_system.sites == []
    assert the_spin_system.abundance == 10.0
    assert the_spin_system.transition_pathways is None

    assert the_spin_system.json() == {"abundance": "10.0 %"}
    assert the_spin_system.reduced_dict() == {
        "sites": [],
        "abundance": 10.0,
    }

    # test-2 # site
    test_site = Site(isotope="29Si", isotropic_chemical_shift=10)

    assert test_site.isotope.symbol == "29Si"
    assert test_site.isotropic_chemical_shift == 10.0
    assert test_site.property_units["isotropic_chemical_shift"] == "ppm"

    assert test_site.json() == {
        "isotope": "29Si",
        "isotropic_chemical_shift": "10.0 ppm",
    }
    assert test_site.reduced_dict() == {
        "isotope": "29Si",
        "isotropic_chemical_shift": 10.0,
    }

    # test-3 # one site spin system
    the_spin_system = SpinSystem(sites=[test_site], abundance=10)
    assert isinstance(the_spin_system.sites[0], Site)
    assert the_spin_system.abundance == 10.0
    assert the_spin_system.json() == {
        "sites": [{"isotope": "29Si", "isotropic_chemical_shift": "10.0 ppm"}],
        "abundance": "10.0 %",
    }
    assert the_spin_system.reduced_dict() == {
        "sites": [{"isotope": "29Si", "isotropic_chemical_shift": 10.0}],
        "abundance": 10,
    }

    # test-4 # two sites spin system
    the_spin_system = SpinSystem(sites=[test_site, test_site], abundance=10)
    assert isinstance(the_spin_system.sites[0], Site)
    assert isinstance(the_spin_system.sites[1], Site)
    assert id(the_spin_system.sites[0]) != id(the_spin_system.sites[1])
    assert the_spin_system.abundance == 10.0
    assert the_spin_system.json() == {
        "sites": [
            {"isotope": "29Si", "isotropic_chemical_shift": "10.0 ppm"},
            {"isotope": "29Si", "isotropic_chemical_shift": "10.0 ppm"},
        ],
        "abundance": "10.0 %",
    }
    assert the_spin_system.reduced_dict() == {
        "sites": [
            {"isotope": "29Si", "isotropic_chemical_shift": 10.0},
            {"isotope": "29Si", "isotropic_chemical_shift": 10.0},
        ],
        "abundance": 10,
    }

    # test-5 # coupling
    test_coupling = Coupling(site_index=[0, 1], isotropic_j=10, dipolar={"D": 100})

    assert test_coupling.site_index == [0, 1]
    assert test_coupling.isotropic_j == 10.0
    assert test_coupling.property_units["isotropic_j"] == "Hz"

    assert test_coupling.dipolar.D == 100.0
    assert test_coupling.dipolar.property_units["D"] == "Hz"

    assert test_coupling.json() == {
        "site_index": [0, 1],
        "isotropic_j": "10.0 Hz",
        "dipolar": {"D": "100.0 Hz"},
    }
    assert test_coupling.reduced_dict() == {
        "site_index": [0, 1],
        "isotropic_j": 10.0,
        "dipolar": {"D": 100.0},
    }

    # test-6 #  two sites and one coupling spin system
    the_spin_system = SpinSystem(
        sites=[test_site, test_site], couplings=[test_coupling], abundance=10
    )
    assert isinstance(the_spin_system.sites[0], Site)
    assert isinstance(the_spin_system.sites[1], Site)
    assert isinstance(the_spin_system.couplings[0], Coupling)
    assert id(the_spin_system.sites[0]) != id(the_spin_system.sites[1])
    assert the_spin_system.abundance == 10.0
    assert the_spin_system.json() == {
        "sites": [
            {"isotope": "29Si", "isotropic_chemical_shift": "10.0 ppm"},
            {"isotope": "29Si", "isotropic_chemical_shift": "10.0 ppm"},
        ],
        "couplings": [
            {
                "site_index": [0, 1],
                "isotropic_j": "10.0 Hz",
                "dipolar": {"D": "100.0 Hz"},
            }
        ],
        "abundance": "10.0 %",
    }
    assert the_spin_system.reduced_dict() == {
        "sites": [
            {"isotope": "29Si", "isotropic_chemical_shift": 10.0},
            {"isotope": "29Si", "isotropic_chemical_shift": 10.0},
        ],
        "couplings": [
            {"site_index": [0, 1], "isotropic_j": 10.0, "dipolar": {"D": 100.0}}
        ],
        "abundance": 10,
    }

    # test-5
    the_spin_system = SpinSystem(
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
        couplings=[{"site_index": [0, 1], "isotropic_j": 34}],
        abundance=4.23,
    )
    assert the_spin_system.name == "Just a test"
    assert the_spin_system.description == "The same"
    assert the_spin_system.sites[0].isotope.symbol == "1H"
    assert the_spin_system.sites[0].isotropic_chemical_shift == 0
    assert the_spin_system.sites[1].isotope.symbol == "17O"
    assert the_spin_system.sites[1].isotropic_chemical_shift == -10
    assert the_spin_system.sites[1].quadrupolar.Cq == 5.1e6
    assert the_spin_system.sites[1].quadrupolar.eta == 0.5
    assert the_spin_system.couplings[0].site_index == [0, 1]
    assert the_spin_system.couplings[0].isotropic_j == 34.0
    assert the_spin_system.abundance == 4.23

    serialize = the_spin_system.json()
    assert serialize == {
        "name": "Just a test",
        "description": "The same",
        "sites": [
            {"isotope": "1H"},
            {
                "isotope": "17O",
                "isotropic_chemical_shift": "-10.0 ppm",
                "quadrupolar": {"Cq": "5100000.0 Hz", "eta": 0.5},
            },
        ],
        "couplings": [{"site_index": [0, 1], "isotropic_j": "34.0 Hz"}],
        "abundance": "4.23 %",
    }
    assert the_spin_system == SpinSystem.parse_dict_with_units(serialize)

    reduced_dict = the_spin_system.reduced_dict()
    assert reduced_dict == {
        "name": "Just a test",
        "description": "The same",
        "sites": [
            {"isotope": "1H", "isotropic_chemical_shift": 0},
            {
                "isotope": "17O",
                "isotropic_chemical_shift": -10.0,
                "quadrupolar": {"Cq": 5100000, "eta": 0.5},
            },
        ],
        "couplings": [{"site_index": [0, 1], "isotropic_j": 34.0}],
        "abundance": 4.23,
    }
    assert the_spin_system == SpinSystem(**reduced_dict)


def test_parse_json_spin_system():
    good_json = {"sites": [], "abundance": "10%"}

    good_json2 = {
        "sites": [{"isotope": "1H", "isotropic_chemical_shift": "0 ppm"}],
        "abundance": "10%",
    }

    # test-1
    iso1 = SpinSystem.parse_dict_with_units(good_json)
    assert len(iso1.sites) == 0
    assert iso1.abundance == 10
    assert iso1.json() == {"abundance": "10.0 %"}
    assert iso1.reduced_dict() == {
        "sites": [],
        "abundance": 10,
    }

    iso2 = SpinSystem.parse_dict_with_units(good_json2)
    assert len(iso2.sites) == 1
    assert iso2.sites[0].isotope.symbol == "1H"
    assert iso2.sites[0].isotropic_chemical_shift == 0
    assert iso2.abundance == 10
    assert iso2.json() == {"sites": [{"isotope": "1H"}], "abundance": "10.0 %"}
    assert iso2.reduced_dict() == {
        "sites": [{"isotope": "1H", "isotropic_chemical_shift": 0}],
        "abundance": 10,
    }

    bad_json = {"sites": [], "abundance": "10 Hz"}
    with pytest.raises(Exception):
        SpinSystem.parse_dict_with_units(bad_json)


def test_spin_system_methods():
    good_json2 = {
        "sites": [{"isotope": "1H", "isotropic_chemical_shift": "2 ppm"}],
        "abundance": "10%",
    }

    # Deprecated `to_freq_dict`
    # to_freq_dict()
    # iso1 = SpinSystem.parse_dict_with_units(good_json2).to_freq_dict(9.4)
    # result = {
    #     "name": None,
    #     "label": None,
    #     "description": None,
    #     "sites": [
    #         {
    #             "isotope": "1H",
    #             "isotropic_chemical_shift": -2 * 42.57748 * 9.4,  # -gamma * B0 * iso
    #             "name": None,
    #             "label": None,
    #             "description": None,
    #             "quadrupolar": None,
    #             "shielding_antisymmetric": None,
    #             "shielding_symmetric": None,
    #         }
    #     ],
    #     "abundance": 10,
    #     "transition_pathways": None,
    # }
    # assert iso1 == result

    # json()
    iso1 = SpinSystem.parse_dict_with_units(good_json2).json()
    result = {
        "sites": [{"isotope": "1H", "isotropic_chemical_shift": "2.0 ppm"}],
        "abundance": "10.0 %",
    }
    assert iso1 == result

    # reduced_dict()
    assert SpinSystem.parse_dict_with_units(good_json2).reduced_dict() == {
        "sites": [{"isotope": "1H", "isotropic_chemical_shift": 2.0}],
        "abundance": 10,
    }


def get_spin_system_list():
    isotopes = ["19F", "31P", "2H", "6Li", "14N", "27Al", "25Mg", "45Sc", "87Sr"]
    return SpinSystem(sites=[Site(isotope=item) for item in isotopes])


def generate_isotopes(isotopes):
    return [Isotope(symbol=item) for item in isotopes]


def test_get_isotopes():
    isotopes = ["19F", "31P", "2H", "6Li", "14N", "27Al", "25Mg", "45Sc", "87Sr"]
    spin_system = get_spin_system_list()
    assert spin_system.get_isotopes() == generate_isotopes(isotopes)
    assert spin_system.get_isotopes(spin_I=0.5) == generate_isotopes(["19F", "31P"])
    assert spin_system.get_isotopes(spin_I=1) == generate_isotopes(["2H", "6Li", "14N"])
    assert spin_system.get_isotopes(spin_I=1.5) == []
    assert spin_system.get_isotopes(spin_I=2.5) == generate_isotopes(["27Al", "25Mg"])
    assert spin_system.get_isotopes(spin_I=3.5) == generate_isotopes(["45Sc"])
    assert spin_system.get_isotopes(spin_I=4.5) == generate_isotopes(["87Sr"])

    assert spin_system.get_isotopes(symbol=True) == isotopes
    assert spin_system.get_isotopes(spin_I=0.5, symbol=True) == ["19F", "31P"]
    assert spin_system.get_isotopes(spin_I=1, symbol=True) == ["2H", "6Li", "14N"]
    assert spin_system.get_isotopes(spin_I=1.5, symbol=True) == []
    assert spin_system.get_isotopes(spin_I=2.5, symbol=True) == ["27Al", "25Mg"]
    assert spin_system.get_isotopes(spin_I=3.5, symbol=True) == ["45Sc"]
    assert spin_system.get_isotopes(spin_I=4.5, symbol=True) == ["87Sr"]


def test_allowed_isotopes():
    assert {"19F", "31P", "129Xe", "1H", "57Fe", "13C", "15N", "29Si"}.issubset(
        set(allowed_isotopes(spin_I=0.5))
    )
    assert {"2H", "6Li", "14N"}.issubset(set(allowed_isotopes(spin_I=1)))
    assert {"7Li", "9Be", "11B", "21Ne", "23Na", "33S", "37Cl", "41K"}.issubset(
        set(allowed_isotopes(spin_I=1.5))
    )
    assert {"17O", "25Mg", "27Al", "47Ti", "55Mn", "67Zn"}.issubset(
        set(allowed_isotopes(spin_I=2.5))
    )
    assert {"43Ca", "45Sc", "49Ti", "51V", "59Co", "123Sb", "133Cs"}.issubset(
        set(allowed_isotopes(spin_I=3.5))
    )
    assert {"73Ge", "83Kr", "87Sr", "93Nb", "113In"}.issubset(
        set(allowed_isotopes(spin_I=4.5))
    )


def test_bad_assignments():
    error = "value is not a valid list"
    with pytest.raises(ValidationError, match=f".*{error}.*"):
        SpinSystem(sites=Site())


def test_equality():
    a = SpinSystem(sites=[Site(isotope="1H")])
    b = SpinSystem(sites=[Site(isotope="1H")])
    assert a == b

    c = SpinSystem(sites=[Site(isotope="1H", isotropic_chemical_shift=16)])
    assert a != c
