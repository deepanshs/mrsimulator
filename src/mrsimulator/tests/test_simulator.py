# -*- coding: utf-8 -*-
"""Test for the base Simulator class."""
from random import randint

import pytest
from mrsimulator import Isotopomer
from mrsimulator import Simulator
from mrsimulator import Site


def test_simulator_assignments():
    a = Simulator()
    assert a.isotopomers == []

    error = "value is not a valid list"
    with pytest.raises(Exception, match=".*{0}.*".format(error)):
        a.isotopomers = ""

    with pytest.raises(Exception, match=".*{0}.*".format(error)):
        a.methods = ""


def test_equality():
    a = Simulator()
    b = Simulator()
    assert a == b

    assert a is not {}

    c = Simulator(isotopomers=[Isotopomer()])
    assert a is not c

    result = {
        "name": "",
        "description": "",
        "isotopomers": [{"abundance": "100 %", "sites": []}],
        "config": {
            "decompose": False,
            "integration_density": 70,
            "integration_volume": "octant",
            "number_of_sidebands": 64,
        },
        "indexes": [],
    }
    assert c.to_dict_with_units(include_methods=True) == result


def get_simulator():
    isotopes = ["19F", "31P", "2H", "6Li", "14N", "27Al", "25Mg", "45Sc", "87Sr"]
    sites = []
    for isotope in isotopes:
        for _ in range(randint(1, 3)):
            sites.append(Site(isotope=isotope))
    sim = Simulator()
    sim.isotopomers.append(Isotopomer(sites=sites))
    return sim


def test_get_isotopes():
    isotopes = {"19F", "31P", "2H", "6Li", "14N", "27Al", "25Mg", "45Sc", "87Sr"}
    sim = get_simulator()
    assert sim.get_isotopes() == isotopes
    assert sim.get_isotopes(I=0.5) == {"19F", "31P"}
    assert sim.get_isotopes(I=1) == {"2H", "6Li", "14N"}
    assert sim.get_isotopes(I=1.5) == set()
    assert sim.get_isotopes(I=2.5) == {"27Al", "25Mg"}
    assert sim.get_isotopes(I=3.5) == {"45Sc"}
    assert sim.get_isotopes(I=4.5) == {"87Sr"}
