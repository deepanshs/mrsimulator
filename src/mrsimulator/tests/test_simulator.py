# -*- coding: utf-8 -*-
"""
Tests for the base Parseable pattern
"""
# import os.path
import pytest
from mrsimulator import Isotopomer
from mrsimulator import Simulator
from mrsimulator import Site
from random import randint
from pydantic import ValidationError


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


def test_config():
    a = Simulator()

    # number of sidebands
    assert a.config.number_of_sidebands == 64
    a.config.number_of_sidebands = 10
    assert a.config.number_of_sidebands == 10

    # integration density
    assert a.config.integration_density == 70
    a.config.integration_density = 20
    assert a.config.integration_density == 20

    # integration volume
    assert a.config.integration_volume == "octant"
    a.config.integration_volume = "hemisphere"
    assert a.config.integration_volume == "hemisphere"

    error = "value is not a valid enumeration member; permitted: 'octant', 'hemisphere'"
    with pytest.raises(ValidationError, match=".*{0}.*".format(error)):
        a.config.integration_volume = "sphere"

    assert a.config.dict() == {
        "number_of_sidebands": 10,
        "integration_volume": "hemisphere",
        "integration_density": 20,
    }
