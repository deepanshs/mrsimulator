# -*- coding: utf-8 -*-
"""Test for the base Simulator class."""
from random import randint

import pytest
from mrsimulator import Simulator
from mrsimulator import Site
from mrsimulator import SpinSystem
from mrsimulator.method.frequency_contrib import freq_default
from mrsimulator.methods import BlochDecaySpectrum
from mrsimulator.spin_system.tests.test_spin_systems import generate_isotopes

__author__ = "Deepansh Srivastava"
__email__ = "srivastava.89@osu.edu"


def test_simulator_assignments():
    a = Simulator()
    assert a.spin_systems == []

    error = "value is not a valid list"
    with pytest.raises(Exception, match=f".*{error}.*"):
        a.spin_systems = ""

    with pytest.raises(Exception, match=f".*{error}.*"):
        a.methods = ""


def test_equality():
    a = Simulator()
    b = Simulator()
    assert a == b

    assert a != {}

    c = Simulator(spin_systems=[SpinSystem()], label="test")
    assert a != c

    result = {
        "label": "test",
        "spin_systems": [{"abundance": "100 %", "sites": []}],
        "config": {
            "decompose_spectrum": "none",
            "integration_density": 70,
            "integration_volume": "octant",
            "number_of_sidebands": 64,
        },
    }
    assert c.json(include_methods=True) == result

    assert c.reduced_dict() == {
        "label": "test",
        "spin_systems": [{"abundance": 100, "sites": []}],
        "methods": [],
        "config": {
            "number_of_sidebands": 64,
            "integration_volume": "octant",
            "integration_density": 70,
            "decompose_spectrum": "none",
        },
    }


def get_simulator():
    isotopes = ["19F", "31P", "2H", "6Li", "14N", "27Al", "25Mg", "45Sc", "87Sr"]
    sites = []
    for isotope in isotopes:
        for _ in range(randint(1, 3)):
            sites.append(Site(isotope=isotope))
    sim = Simulator()
    sim.spin_systems.append(SpinSystem(sites=sites))
    return sim


def test_get_isotopes():
    isotopes = ["14N", "19F", "25Mg", "27Al", "2H", "31P", "45Sc", "6Li", "87Sr"]
    sim = get_simulator()
    assert sim.get_isotopes() == generate_isotopes(isotopes)
    assert sim.get_isotopes(spin_I=0.5) == generate_isotopes(["19F", "31P"])
    assert sim.get_isotopes(spin_I=1) == generate_isotopes(["14N", "2H", "6Li"])
    assert sim.get_isotopes(spin_I=1.5) == []
    assert sim.get_isotopes(spin_I=2.5) == generate_isotopes(["25Mg", "27Al"])
    assert sim.get_isotopes(spin_I=3.5) == generate_isotopes(["45Sc"])
    assert sim.get_isotopes(spin_I=4.5) == generate_isotopes(["87Sr"])

    assert sim.get_isotopes(symbol=True) == isotopes
    assert sim.get_isotopes(spin_I=0.5, symbol=True) == ["19F", "31P"]
    assert sim.get_isotopes(spin_I=1, symbol=True) == ["14N", "2H", "6Li"]
    assert sim.get_isotopes(spin_I=1.5, symbol=True) == []
    assert sim.get_isotopes(spin_I=2.5, symbol=True) == ["25Mg", "27Al"]
    assert sim.get_isotopes(spin_I=3.5, symbol=True) == ["45Sc"]
    assert sim.get_isotopes(spin_I=4.5, symbol=True) == ["87Sr"]


def test_simulator_1():
    sim = Simulator()
    sim.spin_systems = [SpinSystem(sites=[Site(isotope="1H"), Site(isotope="23Na")])]
    sim.methods = [BlochDecaySpectrum()]
    sim.name = "test"
    sim.label = "test0"
    sim.description = "testing-testing 1.2.3"

    red_dict = sim.reduced_dict()
    _ = [item.pop("description") for item in red_dict["methods"]]

    assert red_dict == {
        "name": "test",
        "label": "test0",
        "description": "testing-testing 1.2.3",
        "spin_systems": [
            {
                "sites": [
                    {"isotope": "1H", "isotropic_chemical_shift": 0},
                    {"isotope": "23Na", "isotropic_chemical_shift": 0},
                ],
                "abundance": 100,
            }
        ],
        "methods": [
            {
                "channels": ["1H"],
                "name": "BlochDecaySpectrum",
                "spectral_dimensions": [
                    {
                        "count": 1024,
                        "events": [
                            {
                                "fraction": 1.0,
                                "freq_contrib": freq_default,
                                "magnetic_flux_density": 9.4,
                                "rotor_angle": 0.955316618,
                                "rotor_frequency": 0.0,
                                "transition_query": {"P": {"channel-1": [[-1]]}},
                            }
                        ],
                        "reference_offset": 0.0,
                        "spectral_width": 25000.0,
                    }
                ],
            }
        ],
        "config": {
            "decompose_spectrum": "none",
            "integration_density": 70,
            "integration_volume": "octant",
            "number_of_sidebands": 64,
        },
    }

    # save
    sim.save("test_sim_save.temp")
    sim_load = sim.load("test_sim_save.temp")

    assert sim_load.spin_systems == sim.spin_systems
    assert sim_load.methods == sim.methods
    assert sim_load.name == sim.name
    assert sim_load.description == sim.description
    assert sim_load == sim

    # without units
    sim.save("test_sim_save_no_unit.temp", with_units=False)
    sim_load = sim.load("test_sim_save_no_unit.temp", parse_units=False)
    assert sim_load == sim
