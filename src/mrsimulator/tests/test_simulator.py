# -*- coding: utf-8 -*-
"""Test for the base Simulator class."""
from random import randint

import numpy as np
import pytest
from mrsimulator import Dimension
from mrsimulator import Isotopomer
from mrsimulator import Simulator
from mrsimulator import Site
from mrsimulator.methods import one_d_spectrum


def test_simulator_assignments():
    a = Simulator()
    assert a.isotopomers == []

    error = "value is not a valid list"
    with pytest.raises(Exception, match=".*{0}.*".format(error)):
        a.isotopomers = ""

    with pytest.raises(Exception, match=".*{0}.*".format(error)):
        a.dimensions = ""


def test_equality():
    a = Simulator()
    b = Simulator()
    assert a == b

    assert a is not {}

    c = Simulator(isotopomers=[Isotopomer()])
    assert a is not c

    result = {
        "isotopomers": [
            {"abundance": "100%", "description": "", "name": "", "sites": []}
        ]
    }
    assert c.to_dict_with_units(include_dimensions=True) == result

    result["dimensions"] = [
        {
            "label": "",
            "magnetic_flux_density": "9.4 T",
            "number_of_points": 1024,
            "reference_offset": "0 Hz",
            "rotor_angle": "0.9553166 rad",
            "rotor_frequency": "0 Hz",
            "spectral_width": "10.0 Hz",
        }
    ]
    c.dimensions = [Dimension(spectral_width=10)]
    assert c.to_dict_with_units(include_dimensions=True) == result


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


def test_csdm_object():
    sim = Simulator()
    site = Site(
        isotope="27Al",
        isotropic_chemical_shift=120,
        shielding_symmetric={"zeta": 2.1, "eta": 0.1},
        quadrupole={"Cq": 5.1e6, "eta": 0.5},
    )
    sim.isotopomers = [{"name": "test", "description": "awesome", "sites": [site]}]
    sim.dimensions = [
        {
            "number_of_points": 1024,
            "spectral_width": 100,
            "reference_offset": 0,
            "magnetic_flux_density": 9.4,
            "rotor_frequency": 0,
            "rotor_angle": 0.9553166,
            "isotope": "27Al",
        }
    ]
    x, y = sim.run(method=one_d_spectrum)
    csdm_obj = sim.as_csdm_object()

    assert np.allclose(csdm_obj.dependent_variables[0].components[0], y)
    assert csdm_obj.dimensions[0].count == x.size
