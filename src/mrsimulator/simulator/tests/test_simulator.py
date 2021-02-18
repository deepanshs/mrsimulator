# -*- coding: utf-8 -*-
"""Test for the base Simulator class."""
from random import randint

import numpy as np
import pytest
from mrsimulator import Simulator
from mrsimulator import Site
from mrsimulator import SpinSystem
from mrsimulator.method.frequency_contrib import freq_default
from mrsimulator.methods import BlochDecaySpectrum
from mrsimulator.simulator import __CPU_count__
from mrsimulator.simulator import get_chunks
from mrsimulator.simulator import Sites
from mrsimulator.utils.collection import single_site_system_generator


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
    isotopes = {"19F", "31P", "2H", "6Li", "14N", "27Al", "25Mg", "45Sc", "87Sr"}
    sim = get_simulator()
    assert sim.get_isotopes() == isotopes
    assert sim.get_isotopes(spin_I=0.5) == {"19F", "31P"}
    assert sim.get_isotopes(spin_I=1) == {"2H", "6Li", "14N"}
    assert sim.get_isotopes(spin_I=1.5) == set()
    assert sim.get_isotopes(spin_I=2.5) == {"27Al", "25Mg"}
    assert sim.get_isotopes(spin_I=3.5) == {"45Sc"}
    assert sim.get_isotopes(spin_I=4.5) == {"87Sr"}


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


def test_sites():
    iso = [1.02, 2.12, 13.2, 5.2, 2.1, 1.2]
    zeta = [1.02, 2.12, 13.2, 5.2, 2.1, 1.2]
    eta = [0.1, 0.4, 0.3, 0.6, 0.9, 1.0]
    sites = [
        Site(
            isotope="13C",
            isotropic_chemical_shift=i,
            shielding_symmetric={"zeta": z, "eta": e},
        )
        for i, z, e in zip(iso, zeta, eta)
    ]
    sim = Simulator()
    sim.spin_systems = [SpinSystem(sites=[s]) for s in sites]
    r_sites = sim.sites()
    for i in range(len(sites)):
        assert r_sites[i] == sites[i]

    # test sites to pd
    sites_table = sim.sites().to_pd()

    assert list(sites_table["isotope"]) == ["13C"] * len(iso)
    assert list(sites_table["isotropic_chemical_shift"]) == [
        f"{i} ppm" if i is not None else None for i in iso
    ]
    assert list(sites_table["shielding_symmetric.zeta"]) == [
        f"{i} ppm" if i is not None else None for i in zeta
    ]
    assert list(sites_table["shielding_symmetric.eta"]) == [
        i if i is not None else None for i in eta
    ]

    # test Sites Class
    a = Sites([])

    site = Site(isotope="1H")
    a.append(site)
    assert a[0] == site

    site2 = Site(isotope="17O")
    a[0] = site2
    assert a[0] == site2

    site_dict = {"isotope": "13C"}
    a[0] = site_dict
    assert a[0] == Site(**site_dict)

    with pytest.raises(ValueError, match="Only object of type Site is allowed."):
        a[0] = ""


def test_sites_to_pandas_df():
    isotopes = ["29Si"] * 3 + ["17O"]
    shifts = [-89.0, -89.5, -87.8, 15.0]
    zeta = [59.8, 52.1, 69.4, 12.4]
    eta_n = [0.62, 0.68, 0.6, 0.5]
    Cq = [None, None, None, 5.3e6]
    eta_q = [None, None, None, 0.34]

    spin_systems = single_site_system_generator(
        isotopes=isotopes,
        isotropic_chemical_shifts=shifts,
        shielding_symmetric={"zeta": zeta, "eta": eta_n},
        quadrupolar={"Cq": Cq, "eta": eta_q},
        abundance=1,
    )

    sim = Simulator()
    sim.spin_systems = spin_systems
    pd_o = sim.sites().to_pd()

    assert list(pd_o["isotope"]) == isotopes
    assert list(pd_o["isotropic_chemical_shift"]) == [
        f"{i} ppm" if i is not None else None for i in shifts
    ]
    assert list(pd_o["shielding_symmetric.zeta"]) == [
        f"{i} ppm" if i is not None else None for i in zeta
    ]
    assert list(pd_o["shielding_symmetric.eta"]) == [
        i if i is not None else None for i in eta_n
    ]
    assert list(pd_o["quadrupolar.Cq"]) == [
        f"{i} Hz" if i is not None else None for i in Cq
    ]
    # assert list(pd_o["quadrupolar.eta"]) == [
    #     i if i is not None else None for i in eta_q
    # ]


def test_parallel_chunks():
    def check_chunks(items_list, n_jobs, block):
        chunks = get_chunks(items_list, n_jobs)
        final = get_blocks(block)
        assert chunks == final

    def get_blocks(indexes):
        lst = [np.arange(i) for i in indexes]
        sum_ = 0
        for i, item in enumerate(lst[1:]):
            sum_ += indexes[i]
            item += sum_
        return [item.tolist() for item in lst]

    items_list = np.arange(120).tolist()
    check_chunks(items_list, 3, [40, 40, 40])

    items_list = np.arange(130).tolist()
    check_chunks(items_list, 3, [44, 43, 43])

    items_list = np.arange(185).tolist()
    check_chunks(items_list, 6, [31, 31, 31, 31, 31, 30])

    items_list = np.arange(185).tolist()
    check_chunks(items_list, 8, [24, 23, 23, 23, 23, 23, 23, 23])

    items_list = np.arange(85).tolist()
    check_chunks(items_list, 8, [11, 11, 11, 11, 11, 10, 10, 10])

    div, rem = 85 // __CPU_count__, 85 % __CPU_count__
    lst = [div] * __CPU_count__
    for i in range(rem):
        lst[i] += 1
    items_list = np.arange(85).tolist()
    check_chunks(items_list, -1, lst)
