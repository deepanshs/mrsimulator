"""Test for the base Simulator class."""
import os
from random import randint

import csdmpy as cp
import numpy as np
import pytest
from mrsimulator import Coupling
from mrsimulator import Simulator
from mrsimulator import Site
from mrsimulator import SpinSystem
from mrsimulator.method.lib import BlochDecayCTSpectrum
from mrsimulator.method.lib import BlochDecaySpectrum
from mrsimulator.simulator import __CPU_count__
from mrsimulator.simulator import get_chunks
from mrsimulator.simulator import Sites
from mrsimulator.spin_system.tests.test_spin_systems import generate_isotopes
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
        "spin_systems": [{}],
        "config": {
            "decompose_spectrum": "none",
            "integration_density": 70,
            "integration_volume": "octant",
            "number_of_sidebands": 64,
        },
    }
    assert c.json() == result

    assert c.json(units=False) == {
        "label": "test",
        "spin_systems": [{}],
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
    sim.methods = [BlochDecaySpectrum(channels=["1H"])]
    sim.name = "test"
    sim.label = "test0"
    sim.description = "testing-testing 1.2.3"

    red_dict = sim.json(units=False)
    _ = [item.pop("description") for item in red_dict["methods"]]

    assert red_dict == {
        "name": "test",
        "label": "test0",
        "description": "testing-testing 1.2.3",
        "spin_systems": [
            {
                "sites": [
                    {
                        "isotope": {
                            "spin": 1,
                            "natural_abundance": 99.985,
                            "gyromagnetic_ratio": 42.57747920984721,
                            "quadrupole_moment": 0.0,
                            "atomic_number": 1,
                            "isotope": "1H",
                        },
                        "isotropic_chemical_shift": 0.0,
                    },
                    {
                        "isotope": {
                            "spin": 3,
                            "natural_abundance": 100.0,
                            "gyromagnetic_ratio": 11.268835291595886,
                            "quadrupole_moment": 0.1006,
                            "atomic_number": 11,
                            "isotope": "23Na",
                        },
                        "isotropic_chemical_shift": 0.0,
                    },
                ],
            }
        ],
        "methods": [
            {
                "channels": ["1H"],  # TODO: Adjust line based on serialization decision
                "name": "BlochDecaySpectrum",
                "magnetic_flux_density": 9.4,
                "rotor_angle": 0.9553166181245,
                "rotor_frequency": 0.0,
                "spectral_dimensions": [
                    {
                        "count": 1024,
                        "events": [{"transition_queries": [{"ch1": {"P": [-1]}}]}],
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
    sim_load = Simulator.load("test_sim_save.temp")

    assert sim_load.spin_systems == sim.spin_systems
    assert sim_load.methods == sim.methods
    assert sim_load.name == sim.name
    assert sim_load.description == sim.description
    assert sim_load == sim

    # without units
    sim.save("test_sim_save_no_unit.temp", with_units=False)
    sim_load = Simulator.load("test_sim_save_no_unit.temp", parse_units=False)
    assert sim_load == sim

    os.remove("test_sim_save.temp")
    os.remove("test_sim_save_no_unit.temp")


def test_sim_coesite():
    # coesite
    O17_1 = Site(
        isotope="17O",
        isotropic_chemical_shift=29,
        quadrupolar=dict(Cq=6.05e6, eta=0.000),
    )
    O17_2 = Site(
        isotope="17O",
        isotropic_chemical_shift=41,
        quadrupolar=dict(Cq=5.43e6, eta=0.166),
    )
    O17_3 = Site(
        isotope="17O",
        isotropic_chemical_shift=57,
        quadrupolar=dict(Cq=5.45e6, eta=0.168),
    )
    O17_4 = Site(
        isotope="17O",
        isotropic_chemical_shift=53,
        quadrupolar=dict(Cq=5.52e6, eta=0.169),
    )
    O17_5 = Site(
        isotope="17O",
        isotropic_chemical_shift=58,
        quadrupolar=dict(Cq=5.16e6, eta=0.292),
    )

    sites = [O17_1, O17_2, O17_3, O17_4, O17_5]
    abundance = [0.83, 1.05, 2.16, 2.05, 1.90]  # abundance of each spin system
    spin_systems = [
        SpinSystem(sites=[s], abundance=a) for s, a in zip(sites, abundance)
    ]

    method = BlochDecayCTSpectrum(
        channels=["17O"],
        rotor_frequency=14000,
        spectral_dimensions=[{"count": 2048, "spectral_width": 50000}],
    )

    sim_coesite = Simulator()
    sim_coesite.spin_systems += spin_systems
    sim_coesite.methods += [method]

    sim_coesite.save("sample.mrsim")
    sim_load = Simulator.load("sample.mrsim")
    assert sim_coesite == sim_load

    os.remove("sample.mrsim")


def test_empty_spin_sys_simulator():
    sim = Simulator()
    sim.methods = [
        BlochDecaySpectrum(channels=["1H"], spectral_dimensions=[{"count": 10}])
    ]
    sim.config.decompose_spectrum = "spin_system"
    sim.run()
    assert np.allclose(sim.methods[0].simulation.y[0].components[0], 0)


def test_simulator_2():
    sim = Simulator()
    sim.spin_systems = [
        SpinSystem(
            label="Test",
            sites=[Site(isotope="1H"), Site(isotope="23Na")],
            couplings=[Coupling(site_index=[0, 1], isotropic_j=15)],
        )
    ]
    sim.methods = [
        BlochDecaySpectrum(
            channels=["1H"],
            spectral_dimensions=[{"count": 10}],
            experiment=cp.as_csdm(np.arange(10)),
        )
    ]
    sim.name = "test"
    sim.label = "test0"
    sim.description = "testing-testing 1.2.3"

    sim.config.decompose_spectrum = "spin_system"
    sim.run()

    # save
    sim.save("test_sim_save.temp")
    sim_load = Simulator.load("test_sim_save.temp")

    sim_load_dataset = sim_load.methods[0].simulation
    sim_dataset = sim.methods[0].simulation
    sim_load_dataset._timestamp = ""
    assert sim_load_dataset.dict() == sim_dataset.dict()

    sim_load.methods[0].simulation = None
    sim.methods[0].simulation = None
    assert sim_load.spin_systems == sim.spin_systems
    assert sim_load.methods == sim.methods
    assert sim_load.name == sim.name
    assert sim_load.description == sim.description

    os.remove("test_sim_save.temp")


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
    for i, site in enumerate(sites):
        assert r_sites[i] == site

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
        isotope=isotopes,
        isotropic_chemical_shift=shifts,
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
