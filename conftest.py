# -*- coding: utf-8 -*-
from pprint import pprint

import matplotlib
import matplotlib.pyplot as plt
import numpy as np
import pytest
from mrsimulator import Simulator
from mrsimulator import Site
from mrsimulator import SpinSystem
from mrsimulator.methods import BlochDecayCentralTransitionSpectrum
from mrsimulator.spin_system.isotope import Isotope
from mrsimulator.spin_system.tensors import SymmetricTensor

font = {"weight": "light", "size": 9}
matplotlib.rc("font", **font)


@pytest.fixture(autouse=True)
def add_site(doctest_namespace):

    doctest_namespace["np"] = np
    doctest_namespace["plt"] = plt
    doctest_namespace["SpinSystem"] = SpinSystem
    doctest_namespace["Simulator"] = Simulator
    doctest_namespace["Site"] = Site
    doctest_namespace["SymmetricTensor"] = SymmetricTensor
    doctest_namespace["st"] = SymmetricTensor
    doctest_namespace["pprint"] = pprint
    doctest_namespace["Isotope"] = Isotope

    site1 = Site(
        isotope="13C",
        isotropic_chemical_shift=20,
        shielding_symmetric=SymmetricTensor(zeta=10, eta=0.5),
    )
    doctest_namespace["site1"] = site1

    site2 = Site(
        isotope="1H",
        isotropic_chemical_shift=-4,
        shielding_symmetric=SymmetricTensor(zeta=2.1, eta=0.1),
    )
    doctest_namespace["site2"] = site2

    site3 = Site(
        isotope="27Al",
        isotropic_chemical_shift=120,
        shielding_symmetric=SymmetricTensor(zeta=2.1, eta=0.1),
        quadrupole=SymmetricTensor(Cq=5.1e6, eta=0.5),
    )
    doctest_namespace["site3"] = site3

    spin_system_1H_13C = SpinSystem(sites=[site1, site2])
    doctest_namespace["spin_system_1H_13C"] = spin_system_1H_13C

    spin_system_1 = SpinSystem(sites=[site1])
    doctest_namespace["spin_system_1"] = spin_system_1

    doctest_namespace["spin_systems"] = SpinSystem(sites=[site1, site2, site3])

    spin_systems = [SpinSystem(sites=[site]) for site in [site1, site2, site3]]

    sim = Simulator()
    sim.spin_systems += spin_systems
    doctest_namespace["sim"] = sim

    # coesite
    O17_1 = Site(
        isotope="17O",
        isotropic_chemical_shift=29,
        quadrupolar=SymmetricTensor(Cq=6.05e6, eta=0.000),
    )
    O17_2 = Site(
        isotope="17O",
        isotropic_chemical_shift=41,
        quadrupolar=SymmetricTensor(Cq=5.43e6, eta=0.166),
    )
    O17_3 = Site(
        isotope="17O",
        isotropic_chemical_shift=57,
        quadrupolar=SymmetricTensor(Cq=5.45e6, eta=0.168),
    )
    O17_4 = Site(
        isotope="17O",
        isotropic_chemical_shift=53,
        quadrupolar=SymmetricTensor(Cq=5.52e6, eta=0.169),
    )
    O17_5 = Site(
        isotope="17O",
        isotropic_chemical_shift=58,
        quadrupolar=SymmetricTensor(Cq=5.16e6, eta=0.292),
    )

    sites = [O17_1, O17_2, O17_3, O17_4, O17_5]
    abundance = [0.83, 1.05, 2.16, 2.05, 1.90]  # abundance of each spin system
    spin_systems = [
        SpinSystem(sites=[s], abundance=a) for s, a in zip(sites, abundance)
    ]

    method = BlochDecayCentralTransitionSpectrum(
        channels=["17O"],
        rotor_frequency=14000,
        spectral_dimensions=[{"count": 2048, "spectral_width": 50000}],
    )

    sim_coesite = Simulator()
    sim_coesite.spin_systems += spin_systems
    sim_coesite.methods += [method]

    doctest_namespace["sim_coesite"] = sim_coesite
