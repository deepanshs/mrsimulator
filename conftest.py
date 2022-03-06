# -*- coding: utf-8 -*-
from pprint import pprint

import matplotlib
import matplotlib.pyplot as plt
import numpy as np
import pytest
from mrsimulator import Coupling
from mrsimulator import signal_processing as sp
from mrsimulator import Simulator
from mrsimulator import Site
from mrsimulator import SpinSystem
from mrsimulator.methods import Method2D
from mrsimulator.models import CzjzekDistribution
from mrsimulator.models import ExtCzjzekDistribution
from mrsimulator.spin_system.isotope import Isotope
from mrsimulator.spin_system.tensors import SymmetricTensor
from mrsimulator.transition import Transition
from mrsimulator.transition import TransitionPathway
from sybil import Sybil
from sybil.parsers.codeblock import PythonCodeBlockParser
from sybil.parsers.skip import skip

from plot_directive_parser import PythonPlotParser

font = {"weight": "light", "size": 9}
matplotlib.rc("font", **font)


@pytest.fixture(autouse=True)
def test_models(doctest_namespace):
    cz_model = CzjzekDistribution(0.5)
    doctest_namespace["cz_model"] = cz_model

    S0 = {"Cq": 1e6, "eta": 0.3}
    ext_cz_model = ExtCzjzekDistribution(S0, eps=0.35)
    doctest_namespace["ext_cz_model"] = ext_cz_model


@pytest.fixture(autouse=True)
def add_site(doctest_namespace):

    doctest_namespace["np"] = np
    doctest_namespace["plt"] = plt
    doctest_namespace["SpinSystem"] = SpinSystem
    doctest_namespace["Simulator"] = Simulator
    doctest_namespace["Site"] = Site
    doctest_namespace["Coupling"] = Coupling
    doctest_namespace["SymmetricTensor"] = SymmetricTensor
    doctest_namespace["st"] = SymmetricTensor
    doctest_namespace["pprint"] = pprint
    doctest_namespace["Isotope"] = Isotope
    doctest_namespace["sp"] = sp
    doctest_namespace["Method2D"] = Method2D

    site1 = Site(
        isotope="13C",
        isotropic_chemical_shift=20,
        shielding_symmetric=SymmetricTensor(zeta=10, eta=0.5),
    )
    site2 = Site(
        isotope="1H",
        isotropic_chemical_shift=-4,
        shielding_symmetric=SymmetricTensor(zeta=2.1, eta=0.1),
    )
    site3 = Site(
        isotope="27Al",
        isotropic_chemical_shift=120,
        shielding_symmetric=SymmetricTensor(zeta=2.1, eta=0.1),
        quadrupolar=SymmetricTensor(Cq=5.1e6, eta=0.5),
    )

    doctest_namespace["spin_system_1H_13C"] = SpinSystem(sites=[site1, site2])
    doctest_namespace["spin_systems"] = SpinSystem(sites=[site1, site2, site3])

    spin_systems = [SpinSystem(sites=[site]) for site in [site1, site2, site3]]
    sim = Simulator()
    sim.spin_systems += spin_systems
    doctest_namespace["sim"] = sim

    # Transitions
    t1 = Transition(initial=[0.5, 0.5], final=[0.5, -0.5])
    doctest_namespace["t1"] = t1

    t2 = Transition(initial=[0.5, 0.5], final=[-0.5, 0.5])
    doctest_namespace["t2"] = t2

    path = TransitionPathway([t1, t2])
    doctest_namespace["path"] = path


pytest_collect_file = Sybil(
    parsers=[
        PythonCodeBlockParser(),
        PythonPlotParser(),
        skip,
    ],
    patterns=["*.rst"],
    fixtures=["add_site", "test_models"],
).pytest()
