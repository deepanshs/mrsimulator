# -*- coding: utf-8 -*-
import pytest
from mrsimulator import Isotopomer
from mrsimulator import Simulator
from mrsimulator import Site
from mrsimulator import SymmetricTensor

__all__ = []


@pytest.fixture(autouse=True)
def add_site(doctest_namespace):
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

    isotopomer_1 = Isotopomer(sites=[site1])
    doctest_namespace["isotopomer_1"] = isotopomer_1

    isotopomers = Isotopomer(sites=[site1, site2, site3])
    doctest_namespace["isotopomers"] = isotopomers

    sim = Simulator()
    sim.isotopomers.append(isotopomers)
    doctest_namespace["sim"] = sim
