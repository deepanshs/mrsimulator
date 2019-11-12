# -*- coding: utf-8 -*-
import pytest
from mrsimulator import Site
from mrsimulator import SymmetricTensor

__all__ = []


@pytest.fixture(autouse=True)
def add_site(doctest_namespace):
    doctest_namespace["site1"] = Site(
        isotope="13C",
        isotropic_chemical_shift=20,
        shielding_symmetric=SymmetricTensor(zeta=10, eta=0.5),
        quadrupole=SymmetricTensor(Cq=5.1e6, eta=0.5),
    )
