# -*- coding: utf-8 -*-
from pprint import pprint

import pytest
from mrsimulator import Dimension
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

    dim = Dimension(isotope="27Al", spectral_width=50000, rotor_frequency=12000)
    doctest_namespace["dim"] = dim

    dimension_1 = {
        "number_of_points": 1024,
        "spectral_width": "100 Hz",
        "reference_offset": "0 Hz",
        "magnetic_flux_density": "9.4 T",
        "rotor_frequency": "0 Hz",
        "rotor_angle": "54.935 degree",
        "isotope": "29Si",
    }
    dimension_object = Dimension.parse_dict_with_units(dimension_1)
    doctest_namespace["dimension_object"] = dimension_object
    doctest_namespace["pprint"] = pprint

    doctest_namespace["Simulator"] = Simulator
