# -*- coding: utf-8 -*-
from os import path
from pprint import pprint

import matplotlib
import matplotlib.pyplot as plt
import pytest
from mrsimulator import Isotopomer
from mrsimulator import Simulator
from mrsimulator import Site
from mrsimulator import SymmetricTensor

font = {"family": "Helvetica", "weight": "light", "size": 9}
matplotlib.rc("font", **font)


@pytest.fixture(autouse=True)
def add_site(doctest_namespace):

    doctest_namespace["Isotopomer"] = Isotopomer
    doctest_namespace["Simulator"] = Simulator
    doctest_namespace["Site"] = Site
    doctest_namespace["SymmetricTensor"] = SymmetricTensor
    doctest_namespace["st"] = SymmetricTensor
    doctest_namespace["pprint"] = pprint

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

    doctest_namespace["isotopomers"] = Isotopomer(sites=[site1, site2, site3])

    isotopomers = [Isotopomer(sites=[site]) for site in [site1, site2, site3]]

    sim = Simulator()
    sim.isotopomers += isotopomers
    doctest_namespace["sim"] = sim

    # dim = Method(isotope="27Al", spectral_dimensions=[
    #     dict(spectral_width=50000), rotor_frequency=12000)
    # doctest_namespace["dim"] = dim

    # dimension_1 = {
    #     "number_of_points": 1024,
    #     "spectral_width": "100 Hz",
    #     "reference_offset": "0 Hz",
    #     "magnetic_flux_density": "9.4 T",
    #     "rotor_frequency": "0 Hz",
    #     "rotor_angle": "54.935 degree",
    #     "isotope": "29Si",
    # }
    # dimension_object = Dimension.parse_dict_with_units(dimension_1)
    # doctest_namespace["dimension_object"] = dimension_object

    def plot_save(x, y, filename):
        plt.figure(figsize=(4.5, 2.5))
        plt.plot(x, y, linewidth=1)
        plt.xlim([x.value.max(), x.value.min()])
        plt.xlabel(f"frequency ratio / {str(x.unit)}", **font)
        plt.grid(color="gray", linestyle="--", linewidth=0.75, alpha=0.25)
        plt.tight_layout(pad=0.15)

        filename = path.split(filename)[1]
        filepath = "./docs/_images"
        pth = path.join(filepath, filename)
        plt.savefig(pth + ".pdf")
        plt.savefig(pth + ".png", dpi=100)
        plt.close()

    doctest_namespace["plot_save"] = plot_save

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
    abundance = [0.83, 1.05, 2.16, 2.05, 1.90]  # abundance of each isotopomer
    isotopomers = [Isotopomer(sites=[s], abundance=a) for s, a in zip(sites, abundance)]

    # dimension = Dimension(
    #     isotope="17O",
    #     number_of_points=2046,
    #     spectral_width=50000,
    #     rotor_frequency=14000,
    # )

    # sim_coesite = Simulator()
    # sim_coesite.isotopomers += isotopomers
    # sim_coesite.dimensions += [dimension]

    # doctest_namespace["sim_coesite"] = sim_coesite
