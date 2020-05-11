# -*- coding: utf-8 -*-
import numpy as np
from mrsimulator import Isotopomer
from mrsimulator import Simulator
from mrsimulator import Site
from mrsimulator.methods import BlochDecaySpectrum


def test_two_site_no_coupling_test():
    site1 = Site(
        isotope="29Si",
        isotropic_chemical_shift=10,
        shielding_symmetric={"zeta": 5, "eta": 0},
    )
    site2 = Site(
        isotope="29Si",
        isotropic_chemical_shift=-10,
        shielding_symmetric={"zeta": -5, "eta": 0},
    )

    iso_two_site = [Isotopomer(sites=[site1, site2])]
    iso_single_sites = [Isotopomer(sites=[site1]), Isotopomer(sites=[site2])]

    sim1 = Simulator()
    sim1.isotopomers += iso_two_site
    sim1.methods += [
        BlochDecaySpectrum(
            channels=["29Si"], dimensions=[dict(count=2048, spectral_width=25000)]
        )
    ]
    sim1.run()

    sim2 = Simulator()
    sim2.isotopomers += iso_single_sites
    sim2.methods += [
        BlochDecaySpectrum(
            channels=["29Si"], dimensions=[dict(count=2048, spectral_width=25000)]
        )
    ]
    sim2.run()

    data1 = (sim1.methods[0].simulation / 2).dependent_variables[0].components
    data2 = sim2.methods[0].simulation.dependent_variables[0].components
    assert np.allclose(data1, data2)
