# -*- coding: utf-8 -*-
"""Test for shift and reference offset."""
import numpy as np
from mrsimulator import Simulator
from mrsimulator import Site
from mrsimulator import SpinSystem
from mrsimulator.method.lib import BlochDecaySpectrum

# import matplotlib.pyplot as plt


def setup_test(spin_system, volume="octant", sw=25000, n_gamma=500):
    mth_kwargs = {
        "channels": [spin_system.sites[0].isotope.symbol],
        "spectral_dimensions": [{"count": 1024, "spectral_width": sw}],
    }

    data = []
    for angle in [0, 54.735, 90]:
        method = BlochDecaySpectrum(
            rotor_angle=angle * np.pi / 180, rotor_frequency=0, **mth_kwargs
        )
        sim = Simulator(spin_systems=[spin_system], methods=[method])
        sim.config.integration_volume = volume
        sim.config.number_of_gamma_angles = 1 if angle == 0 else n_gamma
        sim.run(auto_switch=False)

        res = sim.methods[0].simulation.y[0].components[0]
        res /= res.max()
        data.append(res)

    # plt.plot(data[0])
    # plt.plot(data[1], "--")
    # plt.plot(data[2], "-.")
    # plt.show()

    np.testing.assert_almost_equal(data[0], data[1], decimal=2)
    np.testing.assert_almost_equal(data[0], data[2], decimal=1.6)


def test_csa_01():
    site = Site(isotope="13C", shielding_symmetric={"zeta": 50, "eta": 0.5})
    spin_system = SpinSystem(sites=[site])
    setup_test(spin_system, volume="octant", sw=2.5e4)


def test_quad_01():
    site = Site(
        isotope="27Al",
        shielding_symmetric={"zeta": 50, "eta": 0.5},
        quadrupolar={"Cq": 1e6, "eta": 0.1, "alpha": 1, "beta": 2, "gamma": 3},
    )
    spin_system = SpinSystem(sites=[site])
    setup_test(spin_system, volume="hemisphere", sw=1e6)
