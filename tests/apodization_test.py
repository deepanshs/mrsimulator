# -*- coding: utf-8 -*-
"""Apodization test"""
import numpy as np
from mrsimulator import Simulator
from mrsimulator import SpinSystem
from mrsimulator.methods import BlochDecaySpectrum
from mrsimulator.post_simulation import PostSimulator


sim = Simulator()
the_site = {"isotope": "1H", "isotropic_chemical_shift": "0 ppm"}
the_isotopomer = {"name": "site A", "sites": [the_site], "abundance": "80%"}
isotopomer_object = SpinSystem.parse_dict_with_units(the_isotopomer)
sim.spin_systems += [isotopomer_object]

method_1 = BlochDecaySpectrum(
    channels=["1H"],
    magnetic_flux_density=9.4,
    rotor_angle=0,
    rotor_frequency=0,
    spectral_dimensions=[
        {"count": 4096, "spectral_width": 25000, "reference_offset": 0}
    ],
)


PS_1 = PostSimulator(
    scale=1, apodization=[{"args": [200], "function": "Lorentzian", "dimension": 0}]
)

PS_2 = PostSimulator(
    scale=1, apodization=[{"args": [20], "function": "Gaussian", "dimension": 0}]
)

sim.methods += [method_1, method_1]
sim.run()


freqHz = sim.methods[0].spectral_dimensions[0].coordinates_Hz


def test_Lorentzian():
    sim.methods[0].post_simulation = PS_1

    sigma = 200
    test = (sigma / 2) / (np.pi * (freqHz ** 2 + (sigma / 2) ** 2))

    x, y = sim.methods[0].simulation.to_list()
    y = sim.methods[0].apodize().real

    assert np.allclose(
        test / test.max(), y / y.max(), atol=1e-04
    ), "Lorentzian apodization amplitude failed"


def test_Gaussian():
    sim.methods[0].post_simulation = PS_2
    sigma = 20
    test = (1 / (sigma * np.sqrt(2 * np.pi))) * np.exp(-((freqHz / sigma) ** 2) / 2)
    x, y = sim.methods[1].simulation.to_list()
    y = sim.methods[0].apodize().real

    assert np.allclose(
        test / test.max(), y / y.max(), atol=1e-04
    ), "Gaussian apodization amplitude failed"
