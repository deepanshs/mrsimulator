# -*- coding: utf-8 -*-
"""Apodization test"""
import mrsimulator.post_simulation as ps
import numpy as np
from mrsimulator import Simulator
from mrsimulator import SpinSystem
from mrsimulator.methods import BlochDecaySpectrum
from mrsimulator.post_simulation import SignalProcessor

sim = Simulator()
the_site = {"isotope": "1H", "isotropic_chemical_shift": "0 ppm"}
the_spin_system = {"name": "site A", "sites": [the_site], "abundance": "80%"}
spin_system_object = SpinSystem.parse_dict_with_units(the_spin_system)
sim.spin_systems += [spin_system_object]

method_1 = BlochDecaySpectrum(
    channels=["1H"],
    magnetic_flux_density=9.4,
    rotor_angle=0,
    rotor_frequency=0,
    spectral_dimensions=[
        {"count": 4096, "spectral_width": 25000, "reference_offset": 0}
    ],
)

PS_0 = {"dependent_variable": 0, "operations": [ps.Scale(factor=10)]}

PS_1 = {
    "dependent_variable": 0,
    "operations": [
        ps.IFFT(dimension=0),
        ps.Exponential(Lambda=200, dimension=0),
        ps.FFT(dimension=0),
    ],
}

PS_2 = {
    "dependent_variable": 0,
    "operations": [
        ps.IFFT(dimension=0),
        ps.Gaussian(sigma=20, dimension=0),
        ps.FFT(dimension=0),
    ],
}

sim.methods += [method_1]
sim.run()


freqHz = sim.methods[0].spectral_dimensions[0].coordinates_Hz()


def test_scale():
    post_sim = SignalProcessor(data=sim.methods[0].simulation, operations=[PS_0])
    post_sim.apply_operations()
    x0, y0 = sim.methods[0].simulation.to_list()
    x, y = post_sim.data.to_list()

    assert np.allclose(sum(y0 / y), 10), "Scaling failed"


def test_Lorentzian():
    sim.methods[0].post_simulation = PS_1
    post_sim = SignalProcessor(data=sim.methods[0].simulation, operations=[PS_1])
    post_sim.apply_operations()
    x, y = post_sim.data.to_list()

    sigma = 200
    test = (sigma / 2) / (np.pi * (freqHz ** 2 + (sigma / 2) ** 2))

    assert np.allclose(
        test / test.max(), y / y.max(), atol=1e-04
    ), "Lorentzian apodization amplitude failed"


def test_Gaussian():
    post_sim = SignalProcessor(data=sim.methods[0].simulation, operations=[PS_2])
    post_sim.apply_operations()
    x, y = post_sim.data.to_list()

    sigma = 20
    test = (1 / (sigma * np.sqrt(2 * np.pi))) * np.exp(-((freqHz / sigma) ** 2) / 2)

    assert np.allclose(
        test / test.max(), y / y.max(), atol=1e-04
    ), "Gaussian apodization amplitude failed"
