# -*- coding: utf-8 -*-
"""Apodization test"""
import csdmpy as cp
import numpy as np
from mrsimulator import Dimension
from mrsimulator import Isotopomer
from mrsimulator import Simulator
from mrsimulator import Site
from mrsimulator.apodization import Apodization
from mrsimulator.methods import one_d_spectrum


sim = Simulator()

the_site = {"isotope": "1H", "isotropic_chemical_shift": "0 ppm"}

the_isotopomer = {"name": "site A", "sites": [the_site], "abundance": "80%"}

isotopomer_object = Isotopomer.parse_dict_with_units(the_isotopomer)

sim.isotopomers += [isotopomer_object]


dimension = {
    "isotope": "1H",
    "magnetic_flux_density": "9.4 T",
    "rotor_angle": "0 deg",
    "rotor_frequency": "0 kHz",
    "number_of_points": 4096,
    "spectral_width": "25 kHz",
    "reference_offset": "0 kHz",
}

dimension_object = Dimension.parse_dict_with_units(dimension)
sim.dimensions += [dimension_object]

freq, amp = sim.run(method=one_d_spectrum)

freqHz = sim.dimensions[0].coordinates_Hz


def test_Lorentzian():
    sigma = 200
    test = (sigma / 2) / (np.pi * (freqHz ** 2 + (sigma / 2) ** 2))
    y = sim.apodize(Apodization.Lorentzian, sigma=sigma)

    assert np.allclose(
        test / test.max(), y / y.max(), atol=1e-04
    ), "Lorentzian appodization amplitude failed"


def test_Gaussian():
    sigma = 20
    test = (1 / (sigma * np.sqrt(2 * np.pi))) * np.exp(-((freqHz / sigma) ** 2) / 2)
    y = sim.apodize(Apodization.Gaussian, sigma=sigma)

    assert np.allclose(
        test / test.max(), y / y.max(), atol=1e-04
    ), "Gaussian appodization amplitude failed"
