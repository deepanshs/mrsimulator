# -*- coding: utf-8 -*-
"""Apodization test"""
import numpy as np
from mrsimulator import Isotopomer
from mrsimulator import Method
from mrsimulator import Simulator
from mrsimulator.apodization import Apodization


sim = Simulator()

the_site = {"isotope": "1H", "isotropic_chemical_shift": "0 ppm"}

the_isotopomer = {"name": "site A", "sites": [the_site], "abundance": "80%"}

isotopomer_object = Isotopomer.parse_dict_with_units(the_isotopomer)

sim.isotopomers += [isotopomer_object]


# dimension = {
#     "isotope": "1H",
#     "magnetic_flux_density": "9.4 T",
#     "rotor_angle": "0 deg",
#     "rotor_frequency": "0 kHz",
#     "number_of_points": 4096,
#     "spectral_width": "25 kHz",
#     "reference_offset": "0 kHz",
# }

method = {
    "isotope": "1H",
    "spectral_dimensions": [
        {
            "count": 4096,
            "spectral_width": "25 kHz",
            "reference_offset": "0 kHz",
            "events": [
                {
                    "magnetic_flux_density": "9.4 T",
                    "rotor_angle": "0 deg",
                    "rotor_frequency": "0 kHz",
                }
            ],
        }
    ],
}


sim.methods += [Method.parse_dict_with_units(method)]

sim.run()

freqHz = sim.methods[0].spectral_dimensions[0].coordinates_Hz


def test_Lorentzian():
    sigma = 200
    test = (sigma / 2) / (np.pi * (freqHz ** 2 + (sigma / 2) ** 2))
    y = sim.apodize(Apodization.Lorentzian, sigma=sigma)

    assert np.allclose(
        test / test.max(), y / y.max(), atol=1e-04
    ), "Lorentzian apodization amplitude failed"


def test_Gaussian():
    sigma = 20
    test = (1 / (sigma * np.sqrt(2 * np.pi))) * np.exp(-((freqHz / sigma) ** 2) / 2)
    y = sim.apodize(Apodization.Gaussian, sigma=sigma)

    assert np.allclose(
        test / test.max(), y / y.max(), atol=1e-04
    ), "Gaussian apodization amplitude failed"
