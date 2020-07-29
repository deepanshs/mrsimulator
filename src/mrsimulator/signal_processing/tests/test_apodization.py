# -*- coding: utf-8 -*-
"""Apodization test"""
import mrsimulator.signal_processing as sp
import mrsimulator.signal_processing.apodization as apo
import numpy as np
from mrsimulator import Simulator
from mrsimulator import SpinSystem
from mrsimulator.methods import BlochDecaySpectrum

sim = Simulator()
the_site = {"isotope": "1H", "isotropic_chemical_shift": "0 ppm"}
the_spin_system = {"name": "site A", "sites": [the_site], "abundance": "80%"}
spin_system_object = SpinSystem.parse_dict_with_units(the_spin_system)
sim.spin_systems += [spin_system_object, spin_system_object, spin_system_object]
sim.config.decompose_spectrum = "spin_system"

method_1 = BlochDecaySpectrum(
    channels=["1H"],
    magnetic_flux_density=9.4,
    rotor_angle=0,
    rotor_frequency=0,
    spectral_dimensions=[
        {"count": 4096, "spectral_width": 25000, "reference_offset": 0}
    ],
)

PS_0 = [sp.Scale(factor=10)]

PS_1 = [
    sp.IFFT(dim_indx=0),
    apo.Exponential(FWHM="200 Hz", dim_indx=0, dv_indx=0),
    sp.FFT(dim_indx=0),
]

sigma = 20 * 2.354820045030949
PS_2 = [
    sp.IFFT(dim_indx=0),
    apo.Gaussian(FWHM=f"{sigma} Hz", dim_indx=0, dv_indx=[0, 1]),
    sp.FFT(dim_indx=0),
]

PS_3 = [
    sp.IFFT(dim_indx=0),
    apo.Gaussian(FWHM=f"{sigma} Hz", dim_indx=0, dv_indx=None),
    sp.FFT(dim_indx=0),
]

sim.methods += [method_1]
sim.run()


freqHz = sim.methods[0].spectral_dimensions[0].coordinates_Hz()


def test_scale():
    post_sim = sp.SignalProcessor(operations=PS_0)
    data = post_sim.apply_operations(data=sim.methods[0].simulation)
    _, y0, y1, y2 = sim.methods[0].simulation.to_list()
    _, y0_, y1_, y2_ = data.to_list()

    assert y0_.max() / y0.max() == 10, "Scaling failed"
    assert y1_.max() / y1.max() == 10, "Scaling failed"
    assert y2_.max() / y2.max() == 10, "Scaling failed"


def test_Lorentzian():
    post_sim = sp.SignalProcessor(operations=PS_1)
    data = post_sim.apply_operations(data=sim.methods[0].simulation)
    _, y0, y1, y2 = data.to_list()

    FWHM = 200
    test = (FWHM / 2) / (np.pi * (freqHz ** 2 + (FWHM / 2) ** 2))

    assert np.allclose(y1, y2)
    assert np.all(y0 != y1)
    assert np.allclose(
        test / test.max(), y0 / y0.max(), atol=1e-04
    ), "Lorentzian apodization amplitude failed"


def test_Gaussian():
    post_sim = sp.SignalProcessor(operations=PS_2)
    data = post_sim.apply_operations(data=sim.methods[0].simulation)
    _, y0, y1, _ = data.to_list()

    sigma = 20
    test = (1 / (sigma * np.sqrt(2 * np.pi))) * np.exp(-((freqHz / sigma) ** 2) / 2)

    assert np.allclose(y0, y1)
    assert np.allclose(
        test / test.max(), y0 / y0.max(), atol=1e-04
    ), "Gaussian apodization amplitude failed"

    # test None for dv_indx
    post_sim = sp.SignalProcessor(operations=PS_3)
    data = post_sim.apply_operations(data=sim.methods[0].simulation)
    _, y0, y1, y2 = data.to_list()

    sigma = 20
    test = (1 / (sigma * np.sqrt(2 * np.pi))) * np.exp(-((freqHz / sigma) ** 2) / 2)

    assert np.allclose(y0, y1)
    assert np.allclose(y0, y2)
    assert np.allclose(
        test / test.max(), y0 / y0.max(), atol=1e-04
    ), "Gaussian apodization amplitude failed"


def test_scale_class():
    # direct initialization
    a = sp.Scale(factor=200)

    assert a.factor == 200
    assert a.property_units == {}

    # class to dict with units
    dict_ = a.to_dict_with_units()

    assert dict_ == {
        "function": "Scale",
        "factor": 200.0,
    }

    # read from dictionary
    b = sp.Scale.parse_dict_with_units(dict_)

    assert a == b


def test_Exponential_class():
    # direct initialization
    a = apo.Exponential(FWHM="200 s", dim_indx=0, dv_indx=0)

    assert a.FWHM == 200
    assert a.property_units == {"FWHM": "s"}
    assert a.dim_indx == 0
    assert a.dv_indx == 0

    # class to dict with units
    dict_ = a.to_dict_with_units()

    assert dict_ == {
        "function": "apodization",
        "type": "Exponential",
        "FWHM": "200.0 s",
        "dim_indx": 0,
        "dv_indx": 0,
    }

    # read from dictionary
    b = apo.Exponential.parse_dict_with_units(dict_)

    assert a == b


def test_Gaussian_class():
    # direct initialization
    a = apo.Gaussian(FWHM="200 km/s", dim_indx=0, dv_indx=0)

    assert a.FWHM == 200
    assert a.property_units == {"FWHM": "km / s"}
    assert a.dim_indx == 0
    assert a.dv_indx == 0

    # class to dict with units
    dict_ = a.to_dict_with_units()

    assert dict_ == {
        "function": "apodization",
        "type": "Gaussian",
        "FWHM": "200.0 km / s",
        "dim_indx": 0,
        "dv_indx": 0,
    }

    # read from dictionary
    b = apo.Gaussian.parse_dict_with_units(dict_)

    assert a == b
