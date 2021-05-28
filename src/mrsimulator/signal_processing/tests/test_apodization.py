# -*- coding: utf-8 -*-
"""Apodization test"""
import csdmpy as cp
import numpy as np
from mrsimulator import signal_processing as sp
from mrsimulator import Simulator
from mrsimulator import SpinSystem
from mrsimulator.methods import BlochDecaySpectrum

from .test_signal_processing import setup_read_write

__author__ = "Maxwell C. Venetos"
__email__ = "maxvenetos@gmail.com"

sim = Simulator()
the_site = {"isotope": "1H", "isotropic_chemical_shift": "0 ppm"}
the_spin_system = {"name": "site A", "sites": [the_site], "abundance": "80%"}
spin_system_object = SpinSystem.parse_dict_with_units(the_spin_system)
sim.spin_systems += [spin_system_object, spin_system_object, spin_system_object]
sim.config.decompose_spectrum = "spin_system"

sim.methods += [
    BlochDecaySpectrum(
        channels=["1H"],
        magnetic_flux_density=9.4,
        spectral_dimensions=[{"count": 65536, "spectral_width": 25000}],
    )
]
sim.run()
freqHz = sim.methods[0].spectral_dimensions[0].coordinates_Hz()


def test_scale():
    PS_0 = [sp.Scale(factor=10)]
    post_sim = sp.SignalProcessor(operations=PS_0)
    data = post_sim.apply_operations(data=sim.methods[0].simulation.copy())
    _, y0, y1, y2 = sim.methods[0].simulation.to_list()
    _, y0_, y1_, y2_ = data.to_list()

    assert y0_.max() / y0.max() == 10, "Scaling failed"
    assert y1_.max() / y1.max() == 10, "Scaling failed"
    assert y2_.max() / y2.max() == 10, "Scaling failed"


def test_Lorentzian():
    PS_1 = [
        sp.IFFT(dim_index=0),
        sp.apodization.Exponential(FWHM="200 Hz", dim_index=0, dv_index=0),
        sp.FFT(dim_index=0),
    ]
    post_sim = sp.SignalProcessor(operations=PS_1)
    data = post_sim.apply_operations(data=sim.methods[0].simulation.copy())
    _, y0, y1, y2 = data.to_list()

    FWHM = 200
    test = (FWHM / 2) / (np.pi * (freqHz ** 2 + (FWHM / 2) ** 2))

    assert np.allclose(y1, y2)
    assert np.all(y0 != y1)
    assert np.allclose(
        test / test.sum(), y0 / y0.sum(), atol=1e-04
    ), "Lorentzian apodization amplitude failed"
    assert np.allclose(
        y1.sum(), y0.sum()
    ), "Area not conserved after Lorentzian apodization"


def test_Gaussian():
    FWHM = 200 * 2.354820045030949
    PS_2 = [
        sp.IFFT(dim_index=0),
        sp.apodization.Gaussian(FWHM=f"{FWHM} Hz", dim_index=0, dv_index=[0, 1]),
        sp.FFT(dim_index=0),
    ]

    PS_3 = [
        sp.IFFT(dim_index=0),
        sp.apodization.Gaussian(FWHM=f"{FWHM} Hz", dim_index=0, dv_index=None),
        sp.FFT(dim_index=0),
    ]

    post_sim = sp.SignalProcessor(operations=PS_2)
    data = post_sim.apply_operations(data=sim.methods[0].simulation.copy())
    _, y0, y1, _ = data.to_list()

    sigma = 200
    test = (1 / (sigma * np.sqrt(2 * np.pi))) * np.exp(-((freqHz / sigma) ** 2) / 2)

    assert np.allclose(y0, y1), "Gaussian apodization on two dv are not equal."

    assert np.allclose(
        test / test.sum(), y0 / y0.sum(), atol=1e-04
    ), "Gaussian apodization amplitude failed"

    # test None for dv_index
    post_sim = sp.SignalProcessor(operations=PS_3)
    data = post_sim.apply_operations(data=sim.methods[0].simulation.copy())
    _, y0, y1, y2 = data.to_list()

    assert np.allclose(y0, y1), "Gaussian apodization on dv at 0 and 1 are unequal."
    assert np.allclose(y0, y2), "Gaussian apodization on dv at 0 and 2 are unequal."
    assert np.allclose(
        test / test.sum(), y0 / y0.sum(), atol=1e-04
    ), "Gaussian apodization amplitude failed"


def test_scale_class():
    # direct initialization
    a = sp.Scale(factor=200)

    assert a.factor == 200
    assert a.property_units == {}

    # class to dict with units
    dict_ = a.json()

    assert dict_ == {
        "function": "Scale",
        "factor": 200.0,
    }

    # read from dictionary
    b = sp.Scale.parse_dict_with_units(dict_)

    assert a == b


def test_Exponential_class():
    # direct initialization
    a = sp.apodization.Exponential(FWHM="200 s", dim_index=0, dv_index=0)

    assert a.FWHM == 200
    assert a.property_units == {"FWHM": "s"}
    assert a.dim_index == 0
    assert a.dv_index == 0

    py_dict = {
        "function": "apodization",
        "type": "Exponential",
        "FWHM": "200.0 s",
        "dim_index": 0,
        "dv_index": 0,
    }
    setup_read_write(a, py_dict, sp.apodization.Exponential)


def test_Gaussian_class():
    # direct initialization
    a = sp.apodization.Gaussian(FWHM="200 km/s", dim_index=0, dv_index=0)

    assert a.FWHM == 200
    assert a.property_units == {"FWHM": "km / s"}
    assert a.dim_index == 0
    assert a.dv_index == 0

    py_dict = {
        "function": "apodization",
        "type": "Gaussian",
        "FWHM": "200.0 km / s",
        "dim_index": 0,
        "dv_index": 0,
    }
    setup_read_write(a, py_dict, sp.apodization.Gaussian)


def test_2D_area():
    data = np.zeros((256, 128), dtype=float)
    data[128, 64] = 1.0
    csdm_obj = cp.as_csdm(data)

    # test00
    PS = [
        sp.IFFT(dim_index=(0, 1)),
        sp.apodization.Gaussian(FWHM="35", dim_index=0),
        sp.apodization.Gaussian(FWHM="55", dim_index=1),
        sp.FFT(dim_index=(0, 1)),
    ]

    post_sim = sp.SignalProcessor(operations=PS)
    data_new = post_sim.apply_operations(data=csdm_obj.copy())
    _, __, y1 = data_new.to_list()

    assert np.allclose(y1.sum(), data.sum())

    # test01
    PS = [
        sp.IFFT(dim_index=(0, 1)),
        sp.apodization.Gaussian(FWHM="35", dim_index=0),
        sp.apodization.Exponential(FWHM="55", dim_index=1),
        sp.FFT(dim_index=(0, 1)),
    ]

    post_sim = sp.SignalProcessor(operations=PS)
    data_new = post_sim.apply_operations(data=csdm_obj.copy())
    _, __, y1 = data_new.to_list()

    assert np.allclose(y1.sum(), data.sum())
