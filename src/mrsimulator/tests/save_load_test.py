# -*- coding: utf-8 -*-
import os

from mrsimulator import load
from mrsimulator import save
from mrsimulator import signal_processing as sp
from mrsimulator import Simulator
from mrsimulator import Site
from mrsimulator import SpinSystem
from mrsimulator.methods import BlochDecaySpectrum

# mrsimulator save and load test


def setup():
    site = Site(isotope="1H", shielding_symmetric={"zeta": -100, "eta": 0.3})
    spin_sys = SpinSystem(sites=[site, site], abundance=45)

    method = BlochDecaySpectrum(channels=["1H"])
    sim = Simulator(spin_systems=[spin_sys, spin_sys], methods=[method, method])

    sim.run(method_index=0)

    processor = sp.SignalProcessor(
        operations=[
            sp.IFFT(),
            sp.apodization.Exponential(FWHM="2500 Hz"),
            sp.FFT(),
            sp.Scale(factor=20),
        ]
    )
    processors = [processor] * 2

    return sim, processors


def test_save():
    sim, processors = setup()

    assert sim.methods[0].simulation is not None
    assert sim.methods[1].simulation is None

    save("test.mrsim", simulator=sim, signal_processors=processors, params=None)


def test_load():
    sim_r, processors_r, report_r = load("test.mrsim")

    sim, processors = setup()
    sim_r.methods[0].simulation = None
    sim.methods[0].simulation = None
    assert sim_r == sim
    assert processors_r == processors
    assert report_r is None

    os.remove("test.mrsim")
