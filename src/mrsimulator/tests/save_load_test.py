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

    kwargs = dict(simulator=sim, signal_processors=processors, params=None)
    save("test.mrsim", **kwargs), "test file"
    save("test.mrsim.gz", **kwargs), "test gz"


def test_load():
    for file_, tag in zip(["test.mrsim", "test.mrsim.gz"], ["file", "gz"]):
        sim_r, processors_r, report_r = load(file_)

        sim, processors = setup()
        sim_r.methods[0].simulation = None
        sim.methods[0].simulation = None
        assert sim_r == sim, f"test {tag}"
        assert processors_r == processors, f"test {tag}"
        assert report_r is None, f"test {tag}"

        os.remove(file_)
