import os

import pytest
from mrsimulator import __version__
from mrsimulator import Mrsimulator
from mrsimulator import signal_processor as sp
from mrsimulator import Simulator
from mrsimulator import Site
from mrsimulator import SpinSystem
from mrsimulator.method.lib import BlochDecaySpectrum
from pydantic import ValidationError


__author__ = "Matthew D. Giammar"
__email__ = "giammar.7@osu.edu"


def setup_sim_and_processor():
    site = Site(isotope="1H", shielding_symmetric={"zeta": -100, "eta": 0.3})
    spin_sys = SpinSystem(sites=[site, site], abundance=45)

    method = BlochDecaySpectrum(channels=["1H"])
    sim = Simulator(spin_systems=[spin_sys, spin_sys], methods=[method, method])

    processor = sp.SignalProcessor(
        operations=[
            sp.IFFT(),
            sp.apodization.Exponential(FWHM="2500 Hz"),
            sp.FFT(),
            sp.Scale(factor=20),
        ]
    )
    processors = [processor] * 2

    application = {
        "com.github.DeepanshS.mrsimulator": {"foo": "This is some metadata"},
        "com.github.DeepanshS.mrsimulator-app": {"params": "The JSON string of params"},
    }

    return sim, processors, application


def setup_mrsimulator_obj():
    sim, processors, application = setup_sim_and_processor()

    return Mrsimulator(
        simulator=sim,
        signal_processors=processors,
        application=application,
    )


def test_json():
    mrsim = setup_mrsimulator_obj()
    mrsim.simulator.methods[0].simulation = None
    mrsim.simulator.methods[1].simulation = None

    py_dict = {
        "simulator": mrsim.simulator.json(),
        "signal_processors": [sp.json() for sp in mrsim.signal_processors],
        "application": mrsim.application,
        "version": __version__,
    }

    assert mrsim.json() == py_dict


def test_parse():
    mrsim = setup_mrsimulator_obj()
    mrsim.simulator.methods[0].simulation = None
    mrsim.simulator.methods[1].simulation = None
    py_dict = mrsim.json()
    mrsim.simulator.name = "a test Simulator object"
    py_dict["simulator"]["name"] = "a test Simulator object"
    parsed_mrsim = Mrsimulator.parse(py_dict=py_dict)

    assert mrsim == parsed_mrsim


def test_version_with_serialization():
    mrsim = setup_mrsimulator_obj()

    # Set to older version will throw error
    with pytest.raises(ValidationError):
        mrsim.version = "an older version"


def test_save():
    mrsim = setup_mrsimulator_obj()

    assert mrsim.simulator.methods[0].simulation is None
    assert mrsim.simulator.methods[1].simulation is None

    mrsim.save("temp.mrsim")


def test_load():
    mrsim = setup_mrsimulator_obj()
    loaded_mrsim = Mrsimulator.load("temp.mrsim")

    assert mrsim.simulator.methods[0].simulation is None
    assert mrsim.simulator.methods[1].simulation is None
    assert mrsim == loaded_mrsim

    os.remove("temp.mrsim")
