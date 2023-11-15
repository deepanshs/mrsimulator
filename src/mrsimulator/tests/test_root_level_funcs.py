import os

import pytest
from monty.serialization import loadfn
from mrsimulator import __version__
from mrsimulator import dict
from mrsimulator import load
from mrsimulator import parse
from mrsimulator import save
from mrsimulator import signal_processor as sp
from mrsimulator import Simulator
from mrsimulator import Site
from mrsimulator import SpinSystem
from mrsimulator.method.lib import BlochDecaySpectrum
from mrsimulator.utils.error import FileConversionError


MODULE_DIR = os.path.dirname(os.path.abspath(__file__))
TEST_DATA = loadfn(os.path.join(MODULE_DIR, "test_data.json"))

# mrsimulator save and load test


def setup():
    site = Site(isotope="1H", shielding_symmetric={"zeta": -100, "eta": 0.3})
    spin_sys = SpinSystem(sites=[site, site], abundance=45)

    method = BlochDecaySpectrum(channels=["1H"])
    sim = Simulator(spin_systems=[spin_sys, spin_sys], methods=[method, method])

    sim.run(method_index=0)
    sim.methods[0].simulation._timestamp = None

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


def test_save():
    sim, processors, application = setup()

    assert sim.methods[0].simulation is not None
    assert sim.methods[1].simulation is None

    save(
        "test.mrsim",
        simulator=sim,
        signal_processors=processors,
        application=application,
    )


def test_dict():
    sim, processors, application = setup()
    sim.methods[0].simulation = None
    sim.methods[1].simulation = None

    py_dict = {
        "simulator": sim.json(),
        "signal_processors": [sp.json() for sp in processors],
        "application": application,
        "version": __version__,
    }

    assert py_dict == dict(sim, processors, application)


def test_load():
    sim_r, processors_r, application_r = load("test.mrsim")

    sim, processors, application = setup()
    sim_r.methods[0].simulation = None
    sim.methods[0].simulation = None
    assert sim_r == sim
    assert processors_r == processors
    assert application_r == application

    # Load from external URL. May break in the future
    load("https://ssnmr.org/sites/default/files/mrsimulator/test.mrsim")

    os.remove("test.mrsim")


def test_parse_old_struct():
    """Ensures error raised when trying to parse old file struct"""
    old_root_level = TEST_DATA["old_root_level"]

    with pytest.raises(FileConversionError):
        parse(old_root_level)


# def test_mrsim_to_v0_7():
#     sim, processors, application = setup()
#     sim.methods[0].simulation = None

#     old_struct = sim.json()
#     old_struct["signal_processors"] = [sp.json() for sp in processors]
#     old_struct["application"] = application
#     old_struct["version"] = __version__
#     old_struct["some_extra_key"] = "An erroneous key which will be removed"

#     with open("temp_2.mrsim", "w", encoding="utf8") as outfile:
#         json.dump(
#             old_struct,
#             outfile,
#             ensure_ascii=False,
#             sort_keys=False,
#             allow_nan=False,
#             separators=(",", ":"),
#         )

#     # Test error handling of loading old structure
#     e = (
#         "An incompatible JSON root-level structure was detected. Use the method "
#         "mrsim_to_v0_7 to convert to a compliant structure."
#     )
#     with pytest.raises(ValueError, match=e):
#         load("temp_2.mrsim")

#     new_struct = {
#         "simulator": sim.json(),
#         "signal_processors": [sp.json() for sp in processors],
#         "application": application,
#         "version": __version__,
#     }

#     py_dict = mrsim_to_v0_7("temp_2.mrsim", overwrite=True)

#     py_dict["simulator"]["methods"][0]["simulation"] = None
#     new_struct["simulator"]["methods"][0]["simulation"] = None

#     assert py_dict == new_struct

#     os.remove("temp_2.mrsim")
