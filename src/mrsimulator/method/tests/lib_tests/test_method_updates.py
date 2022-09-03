import os

import mrsimulator.method.lib as NamedMethods
import pytest
from mrsimulator import Simulator
from mrsimulator import Site
from mrsimulator import SpinSystem
from mrsimulator.method import Method
from mrsimulator.spin_system.isotope import Isotope
from mrsimulator.spin_system.isotope import set_isotope_mappings
from mrsimulator.utils.error import ImmutableEventError

__author__ = "Deepansh J. Srivastava"
__email__ = "srivastava.89@osu.edu"

methods = [val for k, val in NamedMethods.__dict__.items() if isinstance(val, type)]

sys = SpinSystem(sites=[{"isotope": "1H"}])


def check_sim_save(sim1, sim2, message):
    for mth in sim2.methods:
        mth.simulation._timestamp = ""
        _ = [item.to("ppm", "nmr_frequency_ratio") for item in mth.simulation.x]

    dataset1 = sim1.methods[0].simulation.copy()
    sim1.methods[0].simulation = None
    dataset2 = sim2.methods[0].simulation.copy()
    sim2.methods[0].simulation = None

    e = f"dataset saved and loaded is not equal: type {message}."
    assert dataset1 == dataset2, e
    assert sim1 == sim2, f".mrsim saved and loaded is not equal: type {message}."


def check_methods_save_load(sim1, sim2):
    for mth in sim2.methods:
        mth.simulation._timestamp = ""
        _ = [item.to("ppm", "nmr_frequency_ratio") for item in mth.simulation.x]

    dataset1 = sim1.methods[0].simulation.copy()
    sim1.methods[0].simulation = None
    dataset2 = sim2.methods[0].simulation.copy()
    sim2.methods[0].simulation = None

    assert dataset1 == dataset2, "method saved and loaded is not equal."
    assert sim1 == sim2, ".mrmtd saved and loaded is not equal."


def assert_parsing(method, fn1):
    fn2 = method.parse_dict_with_units(fn1.json())
    assert fn1 == fn2, f"Error with {method} parse with units."

    fn3 = method(**fn1.json(units=False))
    assert fn1 == fn3, f"Error with {method} parse with units."

    event_error = "Event objects are immutable for"
    serialize = fn1.json()
    ent = serialize["spectral_dimensions"][0]["events"]
    ent[0]["transition_queries"][0]["ch1"]["P"] = [-100]

    if method not in [Method]:
        with pytest.raises(ImmutableEventError, match=f".*{event_error}.*"):
            method.parse_dict_with_units(serialize)

        with pytest.raises(ImmutableEventError, match=f".*{event_error}.*"):
            method(**serialize)

        with pytest.raises(ImmutableEventError, match=f".*{event_error}.*"):
            ent[0]["fraction"] = 0.5
            ent.append(
                {"fraction": 0.5, "transition_queries": [{"ch1": {"P": [-100]}}]}
            )
            method.parse_dict_with_units(serialize)


def test_read_write_methods():
    for item in methods:
        kwargs = {
            "channels": ["93Nb"],
            "spectral_dimensions": [
                {"count": 1024, "spectral_width": 100} for _ in range(item.ndim)
            ],
        }

        if item.__name__ == "SSB2D":
            kwargs["rotor_frequency"] = 1e3

        fn1 = item(**kwargs)

        # Parse against self
        assert_parsing(item, fn1)

        # Parse against Method class
        assert_parsing(Method, fn1)

        # save test
        sim = Simulator(spin_systems=[sys], methods=[fn1])

        # save/load with units
        sim.run()
        sim.save("test.mrsim")
        sim2 = Simulator.load("test.mrsim")
        check_sim_save(sim, sim2, "with units")

        # save/load with out units
        sim.run()
        sim.save("test.mrsim", with_units=False)
        sim2 = Simulator.load("test.mrsim", parse_units=False)
        check_sim_save(sim, sim2, "without units")

        # save/load methods
        sim.run()
        sim.export_methods("test.mrmtd")
        sim2 = sim.copy()
        sim2.load_methods("test.mrmtd")
        check_methods_save_load(sim, sim2)

        # save/load spin systems
        sim.run()
        sim.export_spin_systems("test.mrsys")
        sim2 = sim.copy()
        sim2.load_spin_systems("test.mrsys")
        assert sim.spin_systems == sim2.spin_systems

        # save/load spin systems with custom isotopes
        Isotope.new(symbol="C13_a", copy_from="13C")
        systems = [SpinSystem(sites=[Site(isotope=sym)]) for sym in ["13C", "C13_a"]]
        sim.spin_systems = systems
        sim.export_spin_systems("test.mrsys")
        set_isotope_mappings({})  # Reset all mappings
        sim.load_spin_systems("test.mrsys")
        assert sim.spin_systems == systems
        set_isotope_mappings({})  # Reset all mappings

    os.remove("test.mrsim")
    os.remove("test.mrsys")
    os.remove("test.mrmtd")
