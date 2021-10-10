# -*- coding: utf-8 -*-
import os

from mrsimulator import methods as NamedMethods
from mrsimulator import Simulator
from mrsimulator import SpinSystem

__author__ = "Deepansh J. Srivastava"
__email__ = "srivastava.89@osu.edu"

methods = [val for k, val in NamedMethods.__dict__.items() if isinstance(val, type)]

sys = SpinSystem(sites=[{"isotope": "1H"}])
sp = {"count": 1024, "spectral_width": 100}


def test_read_write_obj_mrsim():
    for item in methods:
        kwargs = {"channels": ["93Nb"], "spectral_dimensions": [sp] * item.ndim}

        if item.__name__ == "SSB2D":
            kwargs["rotor_frequency"] = 1e3

        fn1 = item(**kwargs)
        fn2 = item.parse_dict_with_units(fn1.json())
        assert fn1 == fn2, f"Error with {item} parse."

        sim = Simulator(spin_systems=[sys], methods=[fn1])
        sim.run()

        # Test simulator object (.mrsim)
        sim.save("test.mrsim")

        sim2 = Simulator.load("test.mrsim")
        for mth in sim2.methods:
            mth.simulation._timestamp = ""
            _ = [item.to("ppm", "nmr_frequency_ratio") for item in mth.simulation.x]

        assert sim == sim2, f"Error with {item} parse from Simulator.load()."

        # Test spin system objects (.mrsys)
        sim.export_spin_systems("test.mrsys")

        sim2 = sim.copy()
        sim2.load_spin_systems("test.mrsys")

        assert (
            sim.spin_systems == sim2.spin_systems
        ), f"Error with {item} sim.load_spin_systems()."

        # Test methods objects (.mrmth)
        sim.export_methods("test.mrmth")

        sim2 = sim.copy()
        sim2.load_methods("test.mrmth")

        # remove timestamps from csdm obj for testing
        for mth in sim2.methods:
            for item in [mth.simulation, mth.experiment]:
                if item is not None:
                    item._timestamp = ""
            _ = [item.to("ppm", "nmr_frequency_ratio") for item in mth.simulation.x]

        assert sim.methods == sim2.methods, f"Error with {item} sim.load_methods()."

    os.remove("test.mrsim")
    os.remove("test.mrsys")
    os.remove("test.mrmth")
