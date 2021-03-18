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


def test_read_write_methods():
    for item in methods:
        kwargs = {"channels": ["93Nb"], "spectral_dimensions": [sp] * item.ndim}

        if item.__name__ == "SSB2D":
            kwargs["rotor_frequency"] = 1e3

        fn1 = item(**kwargs)
        fn2 = item.parse_dict_with_units(fn1.json())
        assert fn1 == fn2, f"Error with {item} parse."

        sim = Simulator(spin_systems=[sys], methods=[fn1])
        sim.run()
        sim.save("test.mrsim")

        sim2 = Simulator.load("test.mrsim")
        for mth in sim2.methods:
            mth.simulation._timestamp = ""
            _ = [item.to("ppm", "nmr_frequency_ratio") for item in mth.simulation.x]

        data1 = sim.methods[0].simulation.copy()
        sim.methods[0].simulation = None
        data2 = sim2.methods[0].simulation.copy()
        sim2.methods[0].simulation = None

        assert data1 == data2, f"data saved and loaded is not equal: method {item}."
        assert sim == sim2, f"Error in method {item} when using Simulator.load()."

    os.remove("test.mrsim")
