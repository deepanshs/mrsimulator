# -*- coding: utf-8 -*-
"""
Solid-state NMR spectra simulation module for Python
====================================================

mrsimulator is an incredibly fast solid-state NMR spectra simulation package,
capable of simulating spectra from both crystalline and amorphous-like materials.


It aims to provide simple and efficient solutions to the solid-state NMR spectrum
simulation problem. It includes tools for users to create their model spin systems,
simulate and compare the simulation with the measurement, and perform least-squares
minimization, using a collection of pre-defined NMR methods.

See https://mrsimulator.readthedocs.io/en/stable/ for complete documentation.
"""
# version has to be specified at the start.
__author__ = "Deepansh J. Srivastava"
__email__ = "srivastava.89@osu.edu"
__copyright__ = "Copyright 2019-2021, The mrsimulator Project."
__credits__ = ["Deepansh J. Srivastava"]
__license__ = "BSD License"
__maintainer__ = "Deepansh J. Srivastava"
__status__ = "Beta"
__version__ = "0.6.0rc3"

import os

os.environ["OMP_NUM_THREADS"] = "1"
os.environ["OPENBLAS_NUM_THREADS"] = "1"
os.environ["MKL_NUM_THREADS"] = "1"
os.environ["VECLIB_MAXIMUM_THREADS"] = "1"
# os.environ["NUMEXPR_NUM_THREADS"] = "1"

from .spin_system import Site  # lgtm [py/import-own-module] # noqa:F401
from .spin_system import Coupling  # lgtm [py/import-own-module]  # noqa:F401
from .spin_system import SpinSystem  # lgtm [py/import-own-module] # noqa:F401
from .simulator import Simulator
from .method.event import Event  # lgtm [py/import-own-module] # noqa:F401
from .method.spectral_dimension import (  # lgtm [py/import-own-module] # noqa:F401
    SpectralDimension,
)
from .method import Method  # lgtm [py/import-own-module] # noqa:F401
from mrsimulator.utils.importer import import_json
from mrsimulator.signal_processing import SignalProcessor
import json
from lmfit import Parameters


def save(
    filename: str,
    simulator: Simulator,
    signal_processors=None,
    fit_report=None,
    with_units: bool = True,
):
    sim = simulator.json(True, True) if with_units else simulator.reduced_dict()
    if signal_processors is not None:
        sim["signal_processors"] = [item.json() for item in signal_processors]

    sim["params"] = None if fit_report is None else fit_report.params.dumps()

    with open(filename, "w", encoding="utf8") as outfile:
        json.dump(
            sim,
            outfile,
            ensure_ascii=False,
            sort_keys=False,
            allow_nan=False,
            separators=(",", ":"),
        )


def load(filename: str, parse_units: bool = True):
    val = import_json(filename)
    return parse(val, parse_units)


def parse(py_dict, parse_units: bool = True):
    sim = Simulator.parse(py_dict, parse_units)

    signal_processors = (
        [
            SignalProcessor.parse_dict_with_units(item)
            for item in py_dict["signal_processors"]
        ]
        if "signal_processors" in py_dict
        else [SignalProcessor() for _ in sim.methods]
    )

    params = None
    if "params" in py_dict:
        val = py_dict["params"]
        params = None if val is None else Parameters().loads(s=val)

    return sim, signal_processors, params
