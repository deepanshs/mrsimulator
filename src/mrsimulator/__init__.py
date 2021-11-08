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
__version__ = "0.7.0dev1"

import os

os.environ["OMP_NUM_THREADS"] = "1"
os.environ["OPENBLAS_NUM_THREADS"] = "1"
os.environ["MKL_NUM_THREADS"] = "1"
os.environ["VECLIB_MAXIMUM_THREADS"] = "1"
# os.environ["NUMEXPR_NUM_THREADS"] = "1"

from .spin_system import Site  # lgtm [py/import-own-module] # noqa:F401
from .spin_system import Coupling  # lgtm [py/import-own-module]  # noqa:F401
from .spin_system import SpinSystem  # lgtm [py/import-own-module] # noqa:F401
from .simulator import Simulator  # lgtm [py/import-own-module] # noqa:F401
from .method.spectral_dimension import (  # lgtm [py/import-own-module] # noqa:F401
    SpectralDimension,
)
from .method import Method  # lgtm [py/import-own-module] # noqa:F401
from mrsimulator.utils.importer import import_json
from mrsimulator.signal_processing import SignalProcessor
import json
from lmfit import Parameters
from typing import Dict
from typing import List
from copy import deepcopy


class Mrsimulator:
    """The Mrsimulator class.

    Attributes
    ----------

    simulator: A :ref:`simulator_api` object.

    signal_processor: A :ref:`signal_processing_api` object.

    params: An optional parameters object

    application: An optional dict holding metadata.
    """

    # TODO: Improve doc strings (examples, more description)
    # TODO: Decide where params key should live
    # TODO: Create api reference
    # TODO: Make pytests for new class
    # TODO: Create error messages when loading file in wrong place?

    simulator: Simulator = None
    signal_processors: SignalProcessor = None
    params: List[Parameters] = None
    application: Dict = None

    class Config:
        validate_assignment = True
        extra = "forbid"

    @classmethod
    def parse_dict_with_units(cls, py_dict: dict):
        """Parse the physical quantity from a dictionary reperesentation of the
        Mrsimulator object, where the physical quantity is expressed as a string with a
        number and a unit

        Args:
            dict py_dict: A required python dict object.

        Returns:
            A :ref:`mrsimulator_api` object.

        Example
        -------


        """
        py_copy_dict = deepcopy(py_dict)

        if "simulator" in py_copy_dict:
            sim = Simulator.parse_dict_with_units(py_copy_dict["simulator"])
            py_copy_dict["simulator"] = sim

        if "signal_processors" in py_copy_dict:
            sps = [
                SignalProcessor.parse_dict_with_units(sp)
                for sp in py_copy_dict["signal_processors"]
            ]
            py_copy_dict["signal_processors"] = sps

        if "params" in py_copy_dict:
            params = Parameters().loads(py_copy_dict["params"])
            py_copy_dict["params"] = params

        return Mrsimulator(**py_copy_dict)

    def save(self, filename: str, with_units: bool = True):
        """Serialize the Mrsimulator object to a JSON compliant file

        Args:
            bool with_units: If true, the attribute values are serialized as physical
                quantities expressed as a string with a value and a unit. If false, the
                attribute values are serialized as floats. The parameters object is
                serialized to a LMFIT compliant string (add reference).
            str filename: The filename for the serialized file.

        Example
        -------


        """
        py_dict = self.dict(with_units=with_units)

        with open(filename, "w", encoding="utf8") as outfile:
            json.dump(
                py_dict,
                outfile,
                ensure_ascii=False,
                sort_keys=False,
                allow_nan=False,
                separators=(",", ":"),
            )

    def dict(self, with_units: bool = True):
        """Export the Mrsimulator object to a python dictionary.

        Args:
            bool with_units: If true, physical quantities are represented as string with
                units. The default is True.

        Returns:
            Python dictionary

        Example
        -------



        """
        py_dict = {
            "simulator": None,
            "signal_processors": None,
            "params": None,
            "application": self.application,
        }

        if self.simulator is not None:
            py_dict["simulator"] = self.simulator.json(with_units=with_units)

        if self.signal_processors is not None:
            py_dict["signal_processors"] = [sp.json() for sp in self.signal_processors]

        if self.params is not None:
            py_dict["params"] = self.params.dumps()

        return py_dict


def save(
    filename: str,
    simulator: Simulator,
    signal_processors: list = None,
    params: Parameters = None,
    with_units: bool = True,
):
    """Serialize the Simulator, list of SignalProcessor, and lmfit Parameters objects
    to a .mrsim file. Creates a Mrsimulator object and calls save.

    Args:
        str filename: The data is serialized to this file.
        sim: Simulator object.
        signal_processors: A list of PostSimulator objects corresponding to the methods
            in the Simulator object. Default is None.
        params: lmfit Parameters object. Default is None.
        bool with_units: If true, physical quantities are represented as string with
            units. The default is True.
    """
    Mrsimulator(
        simulator=simulator, signal_processors=signal_processors, params=params
    ).save(filename=filename, with_units=with_units)


def dict(
    simulator: Simulator,
    signal_processors: list = None,
    params: Parameters = None,
    with_units: bool = True,
):
    """Export the Simulator, list of SignalProcessor, and lmfit Parameters objects
    to a python dictionary. Creates a Mrsimulator object and calls dict.

    Args:
        sim: Simulator object.
        signal_processors: A list of PostSimulator objects corresponding to the methods
            in the Simulator object. Default is None.
        params: lmfit Parameters object. Default is None.
        bool with_units: If true, physical quantities are represented as string with
            units. The default is True.

    Returns:
        Python dictionary
    """
    return Mrsimulator(
        simulator=simulator, signal_processors=signal_processors, params=params
    ).dict(with_units=with_units)


def load(filename: str, parse_units: bool = True):
    """Load Simulator, list of SignalProcessor and optionally lmfit Parameters objects
    from the .mrsim file.

    Args:
        str filename: The location to the .mrsim file.
        bool parse_units: If true, parse the dictionary for units. The default is True.

    Return:
        Ordered List: Simulator, List[SignalProcessor], Parameters.
    """
    val = import_json(filename)
    return parse(val, parse_units)


def parse(py_dict, parse_units: bool = True):
    """Parse the dictionary object to respective Simulator, SignalProcessor and
    optionally lmfit Parameters object.

    Args:
        dict py_dict: Python dictionary representation of mrsimulator.
        bool parse_units: If true, parse the dictionary for units. Default is True.

    Return:
        Ordered List: Simulator, List[SignalProcessor], Parameters.
    """
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
