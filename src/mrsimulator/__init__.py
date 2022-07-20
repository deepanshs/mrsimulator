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
from datetime import datetime

# version has to be specified at the start.
__author__ = "Deepansh J. Srivastava"
__email__ = "srivastava.89@osu.edu"
__copyright__ = f"Copyright 2019-{datetime.now().year}, The mrsimulator Project."
__credits__ = ["Deepansh J. Srivastava"]
__license__ = "BSD License"
__maintainer__ = "Deepansh J. Srivastava"
__status__ = "Beta"
__version__ = "0.7.0"

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
from mrsimulator.signal_processor import SignalProcessor
from mrsimulator.utils.error import FileConversionError
from mrsimulator.utils.importer import import_json
from mrsimulator.utils.parseable import Parseable
import json
from pydantic import Field
from typing import Dict
from typing import List
from copy import deepcopy


class Mrsimulator(Parseable):
    """The Mrsimulator class.

    Attributes
    ----------

    simulator: A :ref:`simulator_api` object.

    signal_processor: A list of :ref:`signal_processor_api` objects.

    version: The current version of the library represented by a string. This attribute
        is not meant to be modified and serialization will only reflect the current
        version.

    application: An optional dictionary holding metadata.
    """

    simulator: Simulator = None
    signal_processors: List[SignalProcessor] = None
    version: str = Field(__version__, const=True)
    application: Dict = None

    class Config:
        validate_assignment = True
        extra = "forbid"

    @classmethod
    def parse_dict_with_units(cls, py_dict: dict):
        """Parse the physical quantity from a dictionary representation of the
        Mrsimulator object, where the physical quantity is expressed as a string with a
        number and a unit

        Args:
            Dict py_dict: A required python dictionary object.

        Returns:
            A :ref:`mrsimulator_api` object.

        Example
        -------
        >>> mrsim_dict = {
        ...   "simulator": {
        ...     "spin_systems": [
        ...       {
        ...         "sites": [
        ...           {
        ...             "isotope": "1H",
        ...             "isotropic_chemical_shift": "0.0 ppm",
        ...             "shielding_symmetric": { "zeta": "-100.0 ppm", "eta": 0.3 }
        ...           },
        ...           {
        ...             "isotope": "1H",
        ...             "isotropic_chemical_shift": "0.0 ppm",
        ...             "shielding_symmetric": { "zeta": "-100.0 ppm", "eta": 0.3 }
        ...           }
        ...         ],
        ...         "abundance": "45.0 %"
        ...       }
        ...     ],
        ...     "methods": [
        ...       {
        ...         "name": "BlochDecaySpectrum",
        ...         "description": "A one-dimensional Bloch decay spectrum method.",
        ...         "channels": ["1H"],
        ...         "spectral_dimensions": [
        ...           {
        ...             "count": 1024,
        ...             "spectral_width": "25000.0 Hz",
        ...             "events": [{ "transition_queries": [{ "ch1": { "P": [-1] } }] }]
        ...           }
        ...         ],
        ...         "magnetic_flux_density": "9.4 T",
        ...         "rotor_angle": "0.9553166181245 rad",
        ...         "rotor_frequency": "0.0 Hz"
        ...       }
        ...     ],
        ...     "config": {
        ...       "number_of_sidebands": 64,
        ...       "integration_volume": "octant",
        ...       "integration_density": 70,
        ...       "decompose_spectrum": "none"
        ...     }
        ...   },
        ...   "signal_processors": [
        ...     {
        ...       "operations": [
        ...         { "dim_index": 0, "function": "IFFT" },
        ...         {
        ...           "dim_index": 0,
        ...           "FWHM": "2500.0 Hz",
        ...           "function": "apodization",
        ...           "type": "Exponential"
        ...         },
        ...         { "dim_index": 0, "function": "FFT" },
        ...         { "factor": 20.0, "function": "Scale" }
        ...       ]
        ...     }
        ...   ],
        ...   "application": {
        ...     "com.github.DeepanshS.mrsimulator": { "foo": "This is some metadata" }
        ...   }
        ... }
        >>> mrsim = Mrsimulator.parse_dict_with_units(mrsim_dict)
        """
        py_copy_dict = deepcopy(py_dict)

        # Ignore given version by removing key
        _ = py_copy_dict.pop("version", None)

        if "simulator" in py_copy_dict:
            sim = Simulator.parse_dict_with_units(py_copy_dict["simulator"])
            py_copy_dict["simulator"] = sim

        if "signal_processors" in py_copy_dict:
            processors = [
                SignalProcessor.parse_dict_with_units(sp)
                for sp in py_copy_dict["signal_processors"]
            ]
            py_copy_dict["signal_processors"] = processors

        return Mrsimulator(**py_copy_dict)

    @classmethod
    def parse(cls, py_dict: dict, parse_units: bool = True):
        """Parse a dictionary to a Mrsimulator object.

        Args:
            Dict py_dict: Dictionary representation of the Mrsimulator object
            bool parse_units: If true, parse quantity from units string.

        Returns:
            A :ref:`mrsimulator_api` object.
        """
        py_copy_dict = deepcopy(py_dict)

        # Ignore given version by removing key
        _ = py_copy_dict.pop("version", None)

        return (
            Mrsimulator.parse_dict_with_units(py_dict)
            if parse_units
            else Mrsimulator(**py_dict)
        )

    @classmethod
    def load(cls, filename: str, with_units: bool = True):
        """Load the :py:class: `~mrsimulator.Mrsimulator` object from a JSON file by
        parsing a file from a given path

        Args:
            bool parse_units: If true, parse the attribute values from the serialized
                file for physical quantities, expressed as a string with a value and a
                unit.
            str filename: The filename of a JSON serialized Mrsimulator object file.

        Example
        -------

        >>> mrsim = Mrsimulator.load("filename") # doctest: +SKIP
        """
        contents = import_json(filename)
        return Mrsimulator.parse(contents, with_units)

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

        >>> mrsim.save("filename") # doctest: +SKIP
        """
        py_dict = self.json(with_units=with_units)

        with open(filename, "w", encoding="utf8") as outfile:
            json.dump(
                py_dict,
                outfile,
                ensure_ascii=False,
                sort_keys=False,
                allow_nan=False,
                separators=(",", ":"),
            )

    def json(self, with_units: bool = True):
        """Export the Mrsimulator object to a python dictionary.

        Args:
            bool with_units: If true, physical quantities are represented as string with
                units. The default is True.

        Returns:
            Python dictionary
        """
        py_dict = {
            "simulator": None,
            "signal_processors": None,
            "version": __version__,
            "application": self.application,
        }

        if self.simulator is not None:
            py_dict["simulator"] = self.simulator.json(units=with_units)

        if self.signal_processors is not None:
            py_dict["signal_processors"] = [sp.json() for sp in self.signal_processors]

        return py_dict


# ================================= Root Level Methods =================================
def save(
    filename: str,
    simulator: Simulator,
    signal_processors: List = None,
    application: Dict = None,
    with_units: bool = True,
):
    """Serialize the Simulator, list of SignalProcessor, and an application dict
    to a file. Creates a Mrsimulator object and calls save.

    Args:
        str filename: The data is serialized to this file.
        sim: Simulator object.
        signal_processors: A list of PostSimulator objects corresponding to the methods
            in the Simulator object. Default is None.
        application: Dictionary holding metadata to serialize in the file. The
            dictionary will be held in the application key.
        bool with_units: If true, physical quantities are represented as string with
            units. The default is True.
    """
    Mrsimulator(
        simulator=simulator,
        signal_processors=signal_processors,
        application=application,
    ).save(filename=filename, with_units=with_units)


def dict(
    simulator: Simulator,
    signal_processors: List = None,
    application: Dict = None,
    with_units: bool = True,
):
    """Export the Simulator, list of SignalProcessor, and an application dict
    to a python dictionary. Creates a Mrsimulator object with given arguments and calls
    json from the Mrsimulator object.

    Args:
        sim: Simulator object.
        signal_processors: A list of PostSimulator objects corresponding to the methods
            in the Simulator object. Default is None.
        application: Dictionary holding metadata to serialize in the dict. The
            dictionary will be held under the application key.
        bool with_units: If true, physical quantities are represented as string with
            units. The default is True.

    Returns:
        Python dictionary
    """
    return Mrsimulator(
        simulator=simulator,
        signal_processors=signal_processors,
        application=application,
    ).json(with_units=with_units)


def load(filename: str, parse_units: bool = True):
    """Load Simulator object, list of SignalProcessor objects and metadata from a JSON
    serialized file of a :py:class:`~mrsimulator.Mrsimulator` object.

    Args:
        str filename: The location to the .mrsim file.
        bool parse_units: If true, parse the dictionary for units. The default is True.

    Return:
        Ordered List: Simulator, List[SignalProcessor], Dict.
    """
    val = import_json(filename)
    return parse(val, parse_units)


def parse(py_dict, parse_units: bool = True):
    """Parse a dictionary object to the respective Simulator object, list of
    SignalProcessor objects, and the metadata dictionary. If no signal processors are
    provided a list of default SignalProcessor objects with length equal to number of
    methods will be returned.

    Args:
        Dict py_dict: Python dictionary representation of a
            :py:class:`~mrsimulator.Mrsimulator` object.
        bool parse_units: If true, parse the dictionary for units. Default is True.

    Return:
        Ordered List: Simulator, List[SignalProcessor], Dict.
    """
    # Check for difference in keys
    root_keys = set(Mrsimulator().dict().keys())
    if len(set(py_dict.keys()) - root_keys) != 0:
        raise FileConversionError()

    py_copy_dict = deepcopy(py_dict)
    sim = Simulator.parse(py_copy_dict["simulator"], parse_units)

    signal_processors = (
        [
            SignalProcessor.parse_dict_with_units(item)
            for item in py_copy_dict["signal_processors"]
        ]
        if "signal_processors" in py_copy_dict
        else [SignalProcessor() for _ in sim.methods]
    )

    application = py_copy_dict["application"] if "application" in py_copy_dict else None

    return sim, signal_processors, application
