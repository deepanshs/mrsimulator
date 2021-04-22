# -*- coding: utf-8 -*-
"""Base Simulator class."""
import json
from copy import deepcopy
from typing import List

import csdmpy as cp
import numpy as np
import pandas as pd
import psutil
from joblib import delayed
from joblib import Parallel
from mrsimulator import __version__
from mrsimulator import methods as NamedMethods
from mrsimulator import Site
from mrsimulator import SpinSystem
from mrsimulator.base_model import one_d_spectrum
from mrsimulator.method import Method
from mrsimulator.spin_system.isotope import Isotope
from mrsimulator.utils import flatten_dict
from mrsimulator.utils.abstract_list import AbstractList
from mrsimulator.utils.extra import _reduce_dict
from mrsimulator.utils.importer import import_json
from pydantic import BaseModel

from .config import ConfigSimulator

# from IPython.display import JSON

__author__ = "Deepansh Srivastava"
__email__ = "srivastava.89@osu.edu"

__CPU_count__ = psutil.cpu_count()

__named_methods__ = [
    val for k, val in NamedMethods.__dict__.items() if isinstance(val, type)
]
__method_names__ = [item.__name__ for item in __named_methods__]
__sim_methods__ = {k: v for k, v in zip(__method_names__, __named_methods__)}


class Simulator(BaseModel):
    """
    The simulator class.

    Attributes
    ----------

    spin_systems: A list of :ref:`spin_sys_api` or equivalent dict objects (optional).
        A list of :ref:`spin_sys_api` or equivalent dict objects representing a
        collection of isolated NMR spin systems present within the sample. The default
        value is an empty list.

        Example
        -------

        >>> sim = Simulator()
        >>> sim.spin_systems = [
        ...     SpinSystem(sites=[Site(isotope='17O')], abundance=0.015),
        ...     SpinSystem(sites=[Site(isotope='1H')], abundance=1),
        ... ]

        >>> # or equivalently
        >>> sim.spin_systems = [
        ...     {'sites': [{'isotope': '17O'}], 'abundance': 0.015},
        ...     {'sites': [{'isotope': '1H'}], 'abundance': 1},
        ... ]

    methods: A list of :ref:`method_api` or equivalent dict objects (optional).
        A list of :ref:`method_api`  or equivalent dict objects representing an NMR
        methods. The default value is an empty list.

        Example
        -------

        >>> from mrsimulator.methods import BlochDecaySpectrum
        >>> from mrsimulator.methods import BlochDecayCTSpectrum
        >>> sim.methods = [
        ...     BlochDecaySpectrum(channels=['17O'], spectral_width=50000),
        ...     BlochDecayCTSpectrum(channels=['17O'])
        ... ]

    config: :ref:`config_api` object or equivalent dict object (optional).
        The :ref:`config_api` object is used to configure the simulation. The valid
        attributes of the ConfigSimulator object are

        - ``number_of_sidebands``,
        - ``integration_density``,
        - ``integration_volume``, and
        - ``decompose_spectrum``

        Example
        -------

        >>> from mrsimulator.simulator.config import ConfigSimulator
        >>> sim.config = ConfigSimulator(
        ...     number_of_sidebands=32,
        ...     integration_density=64,
        ...     integration_volume='hemisphere',
        ...     decompose_spectrum='spin_system',
        ... )

        >>> # or equivalently
        >>> sim.config = {
        ...     'number_of_sidebands': 32,
        ...     'integration_density': 64,
        ...     'integration_volume': 'hemisphere',
        ...     'decompose_spectrum': 'spin_system',
        ... }

        See :ref:`config_simulator` for details.

    name: str (optional).
        The name or id of the simulation or sample. The default value is None.

        Example
        -------

        >>> sim.name = '1H-17O'
        >>> sim.name
        '1H-17O'

    label: str (optional).
        The label for the simulation or sample. The default value is None.

        Example
        -------

        >>> sim.label = 'Test simulator'
        >>> sim.label
        'Test simulator'

    description: str (optional).
        A description of the simulation or sample. The default value is None.

        Example
        -------

        >>> sim.description = 'Simulation for sample 1'
        >>> sim.description
        'Simulation for sample 1'

    """

    name: str = None
    label: str = None
    description: str = None
    spin_systems: List[SpinSystem] = []
    methods: List[Method] = []
    config: ConfigSimulator = ConfigSimulator()
    indexes = []

    class Config:
        validate_assignment = True

    def __eq__(self, other):
        if not isinstance(other, Simulator):
            return False
        check = [
            self.name == other.name,
            self.label == other.label,
            self.description == other.description,
            self.spin_systems == other.spin_systems,
            self.methods == other.methods,
            self.config == other.config,
            np.all(self.indexes == other.indexes),
        ]
        if np.all(check):
            return True
        return False

    @classmethod
    def parse_dict_with_units(cls, py_dict: dict):
        """
        Parse the physical quantity from a dictionary representation of the Simulator
        object, where the physical quantity is expressed as a string with a number and
        a unit.

        Args:
            dict py_dict: A required python dict object.

        Returns:
            A :ref:`simulator_api` object.

        Example
        -------

        >>> sim_py_dict = {
        ...     'config': {
        ...         'decompose_spectrum': 'none',
        ...         'integration_density': 70,
        ...         'integration_volume': 'octant',
        ...         'number_of_sidebands': 64
        ...     },
        ...     'spin_systems': [
        ...         {
        ...             'abundance': '100 %',
        ...             'sites': [{
        ...                 'isotope': '13C',
        ...                 'isotropic_chemical_shift': '20.0 ppm',
        ...                 'shielding_symmetric': {'eta': 0.5, 'zeta': '10.0 ppm'}
        ...             }]
        ...         },
        ...         {
        ...             'abundance': '100 %',
        ...             'sites': [{
        ...                 'isotope': '1H',
        ...                     'isotropic_chemical_shift': '-4.0 ppm',
        ...                     'shielding_symmetric': {'eta': 0.1, 'zeta': '2.1 ppm'}
        ...             }]
        ...         },
        ...         {
        ...             'abundance': '100 %',
        ...             'sites': [{
        ...                 'isotope': '27Al',
        ...                 'isotropic_chemical_shift': '120.0 ppm',
        ...                 'shielding_symmetric': {'eta': 0.1, 'zeta': '2.1 ppm'}
        ...             }]
        ...         }
        ...     ]
        ... }
        >>> sim = Simulator.parse_dict_with_units(sim_py_dict)
        >>> len(sim.spin_systems)
        3
        """
        py_copy_dict = deepcopy(py_dict)

        if "spin_systems" in py_copy_dict:
            spin_sys = py_copy_dict["spin_systems"]
            spin_sys = [SpinSystem.parse_dict_with_units(obj) for obj in spin_sys]
            py_copy_dict["spin_systems"] = spin_sys

        if "methods" in py_copy_dict:
            methods = py_copy_dict["methods"]
            method_cls = [
                Method
                if obj["name"] not in __method_names__
                else __sim_methods__[obj["name"]]
                for obj in methods
            ]

            methods = [
                fn.parse_dict_with_units(obj) for obj, fn in zip(methods, method_cls)
            ]
            py_copy_dict["methods"] = methods

        return Simulator(**py_copy_dict)

    def get_isotopes(self, spin_I: float = None, symbol: bool = False) -> list:
        """
        List of unique isotopes from the sites within the list of the spin systems
        corresponding to spin quantum number `I`. If `I` is None, a list of all unique
        isotopes is returned instead.

        Args:
            float spin_I: An optional spin quantum number. The valid input are the
                multiples of 0.5.
            bool symbol: If true, return a list of str with isotope symbols.

        Returns:
            A list of :ref:`isotope_api` objects.

        Example
        -------

        >>> sim.get_isotopes()
        [Isotope(symbol='13C'), Isotope(symbol='1H'), Isotope(symbol='27Al')]
        >>> sim.get_isotopes(symbol=True)
        ['13C', '1H', '27Al']

        >>> sim.get_isotopes(spin_I=0.5)
        [Isotope(symbol='13C'), Isotope(symbol='1H')]
        >>> sim.get_isotopes(spin_I=0.5, symbol=True)
        ['13C', '1H']

        >>> sim.get_isotopes(spin_I=1.5)
        []

        >>> sim.get_isotopes(spin_I=2.5)
        [Isotope(symbol='27Al')]
        >>> sim.get_isotopes(spin_I=2.5, symbol=True)
        ['27Al']
        """
        st = []
        for sys in self.spin_systems:
            st += sys.get_isotopes(spin_I, symbol=True)
        st = np.unique(st)
        if not symbol:
            return [Isotope(symbol=item) for item in st]
        return list(st)

    def json(self, include_methods: bool = False, include_version: bool = False):
        """Parse the class object to a JSON compliant python dictionary object, where
        the attribute value with physical quantity is expressed as a string with a
        value and a unit.

        Args:
            bool include_methods: If True, the output dictionary will include the
                serialized method objects. The default value is False.
            bool include_version: If True, add a version key-value pair to the
                serialized output dictionary. The default is False.

        Returns:
            A Dict object.

        Example
        -------

        >>> pprint(sim.json())
        {'config': {'decompose_spectrum': 'none',
                    'integration_density': 70,
                    'integration_volume': 'octant',
                    'number_of_sidebands': 64},
         'spin_systems': [{'abundance': '100.0 %',
                           'sites': [{'isotope': '13C',
                                      'isotropic_chemical_shift': '20.0 ppm',
                                      'shielding_symmetric': {'eta': 0.5,
                                                              'zeta': '10.0 ppm'}}]},
                          {'abundance': '100.0 %',
                           'sites': [{'isotope': '1H',
                                      'isotropic_chemical_shift': '-4.0 ppm',
                                      'shielding_symmetric': {'eta': 0.1,
                                                              'zeta': '2.1 ppm'}}]},
                          {'abundance': '100.0 %',
                           'sites': [{'isotope': '27Al',
                                      'isotropic_chemical_shift': '120.0 ppm',
                                      'shielding_symmetric': {'eta': 0.1,
                                                              'zeta': '2.1 ppm'}}]}]}
        """
        sim = {}
        sim["name"] = self.name
        sim["description"] = self.description
        sim["label"] = self.label
        sim["spin_systems"] = [_.json() for _ in self.spin_systems]
        sim["methods"] = [_.json() for _ in self.methods] if include_methods else None
        sim["config"] = self.config.dict()
        sim["version"] = __version__ if include_version else None

        _ = [sim.pop(k) for k in [k for k in sim.keys() if sim[k] is None]]
        return sim

    def reduced_dict(self, exclude=["property_units", "indexes"]) -> dict:
        """Returns a reduced dictionary representation of the class object by removing
        all key-value pair corresponding to keys listed in the `exclude` argument, and
        keys with value as None.

        Args:
            list exclude: A list of keys to exclude from the dictionary.
        Return: A dict.
        """
        return _reduce_dict(self.dict(), exclude)

    # def pretty(self):
    #     return JSON(self.json(include_methods=True, include_version=True))

    def load_spin_systems(self, filename: str):
        """
        Load a list of spin systems from the given JSON serialized file.

        See an
        `example <https://raw.githubusercontent.com/DeepanshS/mrsimulator-examples/
        master/spin_systems_v0.3.json>`_ of a JSON serialized file. For details, refer
        to the :ref:`load_spin_systems` section of this documentation.

        Args:
            str filename: A local or remote address to a JSON serialized file.

        Example
        -------

        >>> sim.load_spin_systems(filename) # doctest:+SKIP
        """
        contents = import_json(filename)
        # json_data = contents["spin_systems"]
        self.spin_systems = [SpinSystem.parse_dict_with_units(obj) for obj in contents]

    def export_spin_systems(self, filename: str):
        """
        Export a list of spin systems to a JSON serialized file.

        See an
        `example <https://raw.githubusercontent.com/DeepanshS/mrsimulator-examples/
        master/spin_systems_v0.3.json>`_ of a JSON serialized file. For details, refer
        to the :ref:`load_spin_systems` section.

        Args:
            str filename: A filename of the serialized file.

        Example
        -------

        >>> sim.export_spin_systems(filename) # doctest:+SKIP
        """
        spin_sys = [SpinSystem.json(obj) for obj in self.spin_systems]
        with open(filename, "w", encoding="utf8") as outfile:
            json.dump(
                spin_sys, outfile, ensure_ascii=False, sort_keys=False, allow_nan=False
            )

    def run(
        self,
        method_index: list = None,
        n_jobs: int = 1,
        pack_as_csdm: bool = True,
        **kwargs,
    ):
        """Run the simulation and compute spectrum.

        Args:
            method_index: An integer or a list of integers. If provided, only the
                simulations corresponding to the methods at the given index/indexes
                will be computed. The default is None, `i.e.`, the simulation for
                all method will be computed.
            bool pack_as_csdm: If true, the simulation results are stored as a
                `CSDM <https://csdmpy.readthedocs.io/en/stable/api/CSDM.html>`_ object,
                otherwise, as a `ndarray
                <https://numpy.org/doc/1.18/reference/generated/numpy.ndarray.html>`_
                object.
                The simulations are stored as the value of the
                :attr:`~mrsimulator.Method.simulation` attribute of the corresponding
                method.

        Example
        -------

        >>> sim.run() # doctest:+SKIP
        """
        verbose = 0
        if method_index is None:
            method_index = np.arange(len(self.methods))
        elif isinstance(method_index, int):
            method_index = [method_index]
        for index in method_index:
            method = self.methods[index]
            spin_sys = get_chunks(self.spin_systems, n_jobs)
            kwargs_dict = self.config.get_int_dict()
            jobs = (
                delayed(one_d_spectrum)(
                    method=method, spin_systems=sys, **kwargs_dict, **kwargs
                )
                for sys in spin_sys
            )
            amp = Parallel(
                n_jobs=n_jobs,
                verbose=verbose,
                backend="loky",
                # **{
                #     "backend": {
                #         "threads": "threading",
                #         "processes": "multithreading",
                #         None: None,
                #     }["threads"]
                # },
            )(jobs)

            # self.indexes.append(indexes)

            gyromagnetic_ratio = method.channels[0].gyromagnetic_ratio
            B0 = method.spectral_dimensions[0].events[0].magnetic_flux_density
            origin_offset = np.abs(B0 * gyromagnetic_ratio * 1e6)
            for seq in method.spectral_dimensions:
                seq.origin_offset = origin_offset

            if isinstance(amp[0], list):
                simulated_data = []
                for item in amp:
                    simulated_data += item
            if isinstance(amp[0], np.ndarray):
                simulated_data = [np.asarray(amp).sum(axis=0)]

            method.simulation = (
                self._as_csdm_object(simulated_data, method)
                if pack_as_csdm
                else np.asarray(simulated_data)
            )

    def save(self, filename: str, with_units: bool = True):
        """Serialize the simulator object to a JSON file.

        Args:
            bool with_units: If true, the attribute values are serialized as physical
                quantities expressed as a string with a value and a unit. If false,
                the attribute values are serialized as floats.
            str filename: The filename of the serialized file.

        Example
        -------

        >>> sim.save('filename') # doctest: +SKIP
        """
        sim = self.json(True, True) if with_units else self.reduced_dict()
        with open(filename, "w", encoding="utf8") as outfile:
            json.dump(
                sim, outfile, ensure_ascii=False, sort_keys=False, allow_nan=False
            )

    @classmethod
    def load(cls, filename: str, parse_units: bool = True):
        """Load the :class:`~mrsimulator.Simulator` object from a JSON file by parsing.

        Args:
            bool parse_units: If true, parse the attribute values from the serialized
                file for physical quantities, expressed as a string with a value and a
                unit.
            str filename: The filename of a JSON serialized mrsimulator file.

        Returns:
            A :class:`~mrsimulator.Simulator` object.

        Example
        -------

        >>> sim_1 = sim.load('filename') # doctest: +SKIP

        .. seealso::
            :ref:`load_spin_systems`
        """
        val = import_json(filename)
        return Simulator.parse(val, parse_units)

    @classmethod
    def parse(cls, py_dict: dict, parse_units: bool = True):
        """Parse a dictionary for Simulator object.

        Args:
            dict py_dict: Dictionary object.
            bool parse_units: It true, parse quantity from string.
        """
        return (
            Simulator.parse_dict_with_units(py_dict)
            if parse_units
            else Simulator(**py_dict)
        )

    def sites(self):
        """Unique sites within the Simulator object as a list of Site objects.

        Returns:
            A :ref:`sites_api` object.

        Example
        -------

        >>> sites = sim.sites() # doctest: +SKIP
        """
        unique_sites = []
        for sys in self.spin_systems:
            for site in sys.sites:
                if site not in unique_sites:
                    unique_sites.append(site)

        return Sites(unique_sites)

    def _as_csdm_object(self, data: np.ndarray, method: Method) -> cp.CSDM:
        """
        Converts the simulation data from the given method to a CSDM object. Read
        `csdmpy <https://csdmpy.readthedocs.io/en/stable/>`_ for details

        Return:
            A CSDM object.
        """
        new = cp.new()
        for dimension in method.spectral_dimensions[::-1]:
            new.add_dimension(dimension.to_csdm_dimension())
            if new.x[-1].origin_offset != 0:
                new.x[-1].to("ppm", "nmr_frequency_ratio")

        dependent_variable = {
            "type": "internal",
            "quantity_type": "scalar",
            "numeric_type": "float64",
        }
        for index, datum in enumerate(data):
            if len(datum) == 0:
                continue

            dependent_variable["components"] = [datum]
            if self.config.decompose_spectrum == "spin_system":
                self._update_name_description_application(dependent_variable, index)

            new.add_dependent_variable(dependent_variable)
            new.y[-1].encoding = "base64"
        return new

    def _update_name_description_application(self, obj, index):
        """Update the name and description of the dependent variable attributes
        using fields from the spin system."""
        if self.spin_systems == []:
            return
        label = self.spin_systems[index].label
        if label not in ["", None]:
            obj.update({"components_label": [label]})

        name = self.spin_systems[index].name
        name = name if name not in ["", None] else f"spin system {index}"
        obj.update({"name": name})

        description = self.spin_systems[index].description
        if description not in ["", None]:
            obj.update({"description": description})

        obj["application"] = {
            "com.github.DeepanshS.mrsimulator": {
                "spin_systems": [self.spin_systems[index].json()]
            }
        }


class Sites(AbstractList):
    """A list of unique :ref:`site_api` objects within a simulator object."""

    def __init__(self, data=[]):
        super().__init__(data)
        euler = ["alpha", "beta", "gamma"]
        self.site_labels = [
            "name",
            "label",
            "description",
            "isotope",
            "isotropic_chemical_shift",
            *[f"shielding_symmetric.{_}" for _ in ["zeta", "eta", *euler]],
            *[f"quadrupolar.{_}" for _ in ["Cq", "eta", *euler]],
        ]

    def __setitem__(self, index, item):
        """Set an item to the list at index"""
        if isinstance(item, Site):
            self._list[index] = item
        elif isinstance(item, dict):
            self._list[index] = Site(**item)
        else:
            raise ValueError("Only object of type Site is allowed.")

    def to_pd(self):
        """Return sites as a pandas dataframe."""
        row = {item: [] for item in self.site_labels}
        sites = [item.json() for item in self._list]
        for site in sites:
            site_ = flatten_dict(site)
            keys = site_.keys()
            for item in self.site_labels:
                val = None if item not in keys else site_[item]
                row[item].append(val)

        row_len = len(row["name"])
        nones = [None] * row_len
        for item in self.site_labels:
            if row[item] == nones:
                row.pop(item)

        return pd.DataFrame(row)


def get_chunks(items_list, n_jobs):
    """Return the chucks of into list into roughly n_jobs equal chunks

    Args:
        (list) items_list: The input list to divide into n_jobs chunks.
        (int) n_jobs: Number of chunks of input list.
    """
    if n_jobs < 0:
        n_jobs += __CPU_count__ + 1
    list_len = len(items_list)
    n_blocks, n_left = list_len // n_jobs, list_len % n_jobs

    chunks = [0] + [n_blocks] * n_jobs
    for i in range(n_left):
        chunks[i + 1] += 1

    for i in range(1, n_jobs + 1):
        chunks[i] += chunks[i - 1]
    slices = [slice(chunks[i], chunks[i + 1], None) for i in range(n_jobs)]
    return [items_list[item] for item in slices]
