# -*- coding: utf-8 -*-
"""Base Simulator class."""
import json
from copy import deepcopy
from typing import List

import csdmpy as cp
import numpy as np
from mrsimulator import __version__
from mrsimulator import SpinSystem
from mrsimulator.base_model import one_d_spectrum
from mrsimulator.method import Method
from mrsimulator.utils.extra import _reduce_dict
from mrsimulator.utils.importer import import_json
from pydantic import BaseModel

from .config import ConfigSimulator

__author__ = "Deepansh J. Srivastava"
__email__ = "deepansh2012@gmail.com"


class Simulator(BaseModel):
    """
    The simulator class.

    Attributes
    ----------

    spin_systems: A list of :ref:`spin_sys_api` or equivalent dict objects (optional).
        The value is a list of NMR spin systems present within the sample, where each
        spin system is an isolated system. The default value is an empty list.

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

    methods: A list of :ref:`method_api` (optional).
        The value is a list of NMR methods. The default value is an empty list.

        Example
        -------

        >>> from mrsimulator.methods import BlochDecaySpectrum
        >>> from mrsimulator.methods import BlochDecayCentralTransitionSpectrum
        >>> sim.methods = [
        ...     BlochDecaySpectrum(channels=['17O'], spectral_width=50000),
        ...     BlochDecayCentralTransitionSpectrum(channels=['17O'])
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
        The value is the name or id of the simulation or sample. The default value is
        None.

        Example
        -------

        >>> sim.name = '1H-17O'
        >>> sim.name
        '1H-17O'

    label: str (optional).
        The value is a label for the simulation or sample. The default value is None.

        Example
        -------

        >>> sim.label = 'Test simulator'
        >>> sim.label
        'Test simulator'

    description: str (optional).
        The value is a description of the simulation or sample. The default value is
        None.

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

    @classmethod
    def parse_dict_with_units(cls, py_dict):
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
            methods = [Method.parse_dict_with_units(obj) for obj in methods]
            py_copy_dict["methods"] = methods

        return Simulator(**py_copy_dict)

    def get_isotopes(self, spin_I=None) -> set:
        """
        Set of unique isotopes from the sites within the list of the spin systems
        corresponding to spin quantum number `I`. If `I` is None, a set of all unique
        isotopes is returned instead.

        Args:
            float spin_I: An optional spin quantum number. The valid input are the
                multiples of 0.5.

        Returns:
            A Set.

        Example
        -------

        >>> sim.get_isotopes() # doctest:+SKIP
        {'1H', '27Al', '13C'}
        >>> sim.get_isotopes(spin_I=0.5) # doctest:+SKIP
        {'1H', '13C'}
        >>> sim.get_isotopes(spin_I=1.5)
        set()
        >>> sim.get_isotopes(spin_I=2.5)
        {'27Al'}
        """
        st = set()
        for spin_system in self.spin_systems:
            st.update(spin_system.get_isotopes(spin_I))
        return st

    def to_dict_with_units(
        self, include_methods: bool = False, include_version: bool = False
    ):
        """
        Serialize the Simulator object to a JSON compliant python dictionary object
        where physical quantities are represented as string with a value and a unit.

        Args:
            bool include_methods: If True, the output dictionary will include the
                serialized method objects. The default value is False.
            bool include_version: If True, add a version key-value pair to the
                serialized output dictionary. The default is False.

        Returns:
            A Dict object.

        Example
        -------

        >>> pprint(sim.to_dict_with_units())
        {'config': {'decompose_spectrum': 'none',
                    'integration_density': 70,
                    'integration_volume': 'octant',
                    'number_of_sidebands': 64},
         'indexes': [],
         'spin_systems': [{'abundance': '100 %',
                           'sites': [{'isotope': '13C',
                                      'isotropic_chemical_shift': '20.0 ppm',
                                      'shielding_symmetric': {'eta': 0.5,
                                                              'zeta': '10.0 ppm'}}]},
                          {'abundance': '100 %',
                           'sites': [{'isotope': '1H',
                                      'isotropic_chemical_shift': '-4.0 ppm',
                                      'shielding_symmetric': {'eta': 0.1,
                                                              'zeta': '2.1 ppm'}}]},
                          {'abundance': '100 %',
                           'sites': [{'isotope': '27Al',
                                      'isotropic_chemical_shift': '120.0 ppm',
                                      'shielding_symmetric': {'eta': 0.1,
                                                              'zeta': '2.1 ppm'}}]}]}
        """
        sim = {}

        if self.name is not None:
            sim["name"] = self.name

        if self.description is not None:
            sim["description"] = self.description

        if self.label is not None:
            sim["label"] = self.label

        sim["spin_systems"] = [_.to_dict_with_units() for _ in self.spin_systems]

        if include_methods:
            method = [_.to_dict_with_units() for _ in self.methods]
            if len(method) != 0:
                sim["methods"] = method

        sim["config"] = self.config.dict()
        sim["indexes"] = self.indexes
        if include_version:
            sim["version"] = __version__
        return sim

    def reduced_dict(self, exclude=["property_units"]) -> dict:
        """Returns a reduced dictionary representation of the class object by removing
        all key-value pair corresponding to keys listed in the `exclude` argument, and
        keys with value as None.

        Args:
            list exclude: A list of keys to exclude from the dictionary.
        Return: A dict.
         """
        return _reduce_dict(self.dict(), exclude)

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
        spin_systems = [SpinSystem.to_dict_with_units(obj) for obj in self.spin_systems]
        with open(filename, "w", encoding="utf8") as outfile:
            json.dump(
                spin_systems,
                outfile,
                ensure_ascii=False,
                sort_keys=False,
                allow_nan=False,
            )

    def run(self, method_index=None, pack_as_csdm=True, **kwargs):
        """Run the simulation and compute spectrum.

        Args:
            method_index: An integer or a list of integers. If provided, only the
                simulations corresponding to the methods at the given index/indexes
                will be computed. The default is None, `i.e.`, the simulation for
                every method will be computed.
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

        if method_index is None:
            method_index = np.arange(len(self.methods))
        if isinstance(method_index, int):
            method_index = [method_index]
        for index in method_index:
            method = self.methods[index]
            amp, indexes = one_d_spectrum(
                method=method,
                spin_systems=self.spin_systems,
                **self.config.get_int_dict(),
                **kwargs,
            )

            self.indexes.append(indexes)

            if isinstance(amp, list):
                simulated_data = amp
            else:
                simulated_data = [amp]

            if pack_as_csdm:
                method.simulation = self._as_csdm_object(simulated_data, method)
            else:
                method.simulation = np.asarray(simulated_data)

    # """The frequency is in the units of Hz."""
    # gamma = method.isotope.gyromagnetic_ratio
    # B0 = method.spectral_dimensions[0].events[0].magnetic_flux_density
    # larmor_frequency = -gamma * B0
    # reference_offset_in_MHz =method.spectral_dimensions[0].reference_offset / 1e6
    # denom = reference_offset_in_MHz + larmor_frequency
    # freq = method.spectral_dimensions[0].coordinates_Hz / abs(denom)

    # freq *= u.Unit("ppm")
    # return freq, amp

    def save(self, filename: str, with_units=True):
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
        if not with_units:
            with open(filename, "w", encoding="utf8") as outfile:
                json.dump(
                    self.reduced_dict(),
                    outfile,
                    ensure_ascii=False,
                    sort_keys=False,
                    allow_nan=False,
                )
            return

        with open(filename, "w", encoding="utf8") as outfile:
            json.dump(
                self.to_dict_with_units(include_methods=True, include_version=True),
                outfile,
                ensure_ascii=False,
                sort_keys=False,
                allow_nan=False,
            )

    @classmethod
    def load(cls, filename: str, parse_units=True):
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
        contents = import_json(filename)

        if not parse_units:
            return Simulator(**contents)

        return Simulator.parse_dict_with_units(contents)

    def _as_csdm_object(self, data: np.ndarray, method: Method) -> cp.CSDM:
        """
        Converts the simulation data from the given method to a CSDM object. Read
        `csdmpy <https://csdmpy.readthedocs.io/en/latest/>`_ for details

        Return:
            A CSDM object.
        """
        new = cp.new()
        for dimension in method.spectral_dimensions:
            new.add_dimension(dimension.to_csdm_dimension())
            if new.dimensions[-1].origin_offset != 0:
                new.dimensions[-1].to("ppm", "nmr_frequency_ratio")

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
            new.dependent_variables[-1].encoding = "base64"
        return new

    def _update_name_description_application(self, obj, index):
        """Update the name and description of the dependent variable attributes
        using fields from the spin system."""
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
                "spin_systems": [self.spin_systems[index].to_dict_with_units()]
            }
        }
