# -*- coding: utf-8 -*-
"""Base Simulator class."""
import json
from typing import List
from typing import Optional

import csdmpy as cp
import numpy as np
from mrsimulator import __version__
from mrsimulator import Isotopomer
from mrsimulator.apodization import Apodization
from mrsimulator.base_model import one_d_spectrum
from mrsimulator.importer import import_json
from mrsimulator.method import Method
from mrsimulator.simulator_config import ConfigSimulator
from pydantic import BaseModel

from .util import _reduce_dict

# from astropy import units as u

__author__ = "Deepansh J. Srivastava"
__email__ = "deepansh2012@gmail.com"


class Simulator(BaseModel):
    """
    The simulator class.

    Attributes:
        name: An optional string containing the simulation/sample name. The default
            value is an empty string.
        description: An optional string with the simulation/sample description. The
            default value is an empty string.
        isotopomers: A list of the :ref:`isotopomer_api` objects or list of equivalent
            python dictionary object. The default value is an empty list.
        methods: A list of :ref:`method_api` objects or list of equivalent python
            dictionary object. The default value is an empty list.
        config: The :ref:`config_api` object or an equivalent dictionary
            object.
    """

    name: Optional[str] = ""
    description: Optional[str] = ""
    isotopomers: List[Isotopomer] = []
    methods: List[Method] = []
    config: ConfigSimulator = ConfigSimulator()
    indexes = []

    class Config:
        validate_assignment = True
        arbitrary_types_allowed = True

    def __eq__(self, other):
        check = [
            isinstance(other, Simulator),
            self.name == other.name,
            self.description == other.description,
            self.isotopomers == other.isotopomers,
            self.methods == other.methods,
            self.config == other.config,
        ]
        if np.all(check):
            return True
        return False

    def get_isotopes(self, I=None) -> set:
        """
        Set of unique isotopes from sites within the list of isotopomers corresponding
        to spin quantum number `I`. If `I` is unspecified or None, a set of all
        unique isotopes is returned instead.

        Args:
            float I: An optional spin quantum number. The valid input are the multiples
                of 0.5.

        Returns:
            A Set.

        Example:
            >>> sim.get_isotopes() # doctest:+SKIP
            {'1H', '27Al', '13C'}
            >>> sim.get_isotopes(I=0.5) # doctest:+SKIP
            {'1H', '13C'}
            >>> sim.get_isotopes(I=1.5)
            set()
            >>> sim.get_isotopes(I=2.5)
            {'27Al'}
        """
        st = set()
        for isotopomer in self.isotopomers:
            st.update(isotopomer.get_isotopes(I))
        return st

    def dict(self, *args, **kwargs) -> dict:
        """Return the class object as JSON serializable dictionary object."""
        py_dict = super().dict(*args, **kwargs)
        py_dict["config"] = self.config.dict()
        return py_dict

    def to_dict_with_units(
        self, include_methods: bool = False, include_version: bool = False
    ):
        """
        Serialize the Simulator object to a JSON compliant python dictionary object
        where physical quantities are represented as string with a value and a unit.

        Args:
            bool include_methods: If True, the output dictionary will include the
                serialized method objects. The default value is False.
            bool include_version: If True, adds the version key-value pair
                to the serialized output dictionary. The default is False.

        Returns:
            A Dict object.

        Example:
            >>> pprint(sim.to_dict_with_units())
            {'config': {'decompose': False,
                        'integration_density': 70,
                        'integration_volume': 'octant',
                        'number_of_sidebands': 64},
             'description': '',
             'indexes': [],
             'isotopomers': [{'abundance': '100 %',
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
                                                                 'zeta': '2.1 ppm'}}]}],
             'name': ''}
        """
        sim = {}
        sim["name"] = self.name
        sim["description"] = self.description
        sim["isotopomers"] = [_.to_dict_with_units() for _ in self.isotopomers]

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

    def load_isotopomers(self, filename: str):
        """
        Load a list of isotopomers from the given JSON serialized isotopomers file.

        See an
        `example <https://raw.githubusercontent.com/DeepanshS/mrsimulator-test
        /master/isotopomers_ppm.json>`_
        of JSON serialized isotopomers file. For details, refer to the
        :ref:`load_isotopomers` section of this documentation.

        Args:
            str filename: A local or remote address to a JSON serialized isotopomers
                file.

        Example:
            >>> sim.load_isotopomers(filename) # doctest:+SKIP
        """
        contents = import_json(filename)
        json_data = contents["isotopomers"]
        self.isotopomers = [Isotopomer.parse_dict_with_units(obj) for obj in json_data]

    def export_isotopomers(self, filename: str):
        """
        Export a list of isotopomers to a JSON serialized isotopomers file.

        See an
        `example <https://raw.githubusercontent.com/DeepanshS/mrsimulator-test
        /master/isotopomers_ppm.json>`_
        of JSON serialized isotopomers file. For details, refer to the
        :ref:`load_isotopomers` section.

        Args:
            str filename: The list of isotopomers will be serialized to a file with
                the given filename.

        Example:
            >>> sim.export_isotopomers(filename) # doctest:+SKIP
        """
        isotopomers = [Isotopomer.to_dict_with_units(obj) for obj in self.isotopomers]
        with open(filename, "w", encoding="utf8") as outfile:
            json.dump(
                isotopomers,
                outfile,
                ensure_ascii=False,
                sort_keys=False,
                allow_nan=False,
            )

    def run(self, method_index=None, **kwargs):
        """Run the simulation and compute lineshape.

        Args:
            int list method_index: An interger or a list of integers. If provided, only
                the simulations corresponding to the methods at the given index/indexes
                will be computed. The default is None, that is, the simulation for
                every method will computed.

        Example:
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
                isotopomers=self.isotopomers,
                **self.config._dict,
                **kwargs,
            )

            self.indexes.append(indexes)

            if isinstance(amp, list):
                simulated_data = amp
            else:
                simulated_data = [amp]

            method.simulation = self._as_csdm_object(simulated_data, method)

    # """The frequency is in the units of Hz."""
    # gamma = method.isotope.gyromagnetic_ratio
    # B0 = method.spectral_dimensions[0].events[0].magnetic_flux_density
    # larmor_frequency = -gamma * B0
    # reference_offset_in_MHz =method.spectral_dimensions[0].reference_offset / 1e6
    # denom = reference_offset_in_MHz + larmor_frequency
    # freq = method.spectral_dimensions[0].coordinates_Hz / abs(denom)

    # freq *= u.Unit("ppm")
    # return freq, amp

    def save(self, filename: str):
        """Serialize the simulator object to a JSON file.

        Args:
            str filename: A string with the filename of the serialized file.
        """
        with open(filename, "w", encoding="utf8") as outfile:
            json.dump(
                self.to_dict_with_units(include_methods=True, include_version=True),
                outfile,
                ensure_ascii=False,
                sort_keys=False,
                allow_nan=False,
            )

    def load(self, filename: str):
        """Load the :class:`~mrsimulator.Simulator` object from the JSON file.

        Args:
            str filename: A string with the filename of the file holding a mrsimulator
                serialized file.

        Return:
            A :class:`~mrsimulator.Simulator` object.
        """
        sim = Simulator()
        contents = import_json(filename)
        i_data = contents["isotopomers"]
        sim.isotopomers = [Isotopomer.parse_dict_with_units(obj) for obj in i_data]

        m_data = contents["methods"]
        sim.methods = [Method.parse_dict_with_units(obj) for obj in m_data]
        return sim

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
            new.dimensions[-1].to("ppm", "nmr_frequency_ratio")

        dependent_variable = {
            "type": "internal",
            "quantity_type": "scalar",
            "numeric_type": "float64",
        }
        for index, datum in enumerate(data):
            if len(datum) != 0:
                dependent_variable["components"] = [datum]
                name = self.isotopomers[index].name
                if name not in ["", None]:
                    dependent_variable.update({"name": name})

                description = self.isotopomers[index].description
                if description not in ["", None]:
                    dependent_variable.update({"description": description})

                dependent_variable["application"] = {
                    "com.github.DeepanshS.mrsimulator": {
                        "isotopomers": [self.isotopomers[index].to_dict_with_units()]
                    }
                }
                new.add_dependent_variable(dependent_variable)
                new.dependent_variables[-1].encoding = "base64"
        return new

    def apodize(self, fn, dimension=0, **kwargs):
        apodization_filter = Apodization(
            self.methods[0].simulation, dimension=dimension
        )
        return apodization_filter.apodize(fn, **kwargs)
