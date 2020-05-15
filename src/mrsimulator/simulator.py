# -*- coding: utf-8 -*-
"""Base Simulator class."""
import json
from typing import List
from typing import Optional

import csdmpy as cp
import numpy as np
from mrsimulator import Isotopomer
from mrsimulator.apodization import Apodization
from mrsimulator.base_model import one_d_spectrum
from mrsimulator.importer import import_json
from mrsimulator.method import Method
from mrsimulator.simulator_config import ConfigSimulator
from pydantic import BaseModel

# from astropy import units as u

__author__ = "Deepansh J. Srivastava"
__email__ = "deepansh2012@gmail.com"


class Simulator(BaseModel):
    """
    The simulator class.

    Attributes:
        name: An optional string containing the sample name.
        descrition: An optional string containing the sample description.
        isotopomers: List of :ref:`isotopomer_api` objects.
        methods: List of :ref:`method_api` objects.
        config: :ref:`config_api` object.
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

    def get_isotopes(self, I=None):
        """
        Set of unique isotopes from sites in the list of isotopomers corresponding
        to spin quantum number `I`. If `I` is unspecified or None, a set of all
        unique isotopes is returned instead.

        Args:
            I: (optional) The spin quantum number. Valid input are multiples of 0.5.

        Returns:
            A Set

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

    def to_dict_with_units(self, include_methods=False):
        """
        Serialize the Simulator object to a JSON compliant python dictionary object
        with units.

        Args:
            include_methods: A boolean. If True, the output dictionary will include
            a serialized method objects. Default is False.

        Returns:
            Dict object

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
        return sim

    def load_isotopomers(self, filename):
        """
        Load a list of isotopomers from JSON serialized isotopomers file.

        See an
        `example <https://raw.githubusercontent.com/DeepanshS/mrsimulator-test
        /master/isotopomers_ppm.json>`_
        of JSON serialized isotopomers file. For details, refer to the
        :ref:`load_isotopomers` section.

        Args:
            `filename`: A local or remote address to the JSON serialized isotopomers
                        file.
        Example:
            >>> sim.load_isotopomers(filename) # doctest:+SKIP
        """
        contents = import_json(filename)
        json_data = contents["isotopomers"]
        self.isotopomers = [Isotopomer.parse_dict_with_units(obj) for obj in json_data]

    def export_isotopomers(self, filename):
        """
        Export the list of isotopomers to a JSON serialized isotopomers file.

        See an
        `example <https://raw.githubusercontent.com/DeepanshS/mrsimulator-test
        /master/isotopomers_ppm.json>`_
        of JSON serialized isotopomers file. For details, refer to the
        :ref:`load_isotopomers` section.

        Args:
            `filename`: A file name.
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
        """Simulate the lineshape.

        Args:
            method_index: If provided, only update the simulate for the method at
            the given index.

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

            method.simulation = self.as_csdm_object(simulated_data, method)

    # """The frequency is in the units of Hz."""
    # gamma = method.isotope.gyromagnetic_ratio
    # B0 = method.spectral_dimensions[0].events[0].magnetic_flux_density
    # larmor_frequency = -gamma * B0
    # reference_offset_in_MHz =method.spectral_dimensions[0].reference_offset / 1e6
    # denom = reference_offset_in_MHz + larmor_frequency
    # freq = method.spectral_dimensions[0].coordinates_Hz / abs(denom)

    # freq *= u.Unit("ppm")
    # return freq, amp

    def save(self, filename):
        """Serialize the simulator object to a JSON file.

        Args:
            filename: The file name used in serialization.
        """
        with open(filename, "w", encoding="utf8") as outfile:
            json.dump(
                self.to_dict_with_units(include_methods=True),
                outfile,
                ensure_ascii=False,
                sort_keys=False,
                allow_nan=False,
            )

    def load(self, filename):
        """Load the mrsimulator object from the JSON file.

        Args:
            filename: The name of the file holding a mrsimulator serialization.

        Return:
            A Simulator object.
        """
        sim = Simulator()
        contents = import_json(filename)
        i_data = contents["isotopomers"]
        sim.isotopomers = [Isotopomer.parse_dict_with_units(obj) for obj in i_data]

        m_data = contents["methods"]
        sim.methods = [Method.parse_dict_with_units(obj) for obj in m_data]
        return sim

    def as_csdm_object(self, data, method):
        """
        Converts the data to a CSDM object. Read
        `csdmpy <https://csdmpy.readthedocs.io/en/latest/>`_ for details

        Return:
            CSDM object
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
