# -*- coding: utf-8 -*-
"""Base Simulator class."""
from typing import List
from typing import Optional

import csdmpy as cp
import numpy as np
from astropy import units as u
from mrsimulator import Dimension
from mrsimulator import Isotopomer
from mrsimulator.apodization import Apodization
from mrsimulator.base_model import one_d_spectrum
from mrsimulator.importer import import_json
from mrsimulator.method import Method
from mrsimulator.simulator_config import ConfigSimulator
from pydantic import BaseModel

__author__ = "Deepansh J. Srivastava"
__email__ = "deepansh2012@gmail.com"


class Simulator(BaseModel):
    """
    The simulator class.

    Attributes:
        isotopomers: List of :ref:`isotopomer_api` objects.
        dimensions: List of :ref:`dimension_api` objects.
        config: :ref:`config_api` object.
    """

    isotopomers: List[Isotopomer] = []
    dimensions: List[Dimension] = []
    method: Optional[Method] = None
    simulated_data: Optional[List]
    config: ConfigSimulator = ConfigSimulator()

    class Config:
        validate_assignment = True
        arbitrary_types_allowed = True

    def __eq__(self, other):
        check = [
            isinstance(other, Simulator),
            self.isotopomers == other.isotopomers,
            self.method == other.method,
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

    def to_dict_with_units(self, include_dimensions=False):
        """
        Serialize the isotopomers and dimensions attribute from the Simulator object
        to a JSON compliant python dictionary with units.

        Args:
            remove_dimensions: A boolean. If True, the output contains serialized
                dimension objects. Default is False.

        Returns:
            Dict object

        Example:
            >>> pprint(sim.to_dict_with_units())
            {'isotopomers': [{'abundance': '100%',
                              'description': '',
                              'name': '',
                              'sites': [{'isotope': '13C',
                                         'isotropic_chemical_shift': '20.0 ppm',
                                         'shielding_symmetric': {'eta': 0.5,
                                                                 'zeta': '10.0 ppm'}}]},
                             {'abundance': '100%',
                              'description': '',
                              'name': '',
                              'sites': [{'isotope': '1H',
                                         'isotropic_chemical_shift': '-4.0 ppm',
                                         'shielding_symmetric': {'eta': 0.1,
                                                                 'zeta': '2.1 ppm'}}]},
                             {'abundance': '100%',
                              'description': '',
                              'name': '',
                              'sites': [{'isotope': '27Al',
                                         'isotropic_chemical_shift': '120.0 ppm',
                                         'shielding_symmetric': {'eta': 0.1,
                                                                 'zeta': '2.1 ppm'}}]}]}
        """
        sim = {}
        sim["isotopomers"] = [_.to_dict_with_units() for _ in self.isotopomers]

        if include_dimensions:
            dimensions = [_.to_dict_with_units() for _ in self.dimensions]
            if dimensions != []:
                sim["dimensions"] = dimensions

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

    def run(self, **kwargs):
        """Simulate the lineshape.

        Args:
            method: The methods used in simulating line-shapes.

        Example:
            >>> sim.run(method=one_d_spectrum) # doctest:+SKIP
        """
        # isotopomers = [
        #     isotopomer.to_freq_dict(
        #         self.method.sequences[0].events[0].magnetic_flux_density
        #     )
        #     for isotopomer in self.isotopomers
        # ]

        amp = one_d_spectrum(
            method=self.method,
            isotopomers=self.isotopomers,
            **self.config._dict,
            **kwargs,
        )

        if isinstance(amp, list):
            self.simulated_data = amp
        else:
            self.simulated_data = [amp]

        """The frequency is in the units of Hz."""
        gamma = self.method.isotope.gyromagnetic_ratio
        B0 = self.method.sequences[0].events[0].magnetic_flux_density
        larmor_frequency = -gamma * B0

        denom = self.method.sequences[0].reference_offset / 1e6 + larmor_frequency
        freq = self.method.sequences[0].coordinates_Hz / abs(denom)

        freq *= u.Unit("ppm")
        return freq, amp

    def as_csdm_object(self):
        """
        Converts the data to a CSDM object. Read
        `csdmpy <https://csdmpy.readthedocs.io/en/latest/>`_ for details

        Return:
            CSDM object
        """
        new = cp.new()
        for dimension in self.method.sequences:
            new.add_dimension(dimension.to_csdm_dimension())

        dependent_variable = {
            "type": "internal",
            "quantity_type": "scalar",
            "numeric_type": "float64",
        }
        for index, datum in enumerate(self.simulated_data):
            if datum != []:
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
        apodization_filter = Apodization(self.as_csdm_object(), dimension=dimension)
        return apodization_filter.apodize(fn, **kwargs)
