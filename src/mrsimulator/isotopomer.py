# -*- coding: utf-8 -*-
from typing import List, ClassVar, Optional
from mrsimulator.site import Site
from mrsimulator import Parseable
from copy import deepcopy

__author__ = "Deepansh J. Srivastava"
__email__ = ["srivastava.89@osu.edu", "deepansh2012@gmail.com"]


class Isotopomer(Parseable):
    """
    Base isotopmer class representing an isolated spin-system with sites and couplings.

    .. rubric:: Attributes Documentation

    Attributes:
        name: An optional name for the isotopomer. The default is an empty
                string.
        description: An optional description for the isotopomer. The default is
                an empty string.
        sites: A list of Site objects. Default value is an empty list.
        abundance: The abundance of the isotopomer in range [0, 100]. The
                default value is 100. This attribute is useful when
                multiple isotopomers are present.
    """

    name: Optional[str] = ""
    description: Optional[str] = ""
    sites: List[Site] = []
    # couplings: list = [], # TODO: Deepansh what should this look like?
    abundance: float = 100

    property_unit_types: ClassVar = {"abundance": "dimensionless"}

    property_default_units: ClassVar = {"abundance": "pct"}

    @classmethod
    def parse_dict_with_units(cls, py_dict):
        """
        Parse the physical quantities of an isotopomer object when expressed
        as a python dictionary.

        Args:
            py_dict: Python dictionary representation of an isotopomers with
                        physical quantities.
        """
        py_dict_copy = deepcopy(py_dict)
        if "sites" in py_dict_copy:
            py_dict_copy["sites"] = [
                Site.parse_dict_with_units(s) for s in py_dict_copy["sites"]
            ]

        return super().parse_dict_with_units(py_dict_copy)

    def to_freq_dict(self, larmor_frequency):
        """
        Enforces units of Hz by multiplying any ppm values by the Larmor frequency in
        MHz, MHz*ppm -> Hz
        """
        temp_dict = self.dict()

        temp_dict["sites"] = [
            site.to_freq_dict(larmor_frequency) for site in self.sites
        ]

        return temp_dict

    def to_dict_with_units(self):
        """
        Serialize the Isotopomer object to a JSON compliant python dictionary.
        """
        temp_dict = self.dict()
        temp_dict["sites"] = [site.to_dict_with_units() for site in self.sites]
        temp_dict.pop("property_units")

        return temp_dict
