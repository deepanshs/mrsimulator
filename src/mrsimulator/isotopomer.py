# -*- coding: utf-8 -*-
from copy import deepcopy
from typing import ClassVar
from typing import Dict
from typing import List
from typing import Optional

from mrsimulator import Parseable
from mrsimulator.site import Site
from pydantic import Field

__author__ = "Deepansh J. Srivastava"
__email__ = ["srivastava.89@osu.edu", "deepansh2012@gmail.com"]


class Isotopomer(Parseable):
    """
    Base isotopmer class representing an isolated spin-system with sites and couplings.

    Arguments:
        name: An optional name for the isotopomer. The default is an empty
                string.
        description: An optional description for the isotopomer. The default is
                an empty string.
        sites: A list of Site objects or an equivalent python dict object representing
                a nuclear site. Default value is an empty list.
        abundance: The abundance of the isotopomer in unit of %. The default value is
                100. This attribute is useful when multiple isotopomers are present.
    """

    name: Optional[str] = ""
    description: Optional[str] = ""
    sites: List[Site] = []
    # couplings: list = [], # TODO: Deepansh what should this look like?
    abundance: float = Field(default=100, ge=0, le=100)

    property_unit_types: ClassVar = {"abundance": "dimensionless"}
    property_default_units: ClassVar = {"abundance": "pct"}
    property_units: Dict = {"abundance": "pct"}

    @classmethod
    def parse_dict_with_units(cls, py_dict):
        """
        Parse physical quantity from the attributes of an isotopomer object
        expressed as a python dictionary. The physical quantities are expressed
        as string with a number followed by a unit.

        Args:
            py_dict: Python dictionary representation of an isotopomer object with
                    attributes values as string with a physical quantity.
        """
        py_dict_copy = deepcopy(py_dict)
        if "sites" in py_dict_copy:
            py_dict_copy["sites"] = [
                Site.parse_dict_with_units(s) for s in py_dict_copy["sites"]
            ]

        return super().parse_dict_with_units(py_dict_copy)

    def to_freq_dict(self, B0):
        """
        Serialize the Isotopomer object to a JSON compliant python dictionary where the
        attribute values are numbers expressed in default units. The default unit
        for attributes with respective dimensionalities are:
        - frequency: `Hz`
        - angle: `rad`

        Args:
            B0: The macroscopic magnetic flux density

        Return: A python dict
        """
        temp_dict = self.dict()
        temp_dict["sites"] = [site.to_freq_dict(B0) for site in self.sites]
        temp_dict.pop("property_units")
        return temp_dict

    def to_dict_with_units(self):
        """
        Serialize the Isotopomer object to a JSON compliant python dictionary object
        where attribute values are physical quantities expressed as a string with a
        number followed by a unit.

        Return: A python dict
        """
        temp_dict = self.dict()
        temp_dict["sites"] = [site.to_dict_with_units() for site in self.sites]
        temp_dict["abundance"] = f"{self.abundance}%"
        temp_dict.pop("property_units")

        return temp_dict
