# -*- coding: utf-8 -*-
from typing import List, ClassVar
from mrsimulator.site import Site
from mrsimulator import Parseable

__author__ = "Deepansh J. Srivastava"
__email__ = ["srivastava.89@osu.edu", "deepansh2012@gmail.com"]


class Isotopomer(Parseable):
    """
    Base isotopmer class which lists a set of interacting NMR sites.

    .. rubric:: Attributes Documentation

    Attributes:
        sites: A list of Site objects.
        abundance: The fractional abundance of the isotopomer. This attribute
                is useful when multiple isotopomers are present.
    """

    sites: List[Site]
    # couplings: list = [], # TODO: Deepansh what should this look like?
    abundance: float = 100

    property_unit_types: ClassVar = {"abundance": "dimensionless"}

    property_default_units: ClassVar = {"abundance": "pct"}

    @classmethod
    def parse_json_with_units(cls, json_dict):
        """
        Parse the physical quantities of an isotopomer when expressed as a python
        dictionary.

        Args:
            json_dict: Python dictionary representation of an isotopomers containing
                        physical quantities.
        """
        if "sites" in json_dict:
            json_dict["sites"] = [
                Site.parse_json_with_units(s) for s in json_dict["sites"]
            ]

        return super().parse_json_with_units(json_dict)

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
