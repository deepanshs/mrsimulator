# -*- coding: utf-8 -*-
from typing import List, ClassVar
from mrsimulator.site import Site
from mrsimulator import Parseable

__author__ = "Deepansh J. Srivastava"
__email__ = ["srivastava.89@osu.edu", "deepansh2012@gmail.com"]


class Isotopomer(Parseable):
    """
    Base isotopmer class which lists a set of interacting NMR sites
    """

    sites: List[Site]
    # couplings: list = [], # TODO: Deepansh what should this look like?
    abundance: float = 100

    property_unit_types: ClassVar = {"abundance": "dimensionless"}

    property_default_units: ClassVar = {"abundance": "pct"}

    @classmethod
    def parse_json_with_units(cls, json_dict):

        if "sites" in json_dict:
            json_dict["sites"] = [
                Site.parse_json_with_units(s) for s in json_dict["sites"]
            ]

        return super().parse_json_with_units(json_dict)
