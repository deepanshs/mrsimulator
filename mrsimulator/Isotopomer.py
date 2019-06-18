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
    abundance: float

    property_unit_types: ClassVar = {"abundance": "dimensionless"}

    property_default_units: ClassVar = {"abundance": "ppm"}
