# -*- coding: utf-8 -*-
from typings import List
from mrsimulator.site import Site

__author__ = "Deepansh J. Srivastava"
__email__ = ["srivastava.89@osu.edu", "deepansh2012@gmail.com"]


class Isotopomer(BaseModel):
    """
    Base isotopmer class which lists a set of interacting NMR sites
    """

    sites: List(Site)
    # couplings: list = [], # TODO: Deepansh what should this look like?
    abundance: float

    __property_unit_types = {"abundance": "dimensionless"}

    __property_default_units = {"abundance": "ppm"}
