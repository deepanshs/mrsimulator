from typing import ClassVar
from mrsimulator import Parseable

__author__ = "Deepansh J. Srivastava"
__email__ = ["srivastava.89@osu.edu", "deepansh2012@gmail.com"]


class Site(Parseable):

    nucleus: str = "1H"
    isotropic_chemical_shift: float = 0
    anisotropy: float = 0
    asymmetry: float = 0
    alpha: float = 0
    beta: float = 0
    gamma: float = 0

    property_unit_types: ClassVar = {
        "isotropic_chemical_shift": "dimensionless",
        "anisotropy": "dimensionless",
        "alpha": "angle",
        "beta": "angle",
        "gamma": "angle",
    }

    property_default_units: ClassVar = {
        "isotropic_chemical_shift": "ppm",
        "anisotropy": "ppm",
        "alpha": "rad",
        "beta": "rad",
        "gamma": "rad",
    }

    @classmethod
    def parse_json_with_units(cls, json_dict):

        if "shielding_symmetric" in json_dict:
            for k, v in json_dict["shielding_symmetric"].items():
                json_dict[k] = v
            del json_dict["shielding_symmetric"]

        if "orientation" in json_dict:
            for k, v in json_dict["shielding_symmetric"].items():
                json_dict[k] = v
            del json_dict["orientation"]

        return super().parse_json_with_units(json_dict)
