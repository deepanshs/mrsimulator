from typing import ClassVar, Optional
from mrsimulator import Parseable

__author__ = "Deepansh J. Srivastava"
__email__ = ["srivastava.89@osu.edu", "deepansh2012@gmail.com"]


class SymmetricTensor(Parseable):

    anisotropy: Optional[float]
    asymmetry: Optional[float]
    alpha: Optional[float]
    beta: Optional[float]
    gamma: Optional[float]

    property_unit_types: ClassVar = {
        "anisotropy": ["dimensionless", "frequency"],
        "alpha": "angle",
        "beta": "angle",
        "gamma": "angle",
    }

    property_default_units: ClassVar = {
        "anisotropy": ["ppm", "Hz"],
        "alpha": "rad",
        "beta": "rad",
        "gamma": "rad",
    }


class AntisymmetricTensor(Parseable):

    anisotropy: Optional[float]
    alpha: Optional[float]
    beta: Optional[float]

    property_unit_types: ClassVar = {
        "anisotropy": ["dimensionless", "frequency"],
        "alpha": "angle",
        "beta": "angle",
    }

    property_default_units: ClassVar = {
        "anisotropy": ["ppm", "Hz"],
        "alpha": "rad",
        "beta": "rad",
    }
