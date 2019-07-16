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

    def to_freq_dict(self, larmor_frequency):
        """
        Enforces units of Hz by multiplying any ppm values by the Larmor frequency in MHz
        MHz*ppm -> Hz
        """
        temp_dict = self.dict()
        if self.property_units["anisotropy"] is "ppm":
            temp_dict["anisotropy"] *= larmor_frequency

        return temp_dict


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

    def to_freq_dict(self, larmor_frequency):
        """
        Enforces units of Hz by multiplying any ppm values by the Larmor frequency in MHz
        MHz*ppm -> Hz
        """
        temp_dict = self.dict()
        if self.property_units["anisotropy"] is "ppm":
            temp_dict["anisotropy"] *= larmor_frequency

        return temp_dict
