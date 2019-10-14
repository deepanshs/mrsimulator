# -*- coding: utf-8 -*-
from typing import ClassVar, Optional
from mrsimulator import Parseable

__author__ = "Deepansh J. Srivastava"
__email__ = ["srivastava.89@osu.edu", "deepansh2012@gmail.com"]


class SymmetricTensor(Parseable):

    zeta: Optional[float]
    Cq: Optional[float]
    eta: Optional[float]
    alpha: Optional[float]
    beta: Optional[float]
    gamma: Optional[float]

    property_unit_types: ClassVar = {
        "zeta": "dimensionless",
        "Cq": "frequency",
        "alpha": "angle",
        "beta": "angle",
        "gamma": "angle",
    }

    property_default_units: ClassVar = {
        "zeta": "ppm",
        "Cq": "Hz",
        "alpha": "rad",
        "beta": "rad",
        "gamma": "rad",
    }

    def to_freq_dict(self, larmor_frequency):
        """
        Enforces units of Hz by multiplying any ppm values by the Larmor frequency in
        MHz, MHz*ppm -> Hz
        """
        temp_dict = self.dict()
        if "zeta" in self.property_units:
            if self.property_units["zeta"] == "ppm":
                temp_dict["zeta"] *= larmor_frequency

        return temp_dict

    def to_dict_with_units(self):
        """
        Serialize the SymmetricTensor object to a JSON compliant python dictionary
        with units.
        """
        temp_dict = {k: v for k, v in self.dict().items() if v is not None}

        if "zeta" in self.property_units:
            temp_dict["zeta"] = str(temp_dict["zeta"]) + " ppm"

        if "Cq" in self.property_units:
            temp_dict["Cq"] = str(temp_dict["Cq"] / 1.0e6) + " MHz"

        temp_dict.pop("property_units")
        return temp_dict


class AntisymmetricTensor(Parseable):

    zeta: Optional[float]
    alpha: Optional[float]
    beta: Optional[float]

    property_unit_types: ClassVar = {
        "zeta": "dimensionless",
        "alpha": "angle",
        "beta": "angle",
    }

    property_default_units: ClassVar = {"zeta": ["ppm"], "alpha": "rad", "beta": "rad"}

    def to_freq_dict(self, larmor_frequency):
        """
        Enforces units of Hz by multiplying any ppm values by the Larmor frequency in
        MHz, MHz*ppm -> Hz
        """
        temp_dict = self.dict()
        if self.property_units["zeta"] == "ppm":
            temp_dict["zeta"] *= larmor_frequency

        return temp_dict

    def to_dict_with_units(self, larmor_frequency):
        """
        Serialize the AntiSymmetricTensor object to a JSON compliant python dictionary
        with units.
        """
        temp_dict = {k: v for k, v in self.dict().items() if v is not None}
        if "zeta" in self.property_units:
            temp_dict["zeta"] = str(temp_dict["zeta"]) + " ppm"

        return temp_dict
