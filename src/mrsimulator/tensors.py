# -*- coding: utf-8 -*-
"""Base Tensor class."""
from typing import ClassVar
from typing import Dict
from typing import Optional

from mrsimulator import Parseable
from pydantic import Field

__author__ = "Deepansh J. Srivastava"
__email__ = "deepansh2012@gmail.com"


class SymmetricTensor(Parseable):
    """
    Base SymmetricTensor class representing the traceless symmetric part of an
    irreducible second-rank tensor.

    Arguments:
        zeta: The anisotropy parameter of a nuclear shielding tensor expressed using
                Haeberlen convention.
        Cq: The quadrupolar coupling constant derived from an electric field tensor.
        eta: The asymmetry parameter of the SymmetricTensor expressed using
                Haeberlen convention.
        alpha: Euler angle, alpha, given in radian.
        beta: Euler angle, beta, given in radian.
        gamma: Euler angle, gamma, given in radian.
    """

    zeta: Optional[float]
    Cq: Optional[float]
    eta: Optional[float] = Field(default=0, ge=0, le=1)
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
    property_units: Dict = {
        "zeta": "ppm",
        "Cq": "Hz",
        "alpha": "rad",
        "beta": "rad",
        "gamma": "rad",
    }

    def to_freq_dict(self, larmor_frequency):
        """
        Serialize the SymmetricTensor object to a JSON compliant python dictionary
        where the attribute values are numbers expressed in default units. The default
        unit for attributes with respective dimensionalities are:
        - frequency: `Hz`
        - angle: `rad`

        Args:
            larmor_frequency: The larmor frequency in MHz.

        Return:
            A python dict
        """
        temp_dict = self.dict()
        if temp_dict["zeta"] is not None:
            temp_dict["zeta"] *= larmor_frequency
        temp_dict.pop("property_units")
        return temp_dict

    def to_dict_with_units(self):
        """
        Serialize the SymmetricTensor object to a JSON compliant python dictionary
        with units.

        Return:
            Python dict
        """
        temp_dict = {
            k: f"{v} {self.property_units[k]}"
            for k, v in self.dict().items()
            if v is not None and k not in ["property_units", "eta"]
        }
        if self.eta is not None:
            temp_dict["eta"] = self.eta
        return temp_dict


class AntisymmetricTensor(Parseable):
    """
    Base SymmetricTensor class representing the traceless symmetric part of an
    irreducible second-rank tensor.

    Arguments:
        zeta: The anisotropy parameter of the AntiSymmetricTensor expressed using
                Haeberlen convention.
        alpha: Euler angle, alpha, given in radian.
        beta: Euler angle, beta, given in radian.
    """

    zeta: Optional[float]
    alpha: Optional[float]
    beta: Optional[float]

    property_unit_types: ClassVar = {
        "zeta": "dimensionless",
        "alpha": "angle",
        "beta": "angle",
    }
    property_default_units: ClassVar = {"zeta": "ppm", "alpha": "rad", "beta": "rad"}
    property_units: Dict = {"zeta": "ppm", "alpha": "rad", "beta": "rad"}

    def to_freq_dict(self, larmor_frequency):
        """
        Serialize the AntisymmetricTensor object to a JSON compliant python dictionary
        where the attribute values are numbers expressed in default units. The default
        unit for attributes with respective dimensionalities are:
        - frequency: `Hz`
        - angle: `rad`

        Args:
            larmor_frequency: The larmor frequency in MHz.

        Return:
            Python dict
        """
        temp_dict = self.dict()
        if temp_dict["zeta"] is not None:
            temp_dict["zeta"] *= larmor_frequency
        temp_dict.pop("property_units")
        return temp_dict

    def to_dict_with_units(self):
        """
        Serialize the AntiSymmetricTensor object to a JSON compliant python dictionary
        with units.

        Return:
            Python dict
        """
        temp_dict = {
            k: f"{v} {self.property_units[k]}"
            for k, v in self.dict().items()
            if v is not None and k != "property_units"
        }
        return temp_dict
