# -*- coding: utf-8 -*-
"""Base Tensor class."""
from typing import ClassVar
from typing import Dict
from typing import Optional

from pydantic import Field

from .parseable import Parseable

__author__ = "Deepansh J. Srivastava"
__email__ = "deepansh2012@gmail.com"


class SymmetricTensor(Parseable):
    """
    Base SymmetricTensor class representing the traceless symmetric part of an
    irreducible second-rank tensor.

    Attributes:
        zeta: The anisotropy parameter of the nuclear shielding tensor expressed
            using the Haeberlen convention. The default value is None.
        Cq: The quadrupolar coupling constant derived from the electric field
            gradient tensor. The default value is None.
        eta: The asymmetry parameter of the SymmetricTensor expressed using
            the Haeberlen convention. The default value is None.
        alpha: Euler angle, alpha, given in radian. The default value is None.
        beta: Euler angle, beta, given in radian. The default value is None.
        gamma: Euler angle, gamma, given in radian. The default value is None.
    """

    zeta: Optional[float]
    Cq: Optional[float]
    eta: Optional[float] = Field(default=None, ge=0, le=1)
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

    def to_freq_dict(self, larmor_frequency: float) -> dict:
        """
        Serialize the SymmetricTensor object to a JSON compliant python dictionary
        where the attribute values are numbers expressed in default units. The default
        unit for attributes with respective dimensionalities are:
        - frequency: `Hz`
        - angle: `rad`

        Args:
            float larmor_frequency: The larmor frequency in MHz.

        Return:
            A python dict
        """
        temp_dict = self.dict()
        if temp_dict["zeta"] is not None:
            temp_dict["zeta"] *= larmor_frequency
        temp_dict.pop("property_units")
        return temp_dict


class AntisymmetricTensor(Parseable):
    """
    Base SymmetricTensor class representing the traceless symmetric part of an
    irreducible second-rank tensor.

    Attributes:
        zeta: The anisotropy parameter of the AntiSymmetricTensor expressed using
            the Haeberlen convention. The default value is None.
        alpha: Euler angle, alpha, given in radian. The default value is None.
        beta: Euler angle, beta, given in radian. The default value is None.
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

    def to_freq_dict(self, larmor_frequency: float) -> dict:
        """
        Serialize the AntisymmetricTensor object to a JSON compliant python dictionary
        where the attribute values are numbers expressed in default units. The default
        unit for attributes with respective dimensionalities are:
        - frequency: `Hz`
        - angle: `rad`

        Args:
            float larmor_frequency: The larmor frequency in MHz.

        Return:
            Python dict
        """
        temp_dict = self.dict()
        if temp_dict["zeta"] is not None:
            temp_dict["zeta"] *= larmor_frequency
        temp_dict.pop("property_units")
        return temp_dict
