# -*- coding: utf-8 -*-
"""Base Tensor class."""
from typing import ClassVar
from typing import Dict
from typing import Optional

from mrsimulator.utils.parseable import Parseable
from pydantic import Field

__author__ = "Deepansh Srivastava"
__email__ = "srivastava.89@osu.edu"


class SymmetricTensor(Parseable):
    r"""
    Base SymmetricTensor class representing the traceless symmetric part of an
    irreducible second-rank tensor.

    Attributes
    ----------

    zeta: float (optional).
        The anisotropy parameter of the nuclear shielding tensor, in ppm, expressed
        using the Haeberlen convention. The default value is None.

        Example
        -------

        >>> shielding = SymmetricTensor()
        >>> shielding.zeta = 10

    Cq: float (optional).
        The quadrupolar coupling constant, in Hz, derived from the electric field
        gradient tensor. The default value is None.

        Example
        -------

        >>> efg = SymmetricTensor()
        >>> efg.Cq = 10e6

    eta: float (optional).
        The asymmetry parameter of the SymmetricTensor expressed using the Haeberlen
        convention. The default value is None.

        Example
        -------

        >>> shielding.eta = 0.1
        >>> efg.eta = 0.5

    alpha: float (optional).
        Euler angle, :math:`\alpha`, in radians. The default value is None.

        Example
        -------

        >>> shielding.alpha = 0.15
        >>> efg.alpha = 1.5

    beta: float (optional).
        Euler angle, :math:`\beta`, in radians. The default value is None.

        Example
        -------

        >>> shielding.beta = 3.1415
        >>> efg.beta = 1.1451

    gamma: float (optional).
        Euler angle, :math:`\gamma`, in radians. The default value is None.

        Example
        -------

        >>> shielding.gamma = 2.1
        >>> efg.gamma = 0

    Example
    -------

    >>> shielding = SymmetricTensor(zeta=10, eta=0.1, alpha=0.15, beta=3.14, gamma=2.1)
    >>> efg = SymmetricTensor(Cq=10e6, eta=0.5, alpha=1.5, beta=1.1451, gamma=0)
    """

    zeta: float = None
    Cq: float = None
    D: float = None
    eta: float = Field(default=None, ge=0.0, le=1.0)
    alpha: float = None
    beta: float = None
    gamma: float = None

    property_unit_types: ClassVar = {
        "zeta": ["dimensionless", "frequency"],
        "Cq": "frequency",
        "D": "frequency",
        "alpha": "angle",
        "beta": "angle",
        "gamma": "angle",
    }
    property_default_units: ClassVar = {
        "zeta": ["ppm", "Hz"],
        "Cq": "Hz",
        "D": "Hz",
        "alpha": "rad",
        "beta": "rad",
        "gamma": "rad",
    }
    property_units: Dict = {
        "zeta": "ppm",
        "Cq": "Hz",
        "D": "Hz",
        "alpha": "rad",
        "beta": "rad",
        "gamma": "rad",
    }

    # Deprecated
    # def to_freq_dict(self, larmor_frequency: float) -> dict:
    #     """
    #     Serialize the SymmetricTensor object to a JSON compliant python dictionary
    #     where the attribute values are numbers expressed in default units. The default
    #     unit for attributes with respective dimensionalities are:
    #     - frequency: `Hz`
    #     - angle: `rad`

    #     Args:
    #         float larmor_frequency: The larmor frequency in MHz.

    #     Return:
    #         A python dict
    #     """
    #     temp_dict = self.dict()
    #     if temp_dict["zeta"] is not None:
    #         temp_dict["zeta"] *= larmor_frequency
    #     temp_dict.pop("property_units")
    #     return temp_dict


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

    # Deprecated
    # def to_freq_dict(self, larmor_frequency: float) -> dict:
    #     """
    #     Serialize the AntisymmetricTensor object to a JSON compliant python dictionary
    #     where the attribute values are numbers expressed in default units. The default
    #     unit for attributes with respective dimensionalities are:
    #     - frequency: `Hz`
    #     - angle: `rad`

    #     Args:
    #         float larmor_frequency: The larmor frequency in MHz.

    #     Return:
    #         Python dict
    #     """
    #     temp_dict = self.dict()
    #     if temp_dict["zeta"] is not None:
    #         temp_dict["zeta"] *= larmor_frequency
    #     temp_dict.pop("property_units")
    #     return temp_dict
