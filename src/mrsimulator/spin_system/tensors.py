"""Base Tensor class."""
from typing import ClassVar
from typing import Dict
from typing import Optional

from mrsimulator.utils.euler_angles import combine_euler_angles
from mrsimulator.utils.parseable import Parseable
from pydantic.v1 import Field

__author__ = "Deepansh Srivastava"
__email__ = "srivastava.89@osu.edu"


class SymmetricTensor(Parseable):
    r"""Base SymmetricTensor class representing the traceless symmetric part of an
    irreducible second-rank tensor.

    Attributes
    ----------

    zeta: float (optional).
        The anisotropy parameter of the nuclear shielding tensor, in ppm, is expressed
        using the Haeberlen convention. The default value is None.

        Example
        -------

        >>> shielding = SymmetricTensor()
        >>> shielding.zeta = 10

    Cq: float (optional).
        The quadrupolar coupling constant, in Hz, is derived from the electric field
        gradient tensor. The default value is None.

        Example
        -------

        >>> efg = SymmetricTensor()
        >>> efg.Cq = 10e6

    eta: float (optional).
        The asymmetry parameter of the SymmetricTensor is expressed using the Haeberlen
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

    property_unit_types: ClassVar[Dict] = {
        "zeta": ["dimensionless", "frequency"],
        "Cq": "frequency",
        "D": "frequency",
        "alpha": "angle",
        "beta": "angle",
        "gamma": "angle",
    }
    property_default_units: ClassVar[Dict] = {
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

    class Config:
        extra = "forbid"

    def rotate(self, euler_angles: list) -> None:
        """Rotate the tensor by the given list of Euler angle rotations. Euler angles
        are given as a list of (alpha, beta, gamma) tuples, and rotations happen in the
        Haeberlen (ZYZ) convention.

        Arguments:
            (list) euler_angles: An ordered list of angle tuples (alpha, beta, gamma)
                to rotate through.

        Example
        -------

        >>> tensor = SymmetricTensor(zeta=10, eta=0.3, alpha=1, beta=1, gamma=1)
        >>> angles = [(3.1415, 0, -3.1415), (1.5701, 1.5701, 1.5701)]
        >>> tensor.rotate(angles)
        """
        # If the tensor (alpha, beta, gamma) is all initialized to None, then assume
        # all zero
        initial_angles = (
            self.alpha if self.alpha is not None else 0,
            self.beta if self.beta is not None else 0,
            self.gamma if self.gamma is not None else 0,
        )

        alpha, beta, gamma = combine_euler_angles([initial_angles] + euler_angles)

        self.alpha = alpha
        self.beta = beta
        self.gamma = gamma


class AntisymmetricTensor(Parseable):
    """Base AntiSymmetricTensor class representing the traceless symmetric part of an
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

    property_unit_types: ClassVar[Dict] = {
        "zeta": "dimensionless",
        "alpha": "angle",
        "beta": "angle",
    }
    property_default_units: ClassVar[Dict] = {
        "zeta": "ppm",
        "alpha": "rad",
        "beta": "rad",
    }
    property_units: Dict = {"zeta": "ppm", "alpha": "rad", "beta": "rad"}

    class Config:
        extra = "forbid"
