"""The Event class."""
from typing import ClassVar
from typing import Dict
from typing import Union

import numpy as np
from pydantic import validator

from .apodization import Apodization
from .utils import _str_to_quantity
from .utils import CONST

__author__ = "Deepansh Srivastava"


class Correction(Apodization):
    module_name: ClassVar[str] = __name__

    @property
    def function(self):
        return "correction"


class Phase(Correction):
    r"""Phase modulate the dependent variables of CSDM dataset following

    .. math::
        f(x) = e^{-i \phi_0} * e^{-i 2 \pi \phi_1 * x},

    where :math:`x` are the dimension coordinates, and :math:`\phi_0`,
    :math:`\phi_1` are the constant and linear phase, respectively.

    Arguments
    ---------

    constant:
        The constant phase, :math:`\phi_0` in radians. Default is 0.

    linear:
        The linear phase, :math:`\phi_1`. Default is 0.

    dim_index:
        CSDM dimension index along which the operation applies. The
        default is the dimension at index 0.

    dv_index:
        CSDM dependent variable index where the operation applies. If
        unspecified, the operation applies to all dependent variables.

    Example
    -------

    >>> operation4 = sp.correction.Phase(constant=2.14, linear="25.5 Hz", dim_index=0)
    """

    constant: Union[float, str] = 0
    linear: Union[float, str] = 0
    property_units: Dict = {"constant": CONST, "linear": CONST}

    @validator("constant")
    def str_to_quantity_constant(cls, v, values):
        return _str_to_quantity(v, values, "constant")

    @validator("linear")
    def str_to_quantity_linear(cls, v, values):
        return _str_to_quantity(v, values, "linear")

    def fn(self, x):
        x = self.get_coordinates_in_units(x, unit=1.0 / self.property_units["linear"])
        vec = np.exp(-1j * self.constant) * np.exp(-1j * 2 * np.pi * self.linear * x)
        return vec


class AutoPhaseCorrect(Apodization):
    def operate(self, dataset):
        """Apply the operation function.

        Args:
            dataset: A CSDM object.
        """
        max_index = np.unravel_index(int(np.argmax(dataset)), dataset.shape[::-1])
        max_index = max_index[::-1]

        for dim_idx, dim in enumerate(dataset.x):
            dim.coordinates_offset = -dim.coordinates[max_index[dim_idx]]

        dv_indexes = self._get_dv_indexes(self.dv_index, n=len(dataset.y))
        for index in dv_indexes:
            phase_0 = np.angle(dataset[max_index].y[index].components[0])
            multiplier = np.exp(-1j * phase_0)
            dataset.y[index].components *= multiplier

        return dataset
