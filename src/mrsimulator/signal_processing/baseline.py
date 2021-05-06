# -*- coding: utf-8 -*-
from typing import ClassVar
from typing import Dict
from typing import Union

import numpy as np
from pydantic import validator

from ._base import ModuleOperation
from .utils import _str_to_quantity
from .utils import CONST

__author__ = "Deepansh Srivastava"
__email__ = "srivastava.89@osu.edu"


class Baseline(ModuleOperation):
    module_name: ClassVar = __name__

    @property
    def function(self):
        return "baseline"


class Polynomial(Baseline):
    r"""Add a baseline polynomial to all dependent variables (y) in the CSDM object.

    The baseline function is

    .. math::
            f(x) = \sum_{i=0}^3 x_i \times x^i,

    where :math:`x` is the CSDM dimension at index ``dim_index``.

    Args:
        str x0: Constant zeroth order offset. The default value is 0.
        str x1: First order factor. The default value is 0.
        str x2: Second order factor. The default value is 0.
        str x3: Third order factor. The default value is 0.
        int dim_index: The index of the CSDM dimension along which the operation is
            applied. The default is the dimension at index 0.

    Example
    -------

    >>> from mrsimulator import signal_processing as sp
    >>> operation1 = sp.baseline.Polynomial(x0=20, x1=-10)
    """

    x0: Union[float, str] = 0
    x1: Union[float, str] = 0
    x2: Union[float, str] = 0
    x3: Union[float, str] = 0
    dim_index: int = 0
    property_units: Dict = {"x0": CONST, "x1": CONST, "x2": CONST, "x3": CONST}

    @validator("x0")
    def str_to_quantity0(cls, v, values):
        return _str_to_quantity(v, values, "x0")

    @validator("x1")
    def str_to_quantity1(cls, v, values):
        return _str_to_quantity(v, values, "x1")

    @validator("x2")
    def str_to_quantity2(cls, v, values):
        return _str_to_quantity(v, values, "x2")

    @validator("x3")
    def str_to_quantity3(cls, v, values):
        return _str_to_quantity(v, values, "x3")

    def operate(self, data):
        """Applies the operation.

        Args:
            data: CSDM object
        """
        x = data.dimensions[self.dim_index].coordinates
        if self.x1 == self.x2 == self.x3 == 0:
            data += self.x0
            return data

        if "x1" in self.property_units:
            unit = 1 * self.property_units["x1"]
        elif "x2" in self.property_units:
            unit = np.sqrt(1 * self.property_units["x2"]).unit
        elif "x3" in self.property_units:
            unit = np.cbrt(1 * self.property_units["x3"]).unit

        d1 = self.get_coordinates_in_units(x, unit=1.0 / unit)
        d2 = d1 ** 2
        fn = self.x0 + self.x1 * d1 + d2 * (self.x2 + d1 * self.x3)
        for item in data.y:
            item.components += fn
        return data


class ConstantOffset(Baseline):
    r"""Add an offset to the dependent variables (y) of the CSDM object.

    .. math::
        y += \text{offset}

    where :math:`y` is the CSDM dependent variable.

    Args:
        float offset: The offset factor. The default value is 0.

    Example
    -------

    >>> from mrsimulator import signal_processing as sp
    >>> operation1 = sp.baseline.ConstantOffset(offset=20)
    """

    offset: float = 0

    def operate(self, data):
        """Applies the operation for which the class is named for.

        Args:
            data: CSDM object
        """
        data += self.offset
        return data
