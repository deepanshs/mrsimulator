# -*- coding: utf-8 -*-
"""The Event class."""
from sys import modules
from typing import Dict
from typing import Union

import numpy as np
from pydantic import validator

from ._base import AbstractOperation
from .utils import _get_broadcast_shape
from .utils import _str_to_quantity
from .utils import CONST

__author__ = "Deepansh Srivastava"
__email__ = "srivastava.89@osu.edu"


class AbstractAffineTransformation(AbstractOperation):
    dim_index: int = 0
    dv_index: Union[int, list, tuple] = None  # if none apply to all

    @classmethod
    def parse_dict_with_units(cls, py_dict: dict):
        obj = super().parse_dict_with_units(py_dict)
        return getattr(modules[__name__], py_dict["type"])(**obj.dict())

    @property
    def function(self):
        return "affine"

    @property
    def type(self):
        """The type apodization function."""
        return self.__class__.__name__


def get_coordinates(dim):
    """Return the coordinates of dimension, dim, without equivalence units"""

    def get_coords(dim):
        equivalent_fn, equivalent_unit = dim._equivalencies, dim._equivalent_unit
        dim._equivalencies, dim._equivalent_unit = None, None

        coordinates = dim.coordinates

        dim._equivalencies, dim._equivalent_unit = equivalent_fn, equivalent_unit
        return coordinates

    return get_coords(dim) if not hasattr(dim, "subtype") else get_coords(dim.subtype)


class Shear(AbstractAffineTransformation):
    r"""Apply a shear parallel to dimension at index parallel and normal to dimension
    at index dim_index.

    The shear function is an apodization with the following form

    .. math::
        f(x) = e^{-i 2\pi x_0 x_1 a_0},

    where :math:`x_0` are the coordinates of the dimension at index `parallel`,
    :math:`x_1` are the coordinates of the dimension at index `dim_index`, and
    :math:`a_0` is the shear constant.

    Args:
        str factor: The shear factor is given as a string with a value and a unit. The
            default value is 0.
        int dim_index: The shear is applied normal to the CSDM dimension at this index.
            The default is the dimension at index 0.
        int parallel: The shear is applied parallel to the CSDM dimension at this index.
            The default is the dimension at index 1.
        int dv_index: The index of the CSDM dependent variable where the operation is
            applied. If the value is None, the operation will be applied to every
            dependent variable.

    Example
    -------

    >>> operation = sp.affine.Shear(factor='143.4 Hz', dim_index=0, parallel=1)
    """

    factor: Union[float, str] = 0
    parallel: int = 1
    property_units: Dict = {"factor": CONST}

    @validator("factor")
    def str_to_quantity(cls, v, values):
        return _str_to_quantity(v, values, "factor")

    # class Config:
    #     validate_assignment = True

    def operate(self, data):
        """
        Applies the operation for which the class is named for.

        data: CSDM object.
        """
        dims = data.dimensions
        n_dim = len(dims)

        x, y = dims[self.dim_index], dims[self.parallel]
        x_value, y_value = [get_coordinates(_).value for _ in [x, y]]

        vector_x = _get_broadcast_shape(x_value, self.dim_index, n_dim)
        vector_y = _get_broadcast_shape(y_value, self.parallel, n_dim)

        xy = vector_x * vector_y
        unit = 1 / self.property_units["factor"]

        multiplier = unit.to(x.increment.unit * y.increment.unit).value

        # apodization by exp(-i 2π w t a), where w is freq, t is time, and
        # a is the factor.
        apodization_vector = np.exp(-2j * np.pi * xy * self.factor * multiplier)

        dv_indexes = self._get_dv_indexes(self.dv_index, n=len(data.y))

        for i in dv_indexes:
            datum = data.y[i].components
            datum *= apodization_vector

        return data


class Scale(AbstractAffineTransformation):
    r"""Scale the dimension along the specified dimension index.

    Args:
        str factor: The scaling factor. The default is 1.
        int dim_index: The index of the CSDM dimension to scale. The default is the
            dimension at index 0.

    Example
    -------

    >>> operation = sp.affine.Scale(factor=2.14, dim_index=0)
    """

    factor: Union[float, str] = 1
    property_units: Dict = {"factor": CONST}

    @validator("factor")
    def str_to_quantity(cls, v, values):
        return _str_to_quantity(v, values, "factor")

    def operate(self, data):
        """
        Applies the operation for which the class is named for.

        data: CSDM object.
        """
        data_ref = data.x[self.dim_index]
        data_ref.increment *= self.factor
        data_ref.coordinates_offset *= self.factor
        data_ref.reciprocal.coordinates_offset /= self.factor
        return data


# class Translate(AbstractAffineTransformation):
#     r"""Apodize a dependent variable of the CSDM object with a Gaussian function.

#     The apodization function follows

#     .. math::
#         f(x) = e^{-2 \pi^2 \sigma^2  x^2},

#     where :math:`x` are the coordinates of the dimension, and :math:`\sigma` is the
#     standard deviation. The relationship between the standard deviation,
#     :math:`\sigma`,
#     and the full width at half maximum of the reciprocal domain Gaussian function
#     follows

#     .. math::
#         \sigma = \frac{\text{FWHM}}{2\sqrt{2\ln 2}}.

#     Args:
#         str FWHM: The full width at half maximum, FWHM, of the reciprocal domain
#             Gaussian function, given as a string with a value and a unit. The default
#             value is 0.
#         int dim_index: The index of the CSDM dimension along which the operation is
#             applied. The default is the dimension at index 0.
#         int dv_index: The index of the CSDM dependent variable where the operation is
#             applied. If the value is None, the operation will be applied to every
#             dependent variable.

#     Example
#     -------

#     >>> operation4 = sp.apodization.Gaussian(FWHM='143.4 Hz', dim_index=0, dv_index=0)
#     """

#     factor: Union[float, str] = 1
#     property_units: Dict = {"factor": CONST}

#     @validator("factor")
#     def str_to_quantity(cls, v, values):
#         return _str_to_quantity(v, values, "factor")

#     # class Config:
#     #     validate_assignment = True

#     def operate(self, data):
#         """
#         Applies the operation for which the class is named for.

#         data: CSDM object
#         dep_var: int. The index of the dependent variable to apply operation to.
#         """
#         # data.x[self.dim_index].increment *= self.factor
#         data.x[self.dim_index].coordinates_offset += (
#             self.factor * self.property_units["factor"]
#         )
#         return data
