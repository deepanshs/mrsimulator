# -*- coding: utf-8 -*-
"""The Event class."""
from sys import modules
from typing import Dict
from typing import Union

import numpy as np
from csdmpy.units import string_to_quantity
from pydantic import validator

from ._base import AbstractOperation

__author__ = "Deepansh Srivastava"
__email__ = "deepansh2012@gmail.com"


const = string_to_quantity("1")


class AbstractAffineTransformation(AbstractOperation):
    dim_index: int = 0
    normal: int = 0
    dv_index: Union[int, list, tuple] = None  # if none apply to all

    @classmethod
    def parse_dict_with_units(cls, py_dict):
        obj = super().parse_dict_with_units(py_dict)
        return getattr(modules[__name__], py_dict["type"])(**obj.dict())

    @property
    def function(self):
        return "affine"

    @property
    def type(self):
        """The type apodization function."""
        return self.__class__.__name__

    @staticmethod
    def _get_dv_indexes(indexes, n):
        """Return a list of dependent variable indexes.

        Args:
            indexes: An interger, list of integers, or None indicating the dv indexes.
            n: Total number of dependent variables in the CSDM object.
        """
        if indexes is None:
            return np.arange(n)
        if isinstance(indexes, int):
            return [indexes]
        if isinstance(indexes, (list, tuple)):
            return np.asarray(indexes)


def _get_broadcast_shape(array, dim, ndim):
    """Return the broadcast shape of a vector `array` at dimension `dim` for `ndim`
    total dimensions. """
    none = [np.newaxis for _ in range(ndim + 1)]
    if isinstance(dim, int):
        dim = [dim]
    for dim_ in dim:
        none[-dim_ - 1] = slice(None, None, None)

    print(none)
    return array[tuple(none)]


def _str_to_quantity(v, values, prop_name):
    if isinstance(v, str):
        quantity = string_to_quantity(v)
        values["property_units"] = {prop_name: quantity.unit}
        return quantity.value
    if isinstance(v, float):
        return v


class Shear(AbstractAffineTransformation):
    r"""Apodize a dependent variable of the CSDM object with a Gaussian function.

    The apodization function follows

    .. math::
        f(x) = e^{-2 \pi^2 \sigma^2  x^2},

    where :math:`x` are the coordinates of the dimension, and :math:`\sigma` is the
    standard deviation. The relationship between the standard deviation, :math:`\sigma`,
    and the full width at half maximum of the reciprocal domain Gaussian function
    follows

    .. math::
        \sigma = \frac{\text{FWHM}}{2\sqrt{2\ln 2}}.

    Args:
        str FWHM: The full width at half maximum, FWHM, of the reciprocal domain
            Gaussian function, given as a string with a value and a unit. The default
            value is 0.
        int dim_index: The index of the CSDM dimension along which the operation is
            applied. The default is the dimension at index 0.
        int dv_index: The index of the CSDM dependent variable where the operation is
            applied. If the value is None, the operation will be applied to every
            dependent variable.

    Example
    -------

    >>> import mrsimulator.signal_processing.apodization as apo
    >>> operation4 = apo.Gaussian(FWHM='143.4 Hz', dim_index=0, dv_index=0)
    """

    factor: Union[float, str] = 0
    property_units: Dict = {"factor": const}

    @validator("factor")
    def str_to_quantity(cls, v, values):
        return _str_to_quantity(v, values, "factor")

    # class Config:
    #     validate_assignment = True

    def operate(self, data):
        """
        Applies the operation for which the class is named for.

        data: CSDM object
        dep_var: int. The index of the dependent variable to apply operation to.
        """
        data = data.fft(axis=self.dim_index)

        dims = data.dimensions
        n_dim = len(dims)

        x = dims[self.dim_index]
        y = dims[self.normal]

        vector_x = _get_broadcast_shape(x.coordinates.value, self.dim_index, n_dim)
        vector_y = _get_broadcast_shape(y.coordinates.value, self.normal, n_dim)

        xy = vector_x * vector_y
        unit = 1 / self.property_units["factor"]
        multiplier = unit.to(x.increment.unit * y.increment.unit).value

        # apodization by exp(-i 2π w t a), where w is freq, t is time, and
        # a is the factor.
        apodization_vector = np.exp(-2j * np.pi * xy * self.factor * multiplier)

        n_dv = len(data.dependent_variables)
        dv_indexes = self._get_dv_indexes(self.dv_index, n=n_dv)

        for i in dv_indexes:
            datum = data.dependent_variables[i].components
            datum *= apodization_vector

        data = data.fft(axis=self.dim_index)
        return data


class Scale(AbstractAffineTransformation):
    r"""Apodize a dependent variable of the CSDM object with a Gaussian function.

    The apodization function follows

    .. math::
        f(x) = e^{-2 \pi^2 \sigma^2  x^2},

    where :math:`x` are the coordinates of the dimension, and :math:`\sigma` is the
    standard deviation. The relationship between the standard deviation, :math:`\sigma`,
    and the full width at half maximum of the reciprocal domain Gaussian function
    follows

    .. math::
        \sigma = \frac{\text{FWHM}}{2\sqrt{2\ln 2}}.

    Args:
        str FWHM: The full width at half maximum, FWHM, of the reciprocal domain
            Gaussian function, given as a string with a value and a unit. The default
            value is 0.
        int dim_index: The index of the CSDM dimension along which the operation is
            applied. The default is the dimension at index 0.
        int dv_index: The index of the CSDM dependent variable where the operation is
            applied. If the value is None, the operation will be applied to every
            dependent variable.

    Example
    -------

    >>> import mrsimulator.signal_processing.apodization as apo
    >>> operation4 = apo.Gaussian(FWHM='143.4 Hz', dim_index=0, dv_index=0)
    """

    factor: Union[float, str] = 1
    property_units: Dict = {"factor": const}

    @validator("factor")
    def str_to_quantity(cls, v, values):
        return _str_to_quantity(v, values, "factor")

    # class Config:
    #     validate_assignment = True

    @staticmethod
    def fn(wt, factor):
        return np.exp(-2j * np.pi * wt * factor)

    def operate(self, data):
        """
        Applies the operation for which the class is named for.

        data: CSDM object
        dep_var: int. The index of the dependent variable to apply operation to.
        """
        data.dimensions[self.dim_index].increment *= self.factor
        data.dimensions[self.dim_index].coordinates_offset *= self.factor
        return data
