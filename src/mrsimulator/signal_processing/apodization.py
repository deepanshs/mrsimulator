# -*- coding: utf-8 -*-
"""The Event class."""
from sys import modules
from typing import Dict
from typing import Union

import numpy as np
from csdmpy.units import string_to_quantity
from pydantic import validator

from ._base import AbstractOperation

__author__ = "Maxwell C. Venetos"
__email__ = "maxvenetos@gmail.com"


class AbstractApodization(AbstractOperation):
    dim_index: Union[int, list, tuple] = 0
    dv_index: Union[int, list, tuple] = None  # if none apply to all

    @classmethod
    def parse_dict_with_units(cls, py_dict):
        obj = super().parse_dict_with_units(py_dict)
        return getattr(modules[__name__], py_dict["type"])(**obj.dict())

    @property
    def function(self):
        return "apodization"

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

    def _operate(self, data, fn, prop_name, prop_value):
        """A generic operation function.

        Args:
            data: A CSDM object.
            fn: The apodization function.
            prop_name: The argument name for the function fn.
            prop_value: The argument value for the function fn.
        """
        dims = data.dimensions
        ndim = len(dims)

        dim_index = self.dim_index
        if isinstance(dim_index, int):
            dim_index = [dim_index]

        for dim_index_ in dim_index:
            x = dims[dim_index_].coordinates
            unit = 1 / self.property_units[prop_name]
            x = x.to(unit)
            x_value = x.value

            apodization_vactor = _get_broadcast_shape(
                fn(x_value, prop_value), dim_index_, ndim
            )

            n = len(data.dependent_variables)
            dv_indexes = self._get_dv_indexes(self.dv_index, n=n)

            for i in dv_indexes:
                data.dependent_variables[i].components *= apodization_vactor
        return data


def _str_to_quantity(v, values):
    if isinstance(v, str):
        quantity = string_to_quantity(v)
        values["property_units"] = {"FWHM": quantity.unit}
        return quantity.value
    if isinstance(v, float):
        return v


def _get_broadcast_shape(array, dim, ndim):
    """Return the broadcast shape of a vector `array` at dimension `dim` for `ndim`
    total dimensions. """
    none = [None for _ in range(ndim + 1)]
    if isinstance(dim, int):
        dim = [dim]
    for dim_ in dim:
        none[-dim_ - 1] = slice(None, None, None)
    return array[tuple(none)]


class Gaussian(AbstractApodization):
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

    FWHM: Union[float, str] = 0
    property_units: Dict = {"FWHM": ""}

    @validator("FWHM")
    def str_to_quantity(cls, v, values):
        return _str_to_quantity(v, values)

    # class Config:
    #     validate_assignment = True

    @staticmethod
    def fn(x, arg):
        # arg is FWHM
        sigma = arg / 2.354820045030949
        return np.exp(-2 * ((x * sigma * np.pi) ** 2))

    def operate(self, data):
        """
        Applies the operation for which the class is named for.

        data: CSDM object
        dep_var: int. The index of the dependent variable to apply operation to.
        """

        return self._operate(data, fn=self.fn, prop_name="FWHM", prop_value=self.FWHM)


class Exponential(AbstractApodization):
    r"""Apodize a dependent variable of the CSDM object by an exponential function.

    The apodization function follows

    .. math::
        f(x) = e^{-\Gamma |x| \pi},

    where :math:`x` are the coordinates of the dimension, and :math:`\Gamma` is the
    width parameter. The relationship between the width parameter, :math:`\Gamma`, and
    the full width at half maximum for the reciprocal domain Lorentzian function
    follows

    .. math::
        \text{FWHM} = \Gamma.

    Args:
        str FWHM: The full width at half maximum, FWHM, of the reciprocal domain
            Lorentzian function given as a string with a value and a unit. The default
            value is 0.
        int dim_index: The index of the CSDM dimension along which the operation is
            applied. The default is the dimension at index 0.
        int dv_index: The index of the CSDM dependent variable where the operation is
            applied. If the value is None, the operation will be applied to every
            dependent variable.

    Example
    -------

    >>> operation5 = apo.Exponential(FWHM='143.4 m', dim_index=0, dv_index=0)
    """

    FWHM: Union[float, str] = 0
    property_units: Dict = {"FWHM": ""}

    @validator("FWHM")
    def str_to_quantity(cls, v, values):
        return _str_to_quantity(v, values)

    # class Config:
    #     validate_assignment = True

    @staticmethod
    def fn(x, arg):
        return np.exp(-arg * np.pi * np.abs(x))

    def operate(self, data):
        """
        Applies the operation for which the class is named for.

        data: CSDM object
        dep_var: int. The index of the dependent variable to apply operation to.
        """

        return self._operate(data, fn=self.fn, prop_name="FWHM", prop_value=self.FWHM)


# class ExponentialAbs(AbstractApodization):
#     r"""Apodize a dependent variable of the simulation data object by an exponential
#     function. The function follows

#     .. math::
#         f(x) = e^{-\Tau |x| \pi},

#     where :math:`x` are the coordinates of the data dimension and :math:`\Tau` is
#     the width parameter.

#     Args:
#         int dim_index: Data dimension index to apply the function along.
#         float FWHM: The full width at half maximum parameter, :math:`\Tau`.
#         int dv_index: Data dependent variable index to apply the function to. If
#             the type None, the operation will be applied to every dependent variable.

#     Example
#     -------

#     >>> operation5 = apo.Exponential(FWHM=143.4, dim_index=0, dv_index=0)
#     """

#     FWHM: float = 0

#     property_unit_types: ClassVar = {"FWHM": ["time", "frequency"]}
#     property_default_units: ClassVar = {"FWHM": ["s", "Hz"]}
#     property_units: Dict = {"FWHM": "Hz"}

#     @staticmethod
#     def fn(x, arg):
#         return np.exp(-arg * np.pi * np.abs(x))

#     def operate(self, data):
#         """
#         Applies the operation for which the class is named for.

#         data: CSDM object
#         dep_var: int. The index of the dependent variable to apply operation to
#         """

#         return self._operate(data, fn=self.fn, prop_name="FWHM", prop_value=self.FWHM)
