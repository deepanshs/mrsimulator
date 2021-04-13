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

__author__ = "Maxwell C. Venetos"
__email__ = "maxvenetos@gmail.com"


class AbstractApodization(AbstractOperation):
    dim_index: Union[int, list, tuple] = 0
    dv_index: Union[int, list, tuple] = None  # if none apply to all

    @classmethod
    def parse_dict_with_units(cls, py_dict: dict):
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
    def _get_correct_coordinates(x, unit):
        return x.to(unit).value

    def operate(self, data):
        """Apply the operation function.

        Args:
            data: A CSDM object.
        """
        dims = data.dimensions
        ndim = len(dims)

        dim_index = self.dim_index
        if isinstance(dim_index, int):
            dim_index = [dim_index]

        for dim_index_ in dim_index:
            x = dims[dim_index_].coordinates
            apodization_vactor = _get_broadcast_shape(self.fn(x), dim_index_, ndim)

            n = len(data.y)
            dv_indexes = self._get_dv_indexes(self.dv_index, n=n)

            for i in dv_indexes:
                data.y[i].components *= apodization_vactor
        return data


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

    >>> operation4 = sp.apodization.Gaussian(FWHM='143.4 Hz', dim_index=0, dv_index=0)
    """

    FWHM: Union[float, str] = 0
    property_units: Dict = {"FWHM": CONST}

    @validator("FWHM")
    def str_to_quantity(cls, v, values):
        return _str_to_quantity(v, values, "FWHM")

    def fn(self, x):
        x = self._get_correct_coordinates(x, unit=1.0 / self.property_units["FWHM"])
        sigma = self.FWHM / 2.354820045030949
        return 1.0 if self.FWHM == 0.0 else np.exp(-2.0 * (np.pi * sigma * x) ** 2)


class Exponential(AbstractApodization):
    r"""Apodize a dependent variable of the CSDM object by an exponential function.

    The apodization function follows

    .. math::
        f(x) = e^{-\Gamma \pi |x|},

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

    >>> operation5 = sp.apodization.Exponential(FWHM='143.4 m', dim_index=0, dv_index=0)
    """

    FWHM: Union[float, str] = 0
    property_units: Dict = {"FWHM": CONST}

    @validator("FWHM")
    def str_to_quantity(cls, v, values):
        return _str_to_quantity(v, values, "FWHM")

    def fn(self, x):
        x = self._get_correct_coordinates(x, unit=1.0 / self.property_units["FWHM"])
        return 1.0 if self.FWHM == 0.0 else np.exp(-self.FWHM * np.pi * np.abs(x))
