# -*- coding: utf-8 -*-
"""The Event class."""
from typing import ClassVar
from typing import Dict
from typing import Union

import numpy as np
from pydantic import validator
from scipy.special import erfi

from ._base import ModuleOperation
from .utils import _get_broadcast_shape
from .utils import _str_to_quantity
from .utils import CONST

__author__ = "Maxwell C. Venetos"
__email__ = "maxvenetos@gmail.com"


class Apodization(ModuleOperation):
    dim_index: Union[int, list, tuple] = 0
    dv_index: Union[int, list, tuple] = None  # if none apply to all
    module_name: ClassVar = __name__

    @property
    def function(self):
        return "apodization"

    def operate(self, data):
        """Apply the operation function.

        Args:
            data: A CSDM object.
        """
        dims = data.dimensions
        ndim = len(dims)

        dim_index = self.dim_index
        dim_index = [dim_index] if isinstance(dim_index, int) else dim_index

        for i in dim_index:
            x = self.get_coordinates(dims[i])  # dims[i].coordinates
            apodization_vactor = _get_broadcast_shape(self.fn(x), i, ndim)

            dv_indexes = self._get_dv_indexes(self.dv_index, n=len(data.y))

            for index in dv_indexes:
                data.y[index].components *= apodization_vactor
        return data


class Gaussian(Apodization):
    r"""Apodize dependent variable objects of the CSDM data with a Gaussian function.

    The apodization function follows

    .. math::
        f(x) = e^{-2 \pi^2 \sigma^2  x^2},

    where :math:`x` are the coordinates of the dimension, and :math:`\sigma` is the
    standard deviation. The relationship between the standard deviation, :math:`\sigma`,
    and the full width at half maximum of the reciprocal domain Gaussian function
    follows

    .. math::
        \sigma = \frac{\text{FWHM}}{2\sqrt{2\ln 2}}.

    Arguments
    ---------

    FWHM:
        The full width at half maximum, FWHM, of the reciprocal domain Gaussian
        function given as a string with a value and a unit. The default value is 0.

    dim_index:
        The index of the CSDM dimension along which the operation is applied. The
        default is the dimension at index 0.

    dv_index:
        The index of the CSDM dependent variable, where the operation is applied. If
        not provided, the operation will be applied to every dependent variable.

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
        x = self.get_coordinates_in_units(x, unit=1.0 / self.property_units["FWHM"])
        sigma = self.FWHM / 2.354820045030949
        return 1.0 if self.FWHM == 0.0 else np.exp(-2.0 * (np.pi * sigma * x) ** 2)


class Exponential(Apodization):
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

    Arguments
    ---------

    FWHM:
        The full width at half maximum, FWHM, of the reciprocal domain Lorentzian
        function given as a string with a value and a unit. The default value is 0.

    dim_index:
        The index of the CSDM dimension along which the operation is applied. The
        default is the dimension at index 0.

    dv_index:
        The index of the CSDM dependent variable, where the operation is applied. If
        not provided, the operation will be applied to every dependent variable.

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
        x = self.get_coordinates_in_units(x, unit=1.0 / self.property_units["FWHM"])
        return 1.0 if self.FWHM == 0.0 else np.exp(-self.FWHM * np.pi * np.abs(x))


class SkewedGaussian(Apodization):
    r"""Apodize dependent variable objects of the CSDM data with a skewed Gaussian function.

    The apodization function is derived from the skewed Gaussian distribution

    .. math::
        f(x) = 2*\phi(x)*\Phi(\alpha x),

    where :math:`x` are the coordinates of the dimension, and :math:`\phi` is the
    standard normal probability density function, :math:`\Phi` is the cumulative
    distribution function, and :math:`\alpha` is a skewing parameter. The apodization
    function is the fourier transform of the above function which gives another
    skewed Gaussian function and is given by

    .. math::
        f(x) = e^{-2 * (\pi * x)**2}*{1 + i*Erfi(skew*x/\sqrt(2))},

    where skew is given by

    .. math::
        skew = \alpha/\sqrt(1+\alpha**2)

    Arguments
    ---------

    sigma:
        The full width at half maximum, FWHM, of the reciprocal domain Gaussian
        function given as a string with a value and a unit. The default value is 0.

    skew:
        The skewness defining the asymmetry of the Gaussian distribution

    dim_index:
        The index of the CSDM dimension along which the operation is applied. The
        default is the dimension at index 0.

    dv_index:
        The index of the CSDM dependent variable, where the operation is applied. If
        not provided, the operation will be applied to every dependent variable.

    Example
    -------

    >>> operation4 = sp.apodization.Gaussian(FWHM='143.4 Hz', dim_index=0, dv_index=0)
    """
    skew: Union[float, str] = 0
    property_units: Dict = {"skew": CONST}

    @validator("skew")
    def str_to_quantity(cls, v, values):
        return _str_to_quantity(v, values, "skew")

    def fn(self, x):
        x = self.get_coordinates_in_units(x, unit=1.0 / self.property_units["FWHM"])
        return (
            1.0
            if self.skew == 0.0
            else np.exp(-2.0 * (np.pi * x) ** 2)
            * (1 + erfi(self.skew * x / np.sqrt(2)))
        )


class Step(Apodization):
    r"""Apodize a dependent variable of the CSDM object by a step function.

    The apodization function follows

    .. math::
        f(x) = 1 if rising_edge <= x <= falling_edge , ####(upper/lower_bound???)
        else f(x) = 0

    where :math:`x` are the coordinates of the dimension, :math`rising_edge` is the
    start of the step function window, and :math'falling_edge' is the end of the
    step function window.

    Arguments
    ---------

    rising_edge:
        The lowest value in the time domain from which to start the step function
        window. The default value is None which will take the lowest possible
        value for the supplied data.

    falling_edge:
        The highest value in the time domain from which to end the step function
        window. The default value is None which will take the largest possible
        value for the supplied data.

    dim_index:
        The index of the CSDM dimension along which the operation is applied. The
        default is the dimension at index 0.

    dv_index:
        The index of the CSDM dependent variable, where the operation is applied. If
        not provided, the operation will be applied to every dependent variable.

    Example
    -------

    >>> operation6= sp.apodization.Step(rising_edge = '-1 s', falling_edge = '1 s')
    """

    rising_edge: Union[float, str] = 0
    falling_edge: Union[float, str] = 0
    property_units: Dict = {"rising_edge": CONST, "falling_edge": CONST}

    @validator("rising_edge")
    def str_to_quantity_l(cls, v, values):
        return _str_to_quantity(v, values, "rising_edge")

    @validator("falling_edge")
    def str_to_quantity_u(cls, v, values):
        return _str_to_quantity(v, values, "falling_edge")

    def fn(self, x):
        if "falling_edge" in self.property_units:
            unit = 1 * self.property_units["falling_edge"]
        if "rising_edge" in self.property_units:
            unit = 1 * self.property_units["rising_edge"]
        x = self.get_coordinates_in_units(x, unit=1.0 * unit)

        if self.rising_edge is None:
            self.rising_edge = x.min()

        if self.falling_edge is None:
            self.falling_edge = x.max()

        screen = np.where(x > self.rising_edge, 1, 0)
        screen = screen + np.where(x < self.falling_edge, 0, -1)
        return screen
