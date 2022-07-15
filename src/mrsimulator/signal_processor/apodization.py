"""The Event class."""
from typing import ClassVar
from typing import Dict
from typing import List
from typing import Union

import numpy as np
from pydantic import validator
from scipy.special import erf

from ._base import ModuleOperation
from .utils import _get_broadcast_shape
from .utils import _str_to_quantity
from .utils import CONST

__author__ = "Maxwell C. Venetos"
__email__ = "maxvenetos@gmail.com"


class Apodization(ModuleOperation):
    dim_index: Union[int, list, tuple] = 0
    dv_index: Union[int, list, tuple] = None  # if none apply to all
    module_name: ClassVar[str] = __name__

    @property
    def function(self):
        return "apodization"

    def operate(self, dataset):
        """Apply the operation function.

        Args:
            dataset: A CSDM object.
        """
        dims = dataset.dimensions
        ndim = len(dims)

        dim_index = self.dim_index
        dim_index = [dim_index] if isinstance(dim_index, int) else dim_index

        for i in dim_index:
            x = self.get_coordinates(dims[i])  # dims[i].coordinates
            apodization_vector = _get_broadcast_shape(self.fn(x), i, ndim)

            dv_indexes = self._get_dv_indexes(self.dv_index, n=len(dataset.y))

            for index in dv_indexes:
                dataset.y[index].components *= apodization_vector
        return dataset


class MultiDimensionApodization(ModuleOperation):
    dim_index: Union[int, list, tuple] = 0
    dv_index: Union[int, list, tuple] = None  # if none apply to all
    module_name: ClassVar[str] = __name__

    @property
    def function(self):
        return "apodization"

    def operate(self, dataset):
        """Apply the operation function.

        Args:
            dataset: A CSDM object.
        """
        dims = dataset.dimensions
        ndim = len(dims)

        dim_index = self.dim_index
        dim_index = [dim_index] if isinstance(dim_index, int) else dim_index

        apodization_matrix = _get_broadcast_shape(self.fn(), dim_index, ndim)

        dv_indexes = self._get_dv_indexes(self.dv_index, n=len(dataset.y))

        for index in dv_indexes:
            dataset.y[index].components *= apodization_matrix
        return dataset


class Gaussian(Apodization):
    r"""Apodize dependent variables of CSDM dataset with Gaussian function.

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
    r"""Apodize dependent variables of CSDM by exponential function.

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
    r"""Apodize dependent variables of CSDM dataset with skewed Gaussian function.

    The apodization function is derived from the skewed Gaussian distribution

    .. math::
        f(x) = 2\phi(x)\Phi(\alpha x),

    where :math:`x` are the coordinates of the dimension, and :math:`\phi` is the
    standard normal probability density function, :math:`\Phi` is the cumulative
    distribution function, and :math:`\alpha` is a skewing parameter. The apodization
    function is the fourier transform of the above function which gives another
    skewed Gaussian function and is given by

    .. math::
        f(x) = e^{-2  (\pi x)^2}
        \left(1 + i\text{Erfi}\left(\frac{\text{skew}\cdot x}{\sqrt{2}}\right)\right),

    where skew is given by

    .. math::
        \text{skew} = \frac{\alpha}{\sqrt{1+\alpha^2}}

    See https://en.wikipedia.org/wiki/Skew_normal_distribution

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

    >>> operation6 = sp.apodization.SkewedGaussian(skew=2, dim_index=0, dv_index=0)
    """
    skew: Union[float, str] = 0
    FWHM: Union[float, str] = "2.354820045030949 Hz"
    property_units: Dict = {"FWHM": CONST}

    @validator("FWHM")
    def str_to_quantity(cls, v, values):
        return _str_to_quantity(v, values, "FWHM")

    def fn(self, x):
        sigma = self.FWHM / 2.354820045030949
        x = self.get_coordinates_in_units(x, unit=1.0 / self.property_units["FWHM"])
        prob_func = [np.exp(-(0.5) * (sigma * j) ** 2) for j in x]
        cum_prob_func = [0.5 + 0.5 * erf(self.skew * j / np.sqrt(2)) for j in x]
        sg = np.asarray([a * b for a, b in zip(prob_func, cum_prob_func)])
        return 1.0 if self.skew == 0.0 else sg


class TopHat(Apodization):
    r"""Apodize dependent variables of CSDM object by top hat function.

    The apodization function follows

    .. math::
        f(x) = 1 \text{if rising_edge} <= x <= \text{falling_edge} ,
        \text{else} f(x) = 0

    where :math:`x` are the coordinates of the dimension, rising_edge is the
    start of the function window, and falling_edge is the end of the
    function window.

    Arguments
    ---------

    rising_edge:
        The lowest value in the time domain from which to start the function
        window. The default value is None which will take the lowest possible
        value for the supplied dataset.

    falling_edge:
        The highest value in the time domain from which to end the function
        window. The default value is None which will take the largest possible
        value for the supplied dataset.

    dim_index:
        The index of the CSDM dimension along which the operation is applied. The
        default is the dimension at index 0.

    dv_index:
        The index of the CSDM dependent variable, where the operation is applied. If
        not provided, the operation will be applied to every dependent variable.

    Example
    -------

    >>> operation7= sp.apodization.TopHat(rising_edge = "-1 s", falling_edge = "1 s")
    """

    rising_edge: Union[float, str, None] = -np.Inf
    falling_edge: Union[float, str, None] = np.Inf
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
        screen = np.where(
            np.logical_and(x >= self.rising_edge, x < self.falling_edge), 1, 0
        )
        return screen


class TypedArray(np.ndarray):
    """an array metaclass to allow the use of np.ndarrays in Mask class"""

    @classmethod
    def __get_validators__(cls):
        yield cls.validate_type

    @classmethod
    def validate_type(cls, val):
        return np.array(val, dtype=cls.inner_type)


class ArrayMeta(type):
    def __getitem__(self, t):
        return type("Array", (TypedArray,), {"inner_type": t})


class Array(np.ndarray, metaclass=ArrayMeta):

    pass


class Mask(MultiDimensionApodization):
    r"""Apodize dependent variables of CSDM object by user defined mask.

    The apodization function follows

    .. math::
        f(x) = \text{mask}

    where mask is a user-supplied numpy array containing an apodization mask to apply
    to the dataset.

    Arguments
    ---------

    mask:
        mask.

    dim_index:
        The index of the CSDM dimension along which the operation is applied. The
        default is the dimension at index 0.

    dv_index:
        The index of the CSDM dependent variable, where the operation is applied. If
        not provided, the operation will be applied to every dependent variable.

    Example
    -------
    """

    mask: Union[List[float], Array[float], float, int] = None

    class Config:
        arbitrary_types_allowed = True

    def fn(self):
        return self.mask
