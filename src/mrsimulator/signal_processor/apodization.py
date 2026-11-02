"""The Event class."""
from typing import ClassVar
from typing import Dict
from typing import List
from typing import Union

import numpy as np
from astropy.units import Quantity
from pydantic.v1 import validator
from scipy.special import dawsn

from ._base import ModuleOperation
from .utils import _get_broadcast_shape
from .utils import _str_to_quantity
from .utils import CONST

__author__ = "Maxwell C. Venetos"
__email__ = "maxvenetos@gmail.com"


def _multiply_components(dependent_variable, window):
    """Multiply the components of a CSDM dependent variable by the apodization window
    in place. A complex window applied to real-valued components promotes the
    components to complex."""
    components = dependent_variable.components
    if np.iscomplexobj(window) and not np.iscomplexobj(components):
        dependent_variable.components = components * window
    else:
        components *= window


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
                _multiply_components(dataset.y[index], apodization_vector)
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
            _multiply_components(dataset.y[index], apodization_matrix)
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
        return (
            np.ones_like(x)
            if self.FWHM == 0.0
            else np.exp(-2.0 * (np.pi * sigma * x) ** 2)
        )


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
        return (
            np.ones_like(x)
            if self.FWHM == 0.0
            else np.exp(-self.FWHM * np.pi * np.abs(x))
        )


class SkewedGaussian(Apodization):
    r"""Apodize dependent variables of CSDM dataset with skewed Gaussian function.

    The apodization function is the Fourier transform of the skew normal
    distribution,

    .. math::
        p(\nu) = \frac{2}{\omega}\, \phi\left(\frac{\nu}{\omega}\right)
        \Phi\left(\alpha \frac{\nu}{\omega}\right),

    where :math:`\phi` and :math:`\Phi` are the standard normal probability density
    and cumulative distribution functions, :math:`\omega` is the scale, and
    :math:`\alpha` is the shape (skew) parameter. The apodization function is

    .. math::
        f(x) = e^{-2 \pi^2 \omega^2 x^2}
        \left(1 + i\, \text{erfi}\left(\sqrt{2}\, \pi \delta \omega x\right)
        \right),
        \quad \delta = \frac{\alpha}{\sqrt{1 + \alpha^2}},

    where :math:`x` are the coordinates of the dimension, and :math:`\omega` follows
    from the FWHM of the underlying Gaussian,

    .. math::
        \omega = \frac{\text{FWHM}}{2\sqrt{2\ln 2}}.

    For zero skew, the function reduces to the Gaussian apodization of the same FWHM.
    A positive skew gives a line shape with a tail towards higher frequencies. Note,
    the function is complex and is meant for apodizing a complex (time domain) signal.

    See https://en.wikipedia.org/wiki/Skew_normal_distribution

    Arguments
    ---------

    FWHM:
        The full width at half maximum, FWHM, of the reciprocal domain Gaussian
        function (zero skew) given as a string with a value and a unit. The default
        value is 0.

    skew:
        The shape parameter, :math:`\alpha`, defining the asymmetry of the line shape.
        The default value is 0.

    dim_index:
        The index of the CSDM dimension along which the operation is applied. The
        default is the dimension at index 0.

    dv_index:
        The index of the CSDM dependent variable, where the operation is applied. If
        not provided, the operation will be applied to every dependent variable.

    Example
    -------

    >>> operation6 = sp.apodization.SkewedGaussian(FWHM='143.4 Hz', skew=2, dim_index=0)
    """

    skew: Union[float, str] = 0
    FWHM: Union[float, str] = 0
    property_units: Dict = {"FWHM": CONST}

    @validator("FWHM")
    def str_to_quantity(cls, v, values):
        return _str_to_quantity(v, values, "FWHM")

    def fn(self, x):
        x = self.get_coordinates_in_units(x, unit=1.0 / self.property_units["FWHM"])
        if self.FWHM == 0.0:
            return np.ones_like(x)

        omega = self.FWHM / 2.354820045030949
        delta = self.skew / np.sqrt(1.0 + self.skew**2)
        arg = np.pi * omega * np.asarray(x, dtype=float)
        gauss = np.exp(-2.0 * arg**2)

        # erfi(y) = (2/sqrt(pi)) exp(y^2) D(y), where D is the Dawson function. The
        # product with the Gaussian, exp(-2 arg^2 (1 - delta^2)) D(y), stays finite.
        y = np.sqrt(2.0) * delta * arg
        norm = 2.0 / np.sqrt(np.pi)
        skew_part = norm * np.exp(-2.0 * arg**2 * (1.0 - delta**2)) * dawsn(y)
        return gauss + 1j * skew_part


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

    rising_edge: Union[float, str, None] = -np.inf
    falling_edge: Union[float, str, None] = np.inf
    property_units: Dict = {"rising_edge": CONST, "falling_edge": CONST}

    @validator("rising_edge")
    def str_to_quantity_l(cls, v, values):
        return _str_to_quantity(v, values, "rising_edge")

    @validator("falling_edge")
    def str_to_quantity_u(cls, v, values):
        return _str_to_quantity(v, values, "falling_edge")

    def _coordinates_for_edge(self, x, edge):
        """Coordinates x in the unit of the given edge. An edge without a unit is
        compared with the coordinates in the unit of the dimension."""
        unit = self.property_units.get(edge, CONST)
        if isinstance(unit, Quantity):
            return x.value
        return self.get_coordinates_in_units(x, unit=1.0 * unit)

    def fn(self, x):
        rising = self._coordinates_for_edge(x, "rising_edge") >= self.rising_edge
        falling = self._coordinates_for_edge(x, "falling_edge") < self.falling_edge
        return np.where(np.logical_and(rising, falling), 1, 0)


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
