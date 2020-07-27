# -*- coding: utf-8 -*-
"""The Event class."""
from sys import modules
from typing import ClassVar
from typing import Dict
from typing import Union

import numpy as np

from ._base import AbstractOperation

__author__ = "Maxwell C. Venetos"
__email__ = "maxvenetos@gmail.com"


class AbstractApodization(AbstractOperation):
    dim_indx: int = 0
    dv_indx: Union[int, list, tuple] = None  # if none apply to all

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

    def set_property_units(self, unit, prop):
        """
        Populate the property unit attribute of the class based on the dimension unit.
        """
        if unit.physical_type == "frequency":
            self.property_units = {f"{prop}": "s"}
            return
        self.property_units = {f"{prop}": "Hz"}

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
        x = data.dimensions[self.dim_indx].coordinates
        x_value, x_unit = x.value, x.unit

        self.set_property_units(x_unit, prop_name)

        apodization_vactor = fn(x_value, prop_value)

        n = len(data.dependent_variables)
        dv_indexes = self._get_dv_indexes(self.dv_indx, n=n)

        for i in dv_indexes:
            data.dependent_variables[i].components[0] *= apodization_vactor
        return data


class Gaussian(AbstractApodization):
    r"""Apodize a dependent variable of the simulation data object by a Gaussian
    function. The function follows

    .. math::
        f(x) = e^{-2 \pi^2 \frac{\text{FWHM}{2\sqrt{2\ln 2}}}^2  x^2},

    where :math:`x` are the coordinates of the data dimension and FWHM is the full
    width at half maximum of the Gaussian.

    Args:
        float FWHM: The full width at half maximum of the Gaussian, FWHM. The default
            value is 0 and the default unit is the reciprocal of the unit associated
            with the dimension at index `dim_indx`.
        int dim_indx: The index of the CSDM dimension along which the operation is
            applied. The default is the dimension at index 0.
        int dv_indx: The index of the CSDM dependent variable where the operation is
            applied. If the value is None, the operation will be applied to every
            dependent variable.

    Example
    -------

    >>> import mrsimulator.signal_processing.apodization as apo
    >>> operation4 = apo.Gaussian(FWHM=143.4, dim_indx=0, dv_indx=0)
    """

    FWHM: float = 0

    property_unit_types: ClassVar = {"FWHM": ["time", "frequency"]}
    property_default_units: ClassVar = {"FWHM": ["s", "Hz"]}
    property_units: Dict = {"FWHM": "Hz"}

    @staticmethod
    def fn(x, arg):
        # arg is FWHM
        sigma = arg / 2.354820045030949
        return np.exp(-2 * ((x * sigma * np.pi) ** 2))

    def operate(self, data):
        """
        Applies the operation for which the class is named for.

        data: CSDM object
        dep_var: int. The index of the dependent variable to apply operation to
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
#         int dim_indx: Data dimension index to apply the function along.
#         float FWHM: The full width at half maximum parameter, :math:`\Tau`.
#         int dv_indx: Data dependent variable index to apply the function to. If
#             the type None, the operation will be applied to every dependent variable.

#     Example
#     -------

#     >>> operation5 = apo.Exponential(FWHM=143.4, dim_indx=0, dv_indx=0)
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


class Exponential(AbstractApodization):
    r"""Apodize a dependent variable of the simulation data object by an exponential
    function. The function follows

    .. math::
        f(x) = e^{-\Tau |x| \pi},

    where :math:`x` are the coordinates of the data dimension and :math:`\Tau` is
    the width parameter.

    Args:
        float FWHM: The full width at half maximum parameter, :math:`\Tau`. The default
            value is 0 and the default unit is the reciprocal of the unit associated
            with the dimension at index `dim_indx`.
        int dim_indx: The index of the CSDM dimension along which the operation is
            applied. The default is the dimension at index 0.
        int dv_indx: The index of the CSDM dependent variable where the operation is
            applied. If the value is None, the operation will be applied to every
            dependent variable.

    Example
    -------

    >>> operation5 = apo.Exponential(FWHM=143.4, dim_indx=0, dv_indx=0)
    """

    FWHM: float = 0

    property_unit_types: ClassVar = {"FWHM": ["time", "frequency"]}
    property_default_units: ClassVar = {"FWHM": ["s", "Hz"]}
    property_units: Dict = {"FWHM": "Hz"}

    @staticmethod
    def fn(x, arg):
        return np.exp(-arg * np.pi * np.abs(x))

    def operate(self, data):
        """
        Applies the operation for which the class is named for.

        data: CSDM object
        dep_var: int. The index of the dependent variable to apply operation to
        """

        return self._operate(data, fn=self.fn, prop_name="FWHM", prop_value=self.FWHM)
