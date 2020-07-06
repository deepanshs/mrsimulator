# -*- coding: utf-8 -*-
"""The Event class."""
from copy import deepcopy
from sys import modules
from typing import ClassVar
from typing import Dict
from typing import Union

import numpy as np

from ._base import abstractOperation


class AbstractApodization(abstractOperation):
    dim_indx: int = 0
    dep_var_indx: Union[int, list, tuple] = None  # if none apply to all

    @classmethod
    def parse_dict_with_units(cls, py_dict):
        copy_dict = deepcopy(py_dict)
        copy_dict.pop("type")

        obj = super().parse_dict_with_units(copy_dict)
        return getattr(modules[__name__], py_dict["type"])(**obj.dict())

    @property
    def function(self):
        return "apodization"

    @property
    def type(self):
        return self.__class__.__name__

    def set_property_units(self, unit, prop):
        """
        Populate the property unit attribute of the class based on the dimension unit.
        """
        if unit.physical_type == "frequency":
            self.property_units = dict(prop="s")
            return
        self.property_units = dict(prop="Hz")

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

    def operate(self, data, fn, prop_name, prop_value):
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
        dv_indexes = self._get_dv_indexes(self.dep_var_indx, n=n)

        for i in dv_indexes:
            data.dependent_variables[i].components[0] *= apodization_vactor
        return data


class Gaussian(AbstractApodization):
    r"""
    Class for applying a Gaussian function
    to a dependent variable of simulation data

    .. math::
        f(\vec{x}) = \vec{x}*e^{-2*(\vec{x} * \sigma * \pi)^2}

    Args:
        dim_indx: int. Data dimension to apply the function along.
        sigma: float. Standard deviation of Gaussian function
        dep_var_indx: int. Data dependent variable index to apply the function to.
            If type None, will be applied to every dependent variable.

    """

    sigma: float = 0

    property_unit_types: ClassVar = {"sigma": ["time", "frequency"]}
    property_default_units: ClassVar = {"sigma": ["s", "Hz"]}
    property_units: Dict = {"sigma": "Hz"}

    @staticmethod
    def fn(x, arg):
        return np.exp(-((x * arg * np.pi) ** 2) * 2)

    def operate(self, data):
        """
        Applies the operation for which the class is named for.

        data: CSDM object
        dep_var: int. The index of the dependent variable to apply operation to
        """

        return super().operate(
            data, fn=self.fn, prop_name="sigma", prop_value=self.sigma
        )


class Exponential(AbstractApodization):
    r"""
    Class for applying an exponential
    Lorentzian function to a dependent variable of simulation data

    .. math::
        f(\vec{x}) = \vec{x}*e^{-\Lambda * \abs{\vec{x}} * \pi)}

    Args:
        dim_indx: int. Data dimension to apply the function along.
        Lambda: float. Width parameter
        dep_var_indx: int. Data dependent variable index to apply the function to.
            If type None, will be applied to every dependent variable.

    """

    Lambda: float = 0

    property_unit_types: ClassVar = {"Lambda": ["time", "frequency"]}
    property_default_units: ClassVar = {"Lambda": ["s", "Hz"]}
    property_units: Dict = {"Lambda": "Hz"}

    @staticmethod
    def fn(x, arg):
        return np.exp(-arg * np.pi * np.abs(x))

    def operate(self, data):
        """
        Applies the operation for which the class is named for.

        data: CSDM object
        dep_var: int. The index of the dependent variable to apply operation to
        """

        return super().operate(
            data, fn=self.fn, prop_name="Lambda", prop_value=self.Lambda
        )
