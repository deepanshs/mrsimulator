# -*- coding: utf-8 -*-
"""The Event class."""
from typing import ClassVar
from typing import Dict

import numpy as np

from ._base import abstractOperation


class AbstractApodization(abstractOperation):
    @property
    def function(self):
        return "apodization"

    @property
    def type(self):
        return self.__class__.__name__


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

    dim_indx: int = 0
    sigma: float = 0
    dep_var_indx: int = None  # if none apply to all

    property_unit_types: ClassVar = {"sigma": ["time", "frequency"]}
    property_default_units: ClassVar = {"sigma": ["s", "Hz"]}
    property_units: Dict = {"sigma": "Hz"}

    def operate(self, data):
        """
        Applies the operation for which the class is named for.

        data: CSDM object
        dep_var: int. The index of the dependent variable to apply operation to
        """
        # unit = self.property_units["sigma"]
        # data_copy.dimensions[self.dimension].coordinates.to(unit).value
        x = data.dimensions[self.dim_indx].coordinates
        x_value, x_unit = x.value, x.unit
        if x_unit.physical_type == "frequency":
            self.property_units = dict(sigma="s")
        else:
            self.property_units = dict(sigma="Hz")

        if self.dep_var_indx is None:
            for i in range(len(data.dependent_variables)):
                data.dependent_variables[i].components[0] *= np.exp(
                    -((x_value * self.sigma * np.pi) ** 2) * 2
                )

        else:
            data.dependent_variables[self.dep_var_indx].components[0] *= np.exp(
                -((x_value * self.sigma * np.pi) ** 2) * 2
            )
        return data


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

    dim_indx: int = 0
    Lambda: float = 0
    dep_var_indx: int = None  # if none apply to all

    property_unit_types: ClassVar = {"Lambda": ["time", "frequency"]}
    property_default_units: ClassVar = {"Lambda": ["s", "Hz"]}
    property_units: Dict = {"Lambda": "Hz"}

    def operate(self, data):
        """
        Applies the operation for which the class is named for.

        data: CSDM object
        dep_var: int. The index of the dependent variable to apply operation to
        """
        # unit = self.property_units["Lambda"]
        # x = data.dimensions[self.dimension].coordinates.to(unit).value
        x = data.dimensions[self.dim_indx].coordinates
        x_value, x_unit = x.value, x.unit
        if x_unit.physical_type == "frequency":
            self.property_units = dict(Lambda="s")
        else:
            self.property_units = dict(Lambda="Hz")

        if self.dep_var_indx is None:
            for i in range(len(data.dependent_variables)):
                data.dependent_variables[i].components[0] *= np.exp(
                    -self.Lambda * np.pi * np.abs(x_value)
                )
        else:
            data.dependent_variables[self.dep_var_indx].components[0] *= np.exp(
                -self.Lambda * np.pi * np.abs(x_value)
            )
        return data
