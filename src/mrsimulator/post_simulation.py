# -*- coding: utf-8 -*-
"""The Event class."""
from typing import ClassVar
from typing import Dict
from typing import List
from typing import Union

import csdmpy as cp
import numpy as np
from mrsimulator.util.parseable import Parseable
from pydantic import BaseModel


class abstractOperation(Parseable):
    """
    A base class for signal processing operations
    """

    @property
    def function(self):
        return self.__class__.__name__

    def to_dict_with_units(self):
        my_dict = super().to_dict_with_units()
        my_dict["function"] = self.__class__.__name__
        return my_dict

    @classmethod
    def parse_dict_with_units(cls, py_dict):
        my_dict_copy = py_dict.copy()
        if "function" in my_dict_copy.keys():
            my_dict_copy.pop("function")
        return super().parse_dict_with_units(my_dict_copy)

    def operate(self, data):
        return data


class Scale(abstractOperation):
    """
    Class for applying a scaling
    factor to a dependent variable of simulation data
    """

    factor: float = 1

    def operate(self, data, dep_var):
        """
        Applies the operation for which the class is named for.

        $f(\vec(x)) = scale*\vec(x)$

        data: CSDM object
        dep_var: int. The index of the dependent variable to apply operation to
        """
        data.dependent_variables[dep_var].components[0] *= self.factor
        return data


class Gaussian(abstractOperation):
    """
    Class for applying a Gaussian function
    to a dependent variable of simulation data

    $f(\vec{x}) = \vec{x}*e^{-2*(\vec{x} * \sigma * \pi)^2}$

    dimension: int. Data dimension to apply the function along
    sigma: float. Standard deviation of Gaussian function
    """

    dimension: int = 0
    sigma: float = 0

    property_unit_types: ClassVar = {"sigma": ["time", "frequency"]}
    property_default_units: ClassVar = {"sigma": ["s", "Hz"]}
    property_units: Dict = {"sigma": ["s", "Hz"]}

    def operate(self, data, dep_var):
        """
        Applies the operation for which the class is named for.

        data: CSDM object
        dep_var: int. The index of the dependent variable to apply operation to
        """
        # unit = self.property_units["sigma"]
        # data_copy.dimensions[self.dimension].coordinates.to(unit).value
        x = data.dimensions[self.dimension].coordinates.value
        data.dependent_variables[dep_var].components[0] *= np.exp(
            -((x * self.sigma * np.pi) ** 2) * 2
        )
        return data


class Exponential(abstractOperation):
    """
    Class for applying an exponential
    Lorentzian function to a dependent variable of simulation data

    $f(\vec{x}) = \vec{x}*e^{-\Lambda * \abs{\vec{x}} * \pi)}$

    dimension: int. Data dimension to apply the function along
    Lambda: float. Width parameter
    """

    dimension: int = 0
    Lambda: float = 0

    property_unit_types: ClassVar = {"Lambda": ["time", "frequency"]}
    property_default_units: ClassVar = {"Lambda": ["s", "Hz"]}
    property_units: Dict = {"Lambda": ["s", "Hz"]}

    def operate(self, data, dep_var):
        """
        Applies the operation for which the class is named for.

        data: CSDM object
        dep_var: int. The index of the dependent variable to apply operation to
        """
        # unit = self.property_units["Lambda"]
        # x = data.dimensions[self.dimension].coordinates.to(unit).value
        x = data.dimensions[self.dimension].coordinates.value
        data.dependent_variables[dep_var].components[0] *= np.exp(
            -self.Lambda * np.pi * np.abs(x)
        )
        return data


class IFFT(abstractOperation):
    """
    Class for applying an inverse
    Fourier transform to a dependent variable of simulation data

    dimension: int. Data dimension to apply the function along

    """

    dimension: int = 0

    def operate(self, data, dep_var):
        """
        Applies the operation for which the class is named for.

        data: CSDM object
        dep_var: int. The index of the dependent variable to apply operation to
        """
        return data.fft(axis=self.dimension)


class FFT(abstractOperation):
    """
    Class for applying a
    Fourier transform to a dependent variable of simulation data

    dimension: int. Data dimension to apply the function along

    """

    dimension: int = 0

    def operate(self, data, dep_var):
        """
        Applies the operation for which the class is named for.

        data: CSDM object
        dep_var: int. The index of the dependent variable to apply operation to
        """
        return data.fft(axis=self.dimension)


class OperationList(BaseModel):
    """
    Container class to hold a list of operations to be applied to a given
    dependent variable in the data

    dependent_variable: The dependent variable of the data to apply operations to
    operations: List of abstractOperation based operations to sequentially apply to data.
    """

    dependent_variable: int = 0
    operations: List[abstractOperation]


class SignalProcessor(BaseModel):
    """
    Signal processing class to apply lists of various operations to individual dependent variables
    of the data.

    data: CSDM object. From simulation
    operations: List of operation lists
    """

    data: Union[cp.CSDM, np.ndarray] = None
    operations: List[OperationList]

    class Config:
        validate_assignment = True
        arbitrary_types_allowed = True

    def apply_operations(self, **kwargs):
        op_list = self.operations

        copy_data = self.data.copy()

        for item in op_list:
            dep_var = item.dependent_variable
            for filters in item.operations:
                copy_data = filters.operate(copy_data, dep_var=dep_var)

        self.data = copy_data

        return copy_data
