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
    abstractOperation-based class for applying a scaling
    factor to a dependent variable of simulation data
    """

    factor: float = 1

    def operate(self, data, dep_var):
        """
        Applies the operation for which the class is named for.

        data: CSDM object
        dep_var: int. The index of the dependent variable to apply operation to
        """
        data_copy = data.copy()
        data_copy.dependent_variables[dep_var].components[0] *= self.factor
        return data_copy


class Gaussian(abstractOperation):
    """
    abstractOperation-based class for applying a Gaussian function
    to a dependent variable of simulation data

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
        data_copy = data.copy()
        # data_copy.dimensions[self.dimension].coordinates.to(unit).value
        x = data_copy.dimensions[self.dimension].coordinates.value
        data_copy.dependent_variables[dep_var].components[0] *= np.exp(
            -((x * self.sigma * np.pi) ** 2) * 2
        )
        return data_copy


class Exponential(abstractOperation):
    """
    abstractOperation-based class for applying an exponential
    Lorentzian function to a dependent variable of simulation data

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
        data_copy = data.copy()
        # x = data.dimensions[self.dimension].coordinates.to(unit).value
        x = data.dimensions[self.dimension].coordinates.value
        data_copy.dependent_variables[dep_var].components[0] *= np.exp(
            -self.Lambda * np.pi * np.abs(x)
        )
        return data_copy


class IFFT(abstractOperation):
    """
    abstractOperation-based class for applying an inverse
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
        data_copy = data.copy()
        return data_copy.fft(axis=self.dimension)


class FFT(abstractOperation):
    """
    abstractOperation-based class for applying a
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
        data_copy = data.copy()
        return data_copy.fft(axis=self.dimension)


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

        for item in op_list:
            dep_var = item.dependent_variable
            for filters in item.operations:
                self.data = filters.operate(self.data, dep_var=dep_var)

        return self.data
