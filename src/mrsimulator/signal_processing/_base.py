# -*- coding: utf-8 -*-
"""The abstractOperation class."""
from mrsimulator.util.parseable import Parseable


class abstractOperation(Parseable):
    """
    A base class for signal processing operations
    """

    @property
    def function(self):
        return self.__class__.__name__

    def to_dict_with_units(self):
        my_dict = super().to_dict_with_units()
        my_dict["function"] = self.function
        if hasattr(self, "type"):
            my_dict["type"] = self.type
        return my_dict

    @classmethod
    def parse_dict_with_units(cls, py_dict):
        my_dict_copy = py_dict.copy()
        if "function" in my_dict_copy.keys():
            my_dict_copy.pop("function")
        return super().parse_dict_with_units(my_dict_copy)


# class FFT(abstractOperation):
#     """
#     Class for applying a
#     Fourier transform to a dependent variable of simulation data

#     dimension: int. Data dimension to apply the function along

#     """

#     dimension: int = 0

#     def operate(self, data):
#         """
#         Applies the operation for which the class is named for.

#         data: CSDM object
#         dep_var: int. The index of the dependent variable to apply operation to
#         """
#         return data.fft(axis=self.dimension)


# class OperationList(BaseModel):
#     """
#     Container class to hold a list of operations to be applied to a given
#     dependent variable in the data

#     dependent_variable: The dependent variable of the data to apply operations to
#     operations: List of abstractOperation based operations to sequentially apply to data.
#     """

#     dependent_variable: int = 0
#     operations: List[abstractOperation]
