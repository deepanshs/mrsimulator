"""The Operation class."""
from sys import modules
from typing import ClassVar

import numpy as np
from mrsimulator.utils.parseable import Parseable

__author__ = "Maxwell C. Venetos"
__email__ = "maxvenetos@gmail.com"


class Operation(Parseable):
    """A base class for signal processing operations."""

    module_name: ClassVar[str] = None

    @property
    def function(self):
        return self.__class__.__name__

    def json(self) -> dict:
        """Parse the class object to a JSON compliant python dictionary object, where
        the attribute value with physical quantity is expressed as a string with a
        value and a unit."""
        my_dict = super().json()
        my_dict["function"] = self.function
        if hasattr(self, "type"):
            my_dict["type"] = self.type
        return my_dict

    @staticmethod
    def _get_dv_indexes(indexes, n):
        """Return a list of dependent variable indexes.

        Args:
            indexes: An integer, list of integers, or None indicating the dv indexes.
            n: Total number of dependent variables in the CSDM object.
        """
        if indexes is None:
            return np.arange(n)
        if isinstance(indexes, int):
            return [indexes]
        if isinstance(indexes, (list, tuple)):
            return np.asarray(indexes)

    @classmethod
    def parse_dict_with_units(cls, py_dict: dict):
        """Parse dictionary for SignalProcessor

        Args:
            py_dict (dict): A python dictionary representation of the operation.
        """
        my_dict_copy = py_dict.copy()
        if "function" in my_dict_copy.keys():
            my_dict_copy.pop("function")
        return (
            super().parse_dict_with_units(my_dict_copy)
            if "type" not in my_dict_copy
            else getattr(modules[cls.module_name], my_dict_copy["type"])(**my_dict_copy)
        )

    @staticmethod
    def get_coordinates_in_units(x, unit):
        return x.to(unit).value

    @staticmethod
    def get_coordinates(dim):
        """Return the coordinates of dimension, dim, without equivalence units"""

        def get_coords(dim):
            equivalent_fn, equivalent_unit = dim._equivalencies, dim._equivalent_unit
            dim._equivalencies, dim._equivalent_unit = None, None

            coordinates = dim.coordinates

            dim._equivalencies, dim._equivalent_unit = equivalent_fn, equivalent_unit
            return coordinates

        return (
            get_coords(dim) if not hasattr(dim, "subtype") else get_coords(dim.subtype)
        )


class ModuleOperation(Operation):
    @property
    def type(self):
        """The type baseline function."""
        return self.__class__.__name__
