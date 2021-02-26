# -*- coding: utf-8 -*-
"""The AbstractOperation class."""
import numpy as np
from mrsimulator.utils.parseable import Parseable

__author__ = "Maxwell C. Venetos"
__email__ = "maxvenetos@gmail.com"


class AbstractOperation(Parseable):
    """A base class for signal processing operations."""

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
            indexes: An interger, list of integers, or None indicating the dv indexes.
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
        return super().parse_dict_with_units(my_dict_copy)
