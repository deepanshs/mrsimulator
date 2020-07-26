# -*- coding: utf-8 -*-
"""The AbstractOperation class."""
from mrsimulator.utils.parseable import Parseable


__author__ = "Maxwell C. Venetos"
__email__ = "maxvenetos@gmail.com"


class AbstractOperation(Parseable):
    """A base class for signal processing operations."""

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
        """Parse dictionary for SignalProcessor

        Args:
            py_dict (dict): A python dictionary representation of the operation.
        """
        my_dict_copy = py_dict.copy()
        if "function" in my_dict_copy.keys():
            my_dict_copy.pop("function")
        return super().parse_dict_with_units(my_dict_copy)
