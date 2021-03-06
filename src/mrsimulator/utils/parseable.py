# -*- coding: utf-8 -*-
"""Base Parseable class."""
from copy import deepcopy
from enum import Enum
from typing import ClassVar
from typing import Dict

from csdmpy.units import string_to_quantity
from pydantic import BaseModel

from .extra import _reduce_dict

# from IPython.display import JSON

__author__ = "Shyam Dwaraknath"
__email__ = "shyamd@lbl.gov"


class Parseable(BaseModel):
    """
    Base class for all objects that can be parsed easily from JSON with units
    Don't directly use this. Rather inherit from it and implement a data model
    and property units and defaults
    """

    property_unit_types: ClassVar = {}

    property_default_units: ClassVar = {}

    property_units: Dict = {}

    @classmethod
    def parse_dict_with_units(cls, json_dict: dict):
        """Parse the physical quantity from a dictionary representation of the class
        object, where the physical quantity is expressed as a string with a number and
        a unit.

        Args:
            dict json_dict: A required python dict object.
        """
        # Only consider properties with both unit types and default units
        all_props = set(cls.property_unit_types.keys()).intersection(
            set(cls.property_default_units.keys())
        )

        property_units = deepcopy(cls.property_default_units)
        for prop in all_props:
            required_type = cls.property_unit_types[prop]
            default_unit = cls.property_default_units[prop]

            if prop in json_dict:
                # If we have a single option
                if isinstance(required_type, str) and isinstance(default_unit, str):
                    _update_json_dict(
                        prop, json_dict, required_type, default_unit, property_units
                    )

                # If there are multiple type/unit combinations
                elif isinstance(required_type, list) and isinstance(default_unit, list):
                    pos_values = [
                        enforce_units(
                            json_dict[prop], r_type, d_unit, throw_error=False
                        )
                        for r_type, d_unit in zip(required_type, default_unit)
                    ]
                    # If none of the units were enforceable, error! else,
                    # choose the first good one
                    if not ([val is not None for val in pos_values]):
                        raise Exception(f"Could not enforce any units on {prop}")

                    json_dict[prop], property_units[prop] = [
                        d for d in zip(pos_values, default_unit) if d[0] is not None
                    ][0]
        for k, v in property_units.items():
            property_units[k] = v[0] if isinstance(v, list) else v
        return cls(**json_dict, property_units=property_units)

    def reduced_dict(self, exclude=["property_units"]) -> dict:
        """Returns a reduced dictionary representation of the class object by removing
        all key-value pair corresponding to keys listed in the `exclude` argument, and
        keys with value as None.

        Args:
            exclude: A list of keys to exclude from the dictionary.
        Return: A dict.
        """
        return _reduce_dict(self.dict(), exclude)

    # def pretty(self):
    #     return JSON(self.json())

    def json(self) -> dict:
        """Parse the class object to a JSON compliant python dictionary object, where
        the attribute value with physical quantity is expressed as a string with a
        number and a unit."""

        temp_dict = {}
        for k, v in self.dict(exclude={"property_units"}).items():

            # check the dict objects
            if isinstance(v, (dict, Enum)):
                val = getattr(self, k).json()
                if val is not None:
                    temp_dict[k] = val

            # check the list objects
            elif isinstance(v, list):
                temp_dict[k] = get_list(getattr(self, k), v)
            elif v not in [None, ""]:
                temp_dict[k] = v

        temp_keys = temp_dict.keys()
        for key, unit in getattr(self, "property_units").items():
            if key in temp_keys:
                u = unit if unit != "pct" else "%"
                temp_dict[key] = f"{temp_dict[key]} {u}"
        return temp_dict


def enforce_units(value: str, required_type: str, default_unit: str, throw_error=True):
    """ Enforces a required type and default unit on the value. """
    try:
        value = string_to_quantity(value)
        data_type = value.unit.physical_type

        if required_type != data_type:
            raise Exception(
                f"A {required_type} value is required but got a {data_type} instead"
            )

        return value.to(default_unit).value
    except Exception as e:
        if throw_error:
            raise e
        return None


def get_list(member, obj):
    lst = []
    for i, item in enumerate(obj):
        if isinstance(item, list):
            lst.append(get_list(member[i], item))
        elif isinstance(item, (dict, Enum)):
            lst.append(member[i].json())
        elif item not in [None, ""]:
            lst.append(item)
    return lst


def _update_json_dict(prop, json_dict, required_type, default_unit, property_units):
    try:
        json_dict[prop] = enforce_units(json_dict[prop], required_type, default_unit)
        property_units[prop] = default_unit
    except Exception as e:
        raise Exception(f"Error enforcing units for {prop}: {json_dict[prop]}\n{e}")
