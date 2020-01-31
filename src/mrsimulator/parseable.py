# -*- coding: utf-8 -*-
"""Base Parseable class."""
from typing import ClassVar
from typing import Dict

from csdmpy.units import string_to_quantity
from pydantic import BaseModel

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
    def parse_dict_with_units(cls, json_dict):
        # Only consider properties with both unit types and default units
        all_props = set(cls.property_unit_types.keys()).intersection(
            set(cls.property_default_units.keys())
        )

        property_units = {}
        for prop in all_props:
            required_type = cls.property_unit_types[prop]
            default_unit = cls.property_default_units[prop]

            if prop in json_dict:
                # If we have a single option
                if isinstance(required_type, str) and isinstance(default_unit, str):
                    try:
                        json_dict[prop] = enforce_units(
                            json_dict[prop], required_type, default_unit
                        )
                        property_units[prop] = default_unit
                    except Exception as e:
                        raise Exception(
                            f"Error enforcing units for {prop}: {json_dict[prop]}\n{e}"
                        )

                # If there are multiple type/unit combinations
                elif isinstance(required_type, list) and isinstance(default_unit, list):
                    pos_values = [
                        enforce_units(
                            json_dict[prop], r_type, d_unit, throw_error=False
                        )
                        for r_type, d_unit in zip(required_type, default_unit)
                    ]
                    # If none of the units were enforceable, error
                    # else choose the first good one
                    if not ([val is not None for val in pos_values]):
                        raise Exception(f"Could not enforce any units on {prop}")
                    else:
                        json_dict[prop], property_units[prop] = [
                            d for d in zip(pos_values, default_unit) if d[0] is not None
                        ][0]
        return cls(**json_dict, property_units=property_units)


def enforce_units(value: str, required_type: str, default_unit: str, throw_error=True):
    """
    Enforces a required type and default unit on the value

        value
    """
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
        else:
            return None
