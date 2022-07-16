"""Base Parseable class."""
import warnings
from copy import copy
from enum import Enum
from typing import ClassVar
from typing import Dict

from csdmpy.units import string_to_quantity
from pydantic import BaseModel


__author__ = "Shyam Dwaraknath"
__email__ = "shyamd@lbl.gov"

INCLUDE_LIST = [
    "count",
    "magnetic_flux_density",
    "rotor_angle",
    "rotor_frequency",
    "spectral_width",
    "dim_index",
    "dv_index",
    "duration",
    "site_index",
    "isotope",
    "isotropic_chemical_shift",
    "type",
    "function",
    "config",
    "decompose_spectrum",
    "integration_density",
    "integration_volume",
    "number_of_sidebands",
    "name",
    "description",
    "transition_queries",
    "ch1",
    "P",
    "channels",
]

CONST = string_to_quantity("1")


class Parseable(BaseModel):
    """Base class for all objects that can be parsed easily from JSON with units
    Don't directly use this. Rather inherit from it and implement a data model
    and property units and defaults
    """

    name: str = None
    description: str = None
    label: str = None
    property_unit_types: ClassVar[Dict] = {}
    property_default_units: ClassVar[Dict] = {}
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

        property_units = copy(cls.property_default_units)
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
                    if not any([val is not None for val in pos_values]):
                        raise Exception(f"Could not enforce any units on {prop}")

                    json_dict[prop], property_units[prop] = [
                        d for d in zip(pos_values, default_unit) if d[0] is not None
                    ][0]

        for k, v in property_units.items():
            property_units[k] = v[0] if isinstance(v, list) else v
        return cls(**json_dict, property_units=property_units)

    def reduced_dict(self, exclude={}) -> dict:
        """Returns a reduced dictionary representation of the class object by removing
        all key-value pair corresponding to keys listed in the `exclude` argument, and
        keys with value as None.

        Args:
            exclude: A list of keys to exclude from the dictionary.
        Return: A dict.
        """
        warnings.warn(
            "reduced_dict() is deprecated in v0.7, use json(units=True) instead.",
            category=DeprecationWarning,
        )
        return self.json(exclude, units=False)

    def json(self, exclude={}, units=True) -> dict:
        """Parse the class object to a JSON compliant python dictionary object.

        Args:
            exclude: Set of keys that will be excluded from the result.
            units: If true, the attribute value is a physical quantity expressed as a
                string with a number and a unit, else a float.

        Returns: dict
        """
        temp_dict = {}
        for k, v in self.dict(exclude={"property_units", *exclude}).items():
            attr_val = getattr(self, k)
            if k not in INCLUDE_LIST:
                if attr_val == get_default_class_value(self.__class__, k):
                    continue

            # check the dict objects
            if isinstance(v, (dict, Enum)):
                val = attr_val.json(units=units)
                _ = None if val in [None, {}] else temp_dict.update({k: val})

            # check the list objects
            elif isinstance(v, list):
                val = [
                    item if not hasattr(item, "json") else item.json(units=units)
                    for item in attr_val
                ]
                _ = None if val == [] else temp_dict.update({k: val})

            elif v is not None and v != "":
                temp_dict[k] = v

        return temp_dict if not units else self.clear_property_units(temp_dict)

    def clear_property_units(self, temp_dict):
        # if not hasattr(self, "property_units"):
        #     return temp_dict
        temp_keys = temp_dict.keys()
        for key, unit in getattr(self, "property_units").items():
            if key in temp_keys:
                u = unit if unit != "pct" else "%"
                var = f" {u}" if unit != CONST else ""
                temp_dict[key] = f"{temp_dict[key]}{var}"
        return temp_dict


def get_default_class_value(obj, k):
    return (
        getattr(obj(), k)
        if not hasattr(obj, "test_vars")
        else getattr(obj(**obj.test_vars), k)
    )


def enforce_units(value: str, required_type: str, default_unit: str, throw_error=True):
    """Enforces a required type and default unit on the value."""
    try:
        value = string_to_quantity(value)
        data_type = str(value.unit.physical_type).split("/")[0]

        if required_type != data_type:
            raise Exception(
                f"A {required_type} value is required but got a {data_type} instead"
            )

        return value.to(default_unit).value
    except Exception as e:
        if throw_error:
            raise e
        return None


def _update_json_dict(prop, json_dict, required_type, default_unit, property_units):
    try:
        json_dict[prop] = enforce_units(json_dict[prop], required_type, default_unit)
        property_units[prop] = default_unit
    except Exception as e:
        raise Exception(f"Error enforcing units for {prop}: {json_dict[prop]}\n{e}")
