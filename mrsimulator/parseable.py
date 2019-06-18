from pydantic import BaseModel

__author__ = "Shyam Dwaraknath"
__email__ = "shyamd@lbl.gov"


class Parseable(BaseModel):
    """
    Base class for all objects that can be parsed easily from JSON with units
    Don't directly use this. Rather inherit from it and implement a data model
    and property units and defaults
    """

    number_of_points: int = 1024
    spectral_width: float = 100
    reference_offset: float = 0
    magnetic_flux_density: float = 9.4
    rotor_frequency: float = 0
    rotor_angle: float = 0.9553  # 54.935 degrees in radians
    rotor_phase: float = 0
    nucleus: str = "1H"
    spin: int = 1
    natural_abundance: float = 0.04683
    gyromagnetic_ratio: float = -8.465

    __property_unit_types = {}

    __property_default_units = {}

    @classmethod
    def parse_json_with_units(cls, json_dict):

        # Only consider properties with both unit types and default units
        all_props = set(cls.__property_unit_types.keys()).intersection(
            cls.__property_default_units.keys()
        )

        for prop in all_props:
            required_type = cls.__property_unit_types[prop]
            default_unit = cls.__property_default_units[prop]

            if prop in json_dict:
                # If we have a single option
                if isinstance(required_type, str) and isinstance(
                    default_unit, str
                ):
                    json_dict[prop] = enforce_units(
                        json_dict[prop], required_type, default_unit
                    )
                # If there are multiple type/unit combinations
                elif isinstance(required_type, list) and isinstance(
                    default_unit, list
                ):
                    pos_values = [
                        enforce_units(
                            json_dict[prop], required_type, default_unit
                        )
                        for r_type, d_unit in zip(required_type, default_unit)
                    ]
                    # If none of the units were enforceable, error
                    # else choose the first good one
                    if not any(pos_values):
                        raise Exception(
                            f"Could not enforce any units on {prop}"
                        )
                    else:
                        json_dict[prop] = list(
                            filter(None.__ne__, pos_values)
                        )[0]


def enforce_units(value, required_type, default_unit, throw_error=True):
    """
    Enforces a required type and default unit on the value
    """

    value = string_to_quantity(value)
    data_type = value.unit.physical_type

    if required_type != data_type:
        if throw_error:
            raise Exception(
                f"A {required_type} value is required for '{prop}' but got a {data_type} instead"
            )
        else:
            return None

    return value.to(default_unit).value
