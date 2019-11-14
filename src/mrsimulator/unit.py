# -*- coding: utf-8 -*-
"""
Unit library parse extension.

The methods here parse a string containing physical
quantities to numerical value and unit. The unit
manipulation and conversion is done using
Astropy unit library.
"""
import warnings

from astropy import units as u
from astropy.units import cds
from astropy.units import Quantity
from numpy import inf

__author__ = "Deepansh J. Srivastava"
__email__ = "srivastava.89@osu.edu"
__all__ = ["string_to_quantity", "ScalarQuantity"]


tr = u.def_unit(["tr", "turn", "revolution"], 1.0 * u.cycle)
u.add_enabled_units([cds.ppm, tr])

convert = {
    "Å": "Angstrom",
    "°C": "deg_C",
    "°F": "deg_F",
    "°": "deg",
    "µ": "u",
    "ℏ": "/h",
    "Ω": "Ohm",
}


def string_to_quantity(string, dtype=float):
    """
    Parse the string and return a ``Quantity`` object.

    The string must contain a physical quantity. The method parse the string
    for numerical value and unit.

    :returns: ``Quantity`` object
    """
    numeric = "0123456789-+.eE*/j^ ()"
    string = string.strip() + " "

    for i, c in enumerate(string):
        if c not in numeric:
            break

    j = 1
    for j in range(1, i + 1):
        if string[:i][-j] == "(":
            break
        if string[:i][-j] == ")":
            j = j - 1
            break
    if i - j == 0:
        index = i
    else:
        index = i - j

    if index != -1:
        try:
            number = eval(string[:index])
        except ZeroDivisionError:
            number = inf
        except Exception as e:
            raise ValueError(e)
    else:
        index = 0
        number = 1.0

    unit = string[index:].strip()
    if unit != "" and unit != "()":
        if unit[0] == "(" and unit[-1] == ")":
            unit = unit[1:-1]

    for key in convert:
        unit = unit.replace(key, convert[key])
        unit_multiplier = 1
    try:
        unit_qt = u.Unit(unit) * unit_multiplier
        analysis = dtype(number) * unit_qt
        return analysis
    except BaseException as e:
        raise BaseException(e)


def scalar_quantity_format(quantity, numerical_value=True):
    """Convert unit to value object."""
    string = quantity.unit.to_string("fits").strip()
    for key in convert.keys():
        string = string.replace(convert[key], key)
    lst = {"10**-6": "ppm", "10**-2": "%"}
    for key in lst.keys():
        string = string.replace(key, lst[key])

    if numerical_value:
        cat_string = [str(quantity.value)]
    else:
        cat_string = [""]
    subunits = string.split(" ")
    for item in subunits:
        power = False
        if item.find("-") == -1:
            for i, c in enumerate(item):
                if c.isnumeric():
                    if i == 0:
                        break
                    cat_string.append(item[:i] + "^" + item[i:] + " *")
                    power = True
                    break
            if not power:
                cat_string.append(item + " *")
        else:
            l, r = item.split("-")
            cat_string.append(l + "^-" + r + " *")
    string = " ".join(cat_string)[:-2]
    string = string.replace("* / *", "/")
    return string.strip()


# def unit_to_latex(unit):
#     """NotImplemented."""
#     string = unit.to_string("fits").strip()
#     convert_tex = {
#         "Angstrom": "\\AA",
#         "deg_C": "$^\\circ$C",
#         "deg_F": "$^\\circ$C",
#         "deg": "$^\\circ$",
#         "u": "$\\mu$",
#         "/h": "$\\hbar$",
#         "Ohm": "$\\Ohm$",
#         "10**-6": "ppm",
#     }
#     for key in convert_tex:
#         string = string.replace(key, convert_tex[key])
#     lst = {"10**-6": "ppm", "10**-2": "%"}
#     for key in lst.keys():
#         string = string.replace(key, lst[key])

#     cat_string = []
#     subunits = string.split(" ")
#     for item in subunits:
#         power = False
#         if item.find("-") == -1:
#             for i, c in enumerate(item):
#                 if c.isnumeric():
#                     if i == 0:
#                         break
#                     cat_string.append(item[:i] + "$^{" + item[i:] + "}$ ")
#                     power = True
#                     break
#             if not power:
#                 cat_string.append(item + " ")
#         else:
#             l, r = item.split("-")
#             cat_string.append(l + "$^{-" + r + "}$ ")
#     string = " ".join(cat_string)
#     return string


class ScalarQuantity:

    """
    A ScalarQuantity class compliant with the CSDM standards.

    @params: quantity: Return the quantity object from astropy.units
                       library.
    """

    __slots__ = "quantity"

    def __init__(self, quantity_string=None, unit=None):
        self.quantity = self.quantity_object(quantity_string, unit)

    @classmethod
    def quantity_object(cls, quantity_string, unit):
        if isinstance(quantity_string, ScalarQuantity):
            return quantity_string.quantity

        if isinstance(quantity_string, Quantity):
            return quantity_string

        lst = [None, ""]
        if quantity_string in lst and unit in lst:
            return 0.0 * u.Unit("")

        if quantity_string in lst and unit not in lst:
            return 0.0 * unit

        if isinstance(quantity_string, str):
            quantity = string_to_quantity(quantity_string)
            if unit is not None:
                return check_unit_consistency(quantity, unit)
            else:
                return quantity

    def __str__(self):
        return str(self.quantity)

    def format(self, format="quantity"):
        """Format the units according to csdmpy recommendation.

        :ivar: format: The value can either be 'quantity' or 'unit'

        .. doctest::

            >>> a = ScalarQuantity('5 kg m^2 /s')
            >>> print(a.quantity)
            5.0 kg m2 / s
            >>> print(ScalarQuantity(a.quantity).format())
            5.0 kg * m^2 * s^-1
            >>> print(ScalarQuantity(a.quantity).format('unit'))
            kg * m^2 * s^-1
        """

        element = _default_units(self.quantity)
        if format == "unit":
            return scalar_quantity_format(element, numerical_value=False)
        if format == "quantity":
            return scalar_quantity_format(element)


def check_unit_consistency(element, unit):
    if isinstance(unit, str):
        unit = ScalarQuantity(unit).quantity.unit
    if element.unit.physical_type != unit.physical_type:
        options = [
            str(element.unit),
            str(element.unit.physical_type),
            str(unit),
            unit.physical_type,
        ]
        message = (
            "Validation Failed: The unit '{0}' ({1}) is inconsistent "
            "with the unit '{2}' ({3})."
        )
        raise Exception(message.format(*options))
    return element


def _default_units(element):
    if element.unit.physical_type == "frequency":
        element = element.to("Hz")
    return element


def check_quantity_name(element, unit):
    if element is None:
        element = unit.physical_type
        return element

    if unit.physical_type == "unknown":
        warnings.warn(
            (
                "The physical quantity name associated with the unit, {0}, "
                "is not defined in astropy.units package. Continuing with "
                "'{1}' as the physical quantity name."
            ).format(str(unit), element)
        )
        return element

    if element.lower() != unit.physical_type:
        warnings.warn(
            (
                "The physical quantity name, '{0}', is not defined in the "
                "astropy.units package. Continuing with '{0}' as the physical "
                "quantity name for unit {1}."
            ).format(element, str(unit))
        )

    return element.lower()
