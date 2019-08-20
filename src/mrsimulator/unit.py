# -*- coding: utf-8 -*-
"""
Unit library parse extension.

The methods here parse a string containing physical
quantities to numerical value and unit. The unit
manipulation and conversion is done using
Astropy unit library.
"""
from astropy import units as u
from astropy.units import cds
from numpy import inf

# from astropy.units import UnitConversionError
# from astropy.units import cds
# cds.enable()

__author__ = "Deepansh J. Srivastava"
__email__ = ["srivastava.89@osu.edu", "deepansh2012@gmail.com"]

u.add_enabled_units([cds.ppm])
_tr = u.def_unit(["tr", "turn", "cycle", "revolution"], 1 * u.Unit(1))

convert = {
    "Å": "Angstrom",
    "°C": "deg_C",
    "°F": "deg_F",
    "°": "deg",
    "µ": "u",
    "ℏ": "/h",
    "Ω": "Ohm",
    # "ppm" : _ppm,
}


def display_unit(unit):
    """
    Return unit in UTF-8 characters.

    :returns: ``unit`` object
    """
    unit = str(unit)
    for key in convert:
        unit = unit.replace(convert[key], key)
    return unit


# dimensionless_frequency_ratio = [(
#     u.Hz, _ppm, lambda x: 1000.0 * x, lambda x: x / 1000.0
# )]


def string_to_unit(unit):
    """
    Parse the string and return a ``unit`` object.

    The string must only contain a unit.

    :returns: ``unit`` object
    """
    for key in convert:
        unit = unit.replace(key, convert[key])
    # print (unit)
    unit_qt = u.Unit(unit)
    return unit_qt


def is_physical_quantity(string):
    try:
        string_to_quantity(string)
        return True
    except BaseException:
        return False


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


def value_object_format(quantity, numerical_value=True):
    """Convert unit to value object."""
    # mode = 'fits'
    string = quantity.unit.to_string("fits").strip()
    # print('string', string)
    for key in convert.keys():
        string = string.replace(convert[key], key)
    lst = {"10**-6": "ppm", "10**-2": "%"}
    for key in lst.keys():
        string = string.replace(key, lst[key])
    # string = string.replace("10**-6", "ppm")
    # string = string.replace("10**-2", "%")
    # print ('string', string)

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
    # string = string.replace('* ( *', '(')
    # string = string.replace('* ) *', ')')
    return string


def unit_to_latex(unit):
    """NotImplemented."""
    # mode = 'fits'
    string = unit.to_string("fits").strip()
    # print('string', string)
    convert_tex = {
        "Angstrom": "\\AA",
        "deg_C": "$^\\circ$C",
        "deg_F": "$^\\circ$C",
        "deg": "$^\\circ$",
        "u": "$\\mu$",
        "/h": "$\\hbar$",
        "Ohm": "$\\Ohm$",
        "10**-6": "ppm",
    }
    for key in convert_tex:
        string = string.replace(key, convert_tex[key])
    lst = {"10**-6": "ppm", "10**-2": "%"}
    for key in lst.keys():
        string = string.replace(key, lst[key])
    # print (string)

    cat_string = []
    subunits = string.split(" ")
    for item in subunits:
        power = False
        if item.find("-") == -1:
            for i, c in enumerate(item):
                if c.isnumeric():
                    if i == 0:
                        break
                    cat_string.append(item[:i] + "$^{" + item[i:] + "}$ ")
                    power = True
                    break
            if not power:
                cat_string.append(item + " ")
        else:
            l, r = item.split("-")
            cat_string.append(l + "$^{-" + r + "}$ ")
    string = " ".join(cat_string)
    # string = string.replace('* / *', '/')
    # string = string.replace('* ( *', '(')
    # string = string.replace('* ) *', ')')
    return string
