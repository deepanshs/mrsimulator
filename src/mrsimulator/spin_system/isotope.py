# -*- coding: utf-8 -*-
"""Base Isotope class."""
from os import path
from re import match

from monty.serialization import loadfn
from pydantic import BaseModel
from pydantic import validator

__author__ = "Deepansh Srivastava"
__email__ = "srivastava.89@osu.edu"

MODULE_DIR = path.dirname(path.abspath(__file__))
ISOTOPE_DATA = loadfn(path.join(MODULE_DIR, "isotope_data.json"))


class Isotope(BaseModel):
    """The Isotope class.

    Attributes
    ----------

    symbol: str (required)
        The isotope symbol given as the atomic number followed by the atomic symbol.

    Example
    -------

    >>> # 13C isotope information
    >>> carbon = Isotope(symbol='13C')
    >>> carbon.spin
    0.5
    >>> carbon.natural_abundance # in %
    1.11
    >>> carbon.gyromagnetic_ratio # in MHz/T
    10.7084
    >>> carbon.atomic_number
    6
    >>> carbon.quadrupole_moment # in eB
    0.0

    """

    symbol: str

    class Config:
        validate_assignment = True

    @validator("symbol", always=True)
    def get_isotope(cls, v, *, values, **kwargs):
        return format_isotope_string(v)

    def json(self) -> dict:
        """Parse the class object to a JSON compliant python dictionary object, where
        the attribute value with physical quantity is expressed as a string with a
        value and a unit."""
        return self.symbol

    @property
    def spin(self):
        """Spin quantum number, I, of the isotope."""
        isotope_data = get_isotope_data(self.symbol)
        return isotope_data["spin"] / 2.0

    @property
    def natural_abundance(self):
        """Natural abundance of the isotope in units of %."""
        isotope_data = get_isotope_data(self.symbol)
        return isotope_data["natural_abundance"]

    @property
    def gyromagnetic_ratio(self):
        """Reduced gyromagnetic ratio of the nucleus given in units of MHz/T."""
        isotope_data = get_isotope_data(self.symbol)
        return isotope_data["gyromagnetic_ratio"]

    @property
    def quadrupole_moment(self):
        """Quadrupole moment of the nucleus given in units of eB (electron-barn)."""
        isotope_data = get_isotope_data(self.symbol)
        return isotope_data["quadrupole_moment"]

    @property
    def atomic_number(self):
        """Atomic number of the isotope."""
        isotope_data = get_isotope_data(self.symbol)
        return isotope_data["atomic_number"]


def format_isotope_string(isotope_string: str) -> str:
    """Format the isotope string to {A}{symbol}, where A is the isotope number."""
    result = match(r"(\d+)\s*(\w+)", isotope_string)

    if result is None:
        raise Exception(f"Could not parse isotope string {isotope_string}")
    isotope = result.group(2)
    A = result.group(1)

    formatted_string = f"{A}{isotope}"
    if formatted_string not in ISOTOPE_DATA:
        raise Exception(f"Could not parse isotope string {isotope_string}")

    return formatted_string


def get_isotope_data(isotope_string: str) -> dict:
    """Get the isotope's intrinsinc properties from a JSON data file."""
    formatted_isotope_string = format_isotope_string(isotope_string)
    isotope_dict = dict(ISOTOPE_DATA[formatted_isotope_string])
    isotope_dict.update({"isotope": formatted_isotope_string})
    return isotope_dict
