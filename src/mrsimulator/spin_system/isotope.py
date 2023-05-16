"""Base Isotope class."""
from os import path
from re import match
from typing import ClassVar
from typing import Dict

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
    10.708398861439887
    >>> carbon.atomic_number
    6
    >>> carbon.quadrupole_moment # in eB
    0.0
    """

    symbol: str
    # test_vars: ClassVar[Dict] = {"symbol": "1H"}
    custom_isotope_data: ClassVar[Dict] = {}

    class Config:
        extra = "forbid"
        validate_assignment = True

    @validator("symbol", always=True)
    def validate_symbol(cls, v, *, values, **kwargs):
        return format_isotope_string(v)

    def json(self, **kwargs) -> dict:
        return get_isotope_dict(self.symbol)

    def dict(self, **kwargs) -> dict:
        return self.json()

    @classmethod
    def get_isotope(cls, val):
        """Ensuring backwards compatibility with previous serializations and workflows
        means that Isotope objects may need instantiated from string values (isotope
        symbols), dictionary objects (all defining attributes of the Isotope objects),
        or Isotope objects themselves. This utility function parses the type of val and
        returns the corresponding isotope object

        Arguments:
            val: A string, dictionary, or Isotope object representing the isotope
            symbol, isotope data in dictionary form, or
            :py:class:~`mrsimulator.spin_system.isotope.Isotope` object, respectively.

        Returns:
            An instance of the Isotope class
        """
        # Return Isotope object, no further checking needed
        if isinstance(val, Isotope):
            return val

        # Check if string symbol recognized, then return Isotope object
        if isinstance(val, str):
            if val not in get_all_isotope_symbols():
                raise ValueError(f"{val} is an unrecognized Isotope symbol.")

            return Isotope(symbol=val)

        # Value is dictionary, meaning either need to add custom isotope or get symbol
        if isinstance(val, dict):
            if val["isotope"] in get_all_isotope_symbols():  # Already known symbol
                return Isotope(symbol=val["isotope"])

            return Isotope.add_new(**val)

        raise ValueError(f"Type {type(val)} is invalid for this method.")

    @classmethod
    def add_new(
        cls,
        symbol: str,
        spin: float,
        gyromagnetic_ratio: float,
        quadrupole_moment: float = 0,
        natural_abundance: float = 100,
        atomic_number: int = -1,
    ):
        """Add isotope data from a custom Isotope into the stored Isotope data and
        return an instance of the new Isotope. The isotope symbol cannot match an real
        isotope symbol; if the provided symbol matches a known custom isotope symbol,
        then an instance of that isotope is returned.

        Arguments:
            (str) symbol: Required symbol for custom isotope class. String cannot match
                another isotope symbol.
            (float) spin: Required spin number for the isotope. Must be an integer or
                half-integer greater than zero.
            (float) gyromagnetic_ratio: Required gyromagnetic ratio of the isotope given
                in MHz/T.
            (float) quadrupole_moment: Optional quadrupole moment given in eB. Default
                is 0.
            (float) natural_abundance: Optional natural abundance of the isotope given
                as a percentage between 0 and 100. Default is 100.
            (int) atomic_number: Optional atomic number for the custom isotope. Can take
                any integer value and has no bering on simulated spectra.

        Returns:
            An instance of the Isotope class
        """
        # Check for symbol overlap in dictionaries
        if symbol in ISOTOPE_DATA or symbol in Isotope.custom_isotope_data:
            raise ValueError(
                f"Symbol, {symbol}, is already attributed to another Isotope. All "
                "Isotope symbols must be unique; please choose a different symbol."
            )

        # Check for spin integer or half integer
        if spin <= 0 or not float(2 * spin).is_integer():
            raise ValueError(
                f"Isotope spin value must be greater than zero and must be an integer "
                f"or half integer. Got {spin}."
            )

        # Check abundance between 0 and 100, inclusive
        if not 0 <= natural_abundance <= 100:
            raise ValueError(
                "Abundance must be between 0 and 100, inclusive. "
                f"Got {natural_abundance}."
            )

        # Ensure atomic number is an integer value
        if not float(atomic_number).is_integer():
            raise ValueError(
                f"Atomic number must be an integer value. Got {atomic_number}"
            )

        Isotope.custom_isotope_data[symbol] = {
            "spin": int(spin * 2),
            "natural_abundance": natural_abundance,
            "gyromagnetic_ratio": gyromagnetic_ratio,
            "quadrupole_moment": quadrupole_moment,
            "atomic_number": atomic_number,
        }

        return Isotope(symbol=symbol)

    @property
    def spin(self):
        """Spin quantum number, I, of the isotope."""
        isotope_data = get_isotope_dict(self.symbol)
        return isotope_data["spin"] / 2.0

    @property
    def natural_abundance(self):
        """Natural abundance of the isotope in units of %."""
        isotope_data = get_isotope_dict(self.symbol)
        return isotope_data["natural_abundance"]

    @property
    def gyromagnetic_ratio(self):
        """Reduced gyromagnetic ratio of the nucleus given in units of MHz/T."""
        isotope_data = get_isotope_dict(self.symbol)
        return isotope_data["gyromagnetic_ratio"]

    @property
    def quadrupole_moment(self):
        """Quadrupole moment of the nucleus given in units of eB (electron-barn)."""
        isotope_data = get_isotope_dict(self.symbol)
        return isotope_data["quadrupole_moment"]

    @property
    def atomic_number(self):
        """Atomic number of the isotope."""
        isotope_data = get_isotope_dict(self.symbol)
        return isotope_data["atomic_number"]

    def larmor_freq(self, B0=9.4):
        """Return the Larmor frequency of the isotope at a magnetic field strength B0.

        Args:
            float B0: magnetic field strength in T

        Returns:
            float: Larmor frequency in MHz

        Example
        -------

        >>> silicon = Isotope(symbol="29Si")
        >>> freq = silicon.larmor_freq(B0 = 9.4)
        """
        return -self.gyromagnetic_ratio * B0


def format_isotope_string(isotope_string: str) -> str:
    """Format the isotope string to {A}{symbol}, where A is the isotope number."""
    # Skip formatting if the symbol is from a custom isotope
    if isotope_string in Isotope.custom_isotope_data:
        return isotope_string

    result = match(r"(\d+)\s*(\w+)", isotope_string)

    if result is None:
        raise Exception(f"Could not parse isotope string {isotope_string}")

    isotope = result.group(2)
    A = result.group(1)

    formatted_string = f"{A}{isotope}"
    if formatted_string not in ISOTOPE_DATA:
        raise Exception(f"Could not parse isotope string {isotope_string}")

    return formatted_string


def get_isotope_dict(isotope_string: str) -> dict:
    """Get the isotope's intrinsic properties from a JSON data file."""
    formatted_isotope_string = format_isotope_string(isotope_string)

    # isotope_string is always unique, and only exists in one of the dictionaries
    if formatted_isotope_string in ISOTOPE_DATA:
        isotope_dict = dict(ISOTOPE_DATA[formatted_isotope_string])
    else:
        isotope_dict = dict(Isotope.custom_isotope_data[formatted_isotope_string])

    isotope_dict.update({"isotope": formatted_isotope_string})

    return isotope_dict


def get_all_isotope_symbols() -> list:
    """Returns a list of all currently valid isotope symbols"""
    return list(ISOTOPE_DATA.keys()) + list(Isotope.custom_isotope_data.keys())
