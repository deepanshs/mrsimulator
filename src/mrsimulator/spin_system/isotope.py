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
custom_isotope_mapping = {}  # new symbol as key, original symbol as value


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
    >>> carbon.natural_abundance  # in percent
    1.11
    >>> carbon.gyromagnetic_ratio  # in MHz/T
    10.708398861439887
    >>> carbon.atomic_number
    6
    >>> carbon.quadrupole_moment  # in eB
    0.0
    """

    symbol: str
    test_vars: ClassVar[Dict] = {"symbol": "1H"}

    class Config:
        extra = "forbid"
        validate_assignment = True

    @validator("symbol", always=True)
    def validate_symbol(cls, v, *, values, **kwargs):
        if v in custom_isotope_mapping:
            return v
        return format_isotope_string(v)

    @classmethod
    def new(cls, symbol: str, copy_from: str):
        """Create a new isotope with the specified symbol with the same data as the
        isotope specified in copy_from. This way, individual spins in different
        environments can be addressable in queries.

        Arguments:
            (str) symbol: Symbol for the new isotope, cannot match any other isotope
                symbols.
            (str) copy_from: Symbol of real isotope to copy data from

        Returns:
            Instance of new isotope

        Example:
            >>> from mrsimulator.spin_system.isotope import Isotope
            >>> Ha = Isotope.new(symbol="Ha", copy_from="1H")
            >>> Hb = Isotope.new(symbol="Hb", copy_from="1H")
        """
        if symbol in ISOTOPE_DATA or symbol in custom_isotope_mapping:
            raise ValueError(
                "The new symbol for the copied isotope cannot match any previous "
                f"isotope symbols. Got {symbol} which matches a known symbol."
            )

        # Check if there is symbol to map
        if copy_from not in ISOTOPE_DATA:
            raise ValueError(
                "The copy_from symbol must match a known isotope, but "
                f"got {copy_from}, which is unknown."
            )

        custom_isotope_mapping[symbol] = copy_from
        return Isotope(symbol=symbol)

    def json(self, **kwargs) -> str:
        return self.symbol

    @property
    def spin(self):
        """Spin quantum number, I, of the isotope."""
        isotope_data = get_isotope_data(self.symbol)
        return isotope_data["spin"] / 2.0

    @property
    def natural_abundance(self):
        """Natural abundance of the isotope in units of percent."""
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

    def larmor_freq(self, B0=9.4):
        """Return the Larmor frequency of the isotope at a magnetic field strength B0.

        Args:
            float B0: magnetic field strength in T

        Returns:
            float: larmor frequency in MHz

        Example
        -------

        >>> silicon = Isotope(symbol="29Si")
        >>> freq = silicon.larmor_freq(B0 = 9.4)
        """
        return -self.gyromagnetic_ratio * B0


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
    """Get the isotope's intrinsic properties from a JSON data file."""
    if isotope_string in custom_isotope_mapping:
        return get_isotope_data(custom_isotope_mapping[isotope_string])
    else:
        isotope_string = format_isotope_string(isotope_string)
        isotope_dict = dict(ISOTOPE_DATA[isotope_string])
    isotope_dict.update({"isotope": isotope_string})
    return isotope_dict


# @classmethod
# def add_new(
#     cls,
#     symbol: str,
#     spin: float,
#     gyromagnetic_ratio: float,
#     quadrupole_moment: float = 0,
#     natural_abundance: float = 100,
#     atomic_number: int = -1,
# ):
#     """Add isotope data from a custom Isotope into the stored Isotope data and
#     return an instance of the new Isotope. The isotope symbol cannot match an real
#     isotope symbol; if the provided symbol matches a known custom isotope symbol,
#     then an instance of that isotope is returned.

#     Arguments:
#         (str) symbol: Required symbol for custom isotope class. String cannot match
#             another isotope symbol.
#         (float) spin: Required spin number for the isotope. Must be an integer or
#             half-integer greater than zero.
#         (float) gyromagnetic_ratio: Required gyromagnetic ratio of the isotope given
#             in MHz/T.
#         (float) quadrupole_moment: Optional quadrupole moment given in eB. Default
#             is 0.
#         (float) natural_abundance: Optional natural abundance of the isotope given
#             as a percentage between 0 and 100. Default is 100.
#         (int) atomic_number: Optional atomic number for the custom isotope. Can take
#             any integer value and has no bering on simulated spectra.

#     Returns:
#         An Isotope class instance
#     """
#     # Check for symbol overlap in dictionaries
#     if symbol in ISOTOPE_DATA:
#         raise ValueError(
#             "Custom isotope symbol cannot match a real isotope symbol. "
#             f"Got {symbol}."
#         )

#     if symbol in custom_isotope_data:
#         return Isotope(symbol=symbol)

#     # Check for spin integer or half integer
#     if spin <= 0 or not float(2 * spin).is_integer():
#         raise ValueError(
#             f"Isotope spin value must be greater than zero and must be an integer "
#             f"or half integer. Got {spin}."
#         )

#     if not 0 <= natural_abundance <= 100:
#         raise ValueError(
#             "Abundance must be between 0 and 100, inclusive. "
#             f"Got {natural_abundance}."
#         )

#     custom_isotope_data[symbol] = {
#         "spin": int(spin * 2),
#         "natural_abundance": natural_abundance,
#         "gyromagnetic_ratio": gyromagnetic_ratio,
#         "quadrupole_moment": quadrupole_moment,
#         "atomic_number": atomic_number,
#     }
#     return Isotope(symbol=symbol)
