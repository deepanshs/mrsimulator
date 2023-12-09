"""Base Isotope class."""
from os import path
from typing import ClassVar
from typing import Dict

from monty.serialization import loadfn
from pydantic import BaseModel
from pydantic import validator

__author__ = "Deepansh Srivastava"
__email__ = "srivastava.89@osu.edu"

MODULE_DIR = path.dirname(path.abspath(__file__))
ISOTOPE_DATA = loadfn(path.join(MODULE_DIR, "isotope_data.json"))
DEFAULT_ISOTOPE = {
    "spin_multiplicity": 2,
    "gyromagnetic_ratio": 0,
    "quadrupole_moment": 0,
    "natural_abundance": 100,
    "atomic_number": 0,
}


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
    >>> carbon.natural_abundance  # in %
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
    custom_isotope_data: ClassVar[Dict] = {}

    class Config:
        extra = "forbid"
        validate_assignment = True

    @validator("symbol", always=True)
    def validate_symbol(cls, v, *, values, **kwargs):
        if v not in get_all_isotope_symbols():
            raise Exception(f"Isotope symbol `{v}` not recognized.")
        return v

    def json(self, **kwargs) -> dict:
        return (
            self.symbol
            if self.symbol in ISOTOPE_DATA
            else get_isotope_dict(self.symbol)
        )

    # def dict(self, **kwargs) -> dict:
    #     return (
    #         {"symbol": self.symbol} if self.symbol in ISOTOPE_DATA else super().dict()
    #     )

    @classmethod
    def register(cls, symbol: str, copy_from: str = None, **kwargs):
        """Register a new isotope symbol with intrinsic attributes given by kwargs. If
        `copy_from` is not None and matches a known isotope symbol, then the attributes
        are copied from that symbol.

        Arguments:
            (str) symbol: The isotope symbol to register. Must be unique.
            (str) copy_from: An optional isotope symbol to copy attributes from. If
                None, then no attributes are copied.
        """
        if symbol in ISOTOPE_DATA.keys():
            raise KeyError(
                f"`{symbol}` is an immutable registry symbol. Use a different symbol."
            )

        if copy_from is not None and copy_from not in get_all_isotope_symbols():
            raise KeyError(
                f"Cannot copy from unrecognized isotope symbol `{copy_from}`"
            )

        # Grab copied isotope attributes if copy from, otherwise default attributes
        kwargs_new = (
            get_isotope_dict(copy_from)
            if copy_from is not None
            else DEFAULT_ISOTOPE.copy()
        )
        kwargs_new.update(kwargs)

        # Validate the new isotope attributes, then add to lookup dictionary
        _validate_isotope_kwargs(kwargs_new)
        cls.custom_isotope_data[symbol] = kwargs_new

    @classmethod
    def parse(cls, item):
        """Attempt to parse the provided item into an Isotope instance"""
        if isinstance(item, str):
            return Isotope(symbol=item)

        if isinstance(item, dict):
            symbol = item["symbol"]
            if symbol not in ISOTOPE_DATA.keys():
                Isotope.register(**item)
            return Isotope(symbol=symbol)

        if isinstance(item, Isotope):
            return item

        raise ValueError(f"Cannot parse type {type(item)} into an Isotope")

    @property
    def spin(self):
        """Spin quantum number, I, of the isotope."""
        isotope_data = get_isotope_dict(self.symbol)
        return (isotope_data["spin_multiplicity"] - 1) / 2.0

    @property
    def spin_multiplicity(self):
        """Spin multiplicity, (2I+1), of the isotope."""
        isotope_data = get_isotope_dict(self.symbol)
        return isotope_data["spin_multiplicity"]

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
            float: larmor frequency in MHz

        Example
        -------

        >>> silicon = Isotope(symbol="29Si")
        >>> freq = silicon.larmor_freq(B0 = 9.4)
        """
        return -self.gyromagnetic_ratio * B0


def get_isotope_dict(isotope_symbol: str) -> dict:
    """Get the intrinsic properties of isotope with symbol isotope_symbol as a dict"""
    if isotope_symbol not in get_all_isotope_symbols():
        raise Exception(f"Isotope symbol `{isotope_symbol}` not recognized.")

    isotope_dict = get_all_isotope_data()[isotope_symbol]
    isotope_dict.update({"symbol": isotope_symbol})
    return isotope_dict


def get_all_isotope_data() -> dict:
    """Returns an intersected dictionary of the real and custom isotope data"""
    iso_dict = ISOTOPE_DATA.copy()
    iso_dict.update(Isotope.custom_isotope_data)
    return iso_dict


def get_all_isotope_symbols() -> list:
    """Returns a list of all currently valid isotope symbols"""
    return list(ISOTOPE_DATA.keys()) + list(Isotope.custom_isotope_data.keys())


def _validate_isotope_kwargs(iso_dict: dict) -> None:
    """Ensure the attributes in the passed dictionary are within allowed bounds"""
    # Check for spin integer or half integer
    sm = iso_dict["spin_multiplicity"]
    if not isinstance(sm, int) or sm <= 1:
        raise ValueError(
            f"Isotope spin_multiplicity value must an integer value greater than "
            f"one. Got {sm}."
        )

    # Check abundance between 0 and 100, inclusive
    na = iso_dict["natural_abundance"]
    if not 0 <= na <= 100:
        raise ValueError(
            "Abundance must be between 0 and 100, inclusive. " f"Got {na}."
        )

    # Ensure atomic number is an integer value
    an = iso_dict["atomic_number"]
    if not isinstance(an, int):
        raise ValueError(f"Atomic number must be an integer value. Got {an}")
