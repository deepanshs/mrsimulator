"""Base Isotope class."""
from os import path
from re import match
from typing import ClassVar
from typing import Dict

from monty.serialization import loadfn
from pydantic import BaseModel
from pydantic import Field
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
    spin_multiplicity: int = Field(default=2, ge=2)
    gyromagnetic_ratio: float = 0
    quadrupole_moment: float = 0
    natural_abundance: float = Field(default=100, ge=0, le=100)
    atomic_number: int = Field(default=0)

    test_vars: ClassVar[Dict] = {"symbol": "1H"}
    custom_isotope_data: ClassVar[Dict] = {}

    class Config:
        extra = "forbid"
        validate_assignment = True
        allow_mutation = False

    def __init__(self, **kwargs):
        symbol = kwargs["symbol"]
        if symbol not in get_all_isotope_symbols():
            raise Exception(f"Isotope symbol `{symbol}` not recognized.")
        kwargs_new = get_isotope_dict(symbol)
        for k, v in kwargs.items():
            if v != kwargs_new[k]:
                raise ValueError(f"{k} for {symbol} cannot be assigned.")
        super().__init__(**kwargs_new)

    @classmethod
    def register(cls, symbol, copy_from=None, **kwargs):
        if symbol in ISOTOPE_DATA.keys():
            raise KeyError(
                f"`{symbol}` is an immutable registry symbol. Use a different symbol."
            )

        kwargs_new = (
            get_isotope_dict(copy_from)
            if copy_from is not None
            else DEFAULT_ISOTOPE.copy()
        )
        kwargs_new.update(kwargs)
        kwargs_new["symbol"] = symbol
        cls.custom_isotope_data[symbol] = kwargs_new

    @classmethod
    def parse(cls, item):
        if isinstance(item, str):
            return Isotope(symbol=item)
        if isinstance(item, dict):
            return Isotope(**item)
        return item

    @validator("symbol", always=True)
    def validate_symbol(cls, v, *, values, **kwargs):
        return format_isotope_string(v)

    def json(self, **kwargs) -> dict:
        return self.symbol if self.symbol in ISOTOPE_DATA else self.dict()

    def dict(self, **kwargs) -> dict:
        return (
            {"symbol": self.symbol} if self.symbol in ISOTOPE_DATA else super().dict()
        )

    @property
    def spin(self):
        """Spin quantum number, I, of the isotope."""
        isotope_data = get_isotope_dict(self.symbol)
        return (isotope_data["spin_multiplicity"] - 1) / 2.0

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
        raise Exception(f"Could not parse isotope string `{isotope_string}`")

    isotope = result.group(2)
    A = result.group(1)

    formatted_string = f"{A}{isotope}"
    if formatted_string not in ISOTOPE_DATA:
        raise Exception(f"Could not parse isotope string `{isotope_string}`")

    return formatted_string


def get_isotope_dict(isotope_string: str) -> dict:
    """Get the isotope's intrinsic properties from a JSON data file."""
    formatted_isotope_string = format_isotope_string(isotope_string)

    # isotope_string is always unique, and only exists in one of the dictionaries
    if formatted_isotope_string in ISOTOPE_DATA:
        isotope_dict = dict(ISOTOPE_DATA[formatted_isotope_string])
    else:
        isotope_dict = dict(Isotope.custom_isotope_data[formatted_isotope_string])

    isotope_dict.update({"symbol": formatted_isotope_string})

    return isotope_dict


def get_all_isotope_symbols() -> list:
    """Returns a list of all currently valid isotope symbols"""
    return list(ISOTOPE_DATA.keys()) + list(Isotope.custom_isotope_data.keys())


def get_all_isotope_data() -> dict:
    """Returns an intersected dictionary of the real and custom isotope data"""
    iso_dict = ISOTOPE_DATA.copy()
    iso_dict.update(Isotope.custom_isotope_data)
    return iso_dict
