# -*- coding: utf-8 -*-
import os.path
import re
from monty.serialization import loadfn
from mrsimulator import Parseable
from typing import ClassVar

__author__ = "Deepansh J. Srivastava"
__email__ = ["srivastava.89@osu.edu", "deepansh2012@gmail.com"]

MODULE_DIR = os.path.dirname(os.path.abspath(__file__))
ISOTOPE_DATA = loadfn(os.path.join(MODULE_DIR, "isotope_data.json"))


class SpectroscopicDimension(Parseable):
    """
    Base class for an NMR Spectrum

    Args:
        number_of_points: Number of points.
        spectral_width: Width of spectrum region to consider in Hz.
        reference_offset: Offset frequency in Hz.
        magnetic_flux_density: Magnetic flux density in T.
        rotor_frequency: Frequency in Hz.
        rotor_angle: Angle in radians.
        isotope: Isotope in "{A}{Symbol}" notation such as 1H or 29Si.
        spin: nuclear spin quantum number as 2I.
        natural_abundance: fractional natural abundance, IE sums should equal 1.
        gyromagnetic_ratio: In units of (MHz/T).
    """

    number_of_points: int = 1024
    spectral_width: float = 100
    reference_offset: float = 0
    magnetic_flux_density: float = 9.4
    rotor_frequency: float = 0
    rotor_angle: float = 0.9553  # 54.735 degrees in radians
    isotope: str = "1H"
    spin: int = 1
    natural_abundance: float = 0.04683
    gyromagnetic_ratio: float = -8.465

    property_unit_types: ClassVar = {
        "spectral_width": "frequency",
        "reference_offset": "frequency",
        "magnetic_flux_density": "magnetic flux density",
        "rotor_frequency": "frequency",
        "rotor_angle": "angle",
    }

    property_default_units: ClassVar = {
        "spectral_width": "Hz",
        "reference_offset": "Hz",
        "magnetic_flux_density": "T",
        "rotor_frequency": "Hz",
        "rotor_angle": "rad",
    }

    @property
    def larmor_frequency(self):
        """
        Larmor frequency in MHz
        """
        return self.gyromagnetic_ratio * self.magnetic_flux_density

    @classmethod
    def parse_dict_with_units(cls, json_dict):

        # if "direct_dimension" in json_dict:
        #     json_dict.update(json_dict["direct_dimension"])

        if "isotope" in json_dict:
            isotope_data = get_isotope_data(json_dict["isotope"])
            json_dict.update(isotope_data)

        return super().parse_dict_with_units(json_dict)


def get_isotope_data(isotope_string):
    """
    Gets the isotope's intrinsinc NMR properties from a JSON
    data file
    """
    result = re.match(r"(\d+)\s*(\w+)", isotope_string)
    isotope = result.group(2)
    A = result.group(1)

    formatted_isotope_string = f"{A}{isotope}"

    if formatted_isotope_string in ISOTOPE_DATA:
        isotope_dict = dict(ISOTOPE_DATA[formatted_isotope_string])
        isotope_dict.update({"isotope": formatted_isotope_string})
        return isotope_dict
    else:
        raise Exception(f"Could not parse isotope string {formatted_isotope_string}")
