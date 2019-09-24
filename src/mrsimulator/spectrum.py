# -*- coding: utf-8 -*-
import os.path
import re
from monty.serialization import loadfn
from mrsimulator import Parseable
from typing import ClassVar, Optional

__author__ = "Deepansh J. Srivastava"
__email__ = ["srivastava.89@osu.edu", "deepansh2012@gmail.com"]

MODULE_DIR = os.path.dirname(os.path.abspath(__file__))
ISOTOPE_DATA = loadfn(os.path.join(MODULE_DIR, "isotope_data.json"))


class SpectroscopicDimension(Parseable):
    """
    Base class for an NMR spectroscopic dimension.

    .. rubric:: Attributes Documentation

    Attributes:
        number_of_points: A `required` integer with the number of points along
                the spectroscopic dimension.
        spectral_width: A `required` float representing the spectral width of the
                spectroscopic dimension in units of Hz.
        reference_offset: An `optional` float representing the reference offset
                along the spectroscopic dimension in units of Hz. The default
                value is 0.
        magnetic_flux_density: An `optional` float represeting the magnetic flux
                density of the external magnetic field in units of T. The
                default value is 9.4.
        rotor_frequency: An `optional` float representing the sample spinning
                frequency in units of Hz. The default value is 0.
        rotor_angle: An `optional` float representing the angle between the
                sample rotation axis and the external magnetic field in units
                of rad. The default value is 0.9553166.
        isotope: An `optional` isotope string in "{A}{Symbol}" notation such as
                1H or 29Si. The default is None.
    """

    # public attributes
    number_of_points: int
    spectral_width: float
    reference_offset: Optional[float] = 0
    magnetic_flux_density: Optional[float] = 9.4
    rotor_frequency: Optional[float] = 0
    rotor_angle: Optional[float] = 0.9553166  # 54.735 degrees in radians
    isotope: Optional[str] = None

    # private attributes
    spin: int = 1
    natural_abundance: float = 1.0
    gyromagnetic_ratio: float = 0.0

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
        Larmor frequency of the site in MHz.
        """
        return self.gyromagnetic_ratio * self.magnetic_flux_density

    @classmethod
    def parse_dict_with_units(cls, py_dict):
        """
        Parse the physical quantities of a SpectroscopicDimension object
        when expressed as a python dictionary.

        Args:
            py_dict: Python dictionary representation of an isotopomers with
                        physical quantities.
        """
        if "isotope" in py_dict:
            isotope_data = get_isotope_data(py_dict["isotope"])
            py_dict.update(isotope_data)

        return super().parse_dict_with_units(py_dict)


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
