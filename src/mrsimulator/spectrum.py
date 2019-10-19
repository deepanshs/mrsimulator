# -*- coding: utf-8 -*-
import os.path
import re
import numpy as np
import warnings
from monty.serialization import loadfn
from mrsimulator import Parseable
from typing import ClassVar, Optional

__author__ = "Deepansh J. Srivastava"
__email__ = ["srivastava.89@osu.edu", "deepansh2012@gmail.com"]

MODULE_DIR = os.path.dirname(os.path.abspath(__file__))
ISOTOPE_DATA = loadfn(os.path.join(MODULE_DIR, "isotope_data.json"))


class Dimension(Parseable):
    r"""
    Base class for an NMR spectroscopic dimension.

    Args:
        number_of_points: A `required` integer with the number of points, :math:`N`,
                along the dimension.
        spectral_width: A `required` float with the spectral width,
                :math:`\Delta x`, along the dimension in units of Hz.
        reference_offset: An `optional` float with the reference offset, :math:`x_0`
                along the dimension in units of Hz. The default value is 0.
        magnetic_flux_density: An `optional` float with the macroscopic magnetic flux
                density, :math:`B_0`, of the applied external magnetic field in
                units of T. The default value is 9.4.
        rotor_frequency: An `optional` float with the sample spinning frequency
                :math:`\nu_r`, in units of Hz. The default value is 0.
        rotor_angle: An `optional` float with the angle between the
                sample rotation axis and the applied external magnetic field,
                :math:`\theta`, in units of rad. The default value is 0.9553166,
                i.e. the magic angle.
        isotope: An `optional` isotope string in "{A}{Symbol}" notation such as
                1H or 29Si. The default is None.
        label: An `optional`

    .. rubric:: Attribute Documentation

    Attributes:
        number_of_points: Number of points, :math:`N`, along the dimension.
        spectral_width: Spectral width, :math:`\Delta x`, along the
                dimension in units of Hz.
        reference_offset: Reference offset, :math:`x_0`, along the
                dimension in units of Hz.
        magnetic_flux_density: The macroscopic magnetic flux density of
                the applied external magnetic field, :math:`B_0`, in units of T.
        rotor_frequency: Sample spinning frequency, :math:`\nu_r`,
                in units of Hz.
        rotor_angle: The angle between the sample rotation axis and the
                applied macroscopic magnetic field, :math:`\theta` in units of rad.
        isotope: The isotope assigned to the dimension.
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
    _spin: float = 1
    _natural_abundance: float = 1.0
    _gyromagnetic_ratio: float = 1.0

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
    def spin(self):
        if self.isotope is None:
            return self._spin
        isotope_data = get_isotope_data(self.isotope)
        return isotope_data["spin"]

    @spin.setter
    def spin(self, value):
        if self.isotope is None:
            self._spin = value

    @property
    def natural_abundance(self):
        if self.isotope is None:
            return self._natural_abundance
        isotope_data = get_isotope_data(self.isotope)
        return isotope_data["natural_abundance"]

    @natural_abundance.setter
    def natural_abundance(self, value):
        if self.isotope is None:
            self._natural_abundance = value

    @property
    def gyromagnetic_ratio(self):
        if self.isotope is None:
            return 0.0
        isotope_data = get_isotope_data(self.isotope)
        return isotope_data["gyromagnetic_ratio"]

    @gyromagnetic_ratio.setter
    def gyromagnetic_ratio(self, value):
        if self.isotope is None:
            self._gyromagnetic_ratio = value

    @spin.setter
    def spin(self, value):
        if self.isotope is None:
            self._spin = value

    @property
    def larmor_frequency(self):
        r"""
        Signed Larmor frequency, :math:`\omega_0=-\gamma B_0`, of the isotope,
        in MHz, where :math:`\gamma` is the gyromagnetic ratio of the
        isotope and :math:`B_0` is the macroscopic magnetic flux density
        of the applied external field.
        """
        return -self.gyromagnetic_ratio * self.magnetic_flux_density

    @property
    def coordinates_Hz(self):
        r"""
        The grid coordinates along the dimension in units of Hz, evaluated as

        .. math::
            x_\text{Hz} = \left([0, 1, ... N-1] - T\right) \frac{\Delta x}{N} + x_0

        where :math:`T=N/2` and :math:`T=(N-1)/2` for even and odd values of
        :math:`N`, respectively.
        """
        n = self.number_of_points
        Tk = int(n / 2)
        return (np.arange(n) - Tk) * (self.spectral_width / n) + self.reference_offset

    @property
    def coordinates_ppm(self):
        r"""
        The grid coordinates along the dimension as dimension frequency ratio
        in ppm. The coordinates are evaluated as

        .. math::
            x_\text{ppm} = \frac{x_\text{Hz}} {x_0 + \omega_0}
        """
        if self.isotope is None:
            raise warnings.warn(
                (
                    "The coordinates along the dimension without an assigned "
                    "isotope cannot be converted to dimensionless ratio."
                )
            )
        else:

            denominator = np.abs(self.reference_offset * 1e-6 + self.larmor_frequency)
            return self.coordinates_Hz / denominator

    @classmethod
    def parse_dict_with_units(cls, py_dict):
        """
        Parse the physical quantities of a Dimension object
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
