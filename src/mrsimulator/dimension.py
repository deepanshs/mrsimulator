# -*- coding: utf-8 -*-
import os.path
import re
import warnings
from typing import ClassVar
from typing import Dict
from typing import Optional

import numpy as np
from monty.serialization import loadfn
from mrsimulator import Parseable
from pydantic import Field
from pydantic import validator

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
        label: An `optional` string label.

    .. rubric:: Attribute Documentation

    Attributes:
    """

    # public attributes
    number_of_points: int = Field(1024, gt=0)
    spectral_width: float = Field(..., gt=0)
    reference_offset: Optional[float] = Field(default=0)
    magnetic_flux_density: Optional[float] = Field(default=9.4, ge=0)
    rotor_frequency: Optional[float] = Field(default=0, ge=0)
    # 54.735 degrees = 0.9553166 radians
    rotor_angle: Optional[float] = Field(default=0.9553166, ge=0, le=1.5707963268)
    isotope: Optional[str] = None
    label: Optional[str] = ""

    # Immutable attributes
    # spin: Optional[float] = None
    # natural_abundance: Optional[float] = None
    # gyromagnetic_ratio: Optional[float] = None
    # quadrupole_moment: Optional[float] = None
    # atomic_number: Optional[int] = None
    # larmor_frequency: Optional[float] = None  # default is Hz

    property_unit_types: ClassVar = {
        "spectral_width": ["frequency", "dimensionless"],
        "reference_offset": ["frequency", "dimensionless"],
        "magnetic_flux_density": "magnetic flux density",
        "rotor_frequency": "frequency",
        "rotor_angle": "angle",
    }

    property_default_units: ClassVar = {
        "spectral_width": ["Hz", "ppm"],
        "reference_offset": ["Hz", "ppm"],
        "magnetic_flux_density": "T",
        "rotor_frequency": "Hz",
        "rotor_angle": "rad",
    }

    property_units: Dict = {
        "spectral_width": "Hz",
        "reference_offset": "Hz",
        "magnetic_flux_density": "T",
        "rotor_frequency": "Hz",
        "rotor_angle": "rad",
        "natural_abundance": "%",
        "gyromagnetic_ratio": "MHz/T",
        "quadrupole_moment": "eB",
        "larmor_frequency": "Hz",
    }

    class Config:
        validate_assignment = True

    @validator("isotope", always=True)
    def get_isotope(cls, v, *, values, **kwargs):
        if v is None:
            return v
        return format_isotope_string(v)

    # @validator("spin", always=True)
    # def immutable_spin(cls, v, *, values, **kwargs):
    #     if values["isotope"] is None:
    #         return None
    #     return get_isotope_data(values["isotope"])["spin"]

    # @validator("natural_abundance", always=True)
    # def immutable_natural_abundance(cls, v, *, values, **kwargs):
    #     if values["isotope"] is None:
    #         return None
    #     return get_isotope_data(values["isotope"])["natural_abundance"]

    # @validator("gyromagnetic_ratio", always=True)
    # def immutable_gyromagnetic_ratio(cls, v, *, values, **kwargs):
    #     if values["isotope"] is None:
    #         return None
    #     return get_isotope_data(values["isotope"])["gyromagnetic_ratio"]

    # @validator("quadrupole_moment", always=True)
    # def immutable_quadrupole_moment(cls, v, *, values, **kwargs):
    #     if values["isotope"] is None:
    #         return None
    #     return get_isotope_data(values["isotope"])["quadrupole_moment"]

    # @validator("atomic_number", always=True)
    # def immutable_atomic_number(cls, v, *, values, **kwargs):
    #     if values["isotope"] is None:
    #         return None
    #     return get_isotope_data(values["isotope"])["atomic_number"]

    # @validator("larmor_frequency", always=True)
    # def immutable_larmor_frequency(cls, v, *, values, **kwargs):
    #     if values["isotope"] is None:
    #         return None
    #     B0 = values["magnetic_flux_density"]
    #     gamma = get_isotope_data(values["isotope"])["gyromagnetic_ratio"]
    #     return -gamma * B0 * 1e6  # larmor freqency in Hz

    @property
    def spin(self):
        """Spin quantum number, I."""
        if self.isotope is None:
            return None
        isotope_data = get_isotope_data(self.isotope)
        return isotope_data["spin"] / 2.0

    @property
    def natural_abundance(self):
        """Natural abundance of the spin."""
        if self.isotope is None:
            return None
        isotope_data = get_isotope_data(self.isotope)
        return isotope_data["natural_abundance"]

    # @natural_abundance.setter
    # def natural_abundance(self, value):
    #     if self.isotope is None:
    #         self._natural_abundance = value

    @property
    def gyromagnetic_ratio(self):
        """Reduced gyromagnetic ratio of the nucleus given in units of MHz/T."""
        if self.isotope is None:
            return None
        isotope_data = get_isotope_data(self.isotope)
        return isotope_data["gyromagnetic_ratio"]

    # @gyromagnetic_ratio.setter
    # def gyromagnetic_ratio(self, value):
    #     if self.isotope is None:
    #         self._gyromagnetic_ratio = value

    @property
    def quadrupole_moment(self):
        """Quadrupole moment of the nucleus given in units of electron-barn."""
        if self.isotope is None:
            return None
        isotope_data = get_isotope_data(self.isotope)
        return isotope_data["quadrupole_moment"]

    # @quadrupole_moment.setter
    # def quadrupole_moment(self, value):
    #     if self.isotope is None:
    #         self._quadrupole_moment = value

    @property
    def atomic_number(self):
        """Atomic number of the nucleus"""
        if self.isotope is None:
            return None
        isotope_data = get_isotope_data(self.isotope)
        return isotope_data["atomic_number"]

    # @atomic_number.setter
    # def atomic_number(self, value):
    #     if self.isotope is None:
    #         self._atomic_number = value

    @property
    def larmor_frequency(self):
        r"""
        Signed Larmor frequency, :math:`\omega_0=-\gamma B_0`, of the isotope,
        in MHz, where :math:`\gamma` is the gyromagnetic ratio of the
        isotope and :math:`B_0` is the macroscopic magnetic flux density
        of the applied external field.
        """
        if self.isotope is None:
            return None
        return -self.gyromagnetic_ratio * self.magnetic_flux_density * 1e6  # in Hz

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
            warnings.warn(
                (
                    "The coordinates along the dimension without an assigned "
                    "isotope cannot be converted to dimensionless frequency ratio."
                )
            )
        else:

            denominator = (self.reference_offset + self.larmor_frequency) / 1e6
            return self.coordinates_Hz / abs(denominator)

    @classmethod
    def parse_dict_with_units(cls, py_dict):
        """
        Parse the physical quantities of a Dimension object
        when expressed as a python dictionary.

        Args:
            py_dict: Python dictionary representation of an isotopomers with
                        physical quantities.
        """
        return super().parse_dict_with_units(py_dict)

    def to_dict_with_units(self):
        """
        Serialize the Spectrum object to a JSON compliant python dictionary with
        units.
        """
        temp_dict = {
            k: v
            for k, v in self.dict(exclude={"property_units"}).items()
            if v is not None
        }
        temp_keys = temp_dict.keys()
        for key, unit in self.property_units.items():
            if key in temp_keys:
                temp_dict[key] = f"{temp_dict[key]} {unit}"

        return temp_dict


def format_isotope_string(isotope_string):
    """Format isotope string to {A}{symbol}, where A is the isotope number"""
    result = re.match(r"(\d+)\s*(\w+)", isotope_string)

    if result is None:
        raise Exception(f"Could not parse isotope string {isotope_string}")
    isotope = result.group(2)
    A = result.group(1)

    formatted_string = f"{A}{isotope}"
    if formatted_string not in ISOTOPE_DATA:
        raise Exception(f"Could not parse isotope string {isotope_string}")

    return formatted_string


def get_isotope_data(isotope_string):
    """
    Gets the isotope's intrinsinc NMR properties from a JSON
    data file.
    """
    formatted_isotope_string = format_isotope_string(isotope_string)
    isotope_dict = dict(ISOTOPE_DATA[formatted_isotope_string])
    isotope_dict.update({"isotope": formatted_isotope_string})
    return isotope_dict
