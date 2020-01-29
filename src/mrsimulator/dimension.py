# -*- coding: utf-8 -*-
"""Base Dimension class."""
import os.path
import re
import warnings
from typing import ClassVar
from typing import Dict
from typing import Optional

import numpy as np
from csdmpy import Dimension as csdm_dimension
from monty.serialization import loadfn
from mrsimulator import Parseable
from pydantic import Field
from pydantic import validator

__author__ = "Deepansh J. Srivastava"
__email__ = "deepansh2012@gmail.com"

MODULE_DIR = os.path.dirname(os.path.abspath(__file__))
ISOTOPE_DATA = loadfn(os.path.join(MODULE_DIR, "isotope_data.json"))


class Dimension(Parseable):
    r"""
    Base class for an NMR spectroscopic dimension.

    Args:
        number_of_points: An optional integer with the number of points, :math:`N`,
                along the dimension. The default value is ``1024``.
        spectral_width: A `required` float with the spectral width,
                :math:`\Delta x`, along the dimension in units of Hz.
        reference_offset: An `optional` float with the reference offset, :math:`x_0`
                along the dimension in units of Hz. The default value is ``0``.
        magnetic_flux_density: An `optional` float with the macroscopic magnetic flux
                density, :math:`B_0`, of the applied external magnetic field in
                units of T. The default value is ``9.4``.
        rotor_frequency: An `optional` float with the sample spinning frequency
                :math:`\nu_r`, in units of Hz. The default value is ``0``.
        rotor_angle: An `optional` float with the angle between the
                sample rotation axis and the applied external magnetic field,
                :math:`\theta`, in units of rad. The default value is ``0.9553166``,
                i.e. the magic angle.
        isotope: An `optional` isotope string in "{A}{Symbol}" notation such as
                1H or 29Si. The default is`` None``.
        label: An `optional` string label. The default is an empty string.

    Example:
        >>> dim = Dimension(isotope='27Al', spectral_width=50000, rotor_frequency=12000)

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
    }

    class Config:
        validate_assignment = True

    @validator("isotope", always=True)
    def get_isotope(cls, v, *, values, **kwargs):
        if v is None:
            return v
        return format_isotope_string(v)

    @property
    def spin(self):
        """
        Spin quantum number, I.

        Example:
            >>> dim.spin
            2.5
        """
        if self.isotope is None:
            return None
        isotope_data = get_isotope_data(self.isotope)
        return isotope_data["spin"] / 2.0

    @property
    def natural_abundance(self):
        """
        Natural abundance of the spin.

        Example:
            >>> dim.natural_abundance
            100.0
        """
        if self.isotope is None:
            return None
        isotope_data = get_isotope_data(self.isotope)
        return isotope_data["natural_abundance"]

    @property
    def gyromagnetic_ratio(self):
        """
        Reduced gyromagnetic ratio of the nucleus given in units of MHz/T.

        Example:
            >>> dim.gyromagnetic_ratio
            11.10309
        """
        if self.isotope is None:
            return None
        isotope_data = get_isotope_data(self.isotope)
        return isotope_data["gyromagnetic_ratio"]

    @property
    def quadrupole_moment(self):
        """
        Quadrupole moment of the nucleus given in units of electron-barn.

        Example:
            >>> dim.quadrupole_moment
            0.15
        """
        if self.isotope is None:
            return None
        isotope_data = get_isotope_data(self.isotope)
        return isotope_data["quadrupole_moment"]

    @property
    def atomic_number(self):
        """
        Atomic number of the nucleus.

        Example:
            >>> dim.atomic_number
            13
        """
        if self.isotope is None:
            return None
        isotope_data = get_isotope_data(self.isotope)
        return isotope_data["atomic_number"]

    @property
    def larmor_frequency(self):
        r"""
        Signed Larmor frequency, :math:`\omega_0=-\gamma B_0`, of the isotope,
        in units of Hz, where :math:`\gamma` is the gyromagnetic ratio of the
        isotope and :math:`B_0` is the macroscopic magnetic flux density
        of the applied external field.

        Example:
            >>> dim.larmor_frequency
            -104369046.0
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

        Example:
            >>> dim.coordinates_Hz[:5]
            array([-25000.      , -24951.171875, -24902.34375 , -24853.515625,
                   -24804.6875  ])
        """
        n = self.number_of_points
        Tk = int(n / 2)
        return (np.arange(n) - Tk) * (self.spectral_width / n) + self.reference_offset

    @property
    def coordinates_ppm(self):
        r"""
        The grid coordinates along the dimension as dimension frequency ratio
        in units of ppm. The coordinates are evaluated as

        .. math::
            x_\text{ppm} = \frac{x_\text{Hz}} {x_0 + \omega_0}

        where :math:`\omega_0` is the Larmor frequency.

        Example:
            >>> dim.coordinates_ppm[:5]
            array([-239.53462217, -239.06678111, -238.59894005, -238.13109899,
                   -237.66325794])
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
            py_dict: Dict object

        Example:
            >>> dimension_1 = {
            ...     "number_of_points": 1024,
            ...     "spectral_width": "100 Hz",
            ...     "reference_offset": "0 Hz",
            ...     "magnetic_flux_density": "9.4 T",
            ...     "rotor_frequency": "0 Hz",
            ...     "rotor_angle": "54.935 degree",
            ...     "isotope": "29Si",
            ... }
            >>> dimension_object = Dimension.parse_dict_with_units(dimension_1)
        """
        return super().parse_dict_with_units(py_dict)

    def to_dict_with_units(self):
        """
        Serialize the Dimension object to a JSON compliant python dictionary with
        units.

        Example:
            >>> pprint(dimension_object.to_dict_with_units())
            {'isotope': '29Si',
             'label': '',
             'magnetic_flux_density': '9.4 T',
             'number_of_points': 1024,
             'reference_offset': '0.0 Hz',
             'rotor_angle': '0.958796624583085 rad',
             'rotor_frequency': '0.0 Hz',
             'spectral_width': '100.0 Hz'}
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

    def to_csdm_dimension(self):
        count = self.number_of_points
        dim = csdm_dimension(
            type="linear",
            count=count,
            increment=f"{self.spectral_width/count} Hz",
            coordinates_offset=f"{self.reference_offset} Hz",
            origin_offset=f"{np.abs(self.larmor_frequency)} Hz",
            complex_fft=True,
        )
        return dim


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
