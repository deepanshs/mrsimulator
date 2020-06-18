# -*- coding: utf-8 -*-
import warnings
from copy import deepcopy
from typing import ClassVar
from typing import Dict
from typing import List

import csdmpy as cp
import numpy as np
from mrsimulator.utils.parseable import Parseable
from pydantic import Field

from .event import Event

__author__ = "Deepansh J. Srivastava"
__email__ = "srivastava.89@osu.edu"


class SpectralDimension(Parseable):
    r"""Base SpectralDimension class defines a spectroscopic dimension of the method.

    Attributes
    ----------

    count: int (optional).
        The number of points, :math:`N`, along the spectroscopic dimension. The default
        value is 1024.

    spectral_width: float (optional).
        The spectral width, :math:`\Delta x`, of the spectroscopic dimension in units
        of Hz. The default value is 25000.

    reference_offset: float (optional).
        The reference offset, :math:`x_0`, of the spectroscopic dimension in units of
        Hz. The default value is 0.

    origin_offset: float (optional).
        The origin offset (Larmor frequency) along the spectroscopic dimension in units
        of Hz. The default value is None. When the value is None, the origin offset is
        set to the Larmor frequency of the isotope from the
        :attr:`~mrsimulator.Method.channels` attribute of the method.

    label: str (optional).
        The value is a label of the spectroscopic dimension. The default value is None.

    description: str (optional).
        The value is a description of the spectroscopic dimension. The default value is
        None.

    events: A list of :ref:`event_api` or equivalent dict objects (optional).
        The value describes a series of events along the spectroscopic dimension.
    """

    count: int = Field(1024, gt=0)
    spectral_width: float = Field(default=25000.0, gt=0)
    reference_offset: float = Field(default=0.0)
    origin_offset: float = None
    label: str = None
    description: str = None
    events: List[Event] = []

    property_unit_types: ClassVar = {
        "spectral_width": ["frequency", "dimensionless"],
        "reference_offset": ["frequency", "dimensionless"],
        "origin_offset": ["frequency", "dimensionless"],
    }

    property_default_units: ClassVar = {
        "spectral_width": ["Hz", "ppm"],
        "reference_offset": ["Hz", "ppm"],
        "origin_offset": ["Hz", "ppm"],
    }

    property_units: Dict = {
        "spectral_width": "Hz",
        "reference_offset": "Hz",
        "origin_offset": "Hz",
    }

    class Config:
        validate_assignment = True

    @classmethod
    def parse_dict_with_units(cls, py_dict: dict):
        """
        Parse the physical quantities of a SpectralDimension object from a
        python dictionary object.

        Args:
            dict py_dict: Dict object
        """
        py_dict_copy = deepcopy(py_dict)
        if "events" in py_dict_copy:
            py_dict_copy["events"] = [
                Event.parse_dict_with_units(e) for e in py_dict_copy["events"]
            ]

        return super().parse_dict_with_units(py_dict_copy)

    def coordinates_Hz(self) -> np.ndarray:
        r"""
        The grid coordinates along the dimension in units of Hz, evaluated as

        .. math::
            x_\text{Hz} = \left([0, 1, ... N-1] - T\right) \frac{\Delta x}{N} + x_0

        where :math:`T=N/2` and :math:`T=(N-1)/2` for even and odd values of
        :math:`N`, respectively.
        """
        n = self.count
        Tk = int(n / 2)
        increment = self.spectral_width / self.count
        return (np.arange(n) - Tk) * increment + self.reference_offset

    def coordinates_ppm(self) -> np.ndarray:
        r"""
        The grid coordinates along the dimension as dimension frequency ratio
        in units of ppm. The coordinates are evaluated as

        .. math::
            x_\text{ppm} = \frac{x_\text{Hz}} {x_0 + \omega_0}

        where :math:`\omega_0` is the Larmor frequency.
        """
        if self.origin_offset is None:
            warnings.warn(
                UserWarning(
                    "The coordinates along the dimension without an origin offset "
                    "cannot be converted to dimensionless frequency ratio."
                )
            )
        else:
            denominator = (self.reference_offset + self.origin_offset) / 1e6
            return self.coordinates_Hz() / abs(denominator)

    def to_csdm_dimension(self) -> cp.Dimension:
        """Return the spectral dimension as a CSDM dimension object."""
        increment = self.spectral_width / self.count
        label = "" if self.label is None else self.label
        description = "" if self.description is None else self.description
        dim = cp.Dimension(
            type="linear",
            count=self.count,
            increment=f"{increment} Hz",
            coordinates_offset=f"{self.reference_offset} Hz",
            label=label,
            description=description,
            complex_fft=True,
            reciprocal={"coordinates_offset": f"{-1/(2*increment)} s"},
        )
        if self.origin_offset is not None:
            dim.origin_offset = f"{self.origin_offset} Hz"
        return dim
