# -*- coding: utf-8 -*-
import warnings
from copy import deepcopy
from typing import ClassVar
from typing import Dict
from typing import List
from typing import Optional

import csdmpy as cp
import numpy as np
from mrsimulator.parseable import Parseable
from pydantic import Field

from .event import Event

__author__ = "Deepansh J. Srivastava"
__email__ = "srivastava.89@osu.edu"


class SpectralDimension(Parseable):
    r"""Base SpectralDimension class defines the dimensions of the method.

    Attributes:
        count: An optional integer with the number of points, :math:`N`, along the
            dimension. The default value is 1024.
        spectral_width: An optional float with the spectral width, :math:`\Delta x`,
            along the dimension in units of Hz. The default value is 25 kHz.
        reference_offset: An `optional` float with the reference offset, :math:`x_0`
                along the dimension in units of Hz. The default value is 0.
        origin_offset: An `optional` float with the origin offset (Larmor frequency)
                along the dimension in units of Hz. The default value is None.
        label: An `optional` string label. The default is an empty string.
        events: A `required` list of Event object or an equivalent list of dict objects
                describing the series of events along the spectroscopic dimension.
    """

    count: int = Field(1024, gt=0)
    spectral_width: float = Field(default=25000.0, gt=0)
    reference_offset: Optional[float] = Field(default=0.0)
    origin_offset: Optional[float] = None
    label: Optional[str] = ""
    events: List[Event]

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

    @property
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

    @property
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
            return self.coordinates_Hz / abs(denominator)

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

    def to_csdm_dimension(self) -> cp.Dimension:
        """Return the spectral dimension as a CSDM dimension object."""
        increment = self.spectral_width / self.count
        dim = cp.Dimension(
            type="linear",
            count=self.count,
            increment=f"{increment} Hz",
            coordinates_offset=f"{self.reference_offset} Hz",
            complex_fft=True,
            reciprocal={"coordinates_offset": f"{-1/(2*increment)} s"},
        )
        if self.origin_offset is not None:
            dim.origin_offset = f"{self.origin_offset} Hz"
        return dim
