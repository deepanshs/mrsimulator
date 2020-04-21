# -*- coding: utf-8 -*-
"""The Event class."""
from copy import deepcopy
from typing import ClassVar
from typing import Dict
from typing import List
from typing import Optional

import numpy as np
from csdmpy import Dimension as csdm_dimension
from mrsimulator import Parseable
from pydantic import BaseModel
from pydantic import Field
from pydantic import validator

from .isotope import Isotope
from .transition import Transition


class TransitionQuery(BaseModel):
    P: List[float] = [-1.0]
    D: List[float] = Field(default=None)
    f: Optional[float] = Field(default=None)
    transitions: List[Transition] = None

    class Config:
        validate_assignment = True
        arbitrary_types_allowed = True

    def to_dict_with_units(self):
        return {k: v for k, v in self.dict().items() if v is not None}


class Event(Parseable):
    fraction: float = 1
    magnetic_flux_density: Optional[float] = Field(default=9.4, ge=0)
    rotor_frequency: Optional[float] = Field(default=0, ge=0)
    # 54.735 degrees = 0.9553166 radians
    rotor_angle: Optional[float] = Field(default=0.9553166, ge=0, le=1.5707963268)
    transition_query: Optional[TransitionQuery] = None
    # transitions: List[float] = [[-0.5, 0.5]]
    # transitions: List[Transition]

    property_unit_types: ClassVar = {
        "magnetic_flux_density": "magnetic flux density",
        "rotor_frequency": "frequency",
        "rotor_angle": "angle",
    }

    property_default_units: ClassVar = {
        "magnetic_flux_density": "T",
        "rotor_frequency": "Hz",
        "rotor_angle": "rad",
    }

    property_units: Dict = {
        "magnetic_flux_density": "T",
        "rotor_frequency": "Hz",
        "rotor_angle": "rad",
    }

    class Config:
        validate_assignment = True
        # arbitrary_types_allowed = True

    @classmethod
    def parse_dict_with_units(cls, py_dict):
        """
        Parse the physical quantities of a Event object from as a python dictionary.

        Args:
            py_dict: Dict object
        """
        return super().parse_dict_with_units(py_dict)

    # def to_dict_with_units(self):
    #     """
    #     Serialize the Event object to a JSON compliant python dictionary with
    #     units.
    #     """
    #     return self._to_dict_with_units(self)
    #     # temp_dict["transitions"] = self.transitions


class GlobalEvent(Event):
    pass


class Sequence(Parseable):
    count: int = Field(1024, gt=0)
    spectral_width: float = Field(..., gt=0)
    reference_offset: Optional[float] = Field(default=0)
    label: Optional[str] = ""
    events: List[Event] = None

    property_unit_types: ClassVar = {
        "spectral_width": ["frequency", "dimensionless"],
        "reference_offset": ["frequency", "dimensionless"],
    }

    property_default_units: ClassVar = {
        "spectral_width": ["Hz", "ppm"],
        "reference_offset": ["Hz", "ppm"],
    }

    property_units: Dict = {"spectral_width": "Hz", "reference_offset": "Hz"}

    class Config:
        validate_assignment = True

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
        n = self.count
        Tk = int(n / 2)
        increment = self.spectral_width / self.count
        return (np.arange(n) - Tk) * increment + self.reference_offset

    # @property
    # def coordinates_ppm(self):
    #     r"""
    #     The grid coordinates along the dimension as dimension frequency ratio
    #     in units of ppm. The coordinates are evaluated as

    #     .. math::
    #         x_\text{ppm} = \frac{x_\text{Hz}} {x_0 + \omega_0}

    #     where :math:`\omega_0` is the Larmor frequency.

    #     Example:
    #         >>> dim.coordinates_ppm[:5]
    #         array([-239.53462217, -239.06678111, -238.59894005, -238.13109899,
    #                -237.66325794])
    #     """
    #     if self.isotope is None:
    #         warnings.warn(
    #             (
    #                 "The coordinates along the dimension without an assigned "
    #                 "isotope cannot be converted to dimensionless frequency ratio."
    #             )
    #         )
    #     else:

    #         denominator = (self.coordinates_offset + self.larmor_frequency) / 1e6
    #         return self.coordinates_Hz / abs(denominator)

    @classmethod
    def parse_dict_with_units(cls, py_dict):
        """
        Parse the physical quantities of a Dimension object from as a python dictionary.

        Args:
            py_dict: Dict object
        """
        py_dict_copy = deepcopy(py_dict)
        if "events" in py_dict_copy:
            py_dict_copy["events"] = [
                Event.parse_dict_with_units(e) for e in py_dict_copy["events"]
            ]

        return super().parse_dict_with_units(py_dict_copy)

    # def to_dict_with_units(self):
    #     """
    #     Serialize the Dimension object to a JSON compliant python dictionary with
    #     units.

    #     """
    #     return self._to_dict_with_units(self)

    def to_csdm_dimension(self):
        increment = self.spectral_width / self.count
        dim = csdm_dimension(
            type="linear",
            count=self.count,
            increment=f"{increment} Hz",
            coordinates_offset=f"{self.reference_offset} Hz",
            # origin_offset=f"{np.abs(self.larmor_frequency)} Hz",
            complex_fft=True,
            reciprocal={"coordinates_offset": f"{-1/(2*increment)} s"},
        )
        return dim


class Method(Parseable):
    name: Optional[str] = ""
    description: Optional[str] = ""
    isotope: Optional[str] = None
    sequences: List[Sequence]

    class Config:
        validate_assignment = True

    @validator("isotope", always=True)
    def get_isotope(cls, v, *, values, **kwargs):
        if v is None:
            return v
        return Isotope(symbol=v)

    @classmethod
    def parse_dict_with_units(cls, py_dict):
        """
        Parse the physical quantities of a Event object from as a python dictionary.

        Args:
            py_dict: Dict object
        """
        py_dict_copy = deepcopy(py_dict)
        if "sequences" in py_dict_copy:
            py_dict_copy["sequences"] = [
                Sequence.parse_dict_with_units(s) for s in py_dict_copy["sequences"]
            ]
        return super().parse_dict_with_units(py_dict_copy)

    def to_dict_with_units(self):
        temp_dict = self.dict(exclude={"sequences", "isotope"})
        temp_dict["sequences"] = [item.to_dict_with_units() for item in self.sequences]
        temp_dict["isotope"] = self.isotope.to_dict_with_units()
        return temp_dict

    def get_transition_pathways(self, isotopomer):
        transitions = isotopomer.all_transitions
        # print(transitions)
        segments = []
        for seq in self.sequences:
            for ent in seq.events:
                segments.append(
                    np.asarray(
                        transitions.filter(**ent.transition_query.to_dict_with_units())
                    )
                )
        # print(segments)
        return cartesian_product(*segments)
        # pathways = 1
        # for item in segments:
        #     pathways = np.kron(pathways, np.asarray(item))
        # return pathways


def cartesian_product(*arrays):
    la = len(arrays)
    dtype = np.result_type(*arrays)
    arr = np.empty([len(a) for a in arrays] + [la], dtype=dtype)
    for i, a in enumerate(np.ix_(*arrays)):
        arr[..., i] = a
    return arr.reshape(-1, la)
