from copy import deepcopy
from typing import ClassVar
from typing import Dict
from typing import List
from typing import Union

import numpy as np
from mrsimulator.utils.parseable import Parseable
from pydantic import Field
from pydantic import validator

from .frequency_contrib import default_freq_contrib
from .frequency_contrib import freq_list_all
from .frequency_contrib import FrequencyEnum
from .query import MixingEnum
from .query import MixingQuery
from .query import TransitionQuery
from .utils import D_symmetry_indexes
from .utils import P_symmetry_indexes

__author__ = "Deepansh J. Srivastava"
__email__ = "srivastava.89@osu.edu"


class BaseEvent(Parseable):
    """Base BaseEvent class. If the value of the attribute is None, the value of the
    corresponding global attribute will be used instead.

    Attributes
    ----------

    magnetic_flux_density:
        The macroscopic magnetic flux density, :math:`H_0`, of the applied external
        magnetic field during the event in units of T. The default value is ``None``.

    rotor_frequency:
        The sample spinning frequency :math:`\nu_r`, during the event in units of Hz.
        The default value is ``None``.

    rotor_angle:
        The angle between the sample rotation axis and the applied external magnetic
        field vector, :math:`\theta`, during the event in units of rad.
        The default value is ``None``, i.e. the magic angle.

    freq_contrib:
        A list of FrequencyEnum enumeration. The default is all frequency enumerations.

    transition_queries:
        A TransitionQuery or an equivalent dict object listing the queries used in
        selecting the active transitions during the event. Only the active transitions
        from this query will contribute to the net frequency.
    """

    magnetic_flux_density: float = Field(default=None, ge=0.0)
    rotor_frequency: float = Field(default=None, ge=0.0)
    rotor_angle: float = Field(default=None, ge=0.0, le=1.5707963268)
    freq_contrib: List[FrequencyEnum] = default_freq_contrib
    transition_queries: List[TransitionQuery] = [TransitionQuery()]

    property_unit_types: ClassVar[Dict] = {
        "magnetic_flux_density": "magnetic flux density",
        "rotor_frequency": "frequency",
        "rotor_angle": "angle",
    }

    property_default_units: ClassVar[Dict] = {
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

    @validator("rotor_frequency", always=True)
    def validate_rotor_frequency(cls, v, **kwargs):
        if v is not None:
            return 1e12 if np.isinf(v) else v
        return v

    @classmethod
    def parse_dict_with_units(cls, py_dict: dict):
        """Parse the physical quantities of an Event object from a python dictionary
        object.

        Args:
            dict py_dict: Dict object
        """
        py_dict_copy = deepcopy(py_dict)
        return super().parse_dict_with_units(py_dict_copy)

    def _freq_contrib_flags(self) -> np.ndarray:
        array = np.zeros(len(freq_list_all), dtype=int)
        array[[item.index() for item in self.freq_contrib]] = 1
        return array

    def combination(self, isotopes, channels):
        """All possible combinations of the event queries over the given channels and
        list of isotopes.

        Args:
            (list) isotopes: List of isotopes in the spin system.
            (list) channels: List of method channels.
        """
        return [
            item.combination(isotopes, channels) for item in self.transition_queries
        ]

    def filter_transitions(self, all_transitions, isotopes, channels):
        """Filter transitions based on the transition query.

        Args:
            (list)  all_transitions: List of all transitions from the spin system.
            (list) isotopes: List of isotopes in the spin system.
            (list) channels: List of method channels.
        """
        symmetry_combinations = self.combination(isotopes, channels)

        segment = []
        for item in symmetry_combinations:
            st = all_transitions[:]
            st = st[P_symmetry_indexes(st, item["P"])] if item["P"].size > 0 else st
            st = st[D_symmetry_indexes(st, item["D"])] if item["D"].size > 0 else st
            segment += [st]
        return np.vstack(segment)


class SpectralEvent(BaseEvent):
    r"""Base SpectralEvent class defines the spin environment and the transition query
    for a segment of the transition pathway.

    Attributes
    ----------

    fraction:
        The weight of the frequency contribution from the event. The default is 1.

    magnetic_flux_density:
        The macroscopic magnetic flux density, :math:`H_0`, of the applied external
        magnetic field during the event in units of T. The default value is ``9.4``.

    rotor_frequency:
        The sample spinning frequency :math:`\nu_r`, during the event in units of Hz.
        The default value is ``0``.

    rotor_angle:
        The angle between the sample rotation axis and the applied external magnetic
        field vector, :math:`\theta`, during the event in units of rad.
        The default value is ``0.9553166``, i.e. the magic angle.

    freq_contrib:
        A list of FrequencyEnum enumeration. The default is all frequency enumerations.

    transition_queries:
        A TransitionQuery or an equivalent dict object listing the queries used in
        selecting the active transitions during the event. Only the active transitions
        from this query will contribute to the net frequency.
    """
    fraction: float = 1.0

    class Config:
        extra = "forbid"
        validate_assignment = True


class ConstantDurationEvent(BaseEvent):  # TransitionModulationEvent
    r"""Base ConstantDurationEvent class defines the spin environment and the
    transition query for a segment of the transition pathway. The frequency from this
    event contribute to the spectrum as amplitudes.

    Attributes
    ----------

    duration:
        The duration of the event in units of µs. The default is 0.

    magnetic_flux_density:
        The macroscopic magnetic flux density, :math:`H_0`, of the applied external
        magnetic field during the event in units of T. The default value is ``9.4``.

    rotor_frequency:
        The sample spinning frequency :math:`\nu_r`, during the event in units of Hz.
        The default value is ``0``.

    rotor_angle:
        The angle between the sample rotation axis and the applied external magnetic
        field vector, :math:`\theta`, during the event in units of rad.
        The default value is ``0.9553166``, i.e. the magic angle.

    freq_contrib:
        A list of FrequencyEnum enumeration. The default is all frequency enumerations.

    transition_queries:
        A TransitionQuery or an equivalent dict object listing the queries used in
        selecting the active transitions during the event. Only the active transitions
        from this query will contribute to the net frequency.
    """
    duration: float

    property_unit_types: ClassVar[Dict] = {
        "duration": "time",
        **BaseEvent.property_unit_types,
    }
    property_default_units: ClassVar[Dict] = {
        "duration": "µs",
        **BaseEvent.property_default_units,
    }
    property_units: Dict = {
        "duration": "µs",
        **BaseEvent().property_default_units,
    }

    test_vars: ClassVar[Dict] = {"duration": 0.0}

    class Config:
        extra = "forbid"
        validate_assignment = True


class MixingEvent(Parseable):  # TransitionMixingEvent
    """Transition mixing class

    Attributes
    ----------

    query:
        The transition mixing query.
    """

    query: Union[MixingQuery, MixingEnum]

    test_vars: ClassVar[Dict] = {"query": {}}

    class Config:
        extra = "forbid"
        validate_assignment = True

    @validator("query", pre=True, always=True)
    def validate_query(cls, v, **kwargs):
        """Validator which tries to convert query to a MixingEnum if query is string"""
        if isinstance(v, str):
            if v in MixingEnum.allowed_enums():
                v = MixingEnum[v]
            else:
                raise ValueError(
                    f"Unrecognized MixingEnum name '{v}'. "
                    f"The allowed types are {MixingEnum.allowed_enums()}"
                )
        return v

    @classmethod
    def parse_dict_with_units(cls, py_dict):
        """
        Parse the physical quantity from a dictionary representation of the MixingEvent
        object, where the physical quantity is expressed as a string with a number and
        a unit.

        Args:
            dict py_dict: A python dict representation of the MixingEvent object.

        Returns:
            A MixingEvent.
        """
        py_dict_copy = deepcopy(py_dict)
        if isinstance(py_dict_copy["query"], dict):
            py_dict_copy["query"] = MixingQuery.parse_dict_with_units(
                py_dict_copy["query"]
            )
        return super().parse_dict_with_units(py_dict_copy)


class Event(Parseable):
    """Event class Object"""

    event: Union[MixingEvent, ConstantDurationEvent, SpectralEvent]
