# -*- coding: utf-8 -*-
from copy import deepcopy
from typing import ClassVar
from typing import Dict
from typing import List
from typing import Union

import numpy as np
from mrsimulator.utils.parseable import Parseable
from pydantic import Field

from .frequency_contrib import default_freq_contrib
from .frequency_contrib import freq_default
from .frequency_contrib import freq_list_all
from .frequency_contrib import FrequencyEnum
from .transition_query import TransitionQuery
from .utils import D_symmetry_indexes
from .utils import P_symmetry_indexes

__author__ = "Deepansh J. Srivastava"
__email__ = "srivastava.89@osu.edu"


class BaseEvent(Parseable):
    """Base BaseEvent class.

    Attributes
    ----------

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

    transition_query:
        A TransitionQuery or an equivalent dict object listing the queries used in
        selecting the active transitions during the event. Only the active transitions
        from this query will contribute to the net frequency.
    """

    magnetic_flux_density: float = Field(default=9.4, ge=0.0)
    rotor_frequency: float = Field(default=0.0, ge=0.0)
    # 54.735 degrees = 0.9553166 radians
    rotor_angle: float = Field(default=0.955316618, ge=0.0, le=1.5707963268)
    freq_contrib: List[FrequencyEnum] = default_freq_contrib
    transition_query: List[TransitionQuery] = [TransitionQuery()]

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

    @classmethod
    def parse_dict_with_units(cls, py_dict: dict):
        """
        Parse the physical quantities of an Event object from a python dictionary
        object.

        Args:
            dict py_dict: Dict object
        """
        py_dict_copy = deepcopy(py_dict)
        return super().parse_dict_with_units(py_dict_copy)

    def json(self) -> dict:
        """Parse the class object to a JSON compliant python dictionary object, where
        the attribute value with physical quantity is expressed as a string with a
        value and a unit."""
        dict_ = super().json()
        # if "user_variables" in dict_.keys():
        #     dict_.pop("user_variables")
        if "fraction" in dict_:
            dict_.pop("fraction") if dict_["fraction"] == 1.0 else None
        dict_.pop("freq_contrib") if dict_["freq_contrib"] == freq_default else None
        return dict_

    def _freq_contrib_flags(self) -> np.ndarray:
        lst_ = [item.index() for item in self.freq_contrib]
        array = np.zeros(len(freq_list_all), dtype=int)
        array[lst_] = 1
        return array

    def permutation(self, isotopes, channels):
        """Permutate the event queries over the given channels and list of isotopes.

        Args:
            (list) isotopes: List of isotopes in the spin system.
            (list) channels: List of method channels.
        """
        return [item.permutation(isotopes, channels) for item in self.transition_query]

    def filter_transitions(self, all_transitions, isotopes, channels):
        """Filter transitions based on the transition query.

        Args:
            (list)  all_transitions: List of all transitions from the spin system.
            (list) isotopes: List of isotopes in the spin system.
            (list) channels: List of method channels.
        """

        symmetry_permutations = self.permutation(isotopes, channels)

        segment = []
        for item in symmetry_permutations:
            selected_transitions = all_transitions[:]
            if item["P"].size > 0:
                indexes = P_symmetry_indexes(selected_transitions, item["P"])
                selected_transitions = selected_transitions[indexes]

            if item["D"].size > 0:
                indexes = D_symmetry_indexes(selected_transitions, item["D"])
                selected_transitions = selected_transitions[indexes]

            segment += [selected_transitions]
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

    transition_query:
        A TransitionQuery or an equivalent dict object listing the queries used in
        selecting the active transitions during the event. Only the active transitions
        from this query will contribute to the net frequency.
    """
    fraction: float = 1.0

    class Config:
        validate_assignment = True


class ConstantDurationEvent(BaseEvent):
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

    transition_query:
        A TransitionQuery or an equivalent dict object listing the queries used in
        selecting the active transitions during the event. Only the active transitions
        from this query will contribute to the net frequency.
    """
    duration: float

    property_unit_types: ClassVar = {
        "duration": "time",
        **BaseEvent.property_unit_types,
    }
    property_default_units: ClassVar = {
        "duration": "µs",
        **BaseEvent.property_default_units,
    }
    property_units: Dict = {
        "duration": "µs",
        **BaseEvent().property_default_units,
    }

    class Config:
        validate_assignment = True


class Event(Parseable):
    event: Union[ConstantDurationEvent, SpectralEvent]
