from copy import deepcopy
from typing import ClassVar
from typing import Dict
from typing import List
from typing import Optional
from typing import Union

import numpy as np
from mrsimulator.utils.parseable import Parseable
from pydantic.v1 import Field
from pydantic.v1 import validator

from .frequency_contrib import default_freq_contrib
from .frequency_contrib import FREQ_ENUM_SHORTCUT
from .frequency_contrib import FREQ_LIST_ALL
from .frequency_contrib import FrequencyEnum
from .query import Rotation
from .query import TransitionQuery
from .utils import D_symmetry_indexes
from .utils import P_symmetry_indexes

# from .query import MixingEnum

__author__ = "Deepansh J. Srivastava"
__email__ = "srivastava.89@osu.edu"


def parse_dict_to_ev_class(py_dict: dict):
    """Takes in a JSON representation of some Event class and parses it to an object of
    the correct Event class

    Arguments:
        (dict) py_dict: JSON representation of the Event class as a dictionary

    Returns:
        Either a SpectralEvent, DelayEvent, or RotationEvent object
    """
    if "duration" in py_dict:
        return DelayEvent.parse_dict_with_units(py_dict)

    if "ch1" in py_dict or "ch2" in py_dict or "ch3" in py_dict:
        return RotationEvent.parse_dict_with_units(py_dict)

    return SpectralEvent.parse_dict_with_units(py_dict)


class FreeEvent(Parseable):
    """Base FreeEvent class. If the value of the attribute is None, the value of the
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
        A list of TransitionQuery or equivalent dict objects. The queries are used in
        selecting the transitions during the event. Only selected transitions
        from this query will contribute to the net frequency.
    """

    magnetic_flux_density: float = Field(default=None, ge=0.0)
    rotor_frequency: float = Field(default=None, ge=0.0)
    rotor_angle: float = Field(default=None, ge=0.0, le=1.5707963268)
    freq_contrib: List[Union[FrequencyEnum, str]] = default_freq_contrib
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

    @validator("freq_contrib", always=True, pre=True)
    def validate_freq_contrib(cls, v, **kwargs):
        additive = set()
        subtractive = set()
        for item in v:
            if isinstance(item, FrequencyEnum):
                additive.add(item.value)

            function = subtractive if item.startswith("!") else additive
            item_enum = item[1:] if item.startswith("!") else item

            # Item is a string of a recognized freq contrib string
            if item_enum in FREQ_LIST_ALL:
                function.add(item_enum)

            # Item is a string of a recognized freq contrib shortcut
            if item_enum in FREQ_ENUM_SHORTCUT:
                function.update(FREQ_ENUM_SHORTCUT[item_enum])

        # Set additive to all frequency enumerations if none passed (default)
        additive = set(FREQ_LIST_ALL) if additive == set() else additive
        return [FrequencyEnum(item) for item in list(additive - subtractive)]

    @classmethod
    def parse_dict_with_units(cls, py_dict: dict):
        """Parse the physical quantities of an Event object from a Python dictionary
        object.

        Args:
            dict py_dict: Dict object
        """
        py_dict_copy = deepcopy(py_dict)
        return super().parse_dict_with_units(py_dict_copy)

    def dict(self, **kwargs) -> dict:
        """Return a JSON-compliant dictionary of the instance of the event."""
        py_dict = super().dict()
        py_dict["freq_contrib"] = [
            fq.value if isinstance(fq, FrequencyEnum) else fq
            for fq in py_dict["freq_contrib"]
        ]
        return py_dict

    def _freq_contrib_flags(self) -> np.ndarray:
        """Get an array (binary) of frequency contributions flags for the event."""
        array = np.zeros(len(FREQ_LIST_ALL), dtype=int)  # Final array to return
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


class SpectralEvent(FreeEvent):
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


class DelayEvent(FreeEvent):
    r"""Base DelayEvent class defines the spin environment and the
    transition query for a segment of the transition pathway. The frequency from this
    event contributes to the spectrum as complex amplitude modulations.

    Attributes
    ----------

    duration:
        The duration of the event in units of s. The default is 0.

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
        **FreeEvent.property_unit_types,
    }
    property_default_units: ClassVar[Dict] = {
        "duration": "s",
        **FreeEvent.property_default_units,
    }
    property_units: Dict = {
        "duration": "s",
        **FreeEvent().property_default_units,
    }

    test_vars: ClassVar[Dict] = {"duration": 0.0}

    class Config:
        extra = "forbid"
        validate_assignment = True


class MixingEvent(Parseable):
    """MixingEvent class for querying transition mixing between events.

    Attributes
    ----------

    ch1:
        An optional Rotation object for the channel at index 0 of the method's
        channels list."

    ch2:
        An optional Rotation object for the channel at index 1 of the method's
        channels list."

    ch3:
        An optional Rotation object for the channel at index 2 of the method's
        channels list."

    Example
    -------

        >>> query = MixingEvent(ch1={"angle": 1.570796, "phase": 3.141593})

    """

    ch1: Optional[Rotation] = Field(
        title="ch1",
        default=Rotation(),
        description=(
            "An optional Rotation object for imposing a rotation on "
            "channel index 0 of the method's channels array."
        ),
    )
    ch2: Optional[Rotation] = Field(
        title="ch2",
        default=None,
        description=(
            "An optional Rotation object for imposing a rotation on "
            "channel index 0 of the method's channels array."
        ),
    )
    ch3: Optional[Rotation] = Field(
        title="ch3",
        default=None,
        description=(
            "An optional Rotation object for imposing a rotation on "
            "channel index 0 of the method's channels array."
        ),
    )

    class Config:
        validate_assignment = True
        extra = "forbid"

    @classmethod
    def parse_dict_with_units(cls, py_dict):
        """
        Parse the physical quantity from a dictionary representation of the Method
        object, where the physical quantity is expressed as a string with a number and
        a unit.

        Args:
            dict py_dict: A Python dict representation of the MixingEvent object.

        Returns:
            A :ref:`method_api` object.
        """
        py_dict_copy = deepcopy(py_dict)
        obj = {k: Rotation.parse_dict_with_units(v) for k, v in py_dict_copy.items()}
        py_dict_copy.update(obj)
        return super().parse_dict_with_units(py_dict_copy)

    @property
    def channels(self) -> List[Rotation]:
        """Returns an ordered list of all channels"""
        return [self.ch1, self.ch2, self.ch3]


class RotationEvent(MixingEvent):
    """Rotation Event class. Same as mixing event"""

    pass


class Event(Parseable):
    """Event class Object"""

    event: Union[DelayEvent, SpectralEvent, MixingEvent]
