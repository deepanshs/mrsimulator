import warnings
from copy import deepcopy
from typing import ClassVar
from typing import Dict
from typing import List
from typing import Union
from warnings import warn

import csdmpy as cp
import numpy as np
from mrsimulator.utils.error import MissingSpectralEventError
from mrsimulator.utils.parseable import Parseable
from pydantic import Field
from pydantic import validator

from .event import ConstantDurationEvent
from .event import Event
from .event import MixingEvent
from .event import SpectralEvent
from .utils import cartesian_product


__author__ = "Deepansh J. Srivastava"
__email__ = "srivastava.89@osu.edu"

CHANNELS = ["ch1", "ch2", "ch3"]


class Reciprocal(Parseable):
    """Reciprocal dimension from CSDM object."""

    coordinates_offset: float = 0.0

    property_unit_types: ClassVar[Dict] = {"coordinates_offset": "time"}
    property_default_units: ClassVar[Dict] = {"coordinates_offset": "s"}
    property_units: Dict = {"coordinates_offset": "s"}


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
    spectral_width: float = Field(default=25000.0)
    reference_offset: float = Field(default=0.0)
    origin_offset: float = None
    reciprocal: Reciprocal = None
    events: List[Union[MixingEvent, ConstantDurationEvent, SpectralEvent]] = []

    property_unit_types: ClassVar[Dict] = {
        "spectral_width": ["frequency", "dimensionless"],
        "reference_offset": ["frequency", "dimensionless"],
        "origin_offset": ["frequency", "dimensionless"],
    }

    property_default_units: ClassVar[Dict] = {
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
        extra = "forbid"
        validate_assignment = True

    @validator("spectral_width", pre=True, always=True)
    def validate_spectral_width(cls, value):
        """Spectral width cannot be zero."""
        if value != 0:
            return value
        raise ValueError("Spectral width cannot be zero.")

    @validator("events", pre=True, always=True)
    def validate_events(v, **kwargs):
        """Ensure at least one spectralEvent and warn is the sum of fraction in
        SpectralEvents is not 1."""
        if v != []:
            new_v = [
                Event(event=item).event if isinstance(item, dict) else item
                for item in v
            ]
            fractions = [e.fraction for e in new_v if isinstance(e, SpectralEvent)]
            if len(fractions) == 0:  # No SpectralEvent in dimension
                raise MissingSpectralEventError()
            total = sum(fractions)
            if total != 1:
                e = (
                    "The fraction attribute of each SpectralEvent in a "
                    f"SpectralDimension should sum to 1. Sum is {total}."
                    "If this was not intentional, check the fraction attributes."
                )
                warn(e)
        return v

    @classmethod
    def parse_dict_with_units(cls, py_dict: dict):
        """Parse the physical quantities of a SpectralDimension object from a python
        dictionary object.

        Args:
            dict py_dict: Dict object
        """
        py_dict_copy = deepcopy(py_dict)
        if "events" in py_dict_copy:
            py_dict_copy["events"] = [
                ConstantDurationEvent.parse_dict_with_units(e)
                if "duration" in e
                else SpectralEvent.parse_dict_with_units(e)
                for e in py_dict_copy["events"]
            ]

        return super().parse_dict_with_units(py_dict_copy)

    def coordinates_Hz(self) -> np.ndarray:
        r"""The grid coordinates along the dimension in units of Hz, evaluated as

        .. math::
            x_\text{Hz} = \left([0, 1, ... N-1] - T\right) \frac{\Delta x}{N} + x_0

        where :math:`T=N/2` and :math:`T=(N-1)/2` for even and odd values of
        :math:`N`, respectively."""
        n = self.count
        Tk = int(n / 2)
        increment = self.spectral_width / self.count
        return (np.arange(n) - Tk) * increment + self.reference_offset

    def coordinates_ppm(self) -> np.ndarray:
        r"""The grid coordinates along the dimension as dimension frequency ratio
        in units of ppm. The coordinates are evaluated as

        .. math::
            x_\text{ppm} = \frac{x_\text{Hz}} {x_0 + \omega_0}

        where :math:`\omega_0` is the Larmor frequency."""
        if self.origin_offset is None:
            warnings.warn(
                UserWarning(
                    "The coordinates along the dimension without an origin offset "
                    "cannot be converted to dimensionless frequency ratio."
                )
            )
            return

        denominator = (self.origin_offset - self.reference_offset) / 1e6
        return self.coordinates_Hz() / abs(denominator)

    def to_csdm_dimension(self) -> cp.Dimension:
        """Return the spectral dimension as a CSDM dimension object."""
        increment = self.spectral_width / self.count
        label = "" if self.label is None else self.label
        description = "" if self.description is None else self.description

        default_reciprocal = {"coordinates_offset": f"{-1/(2*increment)} s"}
        reciprocal = (
            default_reciprocal if self.reciprocal is None else self.reciprocal.json()
        )

        dim = cp.Dimension(
            type="linear",
            count=self.count,
            increment=f"{increment} Hz",
            coordinates_offset=f"{self.reference_offset} Hz",
            label=label,
            description=description,
            complex_fft=True,
            reciprocal=reciprocal,
        )
        if self.origin_offset is not None:
            dim.origin_offset = f"{self.origin_offset} Hz"
        return dim

    # def events_to_dataframe(self) -> pd.DataFrame:
    #     """Returns events list as DataFrame with event number as columns"""
    #     attributes = list(Event().property_units.keys())
    #     attributes.append("fraction")
    #     rows = attributes.copy()
    #     rows.extend(["p", "d"])
    #     df = pd.DataFrame(index=rows)
    #     for i in range(len(self.events)):
    #         _lst = [getattr(self.events[i], att) for att in attributes]
    #         _lst.append(self.events[i].transition_queries.get_p())
    #         _lst.append(self.events[i].transition_queries.get_d())
    #         df[i] = _lst
    #     return df

    def _get_symmetry_pathways(self, symmetry_element: str) -> list:
        """Generate a list of symmetry pathways for the event.

        The output is as follows
        [
            {"ch1": [ symmetry pathway for ch1], "ch2": ..} 1st symmetry pathway
            {"ch1": [ symmetry pathway for ch1], "ch2": ..} 2nd symmetry pathway
        ]

        Args:
            symmetry_element: Symmetry symbol, "P" or "D"

        Example:
            >>> from mrsimulator.method import SpectralDimension
            >>> sp = SpectralDimension(
            ...     events = [{
            ...         "fraction": 0.5,
            ...         "transition_queries": [
            ...             {"ch1": {"P": [1, 1]}, "ch2": {"P": [1], "D": [2]}},
            ...             {"ch1": {"P": [-1, -1]}},
            ...         ]
            ...     },
            ...     {
            ...         "fraction": 0.5,
            ...         "transition_queries": [
            ...             {"ch1": {"P": [-1]}},
            ...         ]
            ...     }]
            ... )
            >>> pprint(sp._get_symmetry_pathways("P"))
            [{'ch1': [[1, 1], [-1]], 'ch2': [[1], None], 'ch3': [None, None]},
             {'ch1': [[-1, -1], [-1]], 'ch2': [None, None], 'ch3': [None, None]}]
        """
        ha, ga = hasattr, getattr
        tq, de = "transition_queries", np.asarray([0])

        indexes = [np.arange(len(ga(e, tq))) if ha(e, tq) else de for e in self.events]
        products = cartesian_product(*indexes)

        return [
            {
                ch: [
                    ga(ga(e.transition_queries[i], ch), symmetry_element)
                    if ha(e, tq) and ga(e.transition_queries[i], ch) is not None
                    else None
                    for e, i in zip(self.events, item)
                ]
                for ch in CHANNELS
            }
            for item in products
        ]
