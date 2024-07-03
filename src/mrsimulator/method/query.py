# from copy import deepcopy
# from enum import Enum
from itertools import permutations
from typing import ClassVar
from typing import Dict
from typing import List
from typing import Optional

import numpy as np
from mrsimulator.transition import Transition
from mrsimulator.utils.parseable import Parseable
from pydantic.v1 import Field

from .utils import cartesian_product
from .utils import get_iso_dict

__author__ = "Deepansh J. Srivastava"
__email__ = "srivastava.89@osu.edu"

ON_FAIL_MESSAGE = (
    "The length of the transition query symmetry elements cannot exceed the "
    "number of channels."
)

SYMMETRY_ELEMENTS = ["P", "D"]


class SymmetryQuery(Parseable):
    """Base SymmetryQuery class.

    Attributes
    ----------

    P:
        A list of p symmetry functions per site. Here p = Δm = :math:`m_f - m_i` is the
        difference between the spin quantum numbers of the final and initial states.

        Example
        -------

        >>> method = Method(channels=['1H'], spectral_dimensions=[{"events": [
        ...     {"fraction": 1}
        ... ]}])
        >>> method.spectral_dimensions[0].events[0].transition_queries[0].ch1.P = [-1]

    D:
        A list of d symmetry functions per site. Here d = :math:`m_f^2 - m_i^2` is the
        difference between the square of the spin quantum numbers of the final and
        initial states.

        Example
        -------

        >>> method.spectral_dimensions[0].events[0].transition_queries[0].ch1.D = [0]

    """

    P: List[int] = Field(
        default=[0],
        description=(
            "A list of p symmetry functions per site. Here p = Δm = (m_f-m_i) is the "
            "difference between spin quantum numbers of the final and initial states."
        ),
    )
    D: List[int] = Field(
        default=None,
        description=(
            "A list of d symmetry functions per site. Here d = m_f^2 - m_i^2 is the "
            "difference between the square of the spin quantum numbers of the final "
            "and initial states."
        ),
    )
    F: List[float] = Field(default=None)
    transitions: List[Transition] = None

    class Config:
        validate_assignment = True

    def query_combination(self, symmetry, n_site_at_channel_id):
        """Combination of symmetry query based on the number of sites in a channel.

        Args:
            (str) symmetry: The symmetry element, 'P' or 'D'.
            (int) n_site_at_channel: Number of sites for the given channel.

        Example:
            Consider the following
                query = {P: [-1], D: [1]}
                n_isotopes = [3]
                channels = ['A']
        then,
        1. P query will expand to [-1, 0, 0], [0, -1, 0], and [0, 0, -1] combinations
        2. D query will expand to [1, 0, 0], [0, 1, 0], and [0, 0, 1] combinations
        """
        query = getattr(self, symmetry)
        if query is None:
            query = [0 if symmetry == "P" else None] * n_site_at_channel_id
            return list(set(permutations(query)))
        query_size = len(query)
        if query_size > n_site_at_channel_id:
            return []
            # raise ValueError(ON_FAIL_MESSAGE)
        query = query + (n_site_at_channel_id - query_size) * [0]
        return list(set(permutations(query)))


class TransitionQuery(Parseable):
    """TransitionQuery class for querying transition symmetry function.

    Attributes
    ----------

    ch1:
        An optional SymmetryQuery object for querying symmetry functions at the channel
        index 0 of the method's channels array."

    ch2:
        An optional SymmetryQuery object for querying symmetry functions at the channel
        index 1 of the method's channels array."

    ch3:
        An optional SymmetryQuery object for querying symmetry functions at the channel
        index 2 of the method's channels array."

    Example
    -------

        >>> query = TransitionQuery(ch1={'P': [1], 'D': [0]}, ch2={'P': [-1]})

    """

    ch1: Optional[SymmetryQuery] = Field(
        title="ch1",
        default=SymmetryQuery(),
        description=(
            "An optional SymmetryQuery object for querying symmetry functions at "
            "channel index 0 of the method's channels array."
        ),
    )
    ch2: Optional[SymmetryQuery] = Field(
        title="ch2",
        default=None,
        description=(
            "An optional SymmetryQuery object for querying symmetry functions at "
            "channel index 1 of the method's channels array."
        ),
    )
    ch3: Optional[SymmetryQuery] = Field(
        title="ch3",
        default=None,
        description=(
            "An optional SymmetryQuery object for querying symmetry functions at "
            "channel index 2 of the method's channels array."
        ),
    )

    class Config:
        validate_assignment = True
        extra = "forbid"

    @staticmethod
    def cartesian_product_indexing(combinations):
        """Return Cartesian product of indexes"""
        combination_length = [np.arange(len(p)) for p in combinations]
        # if combination_length == []:
        #     return np.asarray([])

        index = cartesian_product(*combination_length)  # cartesian indexes
        return np.asarray(
            [np.hstack([combinations[i][j] for i, j in enumerate(ix)]) for ix in index]
        )

    def combination(self, isotopes, channels):
        """Combinations of TransitionQuery based on the number of sites per channel.

        Args:
            (list) isotopes: List of isotope symbols, ['29Si, '13C', '13C', '1H'].
            (int) channels: List of method channels, ['29Si, '13C'].
        """
        iso_dict = get_iso_dict(channels=channels, isotopes=isotopes)
        sites_per_channel = [
            iso_dict[item].size if item in iso_dict else 0 for item in channels
        ]

        expanded_symmetry = {
            sym: self.expand_elements_for_symmetry(
                sym, isotopes, iso_dict, channels, sites_per_channel
            )
            for sym in SYMMETRY_ELEMENTS
        }

        # remove nan queries
        for key, value in expanded_symmetry.items():
            if np.all(np.isnan(value)):
                expanded_symmetry[key] = np.asarray([])

        return expanded_symmetry

    def _get_missing_channel_isotope(self, isotopes, channels):
        missing_ch = []
        unique_isotopes = list(set(isotopes))
        unique_channels = list(set(channels))
        for isotope in unique_isotopes:
            if isotope not in unique_channels:
                missing_ch.append(isotope)
        return missing_ch

    def expand_elements_for_symmetry(
        self, symmetry, isotopes, iso_dict, channels, sites_per_channel
    ):
        sym_combination = []
        live_ch, live_ch_index = [], []
        live_n_sites = []
        for i, channel_id in enumerate(channels):
            ch_obj = getattr(self, f"ch{i+1}")
            ch_obj = ch_obj if ch_obj is not None else SymmetryQuery()
            if channel_id in isotopes:
                sym_combination += [
                    ch_obj.query_combination(symmetry, sites_per_channel[i])
                ]
                live_ch += [channel_id]
                live_ch_index += [i]
                live_n_sites += [sites_per_channel[i]]

        if sym_combination == [[]]:
            return np.asarray(sym_combination)

        linear_isotopes = [[ch] * ns for ch, ns in zip(live_ch, live_n_sites)]
        linear_iso_dict = get_iso_dict(channels, isotopes=np.hstack(linear_isotopes))

        symmetry_expanded = TransitionQuery.cartesian_product_indexing(sym_combination)

        if symmetry_expanded.size == 0:
            return symmetry_expanded

        all_combinations = np.zeros((symmetry_expanded.shape[0], len(isotopes)))

        # set missing channel isotope query to nan for non-P query
        value = 0 if symmetry == "P" else np.nan
        missing_ch = self._get_missing_channel_isotope(isotopes, channels)
        for ch in missing_ch:
            index = np.where(np.asarray(isotopes) == ch)
            all_combinations[:, index] = value

        _ = [
            all_combinations.__setitem__(
                (slice(None, None, None), iso_dict[channels[live_ch_index[i]]]),
                symmetry_expanded[:, linear_iso_dict[live_ch[i]]],
            )
            for i in range(len(sym_combination))
            if channels[live_ch_index[i]] in iso_dict
        ]

        return all_combinations


class Rotation(Parseable):
    """Base Rotation class.

    Attributes
    ----------

    angle:
        The rf rotation angle in units of radians.

    phase:
        The rf rotation phase in units of radians.
    """

    angle: float = Field(default=0.0, ge=0.0)  # in rads
    phase: float = Field(default=0.0)  # in rads

    property_unit_types: ClassVar[Dict] = {"angle": "angle", "phase": "angle"}
    property_default_units: ClassVar[Dict] = {"angle": "rad", "phase": "rad"}
    property_units: Dict = {"angle": "rad", "phase": "rad"}

    class Config:
        validate_assignment = True
        extra = "forbid"

    def json(self, units=True, **kwargs):
        return (
            {k: f"{getattr(self, k)} {u}" for k, u in self.property_units.items()}
            if units
            else {k: getattr(self, k) for k in self.property_units}
        )


# class MixingEnum(Enum):
#     """Enumeration for defining common mixing queries. The
# enumerations are as follows:

#     Attributes
#     ----------

#     TotalMixing:
#         Setting the query attribute to TotalMixing causes all transitions from the
#         previous spectral/delay event to connect to all transitions in the next
#         spectral/delay event. This is the default behavior when no mixing is defined
#         between events.

#     NoMixing:
#         Defines a query where transition mixing is not allowed between connecting
#         events.

#     Example
#     -------

#     The query attribute of the :py:class:`~mrsimulator.method.event.MixingEvent`
# can be
#     set to the Enum itself or a string representing the Enum.

#     >>> from mrsimulator.method import MixingEvenA
#     >>> from mrsimulator.method.query import MixingEnum
#     >>> # From Enum object
#     >>> total_mix = MixingEvent(ch1=MixingEnum.TotalMixing)
#     >>> no_mix = MixingEvent(ch1=MixingEnum.NoMixing)
#     >>> # From string representing Enum
#     >>> total_mix = MixingEvent(ch1="TotalMixing")
#     >>> no_mix = MixingEvent(ch1="NoMixing")
#     """

#     @classmethod
#     def allowed_enums(cls):
#         """Returns list of str corresponding to all valid enumerations"""
#         return [e.name for e in cls]

#     def json(self, **kwargs):
#         """Return a JSON-compliant serialization of enumeration"""
#         temp = self.value.json(**kwargs)
#         return self.value if isinstance(self.value, str) else temp

#     TotalMixing: str = "TotalMixing"
#     NoMixing: Rotation = Rotation(angle=0, phase=0)
