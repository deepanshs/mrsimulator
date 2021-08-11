# -*- coding: utf-8 -*-
from copy import deepcopy
from itertools import permutations
from typing import ClassVar
from typing import Dict
from typing import List
from typing import Optional

import numpy as np
from mrsimulator.transition import Transition
from mrsimulator.utils.parseable import Parseable
from pydantic import Field

from .utils import cartesian_product
from .utils import get_iso_dict

__author__ = "Deepansh J. Srivastava"
__email__ = "srivastava.89@osu.edu"

ON_FAIL_MESSAGE = (
    "The length of the transition query symmetry elements cannot exceed than the "
    "number of channels."
)


class SymmetryQuery(Parseable):
    """Base SymmetryQuery class.

    Attributes
    ----------

    P:
        A list of p symmetry functions per site. Here p = Δm = :math:`m_f - m_i` is the
        difference between the spin quantum numbers of the final and initial states.

        Example
        -------

        >>> method = Method2D(channels=['1H'])
        >>> method.spectral_dimensions[0].events[0].transition_query[0].ch1.P = [-1]

    D:
        A list of d symmetry functions per site. Here d = :math:`m_f^2 - m_i^2` is the
        difference between the square of the spin quantum numbers of the final and
        initial states.

        Example
        -------

        >>> method.spectral_dimensions[0].events[0].transition_query[0].ch1.D = [0]
    """

    P: List[int] = Field(
        default=[-1],
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

    def permutate_query(self, symmetry, n_site_at_channel_id):
        """Permutation of symmetry query based on the number of sites in given channel.

        Args:
            (str) symmetry: The symmetry element, 'P' or 'D'.
            (int) n_site_at_channel: Number of sites for the given channel.

        Example:
            Consider the following
                query = {P: [-1], D: [1]}
                n_isotopes = [3]
                channels = ['A']
        then,
        1. P query will expand and permutate to [-1, 0, 0], [0, -1, 0], and [0, 0, -1]
        2. D query will expand and permutate to [1, 0, 0], [0, 1, 0], and [0, 0, 1]
        """
        query = getattr(self, symmetry)
        if query is None:
            return []
        query_size = len(query)
        if query_size > n_site_at_channel_id:
            return []
            # raise ValueError(ON_FAIL_MESSAGE)
        query = query + (n_site_at_channel_id - query_size) * [0]
        return list(set(permutations(query)))


class TransitionQuery(Parseable):
    """TransitionQuery class for quering transition symmetry function.

    Attributes
    ----------

    ch1:
        An optional SymmetryQuery object for quering symmetry functions at channel
        index 0 of the method's channels array."

    ch2:
        An optional SymmetryQuery object for quering symmetry functions at channel
        index 1 of the method's channels array."

    ch3:
        An optional SymmetryQuery object for quering symmetry functions at channel
        index 2 of the method's channels array."

    Example
    -------

        >>> query = TransitionQuery(ch1={'P': [1], 'D': [0]}, ch2={'P': [-1]})
    """

    ch1: Optional[SymmetryQuery] = Field(
        title="ch1",
        default=SymmetryQuery(),
        description=(
            "An optional SymmetryQuery object for quering symmetry functions at "
            "channel index 0 of the method's channels array."
        ),
    )
    ch2: Optional[SymmetryQuery] = Field(
        title="ch2",
        default=None,
        description=(
            "An optional SymmetryQuery object for quering symmetry functions at "
            "channel index 1 of the method's channels array."
        ),
    )
    ch3: Optional[SymmetryQuery] = Field(
        title="ch3",
        default=None,
        description=(
            "An optional SymmetryQuery object for quering symmetry functions at "
            "channel index 2 of the method's channels array."
        ),
    )

    class Config:
        validate_assignment = True

    @staticmethod
    def cartesian_product_indexing(permutation):
        """Return Cartesian product of indexes"""
        permutation_length = [np.arange(len(p)) for p in permutation if len(p) != 0]
        if permutation_length == []:
            return np.asarray([])

        indexes = cartesian_product(*permutation_length)  # cartesian indexes
        return np.asarray(
            [np.hstack([permutation[i][j] for i, j in enumerate(ix)]) for ix in indexes]
        )

    def permutation(self, isotopes, channels):
        """Permutation of TransitionQuery based on the number of sites per channel.

        Args:
            (list) isotopes: List of isotope symbols, ['29Si , '13C', '13C', '1H'].
            (int) channels: List of method channels, ['29Si , '13C'].
        """

        iso_dict = get_iso_dict(channels=channels, isotopes=isotopes)
        sites_per_channel = [
            iso_dict[item].size if item in iso_dict else 0 for item in channels
        ]

        # expanded_symmetry = {}
        # for symmetry in ["P", "D"]:
        expanded_symmetry = {
            sym: self.expand_elements_for_symmetry(
                sym, isotopes, iso_dict, channels, sites_per_channel
            )
            for sym in ["P", "D"]
        }
        return expanded_symmetry

    def expand_elements_for_symmetry(
        self, symmetry, isotopes, iso_dict, channels, sites_per_channel
    ):
        P_permutated = []
        live_ch, live_ch_index = [], []
        live_n_sites = []
        for i, channel_id in enumerate(channels):
            ch_obj = getattr(self, f"ch{i+1}")
            if ch_obj is not None:
                P_permutated += [ch_obj.permutate_query(symmetry, sites_per_channel[i])]
                live_ch += [channel_id]
                live_ch_index += [i]
                live_n_sites += [sites_per_channel[i]]

        # P_permutated = [item for item in P_permutated if item != []]

        if P_permutated == [[]]:
            return np.asarray(P_permutated)

        linear_isotopes = [[ch] * ns for ch, ns in zip(live_ch, live_n_sites)]
        linear_iso_dict = get_iso_dict(channels, isotopes=np.hstack(linear_isotopes))

        symmetry_expanded = TransitionQuery.cartesian_product_indexing(P_permutated)

        if symmetry_expanded.size == 0:
            return symmetry_expanded

        P_expanded = np.zeros((symmetry_expanded.shape[0], len(isotopes)))
        _ = [
            P_expanded.__setitem__(
                (slice(None, None, None), iso_dict[channels[live_ch_index[i]]]),
                symmetry_expanded[:, linear_iso_dict[live_ch[i]]],
            )
            for i in range(len(P_permutated))
        ]

        return P_expanded


class RFRotation(Parseable):
    """Base RFRotation class.

    Attributes
    ----------

    tip_angle:
        The rf rotation angle in units of radians.

    phase:
        The rf rotation phase in units of radians.
    """

    tip_angle: float = Field(default=0.0, ge=0.0)  # in rads
    phase: float = Field(default=0.0)  # in rads

    property_unit_types: ClassVar[Dict] = {
        "tip_angle": "angle",
        "phase": "angle",
    }

    property_default_units: ClassVar[Dict] = {
        "tip_angle": "rad",
        "phase": "rad",
    }

    property_units: Dict = {
        "tip_angle": "rad",
        "phase": "rad",
    }

    class Config:
        validate_assignment = True


class MixingQuery(Parseable):
    """MixingQuery class for quering transition mixing between events.

    Attributes
    ----------

    ch1:
        An optional RFRotation object for channel at index 0 of the method's channels."

    ch2:
        An optional RFRotation object for channel at index 1 of the method's channels."

    ch3:
        An optional RFRotation object for channel at index 2 of the method's channels."

    Example
    -------

        >>> query = MixingQuery(ch1={"tip_angle": 1.570796, "phase": 3.141593})
    """

    ch1: RFRotation = None
    ch2: RFRotation = None
    ch3: RFRotation = None

    class Config:
        validate_assignment = True

    @classmethod
    def parse_dict_with_units(cls, py_dict):
        """
        Parse the physical quantity from a dictionary representation of the Method
        object, where the physical quantity is expressed as a string with a number and
        a unit.

        Args:
            dict py_dict: A python dict representation of the Method object.

        Returns:
            A :ref:`method_api` object.
        """
        py_dict_copy = deepcopy(py_dict)
        obj = {k: RFRotation.parse_dict_with_units(v) for k, v in py_dict_copy.items()}
        py_dict_copy.update(obj)
        return super().parse_dict_with_units(py_dict_copy)
