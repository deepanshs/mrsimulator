# -*- coding: utf-8 -*-
from itertools import permutations
from typing import List
from typing import Optional

import numpy as np
from mrsimulator.transition import Transition
from mrsimulator.utils.parseable import Base
from pydantic import Field

from .utils import cartesian_product
from .utils import get_iso_dict

__author__ = "Deepansh J. Srivastava"
__email__ = "srivastava.89@osu.edu"

ON_FAIL_MESSAGE = (
    "The length of the transition query symmetry elements cannot exceed than the "
    "number of channels."
)


class SymmetryQuery(Base):
    """Base SymmetryQuery class.

    Attributes
    ----------

    P:
        A list of p symmetry functions per site. Here p = Δm is the difference between
        spin quantum numbers of the final and initial states.

        Example
        -------

        >>> method = Method2D()
        >>> method.spectral_dimensions[0].events[0].transition_query[0].ch1.P = [-1]

    D:
        A list of d symmetry functions per site. Here p = Δm is the difference between
        spin quantum numbers of the final and initial states.

        Example
        -------

        >>> method.spectral_dimensions[0].events[0].transition_query[0].ch1.D = [0]
    """

    P: List[int] = Field(default=[-1])
    D: List[int] = Field(default=None)
    F: List[float] = Field(default=None)
    transitions: List[Transition] = None

    class Config:
        validate_assignment = True
        arbitrary_types_allowed = True

    def permutate_query(self, symmetry, n_site_at_channel_id):
        """Permutation of symmetry query based on the number of sites in given channel.

        (str) symmetry: The symmetry element, 'P' or 'D'.
        (int) n_site_at_channel: Number of sites for the given channel.

        Example. Consider the following
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


class TransitionQuery(Base):
    """Base TransitionQuery class.

    Attributes
    ----------

    ch1: SymmetryQuery
        An optional SymmetryQuery object for quering isotopes at channel index 0 of the
        method's channels array.


    ch2: SymmetryQuery
        An optional SymmetryQuery object for quering isotopes at channel index 1 of the
        method's channels array.

    Example
    -------

        >>> query = TransitionQuery(ch1={'P': [1], 'D': [0]}, ch2={'P': [-1]})
    """

    ch1: Optional[SymmetryQuery] = Field(default=SymmetryQuery())
    ch2: Optional[SymmetryQuery] = Field(default=None)
    ch3: Optional[SymmetryQuery] = Field(default=None)

    @staticmethod
    def cartesian_product_indexing(P_permutated):
        permutation_length = [
            np.arange(len(item)) for item in P_permutated if len(item) != 0
        ]
        if permutation_length == []:
            return np.asarray([])
        # print("perm length", permutation_length)
        cartesian_index = cartesian_product(*permutation_length)
        # print("cartesian_index", cartesian_index)
        # print(
        #     "parse",
        #     [
        #         [P_permutated[i][j] for i, j in enumerate(item)]
        #         for item in cartesian_index
        #     ],
        # )
        return np.asarray(
            [
                np.hstack([P_permutated[i][j] for i, j in enumerate(item)])
                for item in cartesian_index
            ]
        )

    def permutation(self, isotopes, channels):
        """Permutation of SymmetryQuery based on the number of sites per channel.

        Args:
            (list) isotopes: List of isotope symbols, ['29Si , '13C', '13C', '1H'].
            (int) channels: List of method channels, ['29Si , '13C'].

        # Example:
        #     >>> from mrsimulator.method.transition_query import SymmetryQuery
        #     >>> ss = SymmetryQuery(ch1=[-1], ch2=[1])
        #     >>> symmetry = ss.permutation(
        #     ...     isotopes=['1H', '13C', '1H', '13C'], channels=['1H', '13C']
        #     ... )
        #     >>> pprint(symmetry)
        #     [[-1.0, 1.0, 0.0, 0.0],
        #      [-1.0, 0.0, 0.0, 1.0],
        #      [0.0, 1.0, -1.0, 0.0],
        #      [0.0, 0.0, -1.0, 1.0]]
        """

        iso_dict = get_iso_dict(channels=channels, isotopes=isotopes)
        n_sites_per_channel = [
            iso_dict[item].size if item in iso_dict else 0 for item in channels
        ]

        expanded_symmetry = {}
        for symmetry in ["P", "D"]:
            expanded_symmetry[symmetry] = self.expand_elements_for_symmetry(
                symmetry, isotopes, iso_dict, channels, n_sites_per_channel
            )
        return expanded_symmetry

    def expand_elements_for_symmetry(
        self, symmetry, isotopes, iso_dict, channels, n_sites_per_channel
    ):
        P_permutated = []
        live_channel = []
        live_channel_index = []
        live_n_sites = []
        for i, channel_id in enumerate(channels):
            # if channel_id not in iso_dict:
            #     # warnings.warn(warn_message(channel_id))
            #     return np.asarray([])

            channel_obj = getattr(self, f"ch{i+1}")
            if channel_obj is not None:
                P_permutated += [
                    channel_obj.permutate_query(symmetry, n_sites_per_channel[i])
                ]
                live_channel += [channel_id]
                live_channel_index += [i]
                live_n_sites += [n_sites_per_channel[i]]

        # P_permutated = [item for item in P_permutated if item != []]

        # print("symmetry", symmetry)
        # print("perm", P_permutated)
        if P_permutated == [[]]:
            return np.asarray(P_permutated)

        linear_isotopes = np.hstack(
            [[item] * n_item for item, n_item in zip(live_channel, live_n_sites)]
        )
        linear_iso_dict = get_iso_dict(channels=channels, isotopes=linear_isotopes)

        symmetry_expanded = TransitionQuery.cartesian_product_indexing(P_permutated)

        if symmetry_expanded.size == 0:
            return symmetry_expanded

        # print("symm", symmetry_expanded)
        # print("l")
        P_expanded = np.zeros((symmetry_expanded.shape[0], len(isotopes)))
        for i in range(len(P_permutated)):
            # print("linear_iso_dict", linear_iso_dict[live_channel[i]])
            # print(
            #     "symmetry_expanded at iso",
            #     symmetry_expanded[:, linear_iso_dict[live_channel[i]]],
            # )
            # print("live_channel_index", live_channel_index, live_channel)
            P_expanded[
                :, iso_dict[channels[live_channel_index[i]]]
            ] = symmetry_expanded[:, linear_iso_dict[live_channel[i]]]

        # print("P_expanded", P_expanded)
        return P_expanded
