# -*- coding: utf-8 -*-
import warnings
from functools import reduce
from itertools import permutations

import numpy as np

__author__ = ["Deepansh J. Srivastava", "Maxwell C. Venetos"]
__email__ = ["srivastava.89@osu.edu", "maxvenetos@gmail.com"]


def cartesian_product(*arrays):
    la = len(arrays)
    dtype = np.result_type(*arrays)
    arr = np.empty([len(a) for a in arrays] + [la], dtype=dtype)
    for i, a in enumerate(np.ix_(*arrays)):
        arr[..., i] = a
    return arr.reshape(-1, la)


def P_symmetry_indexes(transitions, list_of_P):
    P = transitions[:, 1, :] - transitions[:, 0, :]
    n = P.shape[1]
    return reduce(
        np.union1d,
        [
            reduce(np.intersect1d, [np.where(P[:, i] == search[i]) for i in range(n)])
            for search in list_of_P
        ],
        np.asarray([], dtype=np.int64),
    )


def D_symmetry_indexes(transitions, list_of_D):
    D = transitions[:, 1, :] ** 2 - transitions[:, 0, :] ** 2
    return reduce(
        np.union1d,
        [
            reduce(
                np.intersect1d,
                [np.where(D[:, i] == search[i]) for i in range(D.shape[1])],
            )
            for search in list_of_D
        ],
        np.asarray([], dtype=np.int64),
    )


def get_iso_dict(channel, isotope):
    """
    Parse the spin system sites to determine indices of each isotope that is part of
    the method channel.

    Args:
        channel: List object
        isotope: List object
    """
    iso_dict = {}

    # determine channels for P
    for i, item in enumerate(isotope):
        if item in channel and item not in iso_dict:
            iso_dict[item] = [i]
        elif item in iso_dict:
            iso_dict[item].append(i)

    return iso_dict


def query_permutations(query, isotope, channel, transition_symmetry="P"):
    """
    Determines the transition symmetries that are involved in a given transition query.

    Args:
        query: Dict object
        channel: List object
        isotope: List object
        transition_symmetry: str object. Derived from a transition query
    """

    P_permutated = []
    iso_dict = get_iso_dict(channel=channel, isotope=isotope)

    # get the query symmetry element.
    query_short = query[transition_symmetry]

    for i, items in enumerate(query_short):
        # Check if method isotope is in the spin system
        if channel[i] not in iso_dict:
            warnings.warn(
                f"Method/channel isotope mismatch. Channel asks for {channel[i]} "
                f"but is not in {isotope}"
            )
            return []

        temp_P = []
        iso_ch_length = len(iso_dict[channel[i]])
        for k in range(len(query_short[items])):
            query_item_len = len(query_short[items][k])
            # Check transition query doesn't require more isotopes than present
            if query_item_len > iso_ch_length:
                raise AttributeError(
                    "Failed: The length of the transition query symmetry elements "
                    "cannot larger than the number of channel"
                )

            query_short[items][k] += (iso_ch_length - query_item_len) * [k]
            temp_P += list(permutations(query_short[items][k]))
        P_permutated += [temp_P]

        # P_permutated += [list(permutations(query_short[items][k]))]

    previous_sets = []
    for i, iso_trans_symmetry in enumerate(P_permutated):
        # creating transition symmetries isotope by isotope
        temp_transitions = []
        iso_ch_i = iso_dict[channel[i]]
        for transition in iso_trans_symmetry:
            P_expanded = np.zeros(len(isotope))

            # fill indices of spin system with the sites transition symmetries
            P_expanded[iso_ch_i] = transition

            if previous_sets == []:
                temp_transitions += [P_expanded]
            else:
                # Each isotope is added to the previous isotope to create the
                # full transition symmetry
                temp_transitions += np.asarray(previous_sets) + P_expanded

        previous_sets = temp_transitions

    return temp_transitions
