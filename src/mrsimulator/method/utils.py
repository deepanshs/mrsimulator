# -*- coding: utf-8 -*-
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


def get_symmetry_indexes(fn, list_of_sym):
    n = fn.shape[1]
    return reduce(
        np.union1d,
        [
            reduce(np.intersect1d, [np.where(fn[:, i] == search[i]) for i in range(n)])
            for search in list_of_sym
        ],
        np.asarray([], dtype=np.int64),
    )


def P_symmetry_indexes(transitions, list_of_P):
    P_fn = transitions[:, 1, :] - transitions[:, 0, :]
    return get_symmetry_indexes(P_fn, list_of_P)


def D_symmetry_indexes(transitions, list_of_D):
    D_fn = transitions[:, 1, :] ** 2 - transitions[:, 0, :] ** 2
    return get_symmetry_indexes(D_fn, list_of_D)


def get_iso_dict(channel, isotope):
    """
    Parse the spin system sites to determine indices of each isotope that is part of
    the method channel.

    Args:
        channel: List object
        isotope: List object
    """
    intersection = set(isotope).intersection(set(channel))
    isotope = np.asarray(isotope)
    return {item: (np.where(isotope == item))[0] for item in intersection}


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

    # def warn_message(id_):
    #     return (
    #         f"Channel asks for isotope `{id_}` but it is not present in the list of "
    #         "spin systems."
    #     )

    on_fail_message = (
        "The length of the transition query symmetry elements cannot exceed than the "
        "number of channels."
    )
    for i, items in enumerate(query_short):
        # Check if method isotope is in the spin system
        if channel[i] not in iso_dict:
            # warnings.warn(warn_message(channel[i]))
            return []

        temp_P = []
        iso_ch_length = len(iso_dict[channel[i]])
        for k in range(len(query_short[items])):
            query_item_len = len(query_short[items][k])
            # Check transition query doesn't require more isotopes than present
            if query_item_len > iso_ch_length:
                raise ValueError(on_fail_message)

            query_short[items][k] += (iso_ch_length - query_item_len) * [k]
            temp_P += list(set(permutations(query_short[items][k])))
        P_permutated += [temp_P]

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
                temp_transitions += list(np.asarray(previous_sets) + P_expanded)

        previous_sets = temp_transitions

    return np.asarray(temp_transitions)
