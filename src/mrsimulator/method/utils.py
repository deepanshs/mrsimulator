# -*- coding: utf-8 -*-
from functools import reduce
from itertools import permutations

import numpy as np

__author__ = ["Deepansh J. Srivastava", "Maxwell C. Venetos"]
__email__ = ["srivastava.89@osu.edu", "maxvenetos@gmail.com"]


def expand_spectral_dimension_object(py_dict):
    glb = {}
    list_g = ["magnetic_flux_density", "rotor_frequency", "rotor_angle"]
    for item in list_g:
        if item in py_dict.keys():
            glb[item] = py_dict[item]
    glb_keys = glb.keys()

    for dim in py_dict["spectral_dimensions"]:
        if "events" not in dim:
            dim["events"] = [{}]
        for ev in dim["events"]:
            intersect = set(ev.keys()).intersection(set(glb_keys))
            for k in glb:
                if k not in intersect:
                    ev[k] = glb[k]

    _ = [py_dict.pop(item) for item in glb_keys]

    return py_dict


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
        query: Transition query dict object
        channel: List of channels as atomic number followed by symbol. Eg. '29Si', '1H'
        isotope: List of isotopes within the spin system. Eg. ['29Si', '29Si', '1H']
        transition_symmetry: str object. Derived from a transition query
    """
    P_permutated = []

    iso_dict = get_iso_dict(channel=channel, isotope=isotope)
    query_symmetry = query[transition_symmetry]  # get the query for symmetry element.

    # def warn_message(id_):
    #     return (
    #         f"Channel asks for isotope `{id_}` but it is not present in the list of "
    #         "spin systems."
    #     )

    on_fail_message = (
        "The length of the transition query symmetry elements cannot exceed than the "
        "number of channels."
    )
    for i, channel_id in enumerate(channel):
        # Check if method's channel isotope is present in the spin system
        if channel_id not in iso_dict:
            # warnings.warn(warn_message(channel_id))
            return np.asarray([])

        n_sites_channel_i = iso_dict[channel_id].size
        channel_query = query_symmetry[f"channel-{i+1}"]

        temp_P = []
        for item in channel_query:
            query_item_len = len(item)
            # Check transition query exceed the number of isotopes present
            if query_item_len > n_sites_channel_i:
                raise ValueError(on_fail_message)

            item += (n_sites_channel_i - query_item_len) * [0]
            temp_P += list(set(permutations(item)))
        P_permutated += [temp_P]

    # Expand the permutation to the number of sites in the spin system
    permutation_length = max(len(item) for item in P_permutated)
    P_expanded = np.zeros((permutation_length, len(isotope)))
    for i, iso_trans_symmetry in enumerate(P_permutated):
        # Update the channel-i indexes with the permuted symmetry.
        P_expanded[:, iso_dict[channel[i]]] = iso_trans_symmetry

    return P_expanded
