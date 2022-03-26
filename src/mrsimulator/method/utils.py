# -*- coding: utf-8 -*-
from functools import reduce

import numpy as np
from mrsimulator.utils.error import MixedSpectralDimensionTypeError

__author__ = ["Deepansh J. Srivastava", "Maxwell C. Venetos", "Matthew D. Giammar"]
__email__ = ["srivastava.89@osu.edu", "maxvenetos@gmail.com", "giammar.7@ous.edu"]


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


def get_iso_dict(channels, isotopes):
    """
    Parse the spin system sites to determine indices of each isotope that is part of
    the method channel.

    Args:
        channels: List of method channels
        isotopes: List of spin system isotopes.
    """
    intersection = set(isotopes).intersection(set(channels))
    isotopes = np.asarray(isotopes)
    return {item: (np.where(isotopes == item))[0] for item in intersection}


def nearest_nonmixing_event(event_name, i):
    """Return the indexes of the nearest non mixing events (SpectralEvent and
    ConstantDurationEvent) about a mixing event at index `i`.

    Args:
        event_name: List of event class names.
        i: Int index of the mixing event.
    """
    options = ["SpectralEvent", "ConstantDurationEvent"]
    low_range = event_name[:i]
    high_range = event_name[i:]
    upper = [high_range.index(item) for item in options if item in high_range]
    lower = [low_range[::-1].index(item) for item in options if item in low_range]
    return [
        i - 1 - (min(lower) if lower != [] else 0),
        i + (min(upper) if upper != [] else 0),
    ]


def tip_angle_and_phase_list(symbol, channels, mixing_query):
    """Return a list of tip_angles and phase of size equal to the number of sites within
    the spin system, corresponding to a mixing_query from a MixingEvent.

    If the site matches the channel, append the tip_angle and phase of the corresponding
    channel to the list, otherwise append 0.

    Args:
        symbols: List of site symbols.
        channels: List of method channel symbols.
        mixing_query: Mixing query object of the MixingEvent.
    """
    angle_mappable = map_mix_query_attr_to_ch(mixing_query)
    tip_angle_ = [
        angle_mappable["tip_angle"][channels.index(sym)] if sym in channels else 0
        for sym in symbol
    ]
    phase_ = [
        angle_mappable["phase"][channels.index(sym)] if sym in channels else 0
        for sym in symbol
    ]
    return tip_angle_, phase_


def get_mixing_query(spectral_dimensions, index):
    """Return the mixing query object corresponding to the event at index `index`. The
    indexing is over flattened list of events from all spectral dimensions.

    Args:
        spectral_dimension: A list SpectralDimension objects.
        index: The index of the event from a flatten event list.
    """
    n_events = len(spectral_dimensions[0].events)
    sp = 0
    while index >= n_events:
        index -= n_events
        sp += 1
        n_events = len(spectral_dimensions[sp].events)
    return spectral_dimensions[sp].events[index].mixing_query


def map_mix_query_attr_to_ch(mixing_query):
    """Map the mixing query attributes (tip_angle and phase) to the channel index.
    If the attribute is defined for the channel use the defined value else set it to 0.

    Args:
        spectral_dimension: A list SpectralDimension objects.
        index: The index of the event from a flatten event list.
    """
    attributes = ["tip_angle", "phase"]
    return {
        item: {
            i: getattr(getattr(mixing_query, f"ch{i+1}"), item) or 0
            if getattr(mixing_query, f"ch{i+1}") is not None
            else 0
            for i in range(3)
        }
        for item in attributes
    }


def mixing_query_connect_map(spectral_dimensions):
    """Return a list of mappables corresponding to each mixing event. The mappable
    corresponds to mixing event and the index of next nearest transition indexes.

    Args:
        spectral_dimensions: A list of SpectralDimension objects."""
    mapping = {}
    event_names = [
        evt.__class__.__name__ for dim in spectral_dimensions for evt in dim.events
    ]
    non_mix_index = [i for i, ev in enumerate(event_names) if ev != "MixingEvent"]
    non_mix_index_map = {index: i for i, index in enumerate(non_mix_index)}

    mapping = [
        {
            "mixing_query": get_mixing_query(spectral_dimensions, i),
            "near_index": [
                non_mix_index_map[k] for k in nearest_nonmixing_event(event_names, i)
            ],
        }
        for i, name in enumerate(event_names)
        if name == "MixingEvent"
    ]

    return mapping


# Helper functions for validating a method object
def check_for_number_of_spectral_dimensions(py_dict, is_named_method=False, n=None):
    """Check number of spectral dimensions passed if method is named method or adds a
    single default spectral dimension if no spectral dimension is present and not a
    named method

    Args:
        (dict) py_dict: dict representation of method object under validation
        (bool) is_named_method: True if from methods library, False otherwise
        (int) n: Number of dimensions for named method

    Raises:
        ValueError if number of passed dimensions does not match required number
    """
    # Named method object
    if is_named_method:
        if "spectral_dimensions" not in py_dict:
            py_dict["spectral_dimensions"] = [{} for _ in range(n)]
        else:
            m = len(py_dict["spectral_dimensions"])
            if m != n:
                raise ValueError(
                    f"Method requires exactly {n} spectral dimensions, given {m}."
                )
        return

    # Generic method object
    if "spectral_dimensions" not in py_dict or py_dict["spectral_dimensions"] == []:
        py_dict["spectral_dimensions"] = [{}]


def check_spectral_dimensions_are_dict(py_dict):
    """Check if type of passed spectral dimension is dict

    Returns:
        True if all dict
        False if all SpectralDimension

    Raises:
        MixedSpectralDimensionTypeError if mixed list provided
    """
    sd_is_dict = [isinstance(sd, dict) for sd in py_dict["spectral_dimensions"]]
    if all(sd_is_dict):  # All items dict objects
        return True
    elif not any(sd_is_dict):  # All items in list are SpectralDimension objects
        return False
    else:  # Mixture of dict and obj provided, raise error
        raise MixedSpectralDimensionTypeError()


def check_for_at_least_one_event(py_dict):
    """Update events to [{}] if not present."""
    _ = [
        item.update({"events": [{}]})
        for item in py_dict["spectral_dimensions"]
        if "events" not in item
    ]


# Deprecated
# def query_permutations(query, isotope, channel, transition_symmetry="P"):
#     """
#     Determines the transition symmetries that are involved in a given transition
#     query.

#     Args:
#         query: Transition query dict object
#         channel: List of channels as atomic number followed by symbol. Eg. '29Si',
#         isotope: List of isotopes within the spin system. Eg. ['29Si', '29Si', '1H']
#         transition_symmetry: str object. Derived from a transition query
#     """
#     P_permutated = []

#     iso_dict = get_iso_dict(channels=channel, isotopes=isotope)
#     query_symmetry = query[transition_symmetry]  # get the query for symmetry element.

#     # def warn_message(id_):
#     #     return (
#     #         f"Channel asks for isotope `{id_}` but it is not present in the list "
#     #         "of spin systems."
#     #     )

#     on_fail_message = (
#         "The length of the transition query symmetry elements cannot exceed than "
#         "the number of channels."
#     )
#     for i, channel_id in enumerate(channel):
#         # Check if method's channel isotope is present in the spin system
#         if channel_id not in iso_dict:
#             # warnings.warn(warn_message(channel_id))
#             return np.asarray([])

#         n_sites_channel_i = iso_dict[channel_id].size
#         channel_query = query_symmetry[f"channel-{i+1}"]

#         temp_P = []
#         for item in channel_query:
#             query_item_len = len(item)
#             # Check transition query exceed the number of isotopes present
#             if query_item_len > n_sites_channel_i:
#                 raise ValueError(on_fail_message)

#             item += (n_sites_channel_i - query_item_len) * [0]
#             temp_P += list(set(permutations(item)))
#         P_permutated += [temp_P]

#     # Expand the permutation to the number of sites in the spin system
#     permutation_length = max(len(item) for item in P_permutated)
#     P_expanded = np.zeros((permutation_length, len(isotope)))
#     for i, iso_trans_symmetry in enumerate(P_permutated):
#         # Update the channel-i indexes with the permuted symmetry.
#         P_expanded[:, iso_dict[channel[i]]] = iso_trans_symmetry

#     return P_expanded
