from functools import reduce

import numpy as np
from mrsimulator.utils.error import MissingSpectralDimensionError
from mrsimulator.utils.error import MixedSpectralDimensionTypeError
from scipy.spatial.transform import Rotation

__author__ = ["Deepansh J. Srivastava", "Maxwell C. Venetos", "Matthew D. Giammar"]
__email__ = ["srivastava.89@osu.edu", "maxvenetos@gmail.com", "giammar.7@ous.edu"]


TWO_PI = np.pi * 2


def wrap_between_pi(a: float):
    """Wraps the provided angle between (-pi and pi]"""
    a %= TWO_PI
    a -= np.sign(a) * TWO_PI if abs(a) > np.pi else 0
    return a


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
            reduce(
                np.intersect1d,
                [
                    np.where(fn[:, i] == search[i])
                    for i in range(n)
                    if not np.isnan(search[i])
                ],
            )
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
    # NOTE: returning index of 0 when non-existent might cause issues when first ev mix?
    return [
        i - 1 - (min(lower) if lower != [] else 0),
        i + (min(upper) if upper != [] else 0),
    ]


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
    query = spectral_dimensions[sp].events[index].query

    # Return the query, if is a MixingQuery, otherwise the value of the MixingEnum
    return query if query.__class__.__name__ == "MixingQuery" else query.value


def map_mix_query_attr_to_ch(mixing_query):
    """Map the mixing query attributes (angle and phase) to the channel index.
    If the attribute is defined for the channel use the defined value else set it to 0.

    Args:
        mixing_query: The mixing query to map
    """
    attributes = ["angle", "phase"]
    return {
        i: {
            item: getattr(getattr(mixing_query, f"ch{i+1}"), item) or 0
            if getattr(mixing_query, f"ch{i+1}") is not None
            else 0
            for item in attributes
        }
        for i in range(3)
    }


# def angle_and_phase_list(symbol, channels, mixing_query):
#     """Return a list of angles and phase of size equal to the number of sites within
#     the spin system, corresponding to a mixing_query from a MixingEvent.

#     If the site matches the channel, append the angle and phase of the corresponding
#     channel to the list, otherwise append 0.

#     Args:
#         symbols: List of site symbols.
#         channels: List of method channel symbols.
#         mixing_query: Mixing query object of the MixingEvent.
#     """
#     angle_phase_mappable = map_mix_query_attr_to_ch(mixing_query)
#     angle_ = [
#         angle_mappable["angle"][channels.index(sym)] if sym in channels else 0
#         for sym in symbol
#     ]
#     phase_ = [
#         angle_mappable["phase"][channels.index(sym)] if sym in channels else 0
#         for sym in symbol
#     ]
#     return angle_, phase_


def to_euler_list(symbol, channels, mixing_queries):
    """Takes a list of symbols, list of isotope symbols per channel, and list of mixing
    queries and converts them into a list of sets of three Euler angles describing the
    total rotation of the combined mixing_queries per channel.

    Args:
        symbols: List of site symbols.
        channels: List of method channel symbols.
        mixing_query: List of mixing query objects from sequential MixingEvents.
    """
    angle_phase_mappable = [map_mix_query_attr_to_ch(query) for query in mixing_queries]

    # angles is a list of sets of Euler angles for each symbol (site)
    # ex. [(3.14, 1.57, -3.14), (-1.57, 2, 1.57)]
    angles = [
        combine_mixing_queries([ap[channels.index(sym)] for ap in angle_phase_mappable])
        if sym in channels
        else [np.pi / 2, 0, -np.pi / 2]
        for sym in symbol
    ]

    return angles


def get_grouped_mixing_queries(spec_dims, event_names):
    """Returns a dictionary where each key is the index of the first MixingEvent in a
    group of sequential MixingEvents and the key is the set of angles and phases for
    those mixing queries in the mixing events.

    Args:
        spec_dims: A list of SpectralDimension objects.
        event_names: A list of all class names
    """
    # dict with index of first mixing event in seq as key and list of queries as values
    mixing_query_sets = {}
    previous_event_mix = False
    for i, name in enumerate(event_names):
        if name != "MixingEvent":
            previous_event_mix = False

        # Skip this event if previous event mixing or if this event TotalMixing
        elif previous_event_mix or get_mixing_query(spec_dims, i) == "TotalMixing":
            continue

        # Add this MixingEvent query and the queries from the next contiguous
        # MixingEvent to a list. Only run for the first MixingEvent in a sequence
        else:
            previous_event_mix = True
            mixing_query_sets[i] = []
            j = i
            while j < len(event_names) and event_names[j] == "MixingEvent":
                # Only add query if the query is not the string TotalMixing
                query = get_mixing_query(spec_dims, j)
                mixing_query_sets[i] += [query] if query != "TotalMixing" else []
                j += 1

    return mixing_query_sets


def mixing_query_connect_map(spectral_dimensions):
    """Return a list of mappables corresponding to each mixing event. The mappable
    corresponds to queries described by adjacent mixing events and the index of next
    and previous nearest transition indexes.

    Args:
        spectral_dimensions: A list of SpectralDimension objects.
    """
    event_names = [
        evt.__class__.__name__ for dim in spectral_dimensions for evt in dim.events
    ]
    grouped_mix_map = get_grouped_mixing_queries(spectral_dimensions, event_names)
    non_mix_index = [i for i, ev in enumerate(event_names) if ev != "MixingEvent"]
    non_mix_index_map = {index: i for i, index in enumerate(non_mix_index)}

    return [
        {
            "mixing_query_list": query_list,
            "near_index": [
                non_mix_index_map[k] for k in nearest_nonmixing_event(event_names, i)
            ],
        }
        for i, query_list in grouped_mix_map.items()
    ]


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
        raise MissingSpectralDimensionError()

    # Currently only support 2 SpectralDimensions in C code, limit to 2 for now
    num_spec_dims = len(py_dict["spectral_dimensions"])
    if num_spec_dims > 2:
        raise NotImplementedError(
            "Mrsimulator currently supports a maximum of two spectral dimensions. "
            f"Found {num_spec_dims} spectral dimensions. "
        )


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


def combine_mixing_queries(queries: list):
    """Takes in a list of mixing queries combining them into a single mixing query

    Args:
        queries: List of dicts each representing a MixingQuery object

    Returns:
        Dictionary with angle and phase of combined MixingQuery objects
    """
    if len(queries) == 0:
        raise ValueError(f"List length must be at least 1. Got length {len(queries)}.")

    # NOTE: Do not need to check for length 1 since queries[1:] == [] when length 1
    alpha, beta, gamma = _angle_phase_to_euler_angles(**queries[0])
    for query in queries[1:]:
        alpha, beta, gamma = _add_two_euler_angles(
            alpha, beta, gamma, *_angle_phase_to_euler_angles(**query)
        )

    return alpha, beta, gamma


def _angle_phase_to_euler_angles(angle: float, phase: float):
    """Takes angle and phase of a mixing query and converts to a set of euler angles in
    the ZYZ convention. The returned angles will be constrained between (-pi, pi]

    Args:
        angle (float): Angle of mixing query between [0, 2pi)
        phase (float): Phase of mixing query between [0, 2pi)

    Returns:
        alpha, beta, gamma: Euler angles of the mixing query
    """
    # Wrap angle and phase between -pi and pi
    angle, phase = wrap_between_pi(angle), wrap_between_pi(phase)
    alpha = (np.pi / 2) - phase
    return wrap_between_pi(alpha), wrap_between_pi(angle), wrap_between_pi(-alpha)


def _euler_angles_to_angle_phase(alpha: float, beta: float, gamma: float):
    """Takes a set of euler angles in the ZYZ convention and converts them to a
    mixing angle and phase. Provided alpha and gamma should be opposite of each other,
    otherwise a ValueError is raised since the rotation vector does not lie in the XY
    plane.

    Args:
        alpha, beta, gamma: Set of euler angles to convert

    Returns:
        angle, phase: Angle and phase of the equivalent mixing query

    Raises:
        ValueError: Raised if alpha and gamma are not opposite of each other
    """
    if not np.isclose(alpha, -gamma):
        raise ValueError(
            "Unable to convert the provided Euler angles to an angle and phase"
        )

    phase = wrap_between_pi(gamma + np.pi / 2)

    return beta, phase


def _add_two_euler_angles(a1, b1, g1, a2, b2, g2):
    """Adds two sets of euler angles -- (a1, b1, g1) and (a2, b2, g2) -- together.
    Also checks for edge cases where gimbal lock would occur.

    If the result is the identity matrix, then beta = 0 and alpha, gamma are unbounded.
    As an arbitrary choice, alpha of pi/2 and gamma of -pi/2 are chosen.

    If the two phases are the same (i.e. a1 == a2 and g1 == g2) and b1 + b2 = pi, then
    scipy is unable to uniquely determine the last angle. In this case, the top left
    element of the rotation matrix equals -cos(2*phase) and alpha and gamma are computed
    from there.

    Returns:
        alpha, beta, gamma: The resulting Euler angles
    """
    rot_1 = Rotation.from_euler("zyz", [a1, b1, g1])
    rot_2 = Rotation.from_euler("zyz", [a2, b2, g2])

    result = rot_1 * rot_2
    result_mat = result.as_matrix()

    # Check for identity matrix by comparing the diagonal to an array of ones
    if np.allclose(result_mat.diagonal(), [1.0, 1.0, 1.0]):
        return np.asarray([np.pi / 2, 0, -np.pi / 2])

    # Check if beta is 180 degrees and phase the same
    if np.isclose(result_mat[2][2], -1.0) and np.allclose([a1, g1], [a2, g2]):
        gamma = np.arccos(result_mat[0][0]) - np.pi
        gamma /= 2
        return np.asarray([-gamma, np.pi, gamma])

    return result.as_euler("zyz")


# def add_euler_angles(angles: list):
#     """Adds an arbitrary number of sequential euler rotations into one rotation

#     Args:
#         angles (list): List of euler angles provided in the form of
#             [(a_1, b_1, g_1), (a_2, b_2, g_2), ...(a_n, b_n, g_n)]

#     Returns:
#         Ordered list of floats (alpha, beta, gamma) representing the summed
#         euler rotations
#     """
#     len_error = "Elements in the angles list must be in the form (alpha, beta, gamma)"
#     # Check if arguments less than two elements
#     if len(angles) < 2:
#         raise ValueError("Please provide at least two sets of angles")

#     if len(list[0]) != 3:
#         raise ValueError(len_error)

#     alpha, beta, gamma = list[0][0], list[0][1], list[0][2]
#     for euler_angles in iter(list[1:]):
#         if len(euler_angles[0]) != 3:
#             raise ValueError(len_error)

#         alpha, beta, gamma = _add_two_euler_angles(alpha, beta, gamma, *euler_angles)

#     return alpha, beta, gamma
