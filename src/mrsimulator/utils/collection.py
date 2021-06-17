# -*- coding: utf-8 -*-
from typing import List

import numpy as np
from mrsimulator import Site
from mrsimulator import SpinSystem

__author__ = ["Deepansh Srivastava", "Matthew D. Giammar"]
__email__ = ["srivastava.89@osu.edu", "giammar.7@buckeyemail.osu.edu"]

SHIELDING_SYM_PARAMS = ["zeta", "eta", "alpha", "beta", "gamma"]
QUADRUPOLAR_PARAMS = ["Cq", "eta", "alpha", "beta", "gamma"]
LIST_LEN_ERROR_MSG = (
    "All arguments passed must be the same size. If one attribute is type list of "
    "length n, then all passed lists must be of length n and all other "
    "attributes must be scalar (singular float, int or str)."
)


def single_site_system_generator(
    isotope,
    isotropic_chemical_shift=0,
    shielding_symmetric=None,
    shielding_antisymmetric=None,
    quadrupolar=None,
    abundance=None,
    site_name=None,
    site_label=None,
    site_description=None,
    rtol=1e-3,
):
    r"""Generate and return a list of single-site spin systems from the input parameters

    Args:

        (list) isotope:
            A required string or a list of site isotope.
        (list) isotropic_chemical_shift:
            A float or a list/ndarray of values. The default value is 0.
        (list) shielding_symmetric:
            A shielding symmetric dict-like object, where the keyword value can either
            be a float or a list/ndarray of floats. The default value is None. The
            allowed keywords are ``zeta``, ``eta``, ``alpha``, ``beta``, and ``gamma``.
        (dict) shielding_antisymmetric:
            A shielding symmetric dict-like object, where the keyword value can either
            be a float or a list/ndarray of floats. The default value is None. The
            allowed keywords are ``zeta``, ``alpha``, and ``beta``.
        (dict) quadrupolar:
            A quadrupolar dict-like object, where the keyword value can either be a
            float or a list/ndarray of floats. The default value is None. The allowed
            keywords are ``Cq``, ``eta``, ``alpha``, ``beta``, and ``gamma``.
        (list) abundance:
            A float or a list/ndarray of floats describing the abundance of each spin
            system.
        (list) site_name:
            A string or a list of strings each with a site name. The default is None.
        (list) site_label:
            A string or a list of strings each with a site label. The default is None.
        (list) site_description:
            A string or a list of strings each with a site description. Default is None.
        (float) rtol:
            The relative tolerance. This value is used in determining the cutoff
            abundance given as
            :math:`\tt{abundance}_{\tt{cutoff}} = \tt{rtol} * \tt{max(abundance)}.`
            The spin systems with abundance below this threshold are ignored.

    Returns:
        list of :ref:`spin_sys_api` objects with a single :ref:`site_api`

    Example:
        >>> # single SpinSystem
        >>> sys1 = single_site_system_generator(
        ...     isotope=["1H"],
        ...     isotropic_chemical_shift=10,
        ...     site_name="Single Proton",
        ... )
        >>> print(len(sys1))
        1

        >>> # multiple SpinSystems
        >>> sys2 = single_site_system_generator(
        ...     isotope="1H",
        ...     isotropic_chemical_shift=[10] * 5,
        ...     site_name="5 Protons",
        ... )
        >>> print(len(sys2))
        5

        >>> # multiple SpinSystems with dict arguments
        >>> Cq = [4.2e6] * 12
        >>> sys3 = single_site_system_generator(
        ...     isotope="17O",
        ...     isotropic_chemical_shift=60.0,  # in ppm,
        ...     quadrupolar={"Cq": Cq, "eta": 0.5},  # Cq in Hz
        ... )
        >>> print(len(sys3))
        12

    .. note::
        The parameter value can either be a float or a list/ndarray. If the parameter
        value is a float, the given value is assigned to the respective parameter in all
        the spin systems. If the parameter value is a list or ndarray, its ith value is
        assigned to the respective parameter of the ith spin system. When multiple
        parameter values are given as lists/ndarrays, the length of all the lists must
        be the same.
    """
    sites = generate_site_list(
        isotope=isotope,
        isotropic_chemical_shift=isotropic_chemical_shift,
        shielding_symmetric=shielding_symmetric,
        shielding_antisymmetric=shielding_antisymmetric,
        quadrupolar=quadrupolar,
        site_name=site_name,
        site_label=site_label,
        site_description=site_description,
    )
    n_sites = len(sites)

    if abundance is None:
        abundance = 1 / n_sites
    abundance = _extend_to_nparray(_fix_item(abundance), n_sites)
    n_abd = abundance.size

    if n_sites == 1:
        sites = np.asarray([sites[0] for _ in range(n_abd)])
        n_sites = sites.size

    if n_sites != n_abd:
        raise ValueError(
            "Number of sites does not mach number of abundances. " + LIST_LEN_ERROR_MSG
        )

    keep_idxs = np.where(abundance > rtol * abundance.max())[0]

    return [
        SpinSystem(sites=[site], abundance=abd)
        for site, abd in zip(sites[keep_idxs], abundance[keep_idxs])
    ]


def generate_site_list(
    isotope,
    isotropic_chemical_shift=0,
    shielding_symmetric=None,
    shielding_antisymmetric=None,
    quadrupolar=None,
    site_name=None,
    site_label=None,
    site_description=None,
) -> List[Site]:
    r"""Takes in lists or list-like objects describing attributes of each site and
    returns a list of Site objects

    Params:
        (list) isotope:
            A string or a list of site isotope.
        (list) isotropic_chemical_shift:
            A float or a list/ndarray of values. The default value is 0.
        (dict) shielding_symmetric:
            A shielding symmetric dict-like object, where the keyword value can either
            be a float or a list/ndarray of floats. The default value is None. The
            allowed keywords are ``zeta``, ``eta``, ``alpha``, ``beta``, and ``gamma``.
        (dict) shielding_antisymmetric:
            A shielding symmetric dict-like object, where the keyword value can either
            be a float or a list/ndarray of floats. The default value is None. The
            allowed keywords are ``zeta``, ``alpha``, and ``beta``.
        (dict) quadrupolar:
            A quadrupolar dict-like object, where the keyword value can either be a
            float or a list/ndarray of floats. The default value is None. The allowed
            keywords are ``Cq``, ``eta``, ``alpha``, ``beta``, and ``gamma``.
        (list) site_name:
            A string or a list of strings each with a site name. The default is None.
        (list) site_label:
            A string or a list of strings each with a site label. The default is None.
        (list) site_description:
            A string or a list of strings each with a site description. Default is None.

    Returns:
        (list) sites: List of :ref:`site_api` objects

    Example:
        >>> # 10 hydrogen sites
        >>> sites1 = generate_site_list(
        ... isotope=["1H"] * 10,
        ... isotropic_chemical_shift=-15,
        ... site_name="10 Protons",
        ... )
        >>> print(len(sites1))
        10

        >>> # 10 hydrogen sites with different shifts
        >>> shifts = np.arange(-10, 10, 2)
        >>> sites2 = generate_site_list(
        ... isotope=["1H"] * 10,
        ... isotropic_chemical_shift=shifts,
        ... site_name="10 Proton",
        ... )
        >>> print(len(sites2))
        10

        >>> # multiple Sites with dict arguments
        >>> Cq = [4.2e6] * 12
        >>> sys3 = generate_site_list(
        ... isotope="17O",
        ... isotropic_chemical_shift=60.0,  # in ppm,
        ... quadrupolar={"Cq": Cq, "eta": 0.5},  # Cq in Hz
        ... )
        >>> print(len(sys3))
        12
    """
    attributes = [
        _fix_item(isotope),
        _fix_item(isotropic_chemical_shift),
        _fix_item(site_name),
        _fix_item(site_label),
        _fix_item(site_description),
    ]

    n_sites = _check_lengths(attributes)

    if shielding_symmetric is not None:
        shld_sym, n_dict = _extend_dict_values(shielding_symmetric, n_sites)
        n_sites = max(n_sites, n_dict)
        attributes.append(shld_sym)
    else:
        attributes.append(None)

    if shielding_antisymmetric is not None:
        shld_antisym, n_dict = _extend_dict_values(shielding_antisymmetric, n_sites)
        n_sites = max(n_sites, n_dict)
        attributes.append(shld_antisym)
    else:
        attributes.append(None)

    if quadrupolar is not None:
        quad, n_dict = _extend_dict_values(quadrupolar, n_sites)
        n_sites = max(n_sites, n_dict)
        attributes.append(quad)
    else:
        attributes.append(None)

    # Attributes order is same as below in list comprehension
    attributes = [_extend_to_nparray(attr, n_sites) for attr in attributes]

    return np.asarray(
        [
            Site(
                isotope=iso,
                isotropic_chemical_shift=shift,
                name=name,
                label=label,
                description=desc,
                shielding_symmetric=symm,
                shielding_antisymmetric=antisymm,
                quadrupolar=quad,
            )
            for iso, shift, name, label, desc, symm, antisymm, quad in zip(*attributes)
        ]
    )


def _fix_item(item):
    """Flattens multidimensional arrays into 1d array"""
    if isinstance(item, (list, np.ndarray)):
        return np.asarray(item).ravel()
    return item


def _extend_to_nparray(item, n):
    """If item is already list/array return np.array, otherwise extend to length n"""
    if isinstance(item, (list, np.ndarray)):
        return np.asarray(item)
    return np.asarray([item for _ in range(n)])


def _extend_dict_values(_dict, n_sites):
    """Checks and extends dict values. Returns dict or list of dicts and max length"""
    _dict = {key: _fix_item(val) for key, val in _dict.items()}
    n_sites_dict = _check_lengths(list(_dict.values()))
    if n_sites != 1 and n_sites_dict != 1 and n_sites != n_sites_dict:
        raise ValueError("A list in a dictionary was misshapen. " + LIST_LEN_ERROR_MSG)

    if n_sites_dict == 1:
        _dict = {
            key: val[0] if isinstance(val, (list, np.ndarray)) else val
            for key, val in _dict.items()
        }
        return _dict, 1

    _dict = {key: _extend_to_nparray(val, n_sites_dict) for key, val in _dict.items()}
    return _zip_dict(_dict), n_sites_dict


def _check_lengths(attributes):
    """Ensures all attribute lengths are 1 or maximum attribute length"""
    lengths = np.array([np.asarray(attr).size for attr in attributes])

    if np.all(lengths == 1):
        return 1

    lengths = lengths[np.where(lengths != 1)]
    if np.unique(lengths).size == 1:
        return lengths[0]

    raise ValueError(
        "An array or list was either too short or too long. " + LIST_LEN_ERROR_MSG
    )


# BUG: doctest fails on example code
def _zip_dict(_dict):
    """Makes list of dicts with the same keys and scalar values from dict of lists.
    Single dictionaries of only None will return None.

    Example:
    >>> foo = {'k1': [1, None, 3, 4], 'k2': [5, None, 7, 8], 'k3': [9, None, 11, 12]}
    >>> pprint(_zip_dict(foo))
    [{'k1': 1, 'k2': 5, 'k3': 9},
     None,
     {'k1': 3, 'k2': 7, 'k3': 11},
     {'k1': 4, 'k2': 8, 'k3': 12}]
    """
    lst = [dict(zip(_dict.keys(), v)) for v in zip(*(_dict[k] for k in _dict.keys()))]
    lst = [None if np.all([item is None for item in d.values()]) else d for d in lst]
    return lst


# def _check_input_list_lengths(attributes):
#     """Ensures all input list lengths are the same"""
#     lengths = np.asarray(
#         [
#             attr.size if isinstance(attr, np.ndarray) else list(attr.values())[0].size
#             for attr in attributes
#         ]
#     )
#     if np.any(lengths != lengths[0]):
#         bad_list = attributes[np.where(lengths != lengths[0])[0][0]]
#         good_len = attributes[0].size
#         raise ValueError(
#             "An array or list was either too short or too long. "
#             + LIST_LEN_ERROR_MSG
#             + f"{bad_list} is size ({len(bad_list)}) should be size ({good_len})"
#         )
#     return


# def _clean_item(item, n):
#     """Cleans passed item to np.array"""
#     # Return flattened np.array if item is already list or array
#     if isinstance(item, (list, np.ndarray)):
#         return np.hstack(np.asarray(item, dtype=object))
#     # Return default value extended to specified length
#     return np.asarray([item for _ in range(n)])


# def _get_shielding_info(shielding_symmetric):
#     n_ss, shield_keys = [], []
#     if shielding_symmetric is not None:
#         shield_keys = shielding_symmetric.keys()
#         shielding_symmetric = {
#             item: _flatten_item(shielding_symmetric[item])
#             for item in SHIELDING_SYM_PARAMS
#             if item in shield_keys
#         }
#         n_ss = [
#             _get_length(shielding_symmetric[item])
#             for item in SHIELDING_SYM_PARAMS
#             if item in shield_keys
#         ]
#     return n_ss, shield_keys, shielding_symmetric


# def _get_quad_info(quadrupolar):
#     n_q, quad_keys = [], []
#     if quadrupolar is not None:
#         quad_keys = quadrupolar.keys()
#         quadrupolar = {
#             item: _flatten_item(quadrupolar[item])
#             for item in QUADRUPOLAR_PARAMS
#             if item in quad_keys
#         }
#         n_q = [
#             _get_length(quadrupolar[item])
#             for item in QUADRUPOLAR_PARAMS
#             if item in quad_keys
#         ]
#     return n_q, quad_keys, quadrupolar


# def _populate_quadrupolar(sys, items):
#     n = len(items[0])
#     for i in range(n):
#         if sys[i].sites[0].isotope.spin > 0.5:
#             sys[i].sites[0].quadrupolar = {
#                 "Cq": items[0][i],
#                 "eta": items[1][i],
#                 "alpha": items[2][i],
#                 "beta": items[3][i],
#                 "gamma": items[4][i],
#             }


# def _populate_shielding(sys, items):
#     n = len(items[0])
#     for i in range(n):
#         sys[i].sites[0].shielding_symmetric = {
#             "zeta": items[0][i],
#             "eta": items[1][i],
#             "alpha": items[2][i],
#             "beta": items[3][i],
#             "gamma": items[4][i],
#         }


# def _extend_defaults_to_list(item, n):
#     """Returns np.array if item is array-like, otherwise n-length list of item"""
#     if isinstance(item, (list, np.ndarray)):
#         return np.asarray(item)
#     return np.asarray([item for _ in range(n)])


# def _get_length(item):
#     """Return length of item if item is array-like, otherwise 0"""
#     if isinstance(item, (list, np.ndarray)):
#         return np.asarray(item).size
#     return 0


# def _check_size(n_list):
#     index = np.where(n_list > 0)
#     n_list_reduced = n_list[index]
#     first_item = n_list_reduced[0]
#     if np.all(n_list_reduced == first_item):
#         return first_item
#     raise ValueError(
#         "Each entry can either be a single item or a list of items. If an entry is a "
#         "list, it's length must be equal to the length of other lists present in the "
#         "system."
#     )
