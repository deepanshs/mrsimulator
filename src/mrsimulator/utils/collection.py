# -*- coding: utf-8 -*-
from typing import Dict
from typing import List
from typing import Union

import numpy as np
from mrsimulator import Site
from mrsimulator import SpinSystem

__author__ = ["Deepansh Srivastava", "Matthew D. Giammar"]
__email__ = ["srivastava.89@osu.edu", "giammar.7@buckeyemail.osu.edu"]

SHIELDING_SYM_PARAMS = ["zeta", "eta", "alpha", "beta", "gamma"]
QUADRUPOLAR_PARAMS = ["Cq", "eta", "alpha", "beta", "gamma"]
LIST_LEN_ERROR_MSG = (
    "All arguments must be the same size. If one attribute is a type list of length n, "
    "then all attributes with list types must also be of length n, and all remaining "
    "attributes must be scalar (singular float, int, or str)."
)


def single_site_system_generator(
    isotope: Union[str, List[str]],
    isotropic_chemical_shift: Union[float, List[float], np.ndarray] = 0,
    shielding_symmetric: Dict = None,
    shielding_antisymmetric: Dict = None,
    quadrupolar: Dict = None,
    abundance: Union[float, List[float], np.ndarray] = None,
    site_name: Union[str, List[str]] = None,
    site_label: Union[str, List[str]] = None,
    site_description: Union[str, List[str]] = None,
    rtol: float = 1e-3,
) -> List[SpinSystem]:
    r"""Generate and return a list of single-site spin systems from the input parameters.

    Args:
        isotope:
            A required string or a list of site isotopes.
        isotropic_chemical_shift:
            A float or a list/ndarray of isotropic chemical shifts per site per spin
            system. The default is 0.
        shielding_symmetric:
            A shielding symmetric dict object, where the keyword value can either
            be a float or a list/ndarray of floats. The default value is None. The
            allowed keywords are ``zeta``, ``eta``, ``alpha``, ``beta``, and ``gamma``.
        shielding_antisymmetric:
            A shielding antisymmetric dict object, where the keyword value can either
            be a float or a list/ndarray of floats. The default value is None. The
            allowed keywords are ``zeta``, ``alpha``, and ``beta``.
        quadrupolar:
            A quadrupolar dict object, where the keyword value can either be a float or
            a list/ndarray of floats. The default value is None. The allowed keywords
            are ``Cq``, ``eta``, ``alpha``, ``beta``, and ``gamma``.
        abundance:
            A float or a list/ndarray of floats describing the abundance of each spin
            system.
        site_name:
            A string or a list of strings with site names per site per spin system. The
            default is None.
        site_label:
            A string or a list of strings with site labels per site per spin system. The
            default is None.
        site_description:
            A string or a list of strings with site descriptions per site per spin
            system. The default is None.
        rtol:
            The relative tolerance used in determining the cutoff abundance, given as,
            :math:`\tt{abundance}_{\tt{cutoff}} = \tt{rtol} * \tt{max(abundance)}.`
            The spin systems with abundance below this threshold are ignored.

    Returns:
        List of :ref:`spin_sys_api` objects with a single :ref:`site_api`

    Example:
        **Single spin system:**

        >>> sys1 = single_site_system_generator(
        ...     isotope=["1H"],
        ...     isotropic_chemical_shift=10,
        ...     site_name="Single Proton",
        ... )
        >>> print(len(sys1))
        1

        **Multiple spin system:**

        >>> sys2 = single_site_system_generator(
        ...     isotope="1H",
        ...     isotropic_chemical_shift=[10] * 5,
        ...     site_name="5 Protons",
        ... )
        >>> print(len(sys2))
        5

        **Multiple spin system with dictionary arguments:**

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
        the spin systems. If the parameter value is a list or ndarray, its `ith` value
        is assigned to the respective parameter of the `ith` spin system. When multiple
        parameter values are given as lists/ndarrays, the length of all the lists must
        be the same.
    """
    sites = site_generator(
        isotope=isotope,
        isotropic_chemical_shift=isotropic_chemical_shift,
        shielding_symmetric=shielding_symmetric,
        shielding_antisymmetric=shielding_antisymmetric,
        quadrupolar=quadrupolar,
        name=site_name,
        label=site_label,
        description=site_description,
    )
    n_sites = len(sites)

    abundance = 1 / n_sites if abundance is None else abundance
    abundance = _extend_to_nparray(_fix_item(abundance), n_sites)
    n_abd = abundance.size

    if n_sites == 1:
        sites = np.asarray([sites[0] for _ in range(n_abd)])
        n_sites = sites.size

    if n_sites != n_abd:
        raise ValueError(
            "Number of sites does not match the number of abundances. "
            f"{LIST_LEN_ERROR_MSG}"
        )

    keep_idxs = np.where(abundance > rtol * abundance.max())[0]

    return [
        SpinSystem(sites=[site], abundance=abd)
        for site, abd in zip(sites[keep_idxs], abundance[keep_idxs])
    ]


def site_generator(
    isotope: Union[str, List[str]],
    isotropic_chemical_shift: Union[float, List[float], np.ndarray] = 0,
    shielding_symmetric: Dict = None,
    shielding_antisymmetric: Dict = None,
    quadrupolar: Dict = None,
    name: Union[str, List[str]] = None,
    label: Union[str, List[str]] = None,
    description: Union[str, List[str]] = None,
) -> List[Site]:
    r"""Generate a list of Site objects from lists of site attributes.

    Args:
        isotope:
            A required string or a list of site isotopes.
        isotropic_chemical_shift:
            A float or a list/ndarray of isotropic chemical shifts per site. The default
            is 0.
        shielding_symmetric:
            A shielding symmetric dict object, where the keyword value can either
            be a float or a list/ndarray of floats. The default value is None. The
            allowed keywords are ``zeta``, ``eta``, ``alpha``, ``beta``, and ``gamma``.
        shielding_antisymmetric:
            A shielding antisymmetric dict object, where the keyword value can either
            be a float or a list/ndarray of floats. The default value is None. The
            allowed keywords are ``zeta``, ``alpha``, and ``beta``.
        quadrupolar:
            A quadrupolar dict object, where the keyword value can either be a float or
            a list/ndarray of floats. The default value is None. The allowed keywords
            are ``Cq``, ``eta``, ``alpha``, ``beta``, and ``gamma``.
        name:
            A string or a list of strings with site names per site. The default is None.
        label:
            A string or a list of strings with site labels per site. The default is
            None.
        description:
            A string or a list of strings with site descriptions per site. The default
            is None.

    Returns:
        sites: List of :ref:`site_api` objects

    Example:
        **Generating 10 hydrogen sites:**

        >>> sites1 = site_generator(
        ...     isotope=["1H"] * 10,
        ...     isotropic_chemical_shift=-15,
        ...     name="10 Protons",
        ... )
        >>> print(len(sites1))
        10

        **Generating 10 hydrogen sites with different shifts:**

        >>> shifts = np.arange(-10, 10, 2)
        >>> sites2 = site_generator(
        ...     isotope=["1H"] * 10,
        ...     isotropic_chemical_shift=shifts,
        ...     name="10 Proton",
        ... )
        >>> print(len(sites2))
        10

        **Generating multiple sites with dictionary arguments:**

        >>> Cq = [4.2e6] * 12
        >>> sys3 = site_generator(
        ...     isotope="17O",
        ...     isotropic_chemical_shift=60.0,  # in ppm,
        ...     quadrupolar={"Cq": Cq, "eta": 0.5},  # Cq in Hz
        ... )
        >>> print(len(sys3))
        12
    """
    lst = [isotope, isotropic_chemical_shift, name, label, description]
    attributes = [_fix_item(item) for item in lst]

    n_sites = _check_lengths(attributes)

    lst_extend = [shielding_symmetric, shielding_antisymmetric, quadrupolar]
    for obj in lst_extend:
        if obj is not None:
            obj_ext, n_dict = _extend_dict_values(obj, n_sites)
            n_sites = max(n_sites, n_dict)
            attributes.append(obj_ext)
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
    """Flattens multidimensional arrays into 1d array."""
    if isinstance(item, (list, np.ndarray)):
        return np.asarray(item).ravel()
    return item


def _extend_to_nparray(item, n):
    """If item is already list/array return np.array, otherwise extend to length n."""
    data = item if isinstance(item, (list, np.ndarray)) else [item for _ in range(n)]
    return np.asarray(data)


def _extend_dict_values(_dict, n_sites):
    """Checks and extends dict values. Returns dict or list of dicts and max length."""
    _dict = {key: _fix_item(val) for key, val in _dict.items()}
    n_sites_dict = _check_lengths(list(_dict.values()))
    if n_sites != 1 and n_sites_dict != 1 and n_sites != n_sites_dict:
        raise ValueError(f"A list in a dictionary was misshapen. {LIST_LEN_ERROR_MSG}")

    if n_sites_dict == 1:
        _dict = {
            key: val[0] if isinstance(val, (list, np.ndarray)) else val
            for key, val in _dict.items()
        }
        return _dict, 1

    _dict = {key: _extend_to_nparray(val, n_sites_dict) for key, val in _dict.items()}
    return _zip_dict(_dict), n_sites_dict


def _check_lengths(attributes):
    """Ensures all attribute lengths are 1 or maximum attribute length."""
    lengths = np.array([np.asarray(attr).size for attr in attributes])

    if np.all(lengths == 1):
        return 1

    lengths = lengths[np.where(lengths != 1)]
    if np.unique(lengths).size == 1:
        return lengths[0]

    raise ValueError(
        f"An array or list was either too short or too long. {LIST_LEN_ERROR_MSG}"
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
