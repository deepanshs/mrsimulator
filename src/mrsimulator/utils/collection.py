# -*- coding: utf-8 -*-
import numpy as np
from mrsimulator import Site
from mrsimulator import SpinSystem

__author__ = "Deepansh Srivastava"
__email__ = "srivastava.89@osu.edu"


SHIELDING_SYM_PARAMS = ["zeta", "eta", "alpha", "beta", "gamma"]
QUADRUPOLAR_PARAMS = ["Cq", "eta", "alpha", "beta", "gamma"]


def single_site_system_generator(
    isotopes,
    isotropic_chemical_shifts=0,
    shielding_symmetric=None,
    quadrupolar=None,
    abundance=None,
    site_labels=None,
    site_names=None,
    site_descriptions=None,
    rtol=1e-3,
):
    r"""Generate and return a list of single-site spin systems from the input parameters.

    Args:

        isotopes:
            A string or a list of site isotopes.
        isotropic_chemical_shifts:
            A float or a list/ndarray of values. The default value is 0.
        shielding_symmetric:
            A shielding symmetric like dict object, where the keyword value can either
            be a float or a list/ndarray of floats. The default value is None. The
            allowed keywords are ``zeta``, ``eta``, ``alpha``, ``beta``, and ``gamma``.
        quadrupolar:
            A quadrupolar like dict object, where the keyword value can either be a
            float or a list/ndarray of floats. The default value is None. The allowed
            keywords are ``Cq``, ``eta``, ``alpha``, ``beta``, and ``gamma``.
        abundance:
            A float or a list/ndarray of floats describing the abundance of each spin
            system.
        site_labels:
            A string or a list of strings each with a site label. The default is None.
        site_names:
            A string or a list of strings each with a site name. The default is None.
        site_descriptions:
            A string or a list of strings each with a site description. Default is None.
        rtol:
            The relative tolerance. This value is used in determining the cutoff
            abundance given as
            :math:`\tt{abundance}_\tt{cutoff} = \tt{rtol} * \tt{max(abundance)}.`
            The spin systems with abundance below this threshold are ignored.

    .. note::
        The parameter value can either be a float or a list/ndarray. If the parameter
        value is a float, the given value is assigned to the respective parameter in all
        the spin systems. If the parameter value is a list or ndarray, its ith value is
        assigned to the respective parameter of the ith spin system. When multiple
        parameter values are given as lists/ndarrays, the length of all the lists must
        be the same.
    """
    isotopes = _fix_item(isotopes)
    n_isotopes = _get_length(isotopes)

    isotropic_chemical_shifts = _fix_item(isotropic_chemical_shifts)
    n_iso = _get_length(isotropic_chemical_shifts)

    abundance = _fix_item(abundance)
    n_abd = _get_length(abundance)

    site_labels = _fix_item(site_labels)
    n_site_labels = _get_length(site_labels)

    site_names = _fix_item(site_names)
    n_site_names = _get_length(site_names)

    site_descriptions = _fix_item(site_descriptions)
    n_site_descriptions = _get_length(site_descriptions)

    n_ss, shield_keys, shielding_symmetric = _get_shielding_info(shielding_symmetric)
    n_q, quad_keys, quadrupolar = _get_quad_info(quadrupolar)

    n_len = _check_size(
        np.asarray(
            [
                n_isotopes,
                n_iso,
                *n_ss,
                *n_q,
                n_abd,
                n_site_labels,
                n_site_names,
                n_site_descriptions,
            ]
        )
    )

    if abundance is None:
        abundance = 1 / n_len

    # system with only isotope and isotropic shifts parameters
    isotopes_ = _get_default_lists(isotopes, n_len)
    iso_chemical_shifts_ = _get_default_lists(isotropic_chemical_shifts, n_len)
    abundance_ = _get_default_lists(abundance, n_len)
    site_labels_ = _get_default_lists(site_labels, n_len)
    site_names_ = _get_default_lists(site_names, n_len)
    site_descriptions_ = _get_default_lists(site_descriptions, n_len)

    index = np.where(abundance_ > rtol * abundance_.max())[0]

    sys = [
        SpinSystem(
            sites=[
                Site(
                    isotope=ist__,
                    isotropic_chemical_shift=iso__,
                    label=lab__,
                    name=name__,
                    description=dis__,
                )
            ],
            abundance=abd__,
        )
        for ist__, iso__, lab__, name__, dis__, abd__ in zip(
            isotopes_[index],
            iso_chemical_shifts_[index],
            site_labels_[index],
            site_names_[index],
            site_descriptions_[index],
            abundance_[index],
        )
    ]

    # system with additional shielding symmetric parameters
    if shielding_symmetric is not None:
        lst = [
            _get_default_lists(shielding_symmetric[item], n_len)[index]
            if item in shield_keys
            else [None for _ in range(index.size)]
            for item in SHIELDING_SYM_PARAMS
        ]
        _populate_shielding(sys, lst)

    # system with additional quadrupolar parameters
    if quadrupolar is not None:
        lst = [
            _get_default_lists(quadrupolar[item], n_len)[index]
            if item in quad_keys
            else [None for _ in range(index.size)]
            for item in QUADRUPOLAR_PARAMS
        ]
        _populate_quadrupolar(sys, lst)

    return sys


def _get_shielding_info(shielding_symmetric):
    n_ss, shield_keys = [], []
    if shielding_symmetric is not None:
        shield_keys = shielding_symmetric.keys()
        shielding_symmetric = {
            item: _fix_item(shielding_symmetric[item])
            for item in SHIELDING_SYM_PARAMS
            if item in shield_keys
        }
        n_ss = [
            _get_length(shielding_symmetric[item])
            for item in SHIELDING_SYM_PARAMS
            if item in shield_keys
        ]
    return n_ss, shield_keys, shielding_symmetric


def _get_quad_info(quadrupolar):
    n_q, quad_keys = [], []
    if quadrupolar is not None:
        quad_keys = quadrupolar.keys()
        quadrupolar = {
            item: _fix_item(quadrupolar[item])
            for item in QUADRUPOLAR_PARAMS
            if item in quad_keys
        }
        n_q = [
            _get_length(quadrupolar[item])
            for item in QUADRUPOLAR_PARAMS
            if item in quad_keys
        ]
    return n_q, quad_keys, quadrupolar


def _populate_quadrupolar(sys, items):
    n = len(items[0])
    for i in range(n):
        if sys[i].sites[0].isotope.spin > 0.5:
            sys[i].sites[0].quadrupolar = {
                "Cq": items[0][i],
                "eta": items[1][i],
                "alpha": items[2][i],
                "beta": items[3][i],
                "gamma": items[4][i],
            }


def _populate_shielding(sys, items):
    n = len(items[0])
    for i in range(n):
        sys[i].sites[0].shielding_symmetric = {
            "zeta": items[0][i],
            "eta": items[1][i],
            "alpha": items[2][i],
            "beta": items[3][i],
            "gamma": items[4][i],
        }


def _get_default_lists(item, n):
    if isinstance(item, (list, np.ndarray)):
        return np.asarray(item)

    return np.asarray([item for _ in range(n)])


def _fix_item(item):
    if isinstance(item, (list, np.ndarray)):
        return np.asarray(item).ravel()
    return item


def _get_length(item):
    """Return the length of item it item is a list, else 0."""
    if isinstance(item, (list, np.ndarray)):
        return np.asarray(item).size
    return 0


def _check_size(n_list):
    index = np.where(n_list > 0)
    n_list_reduced = n_list[index]
    first_item = n_list_reduced[0]
    if np.all(n_list_reduced == first_item):
        return first_item
    raise ValueError(
        "Each entry can either be a single item or a list of items. If an entry is a "
        "list, it's length must be equal to the length of other lists present in the "
        "system."
    )
