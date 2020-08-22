# -*- coding: utf-8 -*-
import numpy as np
from mrsimulator import Site
from mrsimulator import SpinSystem


def single_site_system_generator(
    isotopes,
    isotropic_chemical_shifts=0,
    shielding_symmetric=None,
    quadrupolar=None,
    abundance=None,
):
    n_isotopes = get_length(isotopes)
    n_iso = get_length(isotropic_chemical_shifts)
    n_abd = get_length(abundance)

    n_ss = []
    if shielding_symmetric is not None:
        shield_keys = shielding_symmetric.keys()
        shielding_symmetric_keys = ["zeta", "eta", "alpha", "beta", "gamma"]
        for item in shielding_symmetric_keys:
            if item in shield_keys:
                n_ss.append(get_length(shielding_symmetric[item]))

    n_q = []
    if quadrupolar is not None:
        quad_keys = quadrupolar.keys()
        quadrupolar_keys = ["Cq", "eta", "alpha", "beta", "gamma"]
        n_q = [
            get_length(quadrupolar[item])
            for item in quadrupolar_keys
            if item in quad_keys
        ]

    n = np.asarray([n_isotopes, n_iso, *n_ss, *n_q, n_abd])
    n_len = check_size(n)

    if abundance is None:
        abundance = 1 / n_len

    # system with only isotope and isotropic shifts parameters
    isotopes_ = get_default_lists(isotopes, n_len)
    iso_chemical_shifts_ = get_default_lists(isotropic_chemical_shifts, n_len)
    abundance_ = get_default_lists(abundance, n_len)

    sys = [
        SpinSystem(
            sites=[Site(isotope=ist__, isotropic_chemical_shift=iso__,)],
            abundance=abd__,
        )
        for ist__, iso__, abd__ in zip(isotopes_, iso_chemical_shifts_, abundance_)
    ]

    # system with additional shielding symmetric parameters
    if shielding_symmetric is not None:
        lst = []
        for item in shielding_symmetric_keys:
            if item in shield_keys:
                lst.append(get_default_lists(shielding_symmetric[item], n_len))
            else:
                lst.append([None for _ in range(n_len)])
        _populate_shielding(sys, lst)

    # system with additional quadrupolar parameters
    if quadrupolar is not None:
        lst = []
        for item in quadrupolar_keys:
            if item in quad_keys:
                lst.append(get_default_lists(quadrupolar[item], n_len))
            else:
                lst.append([None for _ in range(n_len)])

        _populate_quadrupolar(sys, lst)

    sys_ = [item for item in sys if item.abundance != 0]
    return sys_


def _populate_quadrupolar(sys, items):
    n = len(items[0])
    for i in range(n):
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


def get_default_lists(item, n):
    if isinstance(item, (list, np.ndarray)):
        return np.asarray(item)

    return np.asarray([item for _ in range(n)])


def get_length(item):
    """Return the length of item it item is a list, else 0."""
    if isinstance(item, (list, np.ndarray)):
        return np.asarray(item).size
    return 0


def check_size(n_list):
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
