# -*- coding: utf-8 -*-
import numpy as np
from scipy.sparse import csr_matrix
from scipy.sparse.csgraph import connected_components

# import copy


def split_spin_system(spin_system):  # 9 lines
    """Simplifies user-defined spin systems into irreducible spin system objects.

    Args:
        spin_system: SpinSystem object to simplify.
    """
    ss = spin_system
    return (
        [ss.__class__(sites=[site], abundance=ss.abundance) for site in ss.sites]
        if not ss.couplings
        else build_new_systems(ss)
    )


def new_systems_needed_np(coupling_indexes, n_sites):
    """Calculates the connect spin system components.

    Args:
        coupling_indexes: A list of coupled site indexes.
        n_sites: The number of sites within the spin system.
    """
    adj_matrix = np.zeros((n_sites, n_sites))
    for (idx1, idx2) in coupling_indexes:
        adj_matrix[idx1, idx2] = 1
        adj_matrix[idx2, idx1] = 1
    graph = csr_matrix(adj_matrix)
    n_components, labels = connected_components(graph, return_labels=True)
    systems_needed = [[] for _ in range(n_components)]
    for i, item in enumerate(labels):
        systems_needed[item].append(i)
    return systems_needed


def build_new_systems(spin_system):  # 10 lines
    """Builds a list of irreducible spin systems from a given reducible spin system.

    Args:
        spin_system: SpinSystem object.
    """
    coupling_list = [coupling.site_index for coupling in spin_system.couplings]
    num_sites = len(spin_system.sites)
    set_of_sets = new_systems_needed_np(coupling_list, num_sites)
    return [build_new_system(one_set, spin_system) for one_set in set_of_sets]


def build_new_system(one_set, spin_system):  # 10-12 lines
    """Builds an irreducible spin system from connected site indices

    Args:
        one_set: A set of indexes corresponding to the connected sites.
        spin_system: SpinSystem object.
    """

    def get_couplings(one_set, couplings):
        new = []
        _ = [
            new.append(_)
            for _ in [cp for i in one_set for cp in couplings if i in cp.site_index]
            if _ not in new
        ]
        return new

    abundance = spin_system.abundance
    new_sites = [spin_system.sites[idx] for idx in list(one_set)]
    new_couplings = get_couplings(one_set, spin_system.couplings)

    options = {} if len(new_couplings) == 0 else {"couplings": new_couplings}
    new_sys = spin_system.__class__(sites=new_sites, abundance=abundance, **options)
    fix_coupling_index(new_sys, list(one_set)) if new_sys.couplings else None
    return new_sys


def fix_coupling_index(new_spin_sys, map):  # 8 lines
    """Uses mapping list to correct site indexes of coupling objects."""
    for coupling in new_spin_sys.couplings:
        coupling.site_index = [map.index(item) for item in coupling.site_index]
