# -*- coding: utf-8 -*-
import numpy as np


def split_spin_system(spin_system):  # 9 lines
    """Simplifies user-defined spin systems into irreducible spin system objects"""
    if spin_system.couplings:
        eigenvectors = build_eigenvectors(len(spin_system.sites), spin_system.couplings)
        systems_needed = new_systems_needed_nosets(eigenvectors)
        # new_systems_needed_matrix(eigenvectors)
        new_spinsystem_list = build_new_systems(spin_system, systems_needed)
    else:
        new_spinsystem_list = []
        for site in spin_system.sites:
            new_spinsystem_list.append(
                spin_system.__class__(sites=[site], abundance=spin_system.abundance)
            )
    return new_spinsystem_list


def build_eigenvectors(num_sites, couplings):  # 8 lines
    """builds a matrix that describes which sites are coupled"""
    spin_matrix = np.identity(num_sites)
    for coupling in couplings:
        idx = tuple(coupling.site_index)
        spin_matrix[idx] = 0.5
        idx2 = tuple((list(idx))[::-1])
        spin_matrix[idx2] = 0.5
    _, eigenvectors = np.linalg.eig(spin_matrix)
    return eigenvectors


def new_systems_needed_matrix(eigenvectors):  # 8 lines
    """Uses eigenvector matrix to build a set of sets describing spin systems
    that need to be made"""
    # x = np.where(eigenvectors != 0)
    set_of_sets = set()
    for i in range(len(eigenvectors)):
        running_set = set()
        for j in range(len(eigenvectors)):
            if not np.allclose(eigenvectors[j, i], 0):
                running_set.add(j)
        set_of_sets.add(frozenset(running_set))
    return set_of_sets


def new_systems_needed_nosets(eigenvectors):  # 10 lines
    """Uses eigenvector matrix to build a list of tuples describing
    spin systems that need to be made"""
    systems_needed = list()
    for i in range(len(eigenvectors)):
        system = list()
        for j in range(len(eigenvectors)):
            if not np.allclose(eigenvectors[j, i], 0):
                if j not in system:
                    system.append(j)
        if tuple(system) not in systems_needed:
            systems_needed.append(tuple(system))
    return systems_needed


def build_new_systems(spin_system, set_of_sets):  # 10 lines
    """builds a list of irreducible spin systems from a given reducible spin
    system and a set of sets that tells what systems we need to build."""
    mapping_dict = dict()
    irr_spin_systems = list()
    for one_set in set_of_sets:
        new_system, mapping_dict = build_new_system(
            one_set, spin_system, mapping_dict, len(irr_spin_systems)
        )
        irr_spin_systems.append(new_system)
    for i, spin_sys in enumerate(irr_spin_systems):
        if spin_sys.couplings:
            irr_spin_systems[i] = fix_coupling_index(spin_sys, mapping_dict)
    return irr_spin_systems


def build_new_system(one_set, spin_system, mapping_dict, sys_idx):  # 10-12 lines
    """Takes in a set of site indices and builds an irreducible spin system
    containing those sites and their couplings"""
    new_sites = list()
    new_couplings = list()
    for site_idx in list(one_set):
        mapping_dict[site_idx] = {"new system": sys_idx, "new site": len(new_sites)}
        new_sites.append(spin_system.sites[site_idx])
        for coupling in spin_system.couplings:
            if site_idx in coupling.site_index and coupling not in new_couplings:
                new_couplings.append(coupling)
    if len(new_couplings) == 0:
        new_sys = spin_system.__class__(
            sites=new_sites, abundance=spin_system.abundance
        )
    else:
        new_sys = spin_system.__class__(
            sites=new_sites, couplings=new_couplings, abundance=spin_system.abundance
        )
    return new_sys, mapping_dict


def fix_coupling_index(new_spin_sys, mapping_dict):  # 8 lines
    """Uses mapping_dict to correct site indexes of coupling
    objects"""
    couplings = new_spin_sys.couplings
    for i, coupling in enumerate(couplings):
        old_idx1, old_idx2 = coupling.site_index
        new_idx1 = mapping_dict[old_idx1]["new site"]
        new_idx2 = mapping_dict[old_idx2]["new site"]
        couplings[i].site_index = [new_idx1, new_idx2]
    new_spin_sys.couplings = couplings
    return new_spin_sys
