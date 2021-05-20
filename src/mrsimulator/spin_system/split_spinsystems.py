# -*- coding: utf-8 -*-
import numpy as np
from mrsimulator import SpinSystem


def split_spin_system(spin_system):  # 9 lines
    """Simplifies user-defined spin systems into irreducible spin system objects"""
    if spin_system.couplings:
        coupling_matrix = build_coupling_matrix(
            len(spin_system.sites), spin_system.couplings
        )
        systems_needed = new_systems_needed_matrix(coupling_matrix)
        new_spinsystem_list = build_new_systems(spin_system, systems_needed)
    else:
        new_spinsystem_list = []
        for site in spin_system.sites:
            new_spinsystem_list.append(
                SpinSystem(sites=[site], abundance=spin_system.abundance)
            )
    return new_spinsystem_list


# old, now using new_systems_needed_matrix ###
# def new_systems_needed(spin_system):  # probably just delete at this point
#    """Takes a coupled spin system and returns a set of frozen sets, each
#    containing the site indexes needed to make an immutable spin system"""
#    coupling_index_list = []
#    for item in coupling_list:
#        coupling_index_list.append(item.site_index)
#
#    set_of_sets = set()
#    for coupling in coupling_index_list:
#        running_set = set()
#        running_set.add(coupling[0])
#        running_set.add(coupling[1])
#        for coupling2 in coupling_index_list:
#            if coupling[0] in coupling2 or coupling[1] in coupling2:
#                running_set.add(coupling2[0])
#                running_set.add(coupling2[1])
#        set_of_sets.add(frozenset(running_set))
#    pre_len = len(set_of_sets)
#    keep_going = True
#    while keep_going:
#        good_set = set()
#        for group1 in set_of_sets:
#            group1 = set(group1)
#            for group2 in set_of_sets:
#                if not group1.isdisjoint(group2):
#                    group1 = group1.union(group2)
#                    good_set.add(frozenset(group1))
#        post_len = len(good_set)
#        if pre_len == post_len:
#            keep_going = False
#        set_of_sets = copy.deepcopy(good_set)
#        pre_len = copy.deepcopy(post_len)
#    return set_of_sets


def build_coupling_matrix(num_sites, couplings):  # 8 lines
    """builds a matrix that describes which sites are coupled"""
    spin_matrix = np.identity(num_sites)
    for coupling in couplings:
        idx = tuple(coupling.site_index)
        spin_matrix[idx] = 1
        idx2 = tuple((list(idx))[::-1])
        spin_matrix[idx2] = 1
    _, coupling_matrix = np.linalg.eig(spin_matrix)
    return coupling_matrix


def new_systems_needed_matrix(coupling_matrix):  # 8 lines
    set_of_sets = set()
    for i in range(len(coupling_matrix)):
        running_set = set()
        for j in range(len(coupling_matrix)):
            if coupling_matrix[j, i]:
                running_set.add(j)
        set_of_sets.add(frozenset(running_set))
    return set_of_sets


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
    new_sys = SpinSystem(
        sites=new_sites, couplings=new_couplings, abundance=spin_system.abundance
    )
    return new_sys, mapping_dict


def fix_coupling_index(new_spin_sys, mapping_dict):  # 8 lines
    couplings = new_spin_sys.couplings
    for i, coupling in enumerate(couplings):
        old_idx1, old_idx2 = coupling.site_index
        new_idx1 = mapping_dict[old_idx1]["new site"]
        new_idx2 = mapping_dict[old_idx2]["new site"]
        couplings[i].site_index = [new_idx1, new_idx2]
    new_spin_sys.couplings = couplings
    return new_spin_sys
