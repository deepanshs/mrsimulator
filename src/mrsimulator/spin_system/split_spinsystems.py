# -*- coding: utf-8 -*-
import numpy as np
from scipy.sparse import csr_matrix
from scipy.sparse.csgraph import connected_components

# import copy


def split_spin_system(spin_system):  # 9 lines
    """Simplifies user-defined spin systems into irreducible spin
    system objects"""
    if spin_system.couplings:
        coupling_list = [i.site_index for i in spin_system.couplings]
        # eigenvectors = build_eigenvectors(
        #   len(spin_system.sites), spin_system.couplings
        # )
        # eigenvectors = laplacian_eigenvectors(coupling_list, len(spin_system.sites))
        # systems_needed = new_systems_needed_matrix(eigenvectors)
        # systems_needed = new_systems_needed_nosets(eigenvectors)
        systems_needed = new_systems_needed_np(coupling_list, len(spin_system.sites))
        new_spinsystem_list = build_new_systems(spin_system, systems_needed)
    else:
        abundance = spin_system.abundance
        new_spinsystem_list = []
        for site in spin_system.sites:
            new_spinsystem_list.append(
                spin_system.__class__(sites=[site], abundance=abundance)
            )
    print(new_spinsystem_list)
    return new_spinsystem_list


# def build_eigenvectors(num_sites, couplings):  # 8 lines
#     """builds a matrix that describes which sites are coupled"""
#     spin_matrix = np.zeros((num_sites, num_sites))
#     for coupling in couplings:
#         idx = tuple(coupling.site_index)
#         spin_matrix[idx] = 2
#         idx2 = tuple((list(idx))[::-1])
#         spin_matrix[idx2] = 2
#     _, eigenvectors = np.linalg.eig(spin_matrix)
#     # print(eigenvectors)
#     return spin_matrix


def new_systems_needed_np(coupling_list, num_sites):
    adj_matrix = np.zeros((num_sites, num_sites))
    for (idx1, idx2) in coupling_list:
        adj_matrix[idx1, idx2] = 1
        adj_matrix[idx2, idx1] = 1
    graph = csr_matrix(adj_matrix)
    n_components, labels = connected_components(graph, return_labels=True)
    systems_needed = [[] for i in range(n_components)]
    print(labels)
    for i, item in enumerate(labels):
        systems_needed[item].append(i)
    print(systems_needed)
    return systems_needed


# def new_systems_needed_matrix(eigenvectors):  # 8 lines
#     """Uses eigenvector matrix to build a set of sets describing spin systems
#     that need to be made"""
#     # x = np.where(eigenvectors != 0)
#     set_of_sets = set()
#     for i in range(len(eigenvectors)):
#         running_set = set()
#         for j in range(len(eigenvectors)):
#             if not np.allclose(eigenvectors[j, i], 0):
#                 running_set.add(j)
#         set_of_sets.add(frozenset(running_set))
#
#     pre_len = len(set_of_sets)
#     loop_count = 0
#     keep_going = True
#     while keep_going:
#         good_set = set()
#         for group1 in list(set_of_sets):
#             for group2 in (list(set_of_sets))[::-1]:
#                 if not group1.isdisjoint(set(group2)):
#                     group1 = set(group1).union(group2)
#                     good_set.add(frozenset(group1))
#         post_len = len(good_set)
#         if loop_count > 0 and pre_len == post_len:
#             keep_going = False
#         set_of_sets = copy.deepcopy(good_set)
#         pre_len = copy.deepcopy(post_len)
#         loop_count += 1
#     return set_of_sets


# def new_systems_needed_nosets(eigenvectors):  # 10 lines
#     """Uses eigenvector matrix to build a list of tuples describing
#     spin systems that need to be made"""
#     systems_needed = list()
#     for i in range(len(eigenvectors)):
#         system = list()
#         for j in range(len(eigenvectors)):
#             if not np.allclose(eigenvectors[j, i], 0):
#                 if j not in system:
#                     system.append(j)
#         if tuple(system) not in systems_needed:
#             systems_needed.append(tuple(system))
#     return systems_needed


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
    abundance = spin_system.abundance
    for site_idx in list(one_set):
        mapping_dict[site_idx] = {"new system": sys_idx, "new site": len(new_sites)}
        new_sites.append(spin_system.sites[site_idx])
        for coupling in spin_system.couplings:
            if site_idx in coupling.site_index:
                if coupling not in new_couplings:
                    new_couplings.append(coupling)
    if len(new_couplings) == 0:
        new_sys = spin_system.__class__(sites=new_sites, abundance=abundance)
    else:
        new_sys = spin_system.__class__(
            sites=new_sites, couplings=new_couplings, abundance=abundance
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
