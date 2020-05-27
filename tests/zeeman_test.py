# -*- coding: utf-8 -*-
"""Zeeman State Tests"""
from itertools import product

from mrsimulator import Isotopomer


def test_zeeman_energy_states():

    # Spin 1/2
    site_H = {"isotope": "1H", "isotropic_chemical_shift": 0}

    # Spin 5/2
    site_O = {"isotope": "17O", "isotropic_chemical_shift": 0}

    iso_H = Isotopomer(sites=[site_H])
    iso_O = Isotopomer(sites=[site_O])
    iso_OH = Isotopomer(sites=[site_O, site_H])

    H_Zeeman = [-0.5, 0.5]
    O_Zeeman = [-2.5, -1.5, -0.5, 0.5, 1.5, 2.5]

    for i, state in enumerate(iso_H.Zeeman_energy_states):
        assert state.tolist() == [H_Zeeman[i]]

    for i, state in enumerate(iso_O.Zeeman_energy_states):
        assert state.tolist() == [O_Zeeman[i]]

    for i, state in enumerate(iso_OH.Zeeman_energy_states):
        assert tuple(state.tolist()) == list(product(O_Zeeman, H_Zeeman))[i]


def test_all_transitions():

    # Spin 1/2
    site_H = {"isotope": "1H", "isotropic_chemical_shift": 0}

    # Spin 5/2
    site_O = {"isotope": "17O", "isotropic_chemical_shift": 0}

    iso_H = Isotopomer(sites=[site_H])
    iso_O = Isotopomer(sites=[site_O])
    iso_OH = Isotopomer(sites=[site_O, site_H])

    H_Zeeman = [-0.5, 0.5]
    O_Zeeman = [-2.5, -1.5, -0.5, 0.5, 1.5, 2.5]

    for i, transition in enumerate(iso_H.all_transitions):
        assert (
            tuple(transition.initial + transition.final) == list(product(H_Zeeman, H_Zeeman))[i]
        )
        # assert transition.tolist() == list(product(H_Zeeman, H_Zeeman))[i]

    for i, transition in enumerate(iso_O.all_transitions):
        assert (
            tuple(transition.initial + transition.final) == list(product(O_Zeeman, O_Zeeman))[i]
        )

    for i, transition in enumerate(iso_OH.all_transitions):
        OH_state = list(product(O_Zeeman, H_Zeeman))
        assert (
            tuple([tuple(transition.initial), tuple(transition.final)]) == list(product(OH_state, OH_state))[i]
        )
