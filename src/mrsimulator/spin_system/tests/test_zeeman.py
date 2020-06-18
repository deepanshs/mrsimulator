# -*- coding: utf-8 -*-
"""Zeeman State Tests"""
from itertools import product

from mrsimulator import SpinSystem


def test_zeeman_energy_states():

    # Spin 1/2
    site_H = {"isotope": "1H", "isotropic_chemical_shift": 0}

    # Spin 5/2
    site_O = {"isotope": "17O", "isotropic_chemical_shift": 0}

    iso_H = SpinSystem(sites=[site_H])
    iso_O = SpinSystem(sites=[site_O])
    iso_OH = SpinSystem(sites=[site_O, site_H])

    H_Zeeman = [-0.5, 0.5]
    O_Zeeman = [-2.5, -1.5, -0.5, 0.5, 1.5, 2.5]

    for i, state in enumerate(iso_H.zeeman_energy_states()):
        assert state.tolist() == [H_Zeeman[i]]

    for i, state in enumerate(iso_O.zeeman_energy_states()):
        assert state.tolist() == [O_Zeeman[i]]

    for i, state in enumerate(iso_OH.zeeman_energy_states()):
        assert tuple(state.tolist()) == list(product(O_Zeeman, H_Zeeman))[i]


def test_all_transitions():

    # Spin 1/2
    site_H = {"isotope": "1H", "isotropic_chemical_shift": 0}

    # Spin 5/2
    site_O = {"isotope": "17O", "isotropic_chemical_shift": 0}

    iso_H = SpinSystem(sites=[site_H])
    iso_O = SpinSystem(sites=[site_O])
    iso_OH = SpinSystem(sites=[site_O, site_H])

    H_Zeeman = [-0.5, 0.5]
    O_Zeeman = [-2.5, -1.5, -0.5, 0.5, 1.5, 2.5]

    iso_H_transitions = iso_H.all_transitions()
    for i, transition in enumerate(iso_H_transitions):
        res = list(product(H_Zeeman, H_Zeeman))[i]
        assert tuple(transition.initial + transition.final) == res
        # assert transition.tolist() == list(product(H_Zeeman, H_Zeeman))[i]

    iso_O_transitions = iso_O.all_transitions()
    for i, transition in enumerate(iso_O_transitions):
        res = list(product(O_Zeeman, O_Zeeman))[i]
        assert tuple(transition.initial + transition.final) == res

    iso_OH_transitions = iso_OH.all_transitions()
    for i, transition in enumerate(iso_OH_transitions):
        OH_state = list(product(O_Zeeman, H_Zeeman))
        res = list(product(OH_state, OH_state))[i]
        assert tuple([tuple(transition.initial), tuple(transition.final)]) == res
