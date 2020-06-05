# -*- coding: utf-8 -*-
# from mrsimulator.abstract_list import TransitionList
from mrsimulator.transition import Transition


def test_transition_1():
    # set up
    a = {"initial": [0.5, 1.5], "final": [-0.5, 1.5]}
    tran = Transition(**a)
    assert tran.initial == [0.5, 1.5]
    assert tran.final == [-0.5, 1.5]

    # to dict with unit
    assert tran.to_dict_with_units() == a


# def test_transition_list_1():
#     a = TransitionList()
#     assert a == list()

#     a.append({})
