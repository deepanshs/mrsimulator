# -*- coding: utf-8 -*-
from mrsimulator import Site
from mrsimulator import SpinSystem
from mrsimulator.methods import Method1D


def test_warnings():
    s = SpinSystem(sites=[Site(isotope="23Na")])
    m = Method1D(channels=["1H"])
    assert m.get_transition_pathways(s) == []
