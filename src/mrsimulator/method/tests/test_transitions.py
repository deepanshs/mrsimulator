# -*- coding: utf-8 -*-
from mrsimulator import SpinSystem
from mrsimulator.method import Method
from mrsimulator.transition import TransitionPathway

__author__ = "Deepansh J. Srivastava"
__email__ = "srivastava.89@osu.edu"

method1 = Method(
    channels=["13C"],
    spectral_dimensions=[{"events": [{"transition_query": [{"ch1": {"P": [-1]}}]}]}],
)

method2 = Method(
    channels=["13C"],
    spectral_dimensions=[{"events": [{"transition_query": [{"ch1": {"P": [-2]}}]}]}],
)


def check_transition_set(got, expected):
    assert len(got) == len(expected), "Inconsisten transition pathway count"
    for item in got:
        assert item in expected, f"Transition pathways not found: {item}"


def test_00():
    system = SpinSystem(sites=[{"isotope": "13C"}])
    tr = method1.get_transition_pathways(system)

    expected = [TransitionPathway([{"final": [-0.5], "initial": [0.5]}])]
    check_transition_set(tr, expected)


def test_01():
    system = SpinSystem(sites=[{"isotope": "13C"}, {"isotope": "1H"}])
    tr = method1.get_transition_pathways(system)
    expected = [
        TransitionPathway([{"final": [-0.5, -0.5], "initial": [0.5, -0.5]}]),
        TransitionPathway([{"final": [-0.5, 0.5], "initial": [0.5, 0.5]}]),
    ]
    check_transition_set(tr, expected)


def test_02():
    system = SpinSystem(sites=[{"isotope": "13C"}, {"isotope": "13C"}])
    tr = method1.get_transition_pathways(system)

    expected = [
        TransitionPathway([{"final": [-0.5, -0.5], "initial": [0.5, -0.5]}]),
        TransitionPathway([{"final": [-0.5, 0.5], "initial": [0.5, 0.5]}]),
        TransitionPathway([{"final": [-0.5, -0.5], "initial": [-0.5, 0.5]}]),
        TransitionPathway([{"final": [0.5, -0.5], "initial": [0.5, 0.5]}]),
    ]
    check_transition_set(tr, expected)


def test_03():
    system = SpinSystem(sites=[{"isotope": "13C"}])
    tr = method2.get_transition_pathways(system)

    expected = []
    check_transition_set(tr, expected)


def test_04():
    system = SpinSystem(
        sites=[{"isotope": "13C"}, {"isotope": "13C"}, {"isotope": "14N"}]
    )
    tr = method1.get_transition_pathways(system)

    expected = [
        TransitionPathway([{"final": [-0.5, -0.5, -1], "initial": [0.5, -0.5, -1]}]),
        TransitionPathway([{"final": [-0.5, 0.5, -1], "initial": [0.5, 0.5, -1]}]),
        TransitionPathway([{"final": [-0.5, -0.5, -1], "initial": [-0.5, 0.5, -1]}]),
        TransitionPathway([{"final": [0.5, -0.5, -1], "initial": [0.5, 0.5, -1]}]),
        TransitionPathway([{"final": [-0.5, -0.5, 0], "initial": [0.5, -0.5, 0]}]),
        TransitionPathway([{"final": [-0.5, 0.5, 0], "initial": [0.5, 0.5, 0]}]),
        TransitionPathway([{"final": [-0.5, -0.5, 0], "initial": [-0.5, 0.5, 0]}]),
        TransitionPathway([{"final": [0.5, -0.5, 0], "initial": [0.5, 0.5, 0]}]),
        TransitionPathway([{"final": [-0.5, -0.5, 1], "initial": [0.5, -0.5, 1]}]),
        TransitionPathway([{"final": [-0.5, 0.5, 1], "initial": [0.5, 0.5, 1]}]),
        TransitionPathway([{"final": [-0.5, -0.5, 1], "initial": [-0.5, 0.5, 1]}]),
        TransitionPathway([{"final": [0.5, -0.5, 1], "initial": [0.5, 0.5, 1]}]),
    ]
    check_transition_set(tr, expected)
