# -*- coding: utf-8 -*-
from mrsimulator.transition.transition_list import TransitionPathway

__author__ = "Deepansh Srivastava"
__email__ = "srivastava.89@osu.edu"


def test_transition_pathway():
    # set up
    a = {"initial": [0.5, 1.5], "final": [-0.5, 1.5]}
    b = {"initial": [0.5, -1.5], "final": [0.5, 1.5]}
    c = {"initial": [0.5, 0.5], "final": [-0.5, -0.5]}
    trans_path = TransitionPathway([a, b, c])

    assert trans_path[0].initial == [0.5, 1.5]
    assert trans_path[0].final == [-0.5, 1.5]

    assert trans_path[1].initial == [0.5, -1.5]
    assert trans_path[1].final == [0.5, 1.5]

    assert trans_path[2].initial == [0.5, 0.5]
    assert trans_path[2].final == [-0.5, -0.5]

    # to dict with unit
    assert trans_path.json() == [a, b, c]

    assert trans_path.tolist() == [
        0.5,
        1.5,
        -0.5,
        1.5,
        0.5,
        -1.5,
        0.5,
        1.5,
        0.5,
        0.5,
        -0.5,
        -0.5,
    ]

    string = " ⟶ ".join(
        ["|-0.5, 1.5⟩⟨0.5, 1.5|", "|0.5, 1.5⟩⟨0.5, -1.5|", "|-0.5, -0.5⟩⟨0.5, 0.5|"]
    )
    assert str(trans_path) == string
