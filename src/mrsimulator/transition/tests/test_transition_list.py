import pytest
from mrsimulator.transition import Transition
from mrsimulator.transition.pathway import TransitionList

__author__ = "Deepansh Srivastava"
__email__ = "srivastava.89@osu.edu"


def test_transition_1():
    # set up
    a = {"initial": [0.5, 1.5], "final": [-0.5, 1.5]}
    tran = Transition(**a)
    assert tran.initial == [0.5, 1.5]
    assert tran.final == [-0.5, 1.5]

    # to dict with unit
    assert tran.json() == a

    assert tran.tolist() == [0.5, 1.5, -0.5, 1.5]

    # p and Δm
    assert tran.p == -1
    assert tran.delta_m == -1

    assert str(tran) == "|-0.5, 1.5⟩⟨0.5, 1.5|"


def test_transition_list_1():
    # create a list.
    tran_list = TransitionList()

    # append to a list.
    a = {"initial": [0.5, 1.5], "final": [-0.5, 1.5]}
    tran_list.append(a)
    assert tran_list[0] == Transition(**a)

    b = Transition(**{"initial": [-0.5, -1.5], "final": [-1.5, -1.5]})
    tran_list.append(b)
    assert tran_list[1] == b

    error = (
        "Expecting a Transition object or an equivalent python dict object, instead "
        "found str."
    )
    with pytest.raises(ValueError, match=f".*{error}.*"):
        tran_list.append("test")

    # assign an item at a list index.
    c = {"initial": [-1.5, -1.5], "final": [1.5, -1.5]}
    tran_list[1] = c
    assert tran_list[1] == Transition(**c)

    with pytest.raises(ValueError, match=f".*{error}.*"):
        tran_list[1] = "test"

    # inset an item to a list.
    tran_list.insert(1, b)
    assert tran_list[1] == b
    assert tran_list[2] == Transition(**c)

    # length of the list.
    assert len(tran_list) == 3

    # delete an item from the list.
    del tran_list[1]
    assert len(tran_list) == 2
    assert tran_list[0] == Transition(**a)
    assert tran_list[1] == Transition(**c)

    # equality test
    assert tran_list != "a"
    assert tran_list != TransitionList([a])
    assert tran_list != TransitionList([a, b])
    assert tran_list == TransitionList([a, c])

    # appending b for further filter tests.
    tran_list.insert(1, b)
    assert tran_list == TransitionList([a, b, c])

    # filter
    tran_filter_1 = tran_list.filter()
    assert tran_filter_1 == tran_list

    # test for P
    tran_filter_1 = tran_list.filter(P=[-1, 0])
    assert tran_filter_1 == TransitionList([a, b])

    tran_filter_1 = tran_list.filter(P=[-1, -1])
    assert tran_filter_1 == TransitionList([])

    # test for D
    tran_filter_1 = tran_list.filter(D=[0, 0])
    assert tran_filter_1 == TransitionList([a, c])

    # test for P and D
    tran_filter_1 = tran_list.filter(P=[3, 0], D=[0, 0])
    assert tran_filter_1 == TransitionList([c])
