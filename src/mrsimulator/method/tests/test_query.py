import numpy as np
from mrsimulator.method.query import MixingEnum
from mrsimulator.method.query import MixingQuery
from mrsimulator.method.query import RotationQuery
from mrsimulator.method.query import SymmetryQuery
from mrsimulator.method.query import TransitionQuery


def test_SymmetryQuery():
    # parse units
    test = {"P": [-1], "D": [0]}
    sym_obj = SymmetryQuery(**test)
    assert sym_obj.name is None
    assert sym_obj.description is None
    assert sym_obj.label is None
    assert sym_obj.P == [-1], "SymmetryQuery P not equal."
    assert sym_obj.D == [0], "SymmetryQuery D not equal."
    assert sym_obj.F is None
    assert sym_obj.transitions is None


def test_TransitionQuery():
    # direct initialization
    test0 = {"ch1": {"P": [-1, -1]}, "ch2": {"P": [-1], "D": [0]}}
    obj1 = TransitionQuery(**test0)
    assert obj1.name is None
    assert obj1.description is None
    assert obj1.label is None
    assert obj1.ch1 == SymmetryQuery(P=[-1, -1]), "TransitionQuery ch1 not equal."
    assert obj1.ch2 == SymmetryQuery(P=[-1], D=[0]), "TransitionQuery ch2 not equal."
    assert obj1.ch3 is None, "TransitionQuery ch3 not equal."

    test1 = {"ch2": {"P": [-2], "D": [2]}}
    obj2 = TransitionQuery(**test1)
    assert obj2.name is None
    assert obj2.description is None
    assert obj2.label is None
    assert obj2.ch1 == SymmetryQuery(P=[0]), "TransitionQuery ch1 not equal."
    assert obj2.ch2 == SymmetryQuery(P=[-2], D=[2]), "TransitionQuery ch2 not equal."
    assert obj2.ch3 is None, "TransitionQuery ch3 not equal."

    test2 = {"ch1": None, "ch2": {"P": [-1, 1], "D": [2]}}
    obj3 = TransitionQuery(**test2)
    assert obj3.name is None
    assert obj3.description is None
    assert obj3.label is None
    assert obj3.ch1 is None, "TransitionQuery ch1 not equal."
    assert obj3.ch2 == SymmetryQuery(P=[-1, 1], D=[2]), "TransitionQuery ch2 not equal."
    assert obj3.ch3 is None, "TransitionQuery ch3 not equal."


def test_RotationQuery():
    # parse units
    test = {"angle": "90 deg", "phase": "5 deg"}
    rf_obj = RotationQuery.parse_dict_with_units(test)
    assert rf_obj.name is None
    assert rf_obj.description is None
    assert rf_obj.label is None
    assert rf_obj.angle == np.pi / 2.0, "RotationQuery angle not equal."
    assert rf_obj.phase == 5 * np.pi / 180, "RotationQuery phase not equal."
    # ensure the default value is rad
    assert rf_obj.property_units["angle"] == "rad"
    assert rf_obj.property_units["phase"] == "rad"

    # direct initialization
    test = {"angle": 3.1415, "phase": 1.212}
    rf_obj = RotationQuery(**test)
    assert rf_obj.angle == 3.1415, "RotationQuery angle not equal."
    assert rf_obj.phase == 1.212, "RotationQuery phase not equal."

    # Always serialize angle and phase
    rf_obj = RotationQuery(angle=0, phase=0)
    assert rf_obj.json(units=False) == {"angle": 0, "phase": 0}
    assert rf_obj.json(units=True) == {"angle": "0.0 rad", "phase": "0.0 rad"}


def test_MixingQuery():
    # parse units
    test1 = {
        "ch1": {"angle": "2.12 rad", "phase": "-176 deg"},
        "ch2": {"angle": "0.12 rad", "phase": "176 deg"},
    }
    obj1 = MixingQuery.parse_dict_with_units(test1)
    rf1 = RotationQuery(angle=2.12, phase=-176 * np.pi / 180)
    rf2 = RotationQuery(angle=0.12, phase=176 * np.pi / 180)
    assert obj1.name is None
    assert obj1.description is None
    assert obj1.label is None
    assert obj1.ch1 == rf1
    assert obj1.ch2 == rf2
    assert obj1.ch3 is None
    assert obj1.channels == [rf1, rf2, None]

    # direct initialization
    test2 = {"ch2": {"angle": 1.101, "phase": 1.61}}
    obj2 = MixingQuery(**test2)
    assert obj2.name is None
    assert obj2.description is None
    assert obj2.label is None
    assert obj2.ch1 is None
    assert obj2.ch2 == RotationQuery(angle=1.101, phase=1.61)
    assert obj2.ch3 is None

    # JSON tests
    test3 = {"ch1": {"angle": 1.23, "phase": 1.23}}
    obj3 = MixingQuery(**test3)
    assert obj3.json(units=False) == test3
    assert obj3.json(units=True) == {"ch1": {"angle": "1.23 rad", "phase": "1.23 rad"}}


def test_MixingEnum():
    total_mix = MixingEnum["TotalMixing"]
    assert total_mix.value == "TotalMixing"
    assert total_mix.json(units=True) == "TotalMixing"
    assert total_mix.json(units=False) == "TotalMixing"

    no_mix = MixingEnum["NoMixing"]
    assert no_mix.value == MixingQuery(
        ch1={"angle": 0, "phase": 0},
        ch2={"angle": 0, "phase": 0},
        ch3={"angle": 0, "phase": 0},
    )
    assert no_mix.json(units=False) == dict(
        ch1={"angle": 0, "phase": 0},
        ch2={"angle": 0, "phase": 0},
        ch3={"angle": 0, "phase": 0},
    )
    assert no_mix.json(units=True) == dict(
        ch1={"angle": "0.0 rad", "phase": "0.0 rad"},
        ch2={"angle": "0.0 rad", "phase": "0.0 rad"},
        ch3={"angle": "0.0 rad", "phase": "0.0 rad"},
    )

    # NOTE: This test will need to be updated as more enumerations are added
    assert set(MixingEnum.allowed_enums()) == {"TotalMixing", "NoMixing"}
