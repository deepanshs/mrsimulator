# -*- coding: utf-8 -*-
import pytest
from mrsimulator import Site
from mrsimulator import SpinSystem
from mrsimulator.method import MixingEvent
from mrsimulator.method import SpectralDimension
from mrsimulator.method.query import MixingQuery
from mrsimulator.method.utils import mixing_query_connect_map
from mrsimulator.method.utils import nearest_nonmixing_event
from mrsimulator.method.utils import tip_angle_and_phase_list
from mrsimulator.methods import Method1D
from mrsimulator.utils.error import MissingSpectralEventError


__author__ = "Deepansh Srivastava"
__email__ = "srivastava.89@osu.edu"


def test_warnings():
    s = SpinSystem(sites=[Site(isotope="23Na")])
    m = Method1D(channels=["1H"])
    assert m.get_transition_pathways(s) == []


ME = "MixingEvent"
SE = "SpectralEvent"
CE = "ConstantDurationEvent"


def test_nearest_mixing_query():
    events = [SE, ME, SE, SE, CE, ME, SE, ME, ME, CE]
    assert nearest_nonmixing_event(events, 1) == [0, 2]
    assert nearest_nonmixing_event(events, 5) == [4, 6]
    assert nearest_nonmixing_event(events, 7) == [6, 9]


def test_tip_angle_and_phase_list():
    symbol = ["29Si", "17O", "29Si", "29Si", "27Al"]
    channels = ["29Si", "27Al"]
    mixing_query = MixingQuery(ch1={"tip_angle": 1.55, "phase": 0.5})
    lst = tip_angle_and_phase_list(symbol, channels, mixing_query)
    assert lst == ([1.55, 0, 1.55, 1.55, 0], [0.5, 0, 0.5, 0.5, 0])

    mixing_query = MixingQuery(
        ch1={"tip_angle": 1.55, "phase": 0.5},
        ch2={"tip_angle": 3.1415, "phase": 1.570},
    )
    lst = tip_angle_and_phase_list(symbol, channels, mixing_query)
    assert lst == ([1.55, 0, 1.55, 1.55, 3.1415], [0.5, 0, 0.5, 0.5, 1.57])


def test_mixing_query_connect_map():
    MX1 = MixingEvent(mixing_query={"ch1": {"tip_angle": 0.12}})
    MX2 = MixingEvent(mixing_query={"ch2": {"tip_angle": 1.12}})
    spectral_dimmensions = [
        SpectralDimension(
            events=[
                {"fraction": 0.5},  # 0
                MX1,
                {"duration": 0.5},  # 1
                {"fraction": 0.5},  # 2
                MX2,
                {"duration": 0.5},  # 3
            ]
        ),
        SpectralDimension(
            events=[
                MX1,
                MX2,
                {"duration": 0.5},  # 4
                MX2,
                {"duration": 0.5},  # 5
                {"fraction": 1},  # 2
            ]
        ),
    ]
    res = mixing_query_connect_map(spectral_dimmensions)
    assert res == [
        {"mixing_query": MX1.mixing_query, "near_index": [0, 1]},
        {"mixing_query": MX2.mixing_query, "near_index": [2, 3]},
        {"mixing_query": MX1.mixing_query, "near_index": [3, 4]},
        {"mixing_query": MX2.mixing_query, "near_index": [3, 4]},
        {"mixing_query": MX2.mixing_query, "near_index": [4, 5]},
    ]

    error = "SpectralDimension requires at least one SpectralEvent"
    with pytest.raises(MissingSpectralEventError, match=f".*{error}.*"):
        SpectralDimension(events=[MX1, MX2, {"duration": 0.5}, MX2, {"duration": 0.5}])
