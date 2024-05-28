import numpy as np
import pytest
from mrsimulator import Site
from mrsimulator import SpinSystem
from mrsimulator.method import Method
from mrsimulator.method import MixingEventA
from mrsimulator.method import SpectralDimension
from mrsimulator.method.utils import combine_mixing_queries
from mrsimulator.method.utils import mixing_query_connect_map
from mrsimulator.method.utils import nearest_nonmixing_event
from mrsimulator.utils.error import MissingSpectralEventError

# from mrsimulator.method.utils import angle_and_phase_list


__author__ = ["Deepansh J. Srivastava", "Matthew D. Giammar"]
__email__ = ["srivastava.89@osu.edu", "giammar.7@osu.edu"]


def test_combine_mixing_queries():
    with pytest.raises(ValueError, match=".*List length must be at least 1.*"):
        combine_mixing_queries([])

    queries = [
        {"angle": 1, "phase": 0},
        {"angle": 1, "phase": 0},
        {"angle": 1, "phase": 0},
    ]
    assert np.allclose(combine_mixing_queries(queries), [np.pi / 2, 3, -np.pi / 2])


def test_warnings():
    s = SpinSystem(sites=[Site(isotope="23Na")])
    m = Method(channels=["1H"], spectral_dimensions=[{}])
    assert m.get_transition_pathways(s) == []


ME = "MixingEventA"
SE = "SpectralEvent"
CE = "DelayEvent"


def test_nearest_mixing_query():
    events = [SE, ME, SE, SE, CE, ME, SE, ME, ME, CE]
    assert nearest_nonmixing_event(events, 1) == [0, 2]
    assert nearest_nonmixing_event(events, 5) == [4, 6]
    assert nearest_nonmixing_event(events, 7) == [6, 9]


def test_mixing_query_connect_map():
    MX1 = MixingEventA(ch1={"angle": 0.12})
    MX2 = MixingEventA(ch2={"angle": 1.12})
    TOTAL_MX = MixingEventA(ch1="TotalMixing")

    # Use MixingEvents with non-enum queries
    spectral_dimensions = [
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
                {"fraction": 1},  # 6
            ]
        ),
    ]
    res = mixing_query_connect_map(spectral_dimensions)
    assert res == [
        {"mixing_query_list": [MX1.query], "near_index": [0, 1]},
        {"mixing_query_list": [MX2.query], "near_index": [2, 3]},
        {"mixing_query_list": [MX1.query, MX2.query], "near_index": [3, 4]},
        {"mixing_query_list": [MX2.query], "near_index": [4, 5]},
    ]

    # Combination of MixingEvents with dict queries and enum queries (total mixing)
    spectral_dimensions = [
        SpectralDimension(
            events=[
                # Connect all, should return no list
                {"fraction": 0.5},  # 0
                TOTAL_MX,
                {"duration": 0.5},  # 1
                # Just MX1
                {"fraction": 0.5},  # 2
                TOTAL_MX,
                MX1,
                {"duration": 0.5},  # 3
            ]
        ),
        SpectralDimension(
            events=[
                # MX1 and MX2
                {"fraction": 0.5},  # 4
                MX1,
                TOTAL_MX,
                MX2,
                {"duration": 0.5},  # 5
                # MX1, MX2, MX1
                {"fraction": 0.5},  # 6
                MX1,
                TOTAL_MX,
                MX2,
                MX1,
                TOTAL_MX,
                {"duration": 0.5},  # 7
            ]
        ),
    ]
    res = mixing_query_connect_map(spectral_dimensions)
    assert res == [
        {"mixing_query_list": [MX1.query], "near_index": [2, 3]},
        {"mixing_query_list": [MX1.query, MX2.query], "near_index": [4, 5]},
        {"mixing_query_list": [MX1.query, MX2.query, MX1.query], "near_index": [6, 7]},
    ]

    error = "SpectralDimension requires at least one SpectralEvent"
    with pytest.raises(MissingSpectralEventError, match=f".*{error}.*"):
        SpectralDimension(events=[MX1, MX2, {"duration": 0.5}, MX2, {"duration": 0.5}])
