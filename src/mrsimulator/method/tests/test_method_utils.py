# -*- coding: utf-8 -*-
import numpy as np
import pytest
from mrsimulator import Site
from mrsimulator import SpinSystem
from mrsimulator.method import Method
from mrsimulator.method import MixingEvent
from mrsimulator.method import SpectralDimension
from mrsimulator.method.utils import _add_two_euler_angles
from mrsimulator.method.utils import _euler_angles_to_angle_phase
from mrsimulator.method.utils import combind_mixing_queries
from mrsimulator.method.utils import mixing_query_connect_map
from mrsimulator.method.utils import nearest_nonmixing_event
from mrsimulator.utils.error import MissingSpectralEventError

# from mrsimulator.method.utils import angle_and_phase_list


__author__ = ["Deepansh J. Srivastava", "Matthew D. Giammar"]
__email__ = ["srivastava.89@osu.edu", "giammar.7@osu.edu"]


def test_add_euler_angles():
    # angle = 2pi/3, phase = 0
    a1 = np.pi / 2
    b1 = np.pi / 3
    g1 = -np.pi / 2
    a2 = np.pi / 2
    b2 = np.pi / 3
    g2 = -np.pi / 2

    assert np.allclose(
        _add_two_euler_angles(a1, b1, g1, a2, b2, g2),
        np.array([np.pi / 2, 2 * np.pi / 3, -np.pi / 2]),
    )

    # angle = 0, phase = undefined
    # Although phase undefined here, we decide to return phase of zero (alpha = pi/2)
    a1 = 0
    b1 = np.pi / 3
    g1 = 0
    a2 = 0
    b2 = -np.pi / 3
    g2 = 0

    assert np.allclose(
        _add_two_euler_angles(a1, b1, g1, a2, b2, g2),
        np.array([np.pi / 2, 0, -np.pi / 2]),
    )

    # angle = 0, phase = undefined
    # Although phase undefined here, we decide to return phase of zero (alpha = pi/2)
    a1 = 0
    b1 = np.pi / 3
    g1 = 0
    a2 = np.pi
    b2 = np.pi / 3
    g2 = -np.pi

    assert np.allclose(
        _add_two_euler_angles(a1, b1, g1, a2, b2, g2),
        np.array([np.pi / 2, 0, -np.pi / 2]),
    )

    # angle = pi, phase = pi / 5
    # When beta = 180 degrees, phase must be analytically derived
    a1 = np.pi / 5
    b1 = np.pi / 2
    g1 = -np.pi / 5
    a2 = np.pi / 5
    b2 = np.pi / 2
    g2 = -np.pi / 5

    assert np.allclose(
        _add_two_euler_angles(a1, b1, g1, a2, b2, g2),
        np.array([np.pi / 5, np.pi, -np.pi / 5]),
    )

    # Out of transverse plane rotation
    a1 = np.pi / 3
    b1 = np.pi / 2
    g1 = -np.pi / 3
    a2 = 0
    b2 = np.pi / 3
    g2 = 0
    assert np.allclose(
        _add_two_euler_angles(a1, b1, g1, a2, b2, g2),
        np.array([1.28976143, 2.01862872, -0.06440383]),
    )


def test_euler_angles_to_angle_phase():
    a = 0
    b = 1
    g = 1
    with pytest.raises(ValueError, match=".*Unable to convert.*"):
        _euler_angles_to_angle_phase(a, b, g)

    a = np.pi
    b = 1
    g = -np.pi
    assert np.allclose(_euler_angles_to_angle_phase(a, b, g), [1, np.pi / 2])


def test_combind_mixing_queries():
    with pytest.raises(ValueError, match=".*List length must be at least 1.*"):
        combind_mixing_queries([])

    queries = [
        {"angle": 1, "phase": 0},
        {"angle": 1, "phase": 0},
        {"angle": 1, "phase": 0},
    ]
    assert np.allclose(combind_mixing_queries(queries), [np.pi / 2, 3, -np.pi / 2])


def test_warnings():
    s = SpinSystem(sites=[Site(isotope="23Na")])
    m = Method(channels=["1H"], spectral_dimensions=[{}])
    assert m.get_transition_pathways(s) == []


ME = "MixingEvent"
SE = "SpectralEvent"
CE = "ConstantDurationEvent"


def test_nearest_mixing_query():
    events = [SE, ME, SE, SE, CE, ME, SE, ME, ME, CE]
    assert nearest_nonmixing_event(events, 1) == [0, 2]
    assert nearest_nonmixing_event(events, 5) == [4, 6]
    assert nearest_nonmixing_event(events, 7) == [6, 9]


def test_mixing_query_connect_map():
    MX1 = MixingEvent(query={"ch1": {"angle": 0.12}})
    MX2 = MixingEvent(query={"ch2": {"angle": 1.12}})
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
        {"mixing_query_list": [MX1.query], "near_index": [0, 1]},
        {"mixing_query_list": [MX2.query], "near_index": [2, 3]},
        {"mixing_query_list": [MX1.query, MX2.query], "near_index": [3, 4]},
        {"mixing_query_list": [MX2.query], "near_index": [4, 5]},
    ]

    error = "SpectralDimension requires at least one SpectralEvent"
    with pytest.raises(MissingSpectralEventError, match=f".*{error}.*"):
        SpectralDimension(events=[MX1, MX2, {"duration": 0.5}, MX2, {"duration": 0.5}])
