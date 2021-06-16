# -*- coding: utf-8 -*-
import numpy as np
from mrsimulator.base_model import transition_connect_factor


A = [0.5, -0.5, -0.5, 0.5]
B = [0.5, 0.5, -0.5, -0.5]
C = [0.5, -0.5, 0.5, -0.5]
D = [1.5, -1.5, 0.5, -0.5]


def test_connection_factor_01():
    """A1* -> A2  |0.5, 0.5 > < -0.5 0.5| --> |-0.5, -0.5> < 0.5, -0.5|"""
    phase = np.asarray([0, 1 / 4, 1 / 2, 1]) * np.pi
    for phi in phase:
        assert np.allclose(transition_connect_factor(0.5, *B, np.pi, phi), 1.0)

    result = np.asarray([1, -1j, -1, 1])
    for phi, res in zip(phase, result):
        assert np.allclose(transition_connect_factor(0.5, *A, np.pi, phi), res)

    for phi, res in zip(phase, result * 0.5):
        assert np.allclose(transition_connect_factor(0.5, *A, np.pi / 2, phi), res)

    for phi, res in zip(phase, result * 0.146447):
        assert np.allclose(transition_connect_factor(0.5, *A, np.pi / 4, phi), res)

    for phi in phase:
        assert np.allclose(transition_connect_factor(0.5, *C, np.pi / 2, phi), 0.5)

    for phi, res in zip(phase, result * 0.375):
        assert np.allclose(transition_connect_factor(1.5, *D, np.pi / 2, phi), res)
