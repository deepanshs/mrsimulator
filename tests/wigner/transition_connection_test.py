import numpy as np
from mrsimulator.base_model import transition_connect_factor


A = [0.5, -0.5, -0.5, 0.5]
B = [0.5, 0.5, -0.5, -0.5]
C = [0.5, -0.5, 0.5, -0.5]
D = [1.5, -1.5, 0.5, -0.5]


def test_connection_factor_01():
    """A1* -> A2  |0.5, 0.5 > < -0.5 0.5| --> |-0.5, -0.5> < 0.5, -0.5|"""
    phase = np.asarray([0, 1 / 4, 1 / 2, 1]) * np.pi
    alpha = np.pi / 2 - phase
    gamma = phase - np.pi / 2
    for a_, g_ in zip(alpha, gamma):
        assert np.allclose(transition_connect_factor(0.5, *B, a_, np.pi, g_), 1.0)
        assert np.allclose(transition_connect_factor(0.5, *C, a_, np.pi / 2, g_), 0.5)

    result = np.asarray([1, -1j, -1, 1])
    for a_, g_, res in zip(alpha, gamma, result):
        assert np.allclose(transition_connect_factor(0.5, *A, a_, np.pi, g_), res)

    for a_, g_, res in zip(alpha, gamma, result * 0.5):
        assert np.allclose(transition_connect_factor(0.5, *A, a_, np.pi / 2, g_), res)

    for a_, g_, res in zip(alpha, gamma, result * 0.146447):
        assert np.allclose(transition_connect_factor(0.5, *A, a_, np.pi / 4, g_), res)

    for a_, g_, res in zip(alpha, gamma, result * 0.375):
        assert np.allclose(transition_connect_factor(1.5, *D, a_, np.pi / 2, g_), res)


# def test_connection_factor_02():
#     """Test transition when phases are out of the XY plane (i.e. alpha != -gamma)"""
#     pass
