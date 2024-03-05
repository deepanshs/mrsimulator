"""Test for c functions."""
import mrsimulator.tests.tests as clib
import numpy as np


def test_exp_Im_alpha_second_rank():
    n = 100
    cos_alpha = (np.random.rand(n) - 0.5) * 2.0
    exp_im_alpha = clib.get_exp_Im_angle(n, cos_alpha, False)

    arr = np.zeros((4, n), dtype=np.complex128)
    alpha = np.arccos(cos_alpha)
    arr[-1] = np.exp(1j * alpha)
    arr[-2] = np.exp(2j * alpha)

    start = 2 * n
    assert np.allclose(exp_im_alpha[start:], arr[2:].ravel())


def test_exp_Im_alpha_fourth_rank():
    n = 100
    cos_alpha = (np.random.rand(n) - 0.5) * 2.0
    exp_im_alpha = clib.get_exp_Im_angle(n, cos_alpha, True)

    arr = np.zeros((4, n), dtype=np.complex128)
    alpha = np.arccos(cos_alpha)
    arr[-1] = np.exp(1j * alpha)
    arr[-2] = np.exp(2j * alpha)
    arr[-3] = np.exp(3j * alpha)
    arr[-4] = np.exp(4j * alpha)

    assert np.allclose(exp_im_alpha, arr.ravel())
