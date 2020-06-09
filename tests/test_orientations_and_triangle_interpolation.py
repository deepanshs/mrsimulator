# -*- coding: utf-8 -*-
"""Test for c functions."""
import mrsimulator.tests.tests as clib
import numpy as np
from mrsimulator.python.orientation import cosine_of_polar_angles_and_amplitudes
from mrsimulator.python.orientation import triangle_interpolation


def test_octahedron_averaging_setup():
    nt = 64
    cos_alpha_py, cos_beta_py, amp_py = cosine_of_polar_angles_and_amplitudes(nt)
    exp_I_alpha_c, exp_I_beta_c, amp_c = clib.cosine_of_polar_angles_and_amplitudes(nt)

    assert np.allclose(cos_alpha_py, exp_I_alpha_c.real, atol=1e-15)
    assert np.allclose(cos_beta_py, exp_I_beta_c.real, atol=1e-15)
    assert np.allclose(amp_py, amp_c, atol=1e-15)


def test_triangle_interpolation():
    f_list = [
        [10.2, 80.3, 80.4],
        [80.2, 80.3, 107.4],
        [-200.1, -200.2, -200.3],
        [-200, -150, -600],
        [-100, -40, 10],
        [-20, 10, 50],
        [10, 30, 50],
        [50.1, 50.4, 50.9],
        [82.3, 100.5, 200],
        [102, 103, 104],
    ]
    for list_ in f_list:
        list_ = np.asarray(list_)
        amp_py = np.zeros(100)
        triangle_interpolation(list_, amp_py)

        amp_c = np.zeros(100)
        clib.triangle_interpolation(list_, amp_c)

        assert np.allclose(amp_py, amp_c, atol=1e-15)
