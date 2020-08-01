# -*- coding: utf-8 -*-
"""Test for c functions."""
import mrsimulator.tests.tests as clib
import numpy as np

from .python_test_for_c_code.orientation import cosine_of_polar_angles_and_amplitudes
from .python_test_for_c_code.orientation import triangle_interpolation


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


def test_triangle_rasterization():

    # triangles within the 2D grids
    f_list = [
        [[6.0, 2.3, 19.0], [15.0, 2.0, 17.9]],
        [[10.0, 10.9, 1.0], [15.0, 2.0, 17.9]],
        [[10.0, 10.9, 1.0], [15.0, 2.0, 15.9]],
        [[10.0, 10.9, 1.0], [15.0, 2.0, 12.0]],
        [[1.5, 10.9, 1.0], [16.0, 1.2, 2.4]],
        [[1.5, 10.9, 1.6], [15.0, 2.0, 1.0]],
        [[1.5, 2.0, 7.4], [1.5, 2.0, 7.3]],
    ]
    for list_ in f_list:
        lst1, lst2 = np.asarray(list_)
        amp1 = np.zeros((20, 20), dtype=np.float64)
        amp2 = np.zeros(20, dtype=np.float64)
        amp3 = np.zeros(20, dtype=np.float64)

        clib.triangle_interpolation2D(lst1, lst2, amp1)
        clib.triangle_interpolation(lst1, amp2)
        clib.triangle_interpolation(lst2, amp3)

        assert np.allclose(amp2, amp1.sum(axis=1), atol=1e-15)
        assert np.allclose(amp3, amp1.sum(axis=0), atol=1e-15)

    # triangles with one or more vertices outside the 2D grids
    f_list2 = [
        [[5.5, -8.0, 17.4], [4.5, 12.0, 17.3]],
        [[5.5, 8.0, 27.4], [4.5, 12.0, 17.3]],
        [[25.5, -8.0, 17.4], [4.5, 12.0, 17.3]],
        [[25.5, -8.0, 27.4], [4.5, 12.0, 17.3]],
        [[17.5, -8.0, -27.4], [4.5, 12.0, 3.0]],
    ]
    for list_ in f_list2:
        lst1, lst2 = np.asarray(list_)
        amp1 = np.zeros((20, 20), dtype=np.float64)
        amp2 = np.zeros(20, dtype=np.float64)

        clib.triangle_interpolation2D(lst1, lst2, amp1)
        clib.triangle_interpolation(lst1, amp2)

        assert np.allclose(amp2, amp1.sum(axis=1), atol=1e-15)
