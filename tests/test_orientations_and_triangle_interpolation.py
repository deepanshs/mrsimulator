# -*- coding: utf-8 -*-
"""Test for c functions."""
# import matplotlib.pyplot as plt
import mrsimulator.tests.tests as clib
import numpy as np

from .python_test_for_c_code.orientation import cosine_of_polar_angles_and_amplitudes
from .python_test_for_c_code.orientation import triangle_interpolation1D


# def two_d_plot_helper(amp, lst1, lst2, title=""):
#     _, ax = plt.subplots(1, 2)
#     ax[0].imshow(amp.real, origin="lower", cmap="gray", aspect="auto")
#     ax[1].imshow(amp.imag, origin="lower", cmap="gray", aspect="auto")
#     ax[0].set_title(f"{title} (real)")
#     ax[1].set_title(f"{title} (imag)")
#     x_t, y_t = np.append(lst2, lst2[0]) - 0.5, np.append(lst1, lst1[0]) - 0.5
#     for i in range(2):
#         ax[i].plot(x_t, y_t, "r", label="vertex")
#         ax[i].scatter(lst2 - 0.5, lst1 - 0.5, s=20, color="yellow", label="vertex")
#     plt.legend()
#     plt.show()


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
        [0.2, 30, 50],
        [-20.2, 0.2, 50],
        [-20.2, 0.7, 50],
        [-20.2, 1.2, 50],
        [50.1, 50.4, 50.9],
        [82.3, 100.5, 200],
        [102, 103, 104],
    ]
    for list_ in f_list:
        list_ = np.asarray(list_)
        amp_py = np.zeros(100)
        triangle_interpolation1D(list_, amp_py)

        amp_c = np.zeros(2 * 100)
        clib.triangle_interpolation1D(list_, amp_c)
        amp_c = amp_c[::2] + 1j * amp_c[1::2]

        # plt.plot(amp_py, "b", label="py")
        # plt.plot(amp_c.real, "r--", label="c")
        # plt.plot(list_ - 0.5, [0, 2 / (list_[2] - list_[0]), 0], "k*--", label="x")
        # plt.legend()
        # plt.show()

        assert np.allclose(amp_py, amp_c.real, atol=1e-15)


def test_delta_interpolation():
    f_list = [
        5.5,
        6.0,
        6.2,
        6.4,
        6.6,
        6.8,
        7.0,
        7.1,
        7.5,
        7.9,
        8.0,
        -0.1,
        -0.9,
        -1.0,
        -5.2,
        9.23,
        10.2,
    ]

    for item in f_list:
        list_ = np.asarray([item + 0.5] * 3)

        # should be
        amp_ = np.zeros(10)
        x1 = int(np.floor(item))
        if x1 >= -1 and x1 <= 8:
            amp_[x1 + 1] = item - x1
        if x1 >= 0 and x1 <= 9:
            amp_[x1] = 1 - (item - x1)

        # from delta interpolarion
        amp_c = np.zeros(2 * 10)
        clib.triangle_interpolation1D(list_, amp_c)
        amp_c = amp_c[::2] + 1j * amp_c[1::2]

        # plt.bar(np.arange(10), amp_c.real, width=1)
        # plt.scatter(item, 1.0)
        # plt.show()

        assert np.allclose(amp_, amp_c.real, atol=1e-15)


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
        [[1.5, 9.0, 9.4], [1.5, 9.0, 9.4]],
        [[1.5, 6.0, 9.4], [1.5, 6.0, 9.4]],
        [[1.5, 6.9, 12.1], [1.5, 5.32, 9.0]],
        [[1.5, 6.9, 12.1], [1.5, 4.8, 9.0]],
        [[1.5, 1.5, 1.5], [1.5, 6.0, 9.4]],
    ]
    for list_ in f_list:
        lst1, lst2 = np.asarray(list_)
        amp1 = np.zeros((20, 2 * 20), dtype=np.float64)
        amp2 = np.zeros(2 * 20, dtype=np.float64)
        amp3 = np.zeros(2 * 20, dtype=np.float64)

        clib.triangle_interpolation2D(lst1, lst2, amp1)
        clib.triangle_interpolation1D(lst1, amp2)
        clib.triangle_interpolation1D(lst2, amp3)

        amp1 = amp1[:, ::2] + 1j * amp1[:, 1::2]
        amp2 = amp2[::2] + 1j * amp2[1::2]
        amp3 = amp3[::2] + 1j * amp3[1::2]

        # two_d_plot_helper(amp1, lst1, lst2, title="Complete triangles")

        assert np.allclose(amp2.real, amp1.real.sum(axis=1), atol=1e-15)
        assert np.allclose(amp3.real, amp1.real.sum(axis=0), atol=1e-15)

    # triangles with one or more vertices outside the 2D grids (top - down)
    f_list2 = [
        [[5.5, -8.0, 17.4], [4.5, 12.0, 17.3]],  # one point down
        [[5.5, 8.0, 27.4], [4.5, 12.0, 17.3]],  # one point up
        [[25.5, -8.0, 17.4], [4.5, 12.0, 17.3]],  # one point up and down each
        [[25.5, -8.0, 27.4], [4.5, 12.0, 17.3]],  # two points up and one point down
        [[-12.5, 28.0, -12.4], [4.5, 12.0, 17.3]],  # two points down and one point up
        [[17.5, -8.0, -27.4], [4.5, 12.0, 3.0]],  # two points down
        [[25.5, 8.0, 27.4], [4.5, 12.0, 18.3]],  # two points up
        [[-12.5, -8.0, -12.4], [4.5, 12.0, 17.3]],  # all out
    ]
    for list_ in f_list2:
        lst1, lst2 = np.asarray(list_)
        amp1 = np.zeros((20, 2 * 20), dtype=np.float64)
        amp2 = np.zeros(2 * 20, dtype=np.float64)

        clib.triangle_interpolation2D(lst1, lst2, amp1)
        clib.triangle_interpolation1D(lst1, amp2)

        amp1 = amp1[:, ::2] + 1j * amp1[:, 1::2]
        amp2 = amp2[::2] + 1j * amp2[1::2]

        # two_d_plot_helper(amp1, lst1, lst2, title="One point clipped triangles")

        assert np.allclose(amp2.real, amp1.real.sum(axis=1), atol=1e-15)

    # triangles with one or more vertices outside the 2D grids (left - right)
    f_list2 = [
        [[5.5, 8.0, 17.4], [21.5, 12.0, 17.3]],  # one point right
        [[5.5, 18.0, 17.4], [21.5, 27.0, 17.3]],  # two point right
        # [[5.5, 8.0, 17.4], [-14.5, 12.0, 17.3]],  # one point left
        # [[25.5, -8.0, 17.4], [4.5, 12.0, 17.3]],  # one point up and down each
        # [[25.5, -8.0, 27.4], [4.5, 12.0, 17.3]],  # two points up and one point down
        # [[-12.5, 28.0, -12.4], [4.5, 12.0, 17.3]],  # two points down and one point up
        # [[17.5, -8.0, -27.4], [4.5, 12.0, 3.0]],  # two points down
        # [[25.5, 8.0, 27.4], [4.5, 12.0, 18.3]],   # two points up
        # [[-12.5, -8.0, -12.4], [4.5, 12.0, 17.3]],  # all out
    ]
    for list_ in f_list2:
        lst1, lst2 = np.asarray(list_)
        amp1 = np.zeros((20, 2 * 20), dtype=np.float64)
        amp2 = np.zeros(2 * 20, dtype=np.float64)

        clib.triangle_interpolation2D(lst1, lst2, amp1)
        clib.triangle_interpolation1D(lst2, amp2)

        amp1 = amp1[:, ::2] + 1j * amp1[:, 1::2]
        amp2 = amp2[::2] + 1j * amp2[1::2]

        # two_d_plot_helper(amp1, lst1, lst2, title="Two points clipped triangles")

        assert np.allclose(amp2.real, amp1.real.sum(axis=0), atol=1e-15)

    # triangles with one or more vertices outside a grid voxel
    f_list2 = [
        [[15.0, 15.5, 15.9], [12.0, 12.5, 12.9]],  # all in
        [[15.0, 15.5, 16.2], [12.0, 12.5, 12.9]],  # one out
        [[15.0, 15.5, 16.2], [12.0, 12.5, 13.9]],  # one out on each dimension
        [[15.0, 15.5, 17.2], [12.0, 12.5, 11.9]],  # one out on each dimension
    ]
    for list_ in f_list2:
        lst1, lst2 = np.asarray(list_)
        amp1 = np.zeros((20, 2 * 20), dtype=np.float64)
        amp2 = np.zeros(2 * 20, dtype=np.float64)

        clib.triangle_interpolation2D(lst1, lst2, amp1)
        clib.triangle_interpolation1D(lst1, amp2)

        amp1 = amp1[:, ::2] + 1j * amp1[:, 1::2]
        amp2 = amp2[::2] + 1j * amp2[1::2]

        # two_d_plot_helper(amp1, lst1, lst2, title="Nearest pixel triangles")

        assert np.allclose(amp2.real, amp1.real.sum(axis=1), atol=1e-15)
