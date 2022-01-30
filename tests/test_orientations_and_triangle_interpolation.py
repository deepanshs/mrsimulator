# -*- coding: utf-8 -*-
"""Test for c functions."""
import matplotlib.pyplot as plt
import mrsimulator.tests.tests as clib
import numpy as np
from mpl_toolkits.axes_grid1 import make_axes_locatable

from .python_test_for_c_code.orientation import cosine_of_polar_angles_and_amplitudes
from .python_test_for_c_code.orientation import triangle_interpolation1D


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
        [0.5, 0.7, 50.1],
        [10.2, 98.6, 99.2],
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
        [99, 100.1, 100.2],
        [98.9, 99.1, 99.5],
        [-0.23, -0.02, 0.1],
        [98.9, 120, 140],
    ]
    for list_ in f_list:
        list_ = np.asarray(list_)
        amp_py = np.zeros(100)
        triangle_interpolation1D(list_, amp_py)

        amp_c = np.zeros(100)
        clib.triangle_interpolation1D(list_, amp_c)

        # plt.plot(amp_py, "b", label="py")
        # plt.plot(amp_c, "r--", label="c")
        # plt.plot(list_ - 0.5, [0, 2 / (list_[2] - list_[0]), 0], "k*--", label="c")
        # plt.legend()
        # plt.show()

        assert np.allclose(amp_py, amp_c, atol=1e-15)


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
        amp_c = np.zeros(10)
        clib.triangle_interpolation1D(list_, amp_c)

        # plt.bar(np.arange(10), amp_c, width=1)
        # plt.scatter([item], [1.0], marker='x', color='k', s=50)
        # plt.xticks(np.arange(15) - 2)
        # plt.grid(axis="x", which='both')
        # plt.show()

        assert np.allclose(amp_, amp_c, atol=1e-15)


def test_triangle_rasterization1():
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
        amp1 = np.zeros((20, 20), dtype=np.float64)
        amp2 = np.zeros(20, dtype=np.float64)
        amp3 = np.zeros(20, dtype=np.float64)

        clib.triangle_interpolation2D(lst1, lst2, amp1)
        clib.triangle_interpolation1D(lst1, amp2)
        clib.triangle_interpolation1D(lst2, amp3)

        # plot_2d_raster(amp1, lst1, lst2, projx=amp3, projy=amp2)

        assert np.allclose(amp2, amp1.sum(axis=1), atol=1e-15)
        assert np.allclose(amp3, amp1.sum(axis=0), atol=1e-15)


def test_triangle_rasterization2():
    # triangles with one or more vertices outside the 2D grids (top - down)
    f_list = [
        [[5.5, -8.0, 17.4], [4.5, 12.0, 17.3]],  # 1 down
        [[5.5, 8.0, 27.4], [4.5, 12.0, 17.3]],  # 1 up
        [[25.5, -8.0, 17.4], [4.5, 12.0, 17.3]],  # 1 up and down each
        [[25.5, -8.0, 27.4], [4.5, 12.0, 17.3]],  # 2 up and 1 down
        [[-12.5, 28.0, -12.4], [4.5, 12.0, 17.3]],  # 2 down and 1 up
        [[17.5, -8.0, -27.4], [4.5, 12.0, 3.0]],  # 2 down
        [[25.5, 8.0, 27.4], [4.5, 12.0, 18.3]],  # 2 up
        [[-12.5, -8.0, -12.4], [4.5, 12.0, 17.3]],  # all out
    ]
    for list_ in f_list:
        lst1, lst2 = np.asarray(list_)
        amp1 = np.zeros((20, 20), dtype=np.float64)
        amp2 = np.zeros(20, dtype=np.float64)

        clib.triangle_interpolation2D(lst1, lst2, amp1)
        clib.triangle_interpolation1D(lst1, amp2)

        # plot_2d_raster(amp1, lst1, lst2, projx=None, projy=amp2)

        assert np.allclose(amp2, amp1.sum(axis=1), atol=1e-15)


def test_triangle_rasterization3():
    # triangles with one or more vertices outside the 2D grids (left - right)
    f_list = [
        [[5.5, 8.0, 17.4], [21.5, 12.0, 17.3], "v"],  # 1 right
        [[5.5, 18.0, 17.4], [21.5, 27.0, 17.3], "v"],  # 2 right
        [[5.5, 8.0, 17.4], [-14.5, 12.0, 17.3], "v"],  # 1 left
        [[5.5, 8.0, 17.4], [-14.5, 12.0, -17.3], "v"],  # 2 left
        # [[25.5, -8.0, 27.4], [4.5, 12.0, 17.3], 'h'],  # 2 up and 1 down
        # [[-12.5, 28.0, -12.4], [4.5, 12.0, 17.3], 'h'],  # 2 down and 1 up
        # [[17.5, -8.0, -27.4], [4.5, 12.0, 3.0], 'h'],  # 2 down
        # [[25.5, 8.0, 27.4], [4.5, 12.0, 18.3], 'h'],   # 2 up
        # [[-12.5, -8.0, -12.4], [4.5, 12.0, 17.3], ''],  # all out
    ]
    for list_ in f_list:
        gp = list_[-1]
        lst1, lst2 = np.asarray(list_[:-1])
        amp1 = np.zeros((20, 20), dtype=np.float64)
        amp_x = np.zeros(20, dtype=np.float64)
        amp_y = np.zeros(20, dtype=np.float64)

        clib.triangle_interpolation2D(lst1, lst2, amp1)
        clib.triangle_interpolation1D(lst2, amp_x)
        clib.triangle_interpolation1D(lst1, amp_y)

        if gp == "v":
            pass
            # plot_2d_raster(amp1, lst1, lst2, projx=amp_x, projy=None)

        if gp == "h":
            # plot_2d_raster(amp1, lst1, lst2, projx=None, projy=amp_y)
            assert np.allclose(amp_y, amp1.sum(axis=1), atol=1e-15)

        if gp == "":
            # plot_2d_raster(amp1, lst1, lst2, projx=None, projy=None)
            assert np.allclose(np.zeros(20), amp1.sum(axis=0), atol=1e-15)
            assert np.allclose(np.zeros(20), amp1.sum(axis=1), atol=1e-15)


def test_triangle_rasterization4():
    # triangles with one or more vertices outside a grid voxel
    f_list = [
        [[15.0, 15.5, 15.9], [12.0, 12.5, 12.9]],  # all in
        [[15.0, 15.5, 16.2], [12.0, 12.5, 12.9]],  # one out
        [[15.0, 15.5, 16.2], [12.0, 12.5, 13.9]],  # one out on each dimension
        [[15.0, 15.5, 17.2], [12.0, 12.5, 11.9]],  # one out on each dimension
    ]
    for list_ in f_list:
        lst1, lst2 = np.asarray(list_)
        amp1 = np.zeros((20, 20), dtype=np.float64)
        amp2 = np.zeros(20, dtype=np.float64)
        amp3 = np.zeros(20, dtype=np.float64)

        clib.triangle_interpolation2D(lst1, lst2, amp1)
        clib.triangle_interpolation1D(lst1, amp2)
        clib.triangle_interpolation1D(lst2, amp3)

        # plot_2d_raster(amp1, lst1, lst2, projx=amp3, projy=amp2)

        assert np.allclose(amp2, amp1.sum(axis=1), atol=1e-15)


def plot_2d_raster(data2d, pts1, pts2, projx=None, projy=None):
    _, ax = plt.subplots()

    ax.imshow(data2d, origin="lower", cmap="gray", aspect="auto")
    x_t, y_t = np.append(pts2, pts2[0]) - 0.5, np.append(pts1, pts1[0]) - 0.5
    ax.plot(x_t, y_t, "r", label="vertex")
    ax.scatter(pts2 - 0.5, pts1 - 0.5, s=40, color="r", edgecolors="k", label="vertex")

    divider = make_axes_locatable(ax)
    axh = divider.append_axes("top", 1.1, pad=0.1, sharex=ax)
    axv = divider.append_axes("right", 1.1, pad=0.1, sharey=ax)

    plt.setp(axh.get_xticklabels() + axv.get_yticklabels(), visible=False)

    kwargs = dict(linewidth=1)
    if projx is not None:
        axh.plot(np.arange(20), projx, "k--", label="1D")
        axh.plot(np.arange(20), data2d.sum(axis=0), "r", label="sum", **kwargs)

    if projy is not None:
        axv.plot(projy, np.arange(20), "k--", label="1D")
        axv.plot(data2d.sum(axis=1), np.arange(20), "r", label="sum", **kwargs)
    plt.legend()
    plt.tight_layout()
    plt.show()
