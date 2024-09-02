"""Test for c functions."""
from copy import deepcopy
from os import mkdir
from os import path

import matplotlib.pyplot as plt
import mrsimulator.tests.tests as clib
import numpy as np
import pytest
from mpl_toolkits.axes_grid1 import make_axes_locatable

from .python_test_for_c_code.orientation import cosine_of_polar_angles_and_amplitudes
from .python_test_for_c_code.orientation import triangle_interpolation1D

SCALE = [1, 10]
__GENERATE_REPORT__ = False

import matplotlib
from matplotlib.backends.backend_pdf import PdfPages

font = {"size": 9}
matplotlib.rc("font", **font)


@pytest.fixture(scope="module")
def report():
    if __GENERATE_REPORT__:  # pragma: no cover
        pdf = PdfPages("reports/interpolation_report.pdf")
    else:
        pdf = None
    return pdf


@pytest.fixture(autouse=True, scope="module")
def report_pdf(report):
    """Generate report"""
    yield
    if __GENERATE_REPORT__:  # pragma: no cover
        report.close()


def plot_x_proj(axh, size, proj_x, task, label, **kwargs):
    """X projection plot"""
    if task == "Low":
        axh.stairs(proj_x, edges=np.arange(size + 1) - 0.5, label=label, **kwargs)
    else:
        axh.plot(np.arange(size), proj_x, label=label, **kwargs)


def plot_y_proj(axv, size, proj_y, task, label, **kwargs):
    """Y projection plot"""
    if task == "Low":
        edges = np.arange(size + 1) - 0.5
        axv.stairs(proj_y, edges=edges, label=label, orientation="horizontal", **kwargs)
    else:
        axv.plot(proj_y, np.arange(size), label=label, **kwargs)


def plot_2d_raster(rep, title=None, report=None):
    """Test report plot"""
    _, axs = plt.subplots(1, 2, figsize=(9, 4))
    len_rep = len(rep)

    for plot_res in range(len_rep):
        data2d, pts1, pts2, proj_x, proj_y = rep[plot_res]
        ax = axs[plot_res]
        pts1 = np.round(pts1 - 0.5, 2)
        pts2 = np.round(pts2 - 0.5, 2)

        ax.imshow(data2d, origin="lower", cmap="gist_ncar_r", aspect="auto")
        x_t, y_t = np.append(pts2, pts2[0]), np.append(pts1, pts1[0])
        ax.plot(x_t, y_t, "r", label="vertex")
        ax.scatter(pts2, pts1, s=40, c="r", edgecolors="k", label="vertex")
        ax.set_xlabel("dimension - 1")
        ax.set_ylabel("dimension - 2")
        ax.grid(linewidth=0.5, alpha=0.5)

        task = "Low" if SCALE[plot_res] == 1 else "High"

        divider = make_axes_locatable(ax)
        axh = divider.append_axes("top", 0.9, pad=0.1, sharex=ax)
        axv = divider.append_axes("right", 1, pad=0.1, sharey=ax)

        plt.setp(axh.get_xticklabels() + axv.get_yticklabels(), visible=False)

        kwargs_proj = dict(linewidth=0.5, color="r")
        kwargs_asg = dict(linewidth=1, color="k", linestyle="--")
        size = data2d.shape[1]
        if proj_x is not None:
            plot_x_proj(axh, size, proj_x, task, label="1D Tent", **kwargs_asg)
        plot_x_proj(
            axh, size, data2d.sum(axis=0), task, label="Proj - 1", **kwargs_proj
        )
        axh.legend()
        axh.set_title(f"{task} resolution grid test\n{title}", fontsize=9)

        size = data2d.shape[0]
        if proj_y is not None:
            plot_y_proj(axv, size, proj_y, task, label="1D Tent", **kwargs_asg)
        plot_y_proj(
            axv, size, data2d.sum(axis=1), task, label="Proj - 2", **kwargs_proj
        )
        axv.legend()

        tri_info = "Triangle\nCoordinates\n"
        tri_info += f"\nA = ({pts2[0]}, {pts1[0]})"
        tri_info += f"\nB = ({pts2[1]}, {pts1[1]})"
        tri_info += f"\nC = ({pts2[2]}, {pts1[2]})\n"
        axv.set_title(tri_info, fontsize=9)
    plt.tight_layout()
    if report is not None:
        report.savefig(dpi=150)
    plt.close()


def test_pdf():
    global __GENERATE_REPORT__
    temp_status = deepcopy(__GENERATE_REPORT__)
    __GENERATE_REPORT__ = True

    is_present = path.isdir("_temp")
    if not is_present:  # pragma: no cover
        mkdir("_temp")
    filename = "_temp/interpolation_report_scrap.pdf"
    report_file = PdfPages(filename)
    amp2d = np.random.rand(100).reshape(10, 10)
    pts1 = np.array([1, 3, 6])
    pts2 = np.array([5, 2, 9])
    proj_x = amp2d.sum(axis=0)
    proj_y = amp2d.sum(axis=1)
    rep = [amp2d, pts1, pts2, proj_x, proj_y]
    plot_2d_raster([rep], title="All points within FOV", report=report_file)
    report_file.close()
    is_file = path.isfile(filename)
    assert is_file

    __GENERATE_REPORT__ = temp_status


def test_octahedron_averaging_setup():
    nt = 64
    cos_alpha_py, cos_beta_py, amp_py = cosine_of_polar_angles_and_amplitudes(nt)
    exp_I_alpha_c, exp_I_beta_c, amp_c = clib.cosine_of_polar_angles_and_amplitudes(nt)

    assert np.allclose(cos_alpha_py, exp_I_alpha_c.real, atol=1e-15)
    assert np.allclose(cos_beta_py, exp_I_beta_c.real, atol=1e-15)
    assert np.allclose(amp_py, amp_c, atol=1e-15)


# def get_height(pts, amp=1):
#     pts = np.sort(pts)
#     return amp * 2.0 / (pts[2] - pts[0])


def test_triangle_interpolation():
    f_list = [
        [-0.91, -4.14, -16.70],
        [-0.91, 2.16, -16.70],
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
        [0.6, 1.6, 45.5],
        [18.35, 18.71, 13.77],
        [11.28, 10.97, 6.36],
        [60, 99.2, 100.2],
        [-0.9, -0.1, 10],
    ]
    for scl in SCALE:
        for i, list_ in enumerate(f_list):
            list_ = np.sort(list_) * scl
            amp_py = np.zeros(100 * scl)
            triangle_interpolation1D(list_, amp_py)

            amp_c = np.zeros(2 * 100 * scl)
            clib.triangle_interpolation1D(list_, amp_c)
            amp_c = amp_c[::2] + 1j * amp_c[1::2]

            # x = np.arange(100 * scl) / scl
            # y = get_height(list_)
            # plt.plot(x, amp_py, "b", label="py")
            # plt.plot(x, amp_c.real, "r--", label="c")
            # plt.plot((list_ - 0.5) / scl, [0, y, 0], "k*--", label="x")
            # plt.title("1D interpolation, span(0, 100)")
            # plt.legend()
            # plt.tight_layout()
            # plt.show()
            # plt.savefig(f"figs/fig_1D_{i}_{scl}.pdf")
            # plt.figure().clear()

            assert np.allclose(amp_py, amp_c.real, atol=1e-15)


def test_delta_interpolation_linear():
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
        -1.8,
        9.23,
        10.2,
    ]
    f_list += list(np.arange(80) / 40 + 5)

    for i, item in enumerate(f_list):
        list_ = np.asarray([item + 0.5] * 3)

        # should be
        amp_ = np.zeros(10)
        x1 = int(np.floor(item))
        if x1 >= -1 and x1 <= 8:
            amp_[x1 + 1] = item - x1
        if x1 >= 0 and x1 <= 9:
            amp_[x1] = 1 - (item - x1)

        # from delta interpolation
        amp_c = np.zeros(2 * 10)
        clib.triangle_interpolation1D(list_, amp_c, type="linear")
        amp_c = amp_c[::2] + 1j * amp_c[1::2]

        # plt.bar(np.arange(10), amp_c.real, width=1)
        # plt.scatter([item], [1.0], marker="x", color="k", s=50)
        # plt.xticks(np.arange(17) - 3)
        # plt.grid(axis="x", which="both")
        # plt.title("1D delta interpolation (linear), span(0, 10)")
        # plt.savefig(f"figs/fig_1D_delta_linear{i}.pdf")
        # plt.figure().clear()

        assert np.allclose(amp_, amp_c.real, atol=1e-15)


def test_gaussian_interpolation():
    sigma = 1.0 / 4.0
    div = 2 * sigma**2
    scale = sigma * np.sqrt(2 * np.pi)
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
        9.23,
        10.2,
        -0.1,
        -0.9,
        -1.0,
        -1.8,
    ]
    f_list += list(np.arange(20) / 10 + 5)

    x = np.arange(20) - 5
    for i, item in enumerate(f_list):
        list_ = np.asarray([item + 0.5] * 3)

        # should be
        gauss = np.exp(-((x - item) ** 2) / div) / scale
        gauss /= gauss.sum()
        gauss = gauss[5:-5]

        # from delta interpolation
        amp_c = np.zeros(20)
        clib.triangle_interpolation1D(list_, amp_c, type="gaussian")
        amp_c = amp_c[::2] + 1j * amp_c[1::2]

        # plt.plot(np.arange(10), gauss, label="gauss")
        # plt.plot(np.arange(10), amp_c, label="calc")
        # plt.legend()
        # plt.show()

        # plt.bar(np.arange(10), amp_c.real, width=1)
        # plt.scatter([item], [1.0], marker="x", color="k", s=50)
        # # plt.xticks(np.arange(17) - 3)
        # plt.grid(axis="x", which="both")
        # plt.title("1D delta interpolation (gaussian), span(0, 10)")
        # plt.savefig(f"figs/fig_1D_delta_gaussian{i}.pdf")
        # plt.figure().clear()

        np.testing.assert_almost_equal(gauss, amp_c, decimal=3)


def test_triangle_rasterization1(report):
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
    for i, list_ in enumerate(f_list):
        rep = []
        for scl in SCALE:
            amp1, amp2, amp3, lst1, lst2 = get_amps_from_interpolation(list_, scl)

            assert np.allclose(amp2, amp1.sum(axis=1), atol=1e-15)
            assert np.allclose(amp3, amp1.sum(axis=0), atol=1e-15)

            rep.append([amp1, lst1, lst2, amp3, amp2])

        plot_2d_raster(rep, title="All points within FOV", report=report)


def test_triangle_rasterization2(report):
    # triangles with one or more vertices outside the 2D grids (top - down)
    f_list = [
        [[5.5, -8.0, 17.4], [4.5, 12.0, 17.3]],  # 1 down
        [[5.5, 8.0, 27.4], [4.5, 12.0, 17.3]],  # 1 up
        [[25.5, -8.0, 17.4], [4.5, 12.0, 17.3]],  # 1 up and down each
        [[17.5, -8.0, -27.4], [4.5, 12.0, 3.0]],  # 2 down
        [[25.5, 8.0, 27.4], [4.5, 12.0, 18.3]],  # 2 up
    ]
    for i, list_ in enumerate(f_list):
        rep = []
        for scl in SCALE:
            lst1, lst2 = np.asarray(list_) * scl
            amp1 = np.zeros((20 * scl, 2 * 20 * scl), dtype=np.float64)
            amp2 = np.zeros(2 * 20 * scl, dtype=np.float64)

            clib.triangle_interpolation2D(lst1, lst2, amp1)
            clib.triangle_interpolation1D(lst1, amp2)

            amp1 = amp1[:, ::2]  # + 1j * amp1[:, 1::2]
            amp2 = amp2[::2]  # + 1j * amp2[1::2]

            assert np.allclose(amp2, amp1.sum(axis=1), atol=1e-15)

            rep.append([amp1, lst1, lst2, None, amp2])

        plot_2d_raster(rep, title="One or more points outside of FOV", report=report)


def test_triangle_rasterization3(report):
    # triangles with one or more vertices outside the 2D grids (left - right)
    f_list = [
        [[5.5, 8.0, 17.4], [21.5, 12.0, 17.3]],  # 1 right
        [[5.5, 18.0, 17.4], [21.5, 27.0, 17.3]],  # 2 right
        [[5.5, 8.0, 17.4], [-14.5, 12.0, 17.3]],  # 1 left
        [[5.5, 8.0, 17.4], [-14.5, 12.0, -17.3]],  # 2 left
        [[5.5, 8.0, 17.4], [-4.5, 12.0, 27.3]],  # 1 left and 1 right
        [[15.5, 8.0, 17.4], [-4.5, 12.0, 27.3]],  # 1 left and 1 right
    ]
    for i, list_ in enumerate(f_list):
        rep = []
        for scl in SCALE:
            amp1, amp_y, amp_x, lst1, lst2 = get_amps_from_interpolation(list_, scl)
            assert np.allclose(amp_x, amp1.sum(axis=0), atol=1e-2)

            rep.append([amp1, lst1, lst2, amp_x, None])
        plot_2d_raster(rep, title="One or more points outside of FOV", report=report)


def test_triangle_rasterization4(report):
    # all points outside but intersecting the view
    f_list = [
        [[-12.5, 18.0, 32.4], [-4.5, 22.0, 27.3]],  # all out
        [[-12.5, 18.0, -22.4], [-4.5, 22.0, 27.3]],  # all out
        [[-12.5, 18.0, 22.4], [-4.5, 22.0, -27.3]],  # all out
        [[-12.5, -8.0, -12.4], [4.5, 12.0, 17.3]],  # all out down
        [[42.5, 28.0, 32.4], [4.5, 12.0, 17.3]],  # all out up
        [[25.5, -8.0, 27.4], [4.5, 12.0, 17.3]],  # 2 up and 1 down
        [[-12.5, 28.0, -12.4], [4.5, 12.0, 17.3]],  # 2 down and 1 up
    ]
    for i, list_ in enumerate(f_list):
        rep = []
        for scl in SCALE:
            amp1, amp_y, amp_x, lst1, lst2 = get_amps_from_interpolation(list_, scl)

            rep.append([amp1, lst1, lst2, amp_x, None])
        plot_2d_raster(rep, title="All points outside of FOV", report=report)
        # assert np.allclose(amp_x.real, amp1.real.sum(axis=0), atol=1e-15)


def test_triangle_rasterization5(report):
    # triangles with one or more vertices outside a grid voxel
    f_list = [
        [[15.0, 15.5, 15.9], [12.0, 12.5, 12.9]],  # all in
        [[15.0, 15.5, 16.2], [12.0, 12.5, 12.9]],  # one out
        [[15.0, 15.5, 16.2], [12.0, 12.5, 13.9]],  # one out on each dimension
        [[15.0, 15.5, 17.2], [12.0, 12.5, 11.9]],  # one out on each dimension
        [[14.5, 15.5, 18.2], [11.5, 12.0, 11.2]],  # one out on each dimension
    ]
    for i, list_ in enumerate(f_list):
        rep = []
        for scl in SCALE:
            amp1, amp2, amp3, lst1, lst2 = get_amps_from_interpolation(list_, scl)
            assert np.allclose(amp2, amp1.sum(axis=1), atol=1e-15)

            rep.append([amp1, lst1, lst2, amp3, amp2])
        plot_2d_raster(rep, title="Small trianges spanning 1-2 bins", report=report)


def get_amps_from_interpolation(list_, scl):
    lst1, lst2 = np.asarray(list_) * scl
    amp1 = np.zeros((20 * scl, 2 * 20 * scl), dtype=np.float64)
    amp2 = np.zeros(2 * 20 * scl, dtype=np.float64)
    amp3 = np.zeros(2 * 20 * scl, dtype=np.float64)

    clib.triangle_interpolation2D(lst1, lst2, amp1)
    clib.triangle_interpolation1D(lst1, amp2)
    clib.triangle_interpolation1D(lst2, amp3)

    amp1 = amp1[:, ::2]  # + 1j * amp1[:, 1::2]
    amp2 = amp2[::2]  # + 1j * amp2[1::2]
    amp3 = amp3[::2]  # + 1j * amp3[1::2]
    return amp1, amp2, amp3, lst1, lst2
