"""Lineshape Test."""
from copy import deepcopy
from os import mkdir
from os import path
from pprint import pformat

import numpy as np
import pytest

from .utils import c_setup
from .utils import c_setup_random_euler_angles

COMMON_PATH = path.join("tests", "spectral_integration_tests")
SIMPSON_TEST_PATH = path.join(COMMON_PATH, "simpson_simulated_lineshapes")
RNMSIM_TEST_PATH = path.join(COMMON_PATH, "rmnsim_lineshapes")
PYTHON_BRUTE_TEST_PATH = path.join(COMMON_PATH, "python_brute_force_lineshapes")
VOLUMES = ["sphere", "hemisphere"]

__GENERATE_REPORT__ = False

import matplotlib
import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages

font = {"size": 9}
matplotlib.rc("font", **font)


@pytest.fixture(scope="module")
def report():
    if __GENERATE_REPORT__:  # pragma: no cover
        pdf = PdfPages("reports/lineshapes_report.pdf")
    else:
        pdf = None
    return pdf


@pytest.fixture(autouse=True, scope="module")
def report_pdf(report):
    """Generate report"""
    yield
    if __GENERATE_REPORT__:  # pragma: no cover
        report.close()


def compile_plots(dim, rep, info, dim2=None, title=None, report=None, label="simpson"):
    """Test report plot"""
    if not __GENERATE_REPORT__:
        return

    if dim2 is None:  # pragma: no cover
        fig, ax = plt.subplots(
            1, 2, figsize=(9, 4), gridspec_kw={"width_ratios": [1, 1]}
        )

        for i, res in enumerate(rep[:1]):
            data_mrsimulator, data_source = res
            ax[i].plot(dim, data_mrsimulator, "k", linewidth=0.75, label="mrsimulator")
            ax[i].plot(dim, data_source, "--r", linewidth=0.75, label=label)
            ax[i].set_xlabel("Frequency / ppm")
            ax[i].legend()

        ax[0].plot(
            dim, data_source - data_mrsimulator, "grey", linewidth=0.75, label="residue"
        )
        ax[0].legend()
        axbig = ax[1]
    else:  # pragma: no cover
        fig, ax = plt.subplots(
            2, 2, sharex=True, figsize=(9, 4), gridspec_kw={"width_ratios": [1, 1]}
        )
        extent = [dim[0], dim[-1], dim2[0], dim2[-1]]
        kwargs = dict(aspect="auto", extent=extent, cmap="gist_ncar_r")
        for i, res in enumerate(rep[:1]):
            data_mrsimulator, data_source = res
            ax[0, 0].imshow(data_mrsimulator, **kwargs)
            ax[0, 0].set_title("mrsimulator")

            ax[1, 0].imshow(data_source, **kwargs)
            ax[1, 0].set_title(label)

            ax[1, 0].set_xlabel("Frequency / ppm")
            ax[0, 0].set_ylabel("Frequency / ppm")
            ax[1, 0].set_ylabel("Frequency / ppm")

        gs = ax[0, 0].get_gridspec()
        for ax in ax[:, -1]:
            ax.remove()
        axbig = fig.add_subplot(gs[:, -1])

    format_kwargs = dict(indent=1, width=100, compact=True, sort_dicts=False)
    str_config = pformat(info["config"], **format_kwargs)
    str_sys = pformat(info["spin_systems"], **format_kwargs)
    str_method = pformat(info["methods"], **format_kwargs)

    str_config = r"$\bf{sim.config}$" + "\n" + str_config + "\n\n"
    str_sys = r"$\bf{sim.spin\_systems}$" + "\n" + str_sys + "\n\n"
    str_method = r"$\bf{sim.methods}$" + "\n" + str_method

    text = str_config + str_sys + str_method
    axbig.set_ylim(0, 1)
    axbig.text(
        0, 1, text, fontsize=7.5, verticalalignment="top", horizontalalignment="left"
    )
    axbig.set_title(f"Test series: {title}")
    axbig.axis("off")
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
    filename = "_temp/lineshapes_report_scrap.pdf"
    report_file = PdfPages(filename)
    dim = np.arange(10)
    res = [np.arange(10), np.arange(10)]
    info = {"config": 1, "spin_systems": 2, "methods": 3}
    compile_plots(dim, [res], info, title="Shielding Sidebands", report=report_file)
    report_file.close()
    is_file = path.isfile(filename)
    assert is_file

    __GENERATE_REPORT__ = temp_status


def check_all_close(res, message, rel_limit):
    """Check if the vectos in res are all close within relative limits"""
    for item in res:
        item0 = item[0] / item[0].sum()
        item1 = item[1] / item[1].sum()
        np.testing.assert_allclose(item0, item1, atol=rel_limit, err_msg=message)


# --------------------------------------------------------------------------- #
# Test against simpson calculations
# --------------------------------------------------------------------------- #


def test_pure_shielding_sideband_simpson(report):
    error_message = (
        "failed to compare shielding sidebands with simpson simulation from file"
    )
    path_ = path.join(SIMPSON_TEST_PATH, "shielding_sidebands")
    total = 8
    for i in range(total):
        filename = path.join(path_, f"test{i:02d}", f"test{i:02d}.json")  # noqa: E231

        res = []
        # euler angle all zero
        data_mrsimulator, data_source, info, dim = c_setup(filename)
        res.append([data_mrsimulator, data_source])

        # random euler angle all zero. Euler angles should not affect the spectrum.
        data_mrsimulator, data_source = c_setup_random_euler_angles(
            filename, "shielding_symmetric"
        )
        res.append([data_mrsimulator, data_source])

        title = f"Shielding Sidebands ({i + 1} / {total})"
        compile_plots(dim, res, info, title=title, report=report)

        message = f"{error_message} test0{i}.json"
        check_all_close(res, message, rel_limit=1e-3)


def test_pure_quadrupolar_sidebands_simpson(report):
    error_message = (
        "failed to compare quadrupolar sidebands with simpson simulation from file"
    )
    path_ = path.join(SIMPSON_TEST_PATH, "quad_sidebands")
    total = 2
    for i in range(total):
        filename = path.join(path_, f"test{i:02d}", f"test{i:02d}.json")  # noqa: E231
        for volume in VOLUMES:
            res = []
            data_mrsimulator, data_source, info, dim = c_setup(
                filename=filename, integration_volume=volume
            )
            res.append([data_mrsimulator, data_source])

            # random euler angle all zero. Euler angles should not affect the spectrum.
            data_mrsimulator, data_source = c_setup_random_euler_angles(
                filename, "quadrupolar"
            )
            res.append([data_mrsimulator, data_source])

        title = f"Quad Sidebands ({i + 1} / {total})"
        compile_plots(dim, res, info, title=title, report=report)

        message = f"{error_message} test0{i:02d}.json"  # noqa: E231
        check_all_close(res, message, rel_limit=1e-3)


def test_csa_plus_quadrupolar_lineshape_simpson(report):
    error_message = (
        "failed to compare quad + csa lineshape with simpson simulation from file"
    )
    path_ = path.join(SIMPSON_TEST_PATH, "csa_quad")
    total = 10
    for i in range(total):
        filename = path.join(path_, f"test{i:02d}", f"test{i:02d}.json")  # noqa: E231
        res = []
        for volume in VOLUMES:
            data_mrsimulator, data_source, info, dim = c_setup(
                filename=filename, integration_volume=volume
            )
            res.append([data_mrsimulator, data_source])

        title = f"Quad + Shielding Sidebands ({i + 1} / {total})"
        compile_plots(dim, res, info, title=title, report=report)

        message = f"{error_message} test0{i:02d}.json"  # noqa: E231
        check_all_close(res, message, rel_limit=5e-2)


def test_1st_order_quadrupolar_lineshape_simpson(report):
    error_message = (
        "failed to compare 1st order quad lineshape with simpson simulation from file"
    )
    path_ = path.join(SIMPSON_TEST_PATH, "quad_1st_order")
    total = 2
    for i in range(total):
        filename = path.join(path_, f"test{i:02d}", f"test{i:02d}.json")  # noqa: E231
        res = []
        for volume in VOLUMES:
            data_mrsimulator, data_source, info, dim = c_setup(
                filename=filename, integration_volume=volume
            )
            res.append([data_mrsimulator, data_source])

        title = f"1st Order Quadrupolar Lineshape ({i + 1} / {total})"
        compile_plots(dim, res, info, title=title, report=report)

        message = f"{error_message} test0{i:02d}.json"  # noqa: E231
        check_all_close(res, message, rel_limit=1e-3)


def test_j_coupling_lineshape_simpson(report):
    error_message = "failed to compare j-coupling with simpson simulation from file"
    path_ = path.join(SIMPSON_TEST_PATH, "j-coupling")
    total = 20
    for i in range(total):
        filename = path.join(path_, f"test{i:02d}", f"test{i:02d}.json")  # noqa: E231
        res = []
        for volume in VOLUMES:
            data_mrsimulator, data_source, info, dim = c_setup(
                filename=filename, integration_volume=volume
            )
            res.append([data_mrsimulator, data_source])

        title = f"J-coupling Spectra ({i + 1} / {total})"
        compile_plots(dim, res, info, title=title, report=report)

        message = f"{error_message} test0{i:02d}.json"  # noqa: E231
        check_all_close(res, message, rel_limit=9e-3)


def test_dipolar_coupling_lineshape_simpson(report):
    error_message = (
        "failed to compare dipolar-coupling with simpson simulation from file"
    )
    path_ = path.join(SIMPSON_TEST_PATH, "dipolar-coupling")
    total = 7
    for i in range(total):
        filename = path.join(path_, f"test{i:02d}", f"test{i:02d}.json")  # noqa: E231
        res = []
        for volume in VOLUMES:
            data_mrsimulator, data_source, info, dim = c_setup(
                filename=filename, integration_volume=volume
            )
            res.append([data_mrsimulator, data_source])

        title = f"Dipolar-coupling Spectra ({i + 1} / {total})"
        compile_plots(dim, res, info, title=title, report=report)

        message = f"{error_message} test0{i:02d}.json"  # noqa: E231
        check_all_close(res, message, rel_limit=4e-3)


def test_2D_sideband_sideband_simpson(report):
    error_message = (
        "failed to compare sideband-sideband with simpson simulation from file"
    )
    path_ = path.join(SIMPSON_TEST_PATH, "sideband_sideband")
    total = 5
    for i in range(total):
        filename = path.join(path_, f"test{i:02d}", f"test{i:02d}.json")  # noqa: E231

        res = []
        # euler angle all zero
        data_mrsimulator, data_source, info, dim = c_setup(filename=filename)
        res.append([data_mrsimulator, data_source])

        title = f"2D sideband-sideband ({i + 1} / {total})"
        compile_plots(
            dim[0], res, info, dim2=dim[1], title=title, report=report, label="Simpson"
        )

        message = f"{error_message} test0{i}.json"
        check_all_close(res, message, rel_limit=1e-4)


# --------------------------------------------------------------------------- #
# Test against rmnsim calculations
# --------------------------------------------------------------------------- #


def test_quad_csa_cross_rmnsim(report):
    error_message = (
        "failed to compare quad-csa cross-term spectra with rmnsim from file"
    )
    total = 7
    for j_, folder in enumerate(["quad_csa_cross1", "quad_csa_cross2"]):
        path_ = path.join(RNMSIM_TEST_PATH, folder)
        for i in range(total):
            filename = path.join(
                path_, f"test{i:02d}", f"test{i:02d}.json"  # noqa: E231
            )

            res = []
            # euler angle all zero
            data_mrsimulator, data_source, info, dim = c_setup(
                filename, number_of_sidebands=1, integration_volume="hemisphere"
            )
            res.append([data_mrsimulator, data_source])

            title = f"Quad-CSA 2nd Order Cross-Term-{j_} ({i + 1} / {total})"
            compile_plots(dim, res, info, title=title, report=report, label="rmnsim")

            message = f"{error_message} test0{i}.json"
            check_all_close(res, message, rel_limit=8e-2)


# --------------------------------------------------------------------------- #
# Test against brute-force NMR calculation where lineshapes
# are averaged over a billion orientations.
# --------------------------------------------------------------------------- #


def test_pure_shielding_static_lineshape_python_brute(report):
    error_message = (
        "failed to compare shielding lineshape with brute force simulation from file"
    )
    path_ = path.join(PYTHON_BRUTE_TEST_PATH, "shielding")
    total = 5
    for i in range(total):
        filename = path.join(path_, f"test{i:02d}", f"test{i:02d}.json")  # noqa: E231

        res = []
        # euler angle all zero
        data_mrsimulator, data_source, info, dim = c_setup(filename=filename)
        res.append([data_mrsimulator, data_source])

        # random euler angle all zero. Euler angles should not affect the spectrum.
        data_mrsimulator, data_source = c_setup_random_euler_angles(
            filename, "shielding_symmetric"
        )
        res.append([data_mrsimulator, data_source])

        title = f"Shielding Static Lineshape (Brute Force) ({i + 1} / {total})"
        compile_plots(dim, res, info, title=title, report=report, label="Brute")

        message = f"{error_message} test0{i}.json"
        check_all_close(res, message, rel_limit=1e-3)


# --------------------------------------------------------------------------- #
# Self-Test
# --------------------------------------------------------------------------- #


def test_pure_quadrupolar_lineshape_python_brute(report):
    error_message = (
        "failed to compare quadrupolar lineshape with brute force simulation from file"
    )
    path_ = path.join(PYTHON_BRUTE_TEST_PATH, "quad")
    total = 19
    for i in range(total):
        filename = path.join(path_, f"test{i:02d}", f"test{i:02d}.json")  # noqa: E231

        res = []
        data_mrsimulator, data_source, info, dim = c_setup(filename=filename)
        res.append([data_mrsimulator, data_source])

        # random euler angle. Euler angles should not affect the spectrum.
        data_mrsimulator, data_source = c_setup_random_euler_angles(
            filename, "quadrupolar"
        )
        res.append([data_mrsimulator, data_source])

        title = f"Quad Lineshape Self-Test ({i + 1} / {total})"
        compile_plots(dim, res, info, title=title, report=report, label="self")

        message = f"{error_message} test0{i:02d}.json"  # noqa: E231
        check_all_close(res, message, rel_limit=1e-3)
