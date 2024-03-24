"""Lineshape Test."""
from os import path
from pprint import pformat

import numpy as np
import pytest

from .utils import c_setup
from .utils import c_setup_random_euler_angles

COMMON_PATH = path.join("tests", "spectral_integration_tests")
SIMPSON_TEST_PATH = path.join(COMMON_PATH, "simpson_simulated_lineshapes")
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
    if __GENERATE_REPORT__:
        pdf = PdfPages("reports/lineshapes.pdf")
    else:
        pdf = None
    return pdf


@pytest.fixture(autouse=True, scope="module")
def report_pdf(report):
    """Generate report"""
    yield
    if __GENERATE_REPORT__:
        report.close()


def compile_plots(dim, rep, info, title=None, report=None, label="simpson"):
    """Test report plot"""
    _, ax = plt.subplots(1, 2, figsize=(9, 4), gridspec_kw={"width_ratios": [1, 1]})

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

    format_kwargs = dict(indent=1, width=100, compact=True, sort_dicts=False)
    str_config = pformat(info["config"], **format_kwargs)
    str_sys = pformat(info["spin_systems"], **format_kwargs)
    str_method = pformat(info["methods"], **format_kwargs)

    str_config = r"$\bf{sim.config}$" + "\n" + str_config + "\n\n"
    str_sys = r"$\bf{sim.spin\_systems}$" + "\n" + str_sys + "\n\n"
    str_method = r"$\bf{sim.methods}$" + "\n" + str_method

    text = str_config + str_sys + str_method
    ax[-1].set_ylim(0, 1)
    ax[-1].text(
        0, 1, text, fontsize=7.5, verticalalignment="top", horizontalalignment="left"
    )
    ax[-1].set_title(f"Test series: {title}")
    ax[-1].axis("off")
    plt.tight_layout()
    if report is not None:
        report.savefig(dpi=150)
    plt.close()


def check_all_close(res, message, rel_limit):
    """Check if the vectos in res are all close within relative limits"""
    for item in res:
        limit = -np.log10(item[1].max()) + rel_limit
        np.testing.assert_almost_equal(item[0], item[1], decimal=limit, err_msg=message)


# --------------------------------------------------------------------------- #
# Test against simpson calculations
# --------------------------------------------------------------------------- #


def test_pure_shielding_sideband_simpson(report):
    error_message = (
        "failed to compare shielding sidebands with simpson simulation from file"
    )
    path_ = path.join(SIMPSON_TEST_PATH, "shielding_sidebands")
    for i in range(8):
        filename = path.join(path_, f"test{i:02d}", f"test{i:02d}.json")

        res = []
        # euler angle all zero
        data_mrsimulator, data_source, info, dim = c_setup(filename)
        res.append([data_mrsimulator, data_source])

        # random euler angle all zero. Euler angles should not affect the spectrum.
        data_mrsimulator, data_source = c_setup_random_euler_angles(
            filename, "shielding_symmetric"
        )
        res.append([data_mrsimulator, data_source])

        if __GENERATE_REPORT__:
            compile_plots(dim, res, info, title="Shielding Sidebands", report=report)

        message = f"{error_message} test0{i}.json"
        check_all_close(res, message, rel_limit=2.5)


def test_pure_quadrupolar_sidebands_simpson(report):
    error_message = (
        "failed to compare quadrupolar sidebands with simpson simulation from file"
    )
    path_ = path.join(SIMPSON_TEST_PATH, "quad_sidebands")
    for i in range(2):
        filename = path.join(path_, f"test{i:02d}", f"test{i:02d}.json")
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

        if __GENERATE_REPORT__:
            compile_plots(dim, res, info, title="Quad Sidebands", report=report)

        message = f"{error_message} test0{i:02d}.json"
        check_all_close(res, message, rel_limit=1.5)


def test_csa_plus_quadrupolar_lineshape_simpson(report):
    error_message = (
        "failed to compare quad + csa lineshape with simpson simulation from file"
    )
    path_ = path.join(SIMPSON_TEST_PATH, "csa_quad")
    for i in range(10):
        print(i)
        filename = path.join(path_, f"test{i:02d}", f"test{i:02d}.json")
        res = []
        for volume in VOLUMES:
            data_mrsimulator, data_source, info, dim = c_setup(
                filename=filename, integration_volume=volume
            )
            res.append([data_mrsimulator, data_source])

        if __GENERATE_REPORT__:
            compile_plots(
                dim,
                res,
                info,
                title="Quad + Shielding Sidebands",
                report=report,
            )

        message = f"{error_message} test0{i:02d}.json"
        check_all_close(res, message, rel_limit=0.9)


def test_1st_order_quadrupolar_lineshape_simpson(report):
    error_message = (
        "failed to compare 1st order quad lineshape with simpson simulation from file"
    )
    path_ = path.join(SIMPSON_TEST_PATH, "quad_1st_order")
    for i in range(2):
        filename = path.join(path_, f"test{i:02d}", f"test{i:02d}.json")
        res = []
        for volume in VOLUMES:
            data_mrsimulator, data_source, info, dim = c_setup(
                filename=filename, integration_volume=volume
            )
            res.append([data_mrsimulator, data_source])

        if __GENERATE_REPORT__:
            compile_plots(
                dim,
                res,
                info,
                title="1st Order Quadrupolar Lineshape",
                report=report,
            )

        message = f"{error_message} test0{i:02d}.json"
        check_all_close(res, message, rel_limit=1.0)


def test_j_coupling_lineshape_simpson(report):
    error_message = "failed to compare j-coupling with simpson simulation from file"
    path_ = path.join(SIMPSON_TEST_PATH, "j-coupling")
    for i in range(20):
        filename = path.join(path_, f"test{i:02d}", f"test{i:02d}.json")
        res = []
        for volume in VOLUMES:
            data_mrsimulator, data_source, info, dim = c_setup(
                filename=filename, integration_volume=volume
            )
            res.append([data_mrsimulator, data_source])

        if __GENERATE_REPORT__:
            compile_plots(dim, res, info, title="J-coupling Spectra", report=report)

        message = f"{error_message} test0{i:02d}.json"
        check_all_close(res, message, rel_limit=1.1)


def test_dipolar_coupling_lineshape_simpson(report):
    error_message = (
        "failed to compare dipolar-coupling with simpson simulation from file"
    )
    path_ = path.join(SIMPSON_TEST_PATH, "dipolar-coupling")
    for i in range(7):
        filename = path.join(path_, f"test{i:02d}", f"test{i:02d}.json")
        res = []
        for volume in VOLUMES:
            data_mrsimulator, data_source, info, dim = c_setup(
                filename=filename, integration_volume=volume
            )
            res.append([data_mrsimulator, data_source])

        if __GENERATE_REPORT__:
            compile_plots(
                dim, res, info, title="Dipolar-coupling Spectra", report=report
            )

        message = f"{error_message} test0{i:02d}.json"
        check_all_close(res, message, rel_limit=1.5)


# --------------------------------------------------------------------------- #
# Test against brute-force NMR calculation where lineshapes
# are averaged over a billion orientations.
# --------------------------------------------------------------------------- #


def test_pure_shielding_static_lineshape_python_brute(report):
    error_message = (
        "failed to compare shielding lineshape with brute force simulation from file"
    )
    path_ = path.join(PYTHON_BRUTE_TEST_PATH, "shielding")
    for i in range(5):
        filename = path.join(path_, f"test{i:02d}", f"test{i:02d}.json")

        res = []
        # euler angle all zero
        data_mrsimulator, data_source, info, dim = c_setup(filename=filename)
        res.append([data_mrsimulator, data_source])

        # random euler angle all zero. Euler angles should not affect the spectrum.
        data_mrsimulator, data_source = c_setup_random_euler_angles(
            filename, "shielding_symmetric"
        )
        res.append([data_mrsimulator, data_source])

        if __GENERATE_REPORT__:
            compile_plots(
                dim,
                res,
                info,
                title="Shielding Static Lineshape (Brute Force)",
                report=report,
                label="Brute",
            )

        message = f"{error_message} test0{i}.json"
        check_all_close(res, message, rel_limit=2)


# --------------------------------------------------------------------------- #
# Self-Test
# --------------------------------------------------------------------------- #


def test_pure_quadrupolar_lineshape_python_brute(report):
    error_message = (
        "failed to compare quadrupolar lineshape with brute force simulation from file"
    )
    path_ = path.join(PYTHON_BRUTE_TEST_PATH, "quad")
    for i in range(19):
        filename = path.join(path_, f"test{i:02d}", f"test{i:02d}.json")

        res = []
        data_mrsimulator, data_source, info, dim = c_setup(filename=filename)
        res.append([data_mrsimulator, data_source])

        # random euler angle. Euler angles should not affect the spectrum.
        data_mrsimulator, data_source = c_setup_random_euler_angles(
            filename, "quadrupolar"
        )
        res.append([data_mrsimulator, data_source])

        if __GENERATE_REPORT__:
            compile_plots(
                dim,
                res,
                info,
                title="Quad Lineshape Self-Test",
                report=report,
                label="self",
            )

        message = f"{error_message} test0{i:02d}.json"
        check_all_close(res, message, rel_limit=1.5)
