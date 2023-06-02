"""Lineshape Test."""
import sys
from os import path

import numpy as np

np.set_printoptions(threshold=sys.maxsize)

from .utils import c_setup
from .utils import c_setup_random_euler_angles

# import matplotlib.pyplot as plt
# SHOW_PLOTS = True

COMMON_PATH = path.join("tests", "spectral_integration_tests")
SIMPSON_TEST_PATH = path.join(COMMON_PATH, "simpson_simulated_lineshapes")
PYTHON_BRUTE_TEST_PATH = path.join(COMMON_PATH, "python_brute_force_lineshapes")

# --------------------------------------------------------------------------- #
# Test against brute-force NMR calculation where lineshapes
# are averaged over a billion orientations.
# --------------------------------------------------------------------------- #


def test_pure_shielding_static_lineshape_python_brute():
    error_message = (
        "failed to compare shielding lineshape with brute force simulation from file"
    )
    path_ = path.join(PYTHON_BRUTE_TEST_PATH, "shielding")
    for i in range(5):
        message = f"{error_message} test0{i}.json"
        filename = path.join(path_, f"test{i:02d}", f"test{i:02d}.json")

        # euler angle all zero
        data_mrsimulator, data_source = c_setup(filename=filename)

        # if SHOW_PLOTS:
        #     plt.plot(data_mrsimulator, "k", label="mrsims")
        #     plt.plot(data_source, "--r", label="py brute force ($10^9$ orientations)")
        #     plt.title("Shielding Static Lineshape")
        #     plt.legend()
        #     plt.show()

        limit = -np.log10(data_source.max()) + 2
        np.testing.assert_almost_equal(
            data_mrsimulator, data_source, decimal=limit, err_msg=message
        )

        # random euler angle all zero. Euler angles should not affect the spectrum.
        data_mrsimulator, data_source = c_setup_random_euler_angles(
            filename, "shielding_symmetric"
        )

        # if SHOW_PLOTS:
        #     plt.plot(data_mrsimulator, "k", label="mrsims")
        #     plt.plot(data_source, "--r", label="py brute force ($10^9$ orientations)")
        #     plt.title("Shielding Static Lineshape")
        #     plt.legend()
        #     plt.show()

        np.testing.assert_almost_equal(
            data_mrsimulator, data_source, decimal=limit, err_msg=message
        )


def test_pure_quadrupolar_lineshape_python_brute():
    error_message = (
        "failed to compare quadrupolar lineshape with brute force simulation from file"
    )
    path_ = path.join(PYTHON_BRUTE_TEST_PATH, "quad")
    for i in range(19):
        message = f"{error_message} test0{i:02d}.json"
        filename = path.join(path_, f"test{i:02d}", f"test{i:02d}.json")
        data_mrsimulator, data_source = c_setup(filename=filename)

        # if SHOW_PLOTS:
        #     plt.plot(data_mrsimulator, "k", label="mrsims")
        #     plt.plot(data_source, "--r", label="py brute force ($10^9$ orientations)")
        #     plt.title("Quad Static Lineshape")
        #     plt.legend()
        #     plt.show()

        limit = -np.log10(data_source.max()) + 1.5
        np.testing.assert_almost_equal(
            data_mrsimulator, data_source, decimal=limit, err_msg=message
        )

        # random euler angle. Euler angles should not affect the spectrum.
        data_mrsimulator, data_source = c_setup_random_euler_angles(
            filename, "quadrupolar"
        )

        # if SHOW_PLOTS:
        #     plt.plot(data_mrsimulator, "k", label="mrsims")
        #     plt.plot(data_source, "--r", label="py brute force ($10^9$ orientations)")
        #     plt.title("Quad Static Lineshape")
        #     plt.legend()
        #     plt.show()

        np.testing.assert_almost_equal(
            data_mrsimulator, data_source, decimal=limit, err_msg=message
        )


# --------------------------------------------------------------------------- #
# Test against simpson calculations
# --------------------------------------------------------------------------- #


def test_pure_shielding_sideband_simpson():
    error_message = (
        "failed to compare shielding sidebands with simpson simulation from file"
    )
    path_ = path.join(SIMPSON_TEST_PATH, "shielding_sidebands")
    for i in range(8):
        message = f"{error_message} test0{i}.json"
        filename = path.join(path_, f"test{i:02d}", f"test{i:02d}.json")

        # euler angle all zero
        data_mrsimulator, data_source = c_setup(filename)

        # if SHOW_PLOTS:
        #     plt.plot(data_mrsimulator, "k", label="mrsims")
        #     plt.plot(data_source, "--r", label="simpson")
        #     plt.title("Shielding Sidebands")
        #     plt.legend()
        #     plt.show()

        limit = -np.log10(data_source.max()) + 2.5
        np.testing.assert_almost_equal(
            data_mrsimulator, data_source, decimal=limit, err_msg=message
        )

        # random euler angle all zero. Euler angles should not affect the spectrum.
        data_mrsimulator, data_source = c_setup_random_euler_angles(
            filename, "shielding_symmetric"
        )

        # if SHOW_PLOTS:
        #     plt.plot(data_mrsimulator, "k", label="mrsims")
        #     plt.plot(data_source, "--r", label="simpson")
        #     plt.title("Shielding Sidebands")
        #     plt.legend()
        #     plt.show()

        np.testing.assert_almost_equal(
            data_mrsimulator, data_source, decimal=limit, err_msg=message
        )


def test_pure_quadrupolar_sidebands_simpson():
    error_message = (
        "failed to compare quadrupolar sidebands with simpson simulation from file"
    )
    path_ = path.join(SIMPSON_TEST_PATH, "quad_sidebands")
    for i in range(2):
        message = f"{error_message} test0{i:02d}.json"
        filename = path.join(path_, f"test{i:02d}", f"test{i:02d}.json")
        data_mrsimulator, data_source = c_setup(
            filename=filename, integration_volume="hemisphere"
        )

        # if SHOW_PLOTS:
        #     plt.plot(data_mrsimulator, "k", label="mrsims")
        #     plt.plot(data_source, "--r", label="simpson")
        #     plt.title("Quad Sidebands")
        #     plt.legend()
        #     plt.show()

        limit = -np.log10(data_source.max()) + 2
        np.testing.assert_almost_equal(
            data_mrsimulator, data_source, decimal=limit, err_msg=message
        )

        # random euler angle all zero. Euler angles should not affect the spectrum.
        data_mrsimulator, data_source = c_setup_random_euler_angles(
            filename, "quadrupolar"
        )

        # if SHOW_PLOTS:
        #     plt.plot(data_mrsimulator, "k", label="mrsims")
        #     plt.plot(data_source, "--r", label="simpson")
        #     plt.title("Quad Sidebands")
        #     plt.legend()
        #     plt.show()

        np.testing.assert_almost_equal(
            data_mrsimulator, data_source, decimal=limit, err_msg=message
        )


def test_csa_plus_quadrupolar_lineshape_simpson():
    error_message = (
        "failed to compare quad + csa lineshape with simpson simulation from file"
    )
    path_ = path.join(SIMPSON_TEST_PATH, "csa_quad")
    for i in range(6):
        message = f"{error_message} test0{i:02d}.json"
        filename = path.join(path_, f"test{i:02d}", f"test{i:02d}.json")
        data_mrsimulator, data_source = c_setup(
            filename=filename, integration_volume="hemisphere"
        )

        # if SHOW_PLOTS:
        #     plt.plot(data_mrsimulator, "k", label="mrsims")
        #     plt.plot(data_source, "--r", label="simpson")
        #     plt.title("Quad + Shielding Sidebands")
        #     plt.legend()
        #     plt.show()

        limit = -np.log10(data_source.max()) + 1
        np.testing.assert_almost_equal(
            data_mrsimulator, data_source, decimal=limit, err_msg=message
        )


def test_j_coupling_lineshape_simpson():
    error_message = "failed to compare j-coupling with simpson simulation from file"
    path_ = path.join(SIMPSON_TEST_PATH, "j-coupling")
    for i in range(20):
        message = f"{error_message} test0{i:02d}.json"
        filename = path.join(path_, f"test{i:02d}", f"test{i:02d}.json")
        data_mrsimulator, data_source = c_setup(
            filename=filename,
            integration_volume="hemisphere",
        )

        # if SHOW_PLOTS:
        #     plt.plot(data_mrsimulator, "k", label="mrsims")
        #     plt.plot(data_source, "--r", label="simpson")
        #     plt.title("J-coupling Spectra")
        #     plt.legend()
        #     plt.show()

        limit = -np.log10(data_source.max()) + 1.1
        np.testing.assert_almost_equal(
            data_mrsimulator, data_source, decimal=limit, err_msg=message
        )


def test_dipolar_coupling_lineshape_simpson():
    error_message = (
        "failed to compare dipolar-coupling with simpson simulation from file"
    )
    path_ = path.join(SIMPSON_TEST_PATH, "dipolar-coupling")
    for i in range(2, 3, 1):
        message = f"{error_message} test0{i:02d}.json"
        filename = path.join(path_, f"test{i:02d}", f"test{i:02d}.json")
        data_mrsimulator, data_source = c_setup(
            filename=filename, integration_volume="hemisphere"
        )

        # print()
        # print(f"{i} mrsim", data_mrsimulator)
        # print(f"{i} simpson", data_source)
        # print()
        # print()
        # print()

        # if SHOW_PLOTS:
        #     plt.plot(data_mrsimulator, "k", label="mrsims")
        #     plt.plot(data_source, "--r", label="simpson")
        #     plt.title("Dipolar-coupling Spectra")
        #     plt.legend()
        #     plt.show()

        limit = -np.log10(data_source.max()) + 1.3
        np.testing.assert_almost_equal(
            data_mrsimulator, data_source, decimal=limit, err_msg=message
        )
