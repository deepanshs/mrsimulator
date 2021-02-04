# -*- coding: utf-8 -*-
"""Lineshape Test."""
from os import path

import numpy as np

from .utils import c_setup
from .utils import c_setup_random_euler_angles

# import matplotlib.pyplot as plt


# --------------------------------------------------------------------------- #
# Test against simpson calculations

COMMON_PATH = path.join("tests", "spectral_integration_tests")
SIMPSON_TEST_PATH = path.join(COMMON_PATH, "simpson_simulated_lineshapes")
PYTHON_BRUTE_TEST_PATH = path.join(COMMON_PATH, "python_brute_force_lineshapes")


def test_pure_shielding_sideband_simulation_against_simpson():
    error_message = (
        "failed to compare shielding lineshape with simpson simulation from file"
    )
    path_ = path.join(SIMPSON_TEST_PATH, "shielding_sidebands")
    for i in range(8):
        message = f"{error_message} test0{i}.json"
        filename = path.join(path_, f"test{i:02d}", f"test{i:02d}.json")

        error_message = (
            f"failed to compare shielding simulation with file test0{i}.json"
        )
        # euler angle all zero
        data_mrsimulator, data_source = c_setup(filename)
        np.testing.assert_almost_equal(
            data_mrsimulator, data_source, decimal=2, err_msg=message
        )

        # random euler angle all zero. Euler angles should not affect the spectrum.
        data_mrsimulator, data_source = c_setup_random_euler_angles(
            filename, "shielding_symmetric"
        )
        np.testing.assert_almost_equal(
            data_mrsimulator, data_source, decimal=2, err_msg=message
        )


# --------------------------------------------------------------------------- #
# Test against brute-force NMR calculation where lineshapes
# are averaged over a billion orientations.


def test_pure_shielding_static_simulation_against_brute_force_lineshapes():
    error_message = (
        "failed to compare shielding lineshape with brute force simulation from file"
    )
    path_ = path.join(PYTHON_BRUTE_TEST_PATH, "shielding")
    for i in range(5):
        message = f"{error_message} test0{i}.json"
        filename = path.join(path_, f"test{i:02d}", f"test{i:02d}.json")

        # euler angle all zero
        data_mrsimulator, data_source = c_setup(filename)
        np.testing.assert_almost_equal(
            data_mrsimulator, data_source, decimal=2, err_msg=message
        )

        # random euler angle all zero. Euler angles should not affect the spectrum.
        data_mrsimulator, data_source = c_setup_random_euler_angles(
            filename, "shielding_symmetric"
        )
        np.testing.assert_almost_equal(
            data_mrsimulator, data_source, decimal=2, err_msg=message
        )


def test_pure_quadrupolar_lineshapes():
    error_message = (
        "failed to compare quadrupolar lineshape simulation with data from file"
    )
    path_ = path.join(PYTHON_BRUTE_TEST_PATH, "quad")
    for i in range(19):
        message = f"{error_message} test0{i:02d}.json"
        filename = path.join(path_, f"test{i:02d}", f"test{i:02d}.json")
        data_mrsimulator, data_source = c_setup(filename)
        np.testing.assert_almost_equal(
            data_mrsimulator, data_source, decimal=2, err_msg=message
        )

        # random euler angle. Euler angles should not affect the spectrum.
        data_mrsimulator, data_source = c_setup_random_euler_angles(
            filename, "quadrupolar"
        )
        np.testing.assert_almost_equal(
            data_mrsimulator, data_source, decimal=1, err_msg=message
        )


def test_pure_quadrupolar_simpson_sidebands():
    error_message = (
        "failed to compare quadrupolar lineshape simulation with data from file"
    )
    path_ = path.join(SIMPSON_TEST_PATH, "quad_sidebands")
    for i in range(2):
        message = f"{error_message} test0{i:02d}.json"
        filename = path.join(path_, f"test{i:02d}", f"test{i:02d}.json")
        data_mrsimulator, data_source = c_setup(
            filename, integration_volume="hemisphere"
        )
        np.testing.assert_almost_equal(
            data_mrsimulator, data_source, decimal=1, err_msg=message
        )


def test_csa_plus_quadrupolar_simpson_lineshapes():
    error_message = (
        "failed to compare quadrupolar lineshape simulation with data from file"
    )
    path_ = path.join(SIMPSON_TEST_PATH, "csa_quad")
    for i in range(6):
        message = f"{error_message} test0{i:02d}.json"
        filename = path.join(path_, f"test{i:02d}", f"test{i:02d}.json")
        data_mrsimulator, data_source = c_setup(
            filename, integration_volume="hemisphere"
        )
        np.testing.assert_almost_equal(
            data_mrsimulator, data_source, decimal=1, err_msg=message
        )


def test_j_coupling_simpson_lineshapes():
    error_message = (
        "failed to compare j-coupling lineshape simulation with data from file"
    )
    path_ = path.join(SIMPSON_TEST_PATH, "j-coupling")
    for i in range(14):
        message = f"{error_message} test0{i:02d}.json"
        filename = path.join(path_, f"test{i:02d}", f"test{i:02d}.json")
        data_mrsimulator, data_source = c_setup(
            filename, integration_volume="hemisphere"
        )
        # plt.plot(data_mrsimulator)
        # plt.plot(data_source)
        # plt.show()
        np.testing.assert_almost_equal(
            data_mrsimulator, data_source, decimal=1, err_msg=message
        )
