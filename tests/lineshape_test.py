# -*- coding: utf-8 -*-
import json
from os import path

import numpy as np
from mrsimulator import Dimension
from mrsimulator import Isotopomer
from mrsimulator import Simulator
from numpy.fft import fft
from numpy.fft import fftshift


def _import_json(filename):
    with open(filename, "rb") as f:
        content = f.read()
        return json.loads(str(content, encoding="UTF-8"))


def _get_header_and_footer(source_file):
    """
    Return the number of rows in the header and footer.
    This assumes that the data is listed in-between lines
    with keywords 'DATA' and 'END' respectively.
    """
    f = open(source_file)
    skip_header = 0
    total_lines = 0
    for line_ in f:
        skip_header += 1
        total_lines += 1
        if "DATA" in line_:
            total_lines += 1
            skip_footer = 1
            break
    f.close()
    if total_lines == skip_header:
        skip_header = 0
        skip_footer = 0
    return skip_header, skip_footer


def get_data(filename):
    """Load a simpson or DMfit output file"""

    # source data
    data_object = _import_json(filename)
    test_data_object = data_object["test_data"]

    source_file = test_data_object["filename"]
    path_ = path.abspath(filename)
    path_ = path.split(path_)[0]
    source_file = path.join(path_, source_file)

    delimiter = None
    periodic = False
    if "delimiter" in test_data_object.keys():
        delimiter = test_data_object["delimiter"]
    if "periodic" in test_data_object.keys():
        periodic = test_data_object["periodic"]

    if test_data_object["type"] == "npz":
        data_source = np.load(source_file)[1]
        data_source /= data_source.max()
        return data_object, data_source

    skip_header, skip_footer = _get_header_and_footer(source_file)

    if test_data_object["type"] == "complex":
        x, y = np.genfromtxt(
            source_file,
            skip_header=skip_header,
            delimiter=delimiter,
            skip_footer=skip_footer,
            unpack=True,
        )
        data_source = x + 1j * y
    if test_data_object["type"] == "real":
        x = np.genfromtxt(
            source_file,
            skip_header=skip_header,
            delimiter=delimiter,
            skip_footer=skip_footer,
            usecols=0,
        )
        data_source = x

    if test_data_object["quantity"] == "time":
        if periodic:
            data_source = fftshift(fft(data_source)).real
        else:
            data_source[0] /= 2.0
            data_source = fftshift(fft(data_source)).real

    data_source /= data_source.max()

    if test_data_object["source"] == "dmfit":
        data_source = data_source[::-1]
        data_source = np.roll(data_source, 1)

    return data_object, data_source


def c_setup(data_object, data_source):
    # mrsimulator
    dimensions = [
        Dimension.parse_dict_with_units(item) for item in data_object["dimensions"]
    ]
    isotopomers = [
        Isotopomer.parse_dict_with_units(item) for item in data_object["isotopomers"]
    ]

    s1 = Simulator(isotopomers=isotopomers, dimensions=dimensions)
    data_mrsimulator = s1.one_d_spectrum(integration_density=120)[1]
    data_mrsimulator /= data_mrsimulator.max()

    return data_mrsimulator, data_source


def c_setup_random_euler_angles(data_object, data_source):
    # mrsimulator
    dimensions = [
        Dimension.parse_dict_with_units(item) for item in data_object["dimensions"]
    ]
    isotopomers = [
        Isotopomer.parse_dict_with_units(isotopomer)
        for isotopomer in data_object["isotopomers"]
    ]
    for isotopomer in isotopomers:
        isotopomer.sites[0].shielding_symmetric.alpha = np.random.rand(1) * 2 * np.pi
        isotopomer.sites[0].shielding_symmetric.beta = np.random.rand(1) * 2 * np.pi
        isotopomer.sites[0].shielding_symmetric.gamma = np.random.rand(1) * 2 * np.pi

    s1 = Simulator(isotopomers=isotopomers, dimensions=dimensions)
    data_mrsimulator = s1.one_d_spectrum(integration_density=120)[1]
    data_mrsimulator /= data_mrsimulator.max()

    return data_mrsimulator, data_source


# --------------------------------------------------------------------------- #
# Test against simpson calculations
path_for_simpson_test_shielding = path.join("tests", "simpson", "shielding")


def test_shielding_simulation_against_simpson():
    error_message = (
        "failed to compare shielding lineshape with simpson simulation from file"
    )
    for i in range(8):
        message = f"{error_message} test0{i}.json"
        filename = path.join(
            path_for_simpson_test_shielding, f"test{i:02d}", f"test{i:02d}.json"
        )

        error_message = (
            f"failed to compare shielding simulation with file test0{i}.json"
        )
        # euler angle all zero
        data_mrsimulator, data_source = c_setup(*get_data(filename))
        np.testing.assert_almost_equal(
            data_mrsimulator, data_source, decimal=2, err_msg=message
        )

        # random euler angle all zero. Euler angles should not affect the spectrum.
        data_mrsimulator, data_source = c_setup_random_euler_angles(*get_data(filename))
        np.testing.assert_almost_equal(
            data_mrsimulator, data_source, decimal=2, err_msg=message
        )


# --------------------------------------------------------------------------- #
# Test against brute-force NMR calculation where lineshapes
# are averaged over a billion orientations.
path_for_python_test_shielding = path.join("tests", "python", "shielding")


def test_shielding_against_brute_force_lineshape_simulation():
    error_message = (
        "failed to compare shielding lineshape with brute force simulation from file"
    )
    for i in range(5):
        message = f"{error_message} test0{i}.json"
        filename = path.join(
            path_for_python_test_shielding, f"test{i:02d}", f"test{i:02d}.json"
        )

        # euler angle all zero
        data_mrsimulator, data_source = c_setup(*get_data(filename))
        np.testing.assert_almost_equal(
            data_mrsimulator, data_source, decimal=2, err_msg=message
        )

        # random euler angle all zero. Euler angles should not affect the spectrum.
        data_mrsimulator, data_source = c_setup_random_euler_angles(*get_data(filename))
        np.testing.assert_almost_equal(
            data_mrsimulator, data_source, decimal=2, err_msg=message
        )


# Test pure quadrupole lineshape simulation
path_for_quad_test_mrsimulator = path.join("tests", "python", "quad")


def test_pure_quadrupolar_lineshapes():
    error_message = (
        "failed to compare quadrupolar lineshape simulation with data from file"
    )
    for i in range(19):
        message = f"{error_message} test0{i:02d}.json"
        filename = path.join(
            path_for_quad_test_mrsimulator, f"test{i:02d}", f"test{i:02d}.json"
        )
        data_mrsimulator, data_source = c_setup(*get_data(filename))
        np.testing.assert_almost_equal(
            data_mrsimulator, data_source, decimal=2, err_msg=message
        )
