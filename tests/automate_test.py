# -*- coding: utf-8 -*-
import json
from os import path

import numpy as np
from numpy.fft import fft
from numpy.fft import fftshift

from mrsimulator import Isotopomer, Spectrum, Simulator

TEST_C = True


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
    spectrum = Spectrum.parse_json_with_units(data_object["spectrum"])
    isotopomer = [
        Isotopomer.parse_json_with_units(item) for item in data_object["isotopomers"]
    ]
    print(spectrum)
    print(isotopomer[0].sites[0])

    s1 = Simulator(isotopomer, spectrum)
    data_mrsimulator = s1.one_d_spectrum(geodesic_polyhedron_frequency=120)[1]
    data_mrsimulator /= data_mrsimulator.max()

    return data_mrsimulator, data_source


# --------------------------------------------------------------------------- #
# Test against simpson calculations


def test00_sim_shielding():
    path_ = path.join("tests", "simpson")
    file_ = path.join(path_, "test00", "test00.json")

    data_mrsimulator, data_source = c_setup(*get_data(file_))
    np.testing.assert_almost_equal(data_mrsimulator, data_source, decimal=2)


def test01_sim_shielding():
    path_ = path.join("tests", "simpson")
    file_ = path.join(path_, "test01", "test01.json")
    data_mrsimulator, data_source = c_setup(*get_data(file_))
    np.testing.assert_almost_equal(data_mrsimulator, data_source, decimal=2)


def test02_sim_shielding():
    path_ = path.join("tests", "simpson")
    file_ = path.join(path_, "test02", "test02.json")
    data_mrsimulator, data_source = c_setup(*get_data(file_))
    np.testing.assert_almost_equal(data_mrsimulator, data_source, decimal=2)


def test03_sim_shielding():
    path_ = path.join("tests", "simpson")
    file_ = path.join(path_, "test03", "test03.json")
    data_mrsimulator, data_source = c_setup(*get_data(file_))
    np.testing.assert_almost_equal(data_mrsimulator, data_source, decimal=2)


def test04_sim_shielding():
    path_ = path.join("tests", "simpson")
    file_ = path.join(path_, "test04", "test04.json")
    data_mrsimulator, data_source = c_setup(*get_data(file_))
    np.testing.assert_almost_equal(data_mrsimulator, data_source, decimal=2)


def test05_sim_shielding():
    path_ = path.join("tests", "simpson")
    file_ = path.join(path_, "test05", "test05.json")
    data_mrsimulator, data_source = c_setup(*get_data(file_))
    np.testing.assert_almost_equal(data_mrsimulator, data_source, decimal=2)


def test06_sim_shielding():
    path_ = path.join("tests", "simpson")
    file_ = path.join(path_, "test06", "test06.json")
    data_mrsimulator, data_source = c_setup(*get_data(file_))
    np.testing.assert_almost_equal(data_mrsimulator, data_source, decimal=2)


def test07_sim_shielding():
    path_ = path.join("tests", "simpson")
    file_ = path.join(path_, "test07", "test07.json")
    data_mrsimulator, data_source = c_setup(*get_data(file_))
    np.testing.assert_almost_equal(data_mrsimulator, data_source, decimal=2)


# --------------------------------------------------------------------------- #
# Test against brute-force NMR calculation where lineshapes
# are averaged over a billion orientations.


def test00_python_shielding():
    path_ = path.join("tests", "python")
    file_ = path.join(path_, "test00", "test00.json")
    data_mrsimulator, data_source = c_setup(*get_data(file_))
    np.testing.assert_almost_equal(data_mrsimulator, data_source, decimal=2)


def test01_python_shielding():
    path_ = path.join("tests", "python")
    file_ = path.join(path_, "test01", "test01.json")
    data_mrsimulator, data_source = c_setup(*get_data(file_))
    np.testing.assert_almost_equal(data_mrsimulator, data_source, decimal=2)


def test02_python_shielding():
    path_ = path.join("tests", "python")
    file_ = path.join(path_, "test02", "test02.json")
    data_mrsimulator, data_source = c_setup(*get_data(file_))
    np.testing.assert_almost_equal(data_mrsimulator, data_source, decimal=2)


def test03_python_shielding():
    path_ = path.join("tests", "python")
    file_ = path.join(path_, "test03", "test03.json")
    data_mrsimulator, data_source = c_setup(*get_data(file_))
    np.testing.assert_almost_equal(data_mrsimulator, data_source, decimal=2)


def test04_python_shielding():
    path_ = path.join("tests", "python")
    file_ = path.join(path_, "test04", "test04.json")
    data_mrsimulator, data_source = c_setup(*get_data(file_))
    np.testing.assert_almost_equal(data_mrsimulator, data_source, decimal=2)


# Test quad


def test_python_quad_array():
    path_ = path.join("tests", "python", "quad")
    for i in range(19):
        file_ = path.join(path_, f"test{i:02d}", f"test{i:02d}.json")
        data_mrsimulator, data_source = c_setup(*get_data(file_))
        np.testing.assert_almost_equal(data_mrsimulator, data_source, decimal=2)
