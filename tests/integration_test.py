# -*- coding: utf-8 -*-
import json
from os import path

import numpy as np
from numpy.fft import fft
from numpy.fft import fftshift

from mrsimulator import Simulator
from mrsimulator.methods import one_d_spectrum
from mrsimulator.python.simulator import simulator

TEST_PY = True


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
    """Load a simpson output file"""

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


def c_test(data_object, data_source):
    # mrsimulator
    spectrum, isotopomer = data_object["spectrum"], data_object["isotopomer"]

    s1 = Simulator(isotopomer)
    s1.spectrum = spectrum
    freq_c, data_mrsimulator_c = s1.run(
        one_d_spectrum, geodesic_polyhedron_frequency=120
    )
    data_mrsimulator_c /= data_mrsimulator_c.max()
    satisfy_c = np.all((data_mrsimulator_c - data_source) < 0.01)
    return satisfy_c


def py_test(data_object, data_source):
    # mrsimulator
    spectrum, isotopomer = data_object["spectrum"], data_object["isotopomer"]

    s1 = Simulator(isotopomer)
    s1.spectrum = spectrum
    freq_py, data_mrsimulator_py = simulator(
        isotopomers=s1._isotopomers_c, spectrum=s1._spectrum_c
    )
    data_mrsimulator_py /= data_mrsimulator_py.max()
    satisfy_py = np.all((data_mrsimulator_py - data_source) < 0.01)
    return satisfy_py


# --------------------------------------------------------------------------- #
# The test pass criterion
# np.all((mrsimulator_vector - test_vector) < 0.01
# --------------------------------------------------------------------------- #

# --------------------------------------------------------------------------- #
# Test against simpson calculations


def test_integration_sim00_c():
    path_ = path.join("tests", "simpson")
    file_ = path.join(path_, "test00", "test00.json")
    assert c_test(*get_data(file_))


def test_integration_sim01_c():
    path_ = path.join("tests", "simpson")
    file_ = path.join(path_, "test01", "test01.json")
    assert c_test(*get_data(file_))


def test_integration_sim02_c():
    path_ = path.join("tests", "simpson")
    file_ = path.join(path_, "test02", "test02.json")
    assert c_test(*get_data(file_))


def test_integration_sim03_c():
    path_ = path.join("tests", "simpson")
    file_ = path.join(path_, "test03", "test03.json")
    assert c_test(*get_data(file_))


def test_integration_sim04_c():
    path_ = path.join("tests", "simpson")
    file_ = path.join(path_, "test04", "test04.json")
    assert c_test(*get_data(file_))


def test_integration_sim05_c():
    path_ = path.join("tests", "simpson")
    file_ = path.join(path_, "test05", "test05.json")
    assert c_test(*get_data(file_))


def test_integration_sim06_c():
    path_ = path.join("tests", "simpson")
    file_ = path.join(path_, "test06", "test06.json")
    assert c_test(*get_data(file_))


def test_integration_sim07_c():
    path_ = path.join("tests", "simpson")
    file_ = path.join(path_, "test07", "test07.json")
    assert c_test(*get_data(file_))


if TEST_PY:

    def test_integration_sim00_py():
        path_ = path.join("tests", "simpson")
        file_ = path.join(path_, "test00", "test00.json")
        assert py_test(*get_data(file_))

    def test_integration_sim01_py():
        path_ = path.join("tests", "simpson")
        file_ = path.join(path_, "test01", "test01.json")
        assert py_test(*get_data(file_))

    def test_integration_sim02_py():
        path_ = path.join("tests", "simpson")
        file_ = path.join(path_, "test02", "test02.json")
        assert py_test(*get_data(file_))

    def test_integration_sim03_py():
        path_ = path.join("tests", "simpson")
        file_ = path.join(path_, "test03", "test03.json")
        assert py_test(*get_data(file_))

    def test_integration_sim04_py():
        path_ = path.join("tests", "simpson")
        file_ = path.join(path_, "test04", "test04.json")
        assert py_test(*get_data(file_))

    def test_integration_sim05_py():
        path_ = path.join("tests", "simpson")
        file_ = path.join(path_, "test05", "test05.json")
        assert py_test(*get_data(file_))

    def test_integration_sim06_py():
        path_ = path.join("tests", "simpson")
        file_ = path.join(path_, "test06", "test06.json")
        assert py_test(*get_data(file_))

    def test_integration_sim07_py():
        path_ = path.join("tests", "simpson")
        file_ = path.join(path_, "test07", "test07.json")
        assert py_test(*get_data(file_))


# --------------------------------------------------------------------------- #
# Test against brute-force NMR calculation where lineshapes
# are averaged over a billion orientations.


def test_integration_python00_c():
    path_ = path.join("tests", "python")
    file_ = path.join(path_, "test00", "test00.json")
    assert c_test(*get_data(file_))


def test_integration_python01_c():
    path_ = path.join("tests", "python")
    file_ = path.join(path_, "test01", "test01.json")
    assert c_test(*get_data(file_))


def test_integration_python02_c():
    path_ = path.join("tests", "python")
    file_ = path.join(path_, "test02", "test02.json")
    assert c_test(*get_data(file_))


def test_integration_python03_c():
    path_ = path.join("tests", "python")
    file_ = path.join(path_, "test03", "test03.json")
    assert c_test(*get_data(file_))


def test_integration_python04_c():
    path_ = path.join("tests", "python")
    file_ = path.join(path_, "test04", "test04.json")
    assert c_test(*get_data(file_))


if TEST_PY:

    def test_integration_python00_py():
        path_ = path.join("tests", "python")
        file_ = path.join(path_, "test00", "test00.json")
        assert py_test(*get_data(file_))

    def test_integration_python01_py():
        path_ = path.join("tests", "python")
        file_ = path.join(path_, "test01", "test01.json")
        assert py_test(*get_data(file_))

    def test_integration_python02_py():
        path_ = path.join("tests", "python")
        file_ = path.join(path_, "test02", "test02.json")
        assert py_test(*get_data(file_))

    def test_integration_python03_py():
        path_ = path.join("tests", "python")
        file_ = path.join(path_, "test03", "test03.json")
        assert py_test(*get_data(file_))

    def test_integration_python04_py():
        path_ = path.join("tests", "python")
        file_ = path.join(path_, "test04", "test04.json")
        assert py_test(*get_data(file_))
