# -*- coding: utf-8 -*-
import json
from os import path

import numpy as np
from mrsimulator import Method
from mrsimulator import signal_processor as sp
from mrsimulator import Simulator
from mrsimulator import SpinSystem
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

    # if test_data_object["source"] == "dmfit":
    #     data_source = data_source[::-1]
    #     data_source = np.roll(data_source, 1)

    return data_object, data_source


def simulator_setup(
    data_object,
    integration_volume="octant",
    integration_density=120,
    number_of_sidebands=90,
):
    methods = [Method.parse_dict_with_units(_) for _ in data_object["methods"]]

    spin_systems = [
        SpinSystem.parse_dict_with_units(item) for item in data_object["spin_systems"]
    ]
    s1 = Simulator(spin_systems=spin_systems, methods=methods)
    s1.config.decompose_spectrum = "spin_system"
    s1.spin_systems[0].name = "test name"
    s1.spin_systems[0].description = "test description"
    s1.config.integration_density = integration_density
    s1.config.number_of_sidebands = number_of_sidebands
    s1.config.integration_volume = integration_volume

    return s1


def simulator_process(sim, data_object):
    sim.run()

    sim_dataset = sim.methods[0].simulation

    if "operations" in data_object:
        processor = sp.SignalProcessor.parse_dict_with_units(
            {"operations": data_object["operations"]}
        )
        sim_dataset = processor.apply_operations(dataset=sim_dataset)

    data_mrsimulator = np.asarray(sim_dataset.to_list()[1:])
    data_mrsimulator = data_mrsimulator.sum(axis=0)
    data_mrsimulator /= data_mrsimulator.sum()

    dv = sim_dataset.y[0]
    assert dv.name == "test name"
    assert dv.description == "test description"

    return data_mrsimulator


def c_setup(
    filename,
    integration_volume="octant",
    integration_density=120,
    number_of_sidebands=90,
):
    # mrsimulator
    data_object, data_source = get_data(filename)
    data_source /= data_source.sum()

    sim = simulator_setup(
        data_object, integration_volume, integration_density, number_of_sidebands
    )
    data_mrsimulator = simulator_process(sim, data_object)
    return data_mrsimulator, data_source


def c_setup_random_euler_angles(filename, group):
    # mrsimulator
    data_object, data_source = get_data(filename)
    data_source /= data_source.sum()

    sim = simulator_setup(data_object, integration_volume="hemisphere")
    pix2 = 2 * np.pi
    if group == "shielding_symmetric":
        for spin_system in sim.spin_systems:
            spin_system.sites[0].shielding_symmetric.alpha = np.random.rand(1) * pix2
            spin_system.sites[0].shielding_symmetric.beta = np.random.rand(1) * pix2
            spin_system.sites[0].shielding_symmetric.gamma = np.random.rand(1) * pix2

    if group == "quadrupolar":
        for spin_system in sim.spin_systems:
            spin_system.sites[0].quadrupolar.alpha = np.random.rand(1) * pix2
            spin_system.sites[0].quadrupolar.beta = np.random.rand(1) * pix2
            spin_system.sites[0].quadrupolar.gamma = np.random.rand(1) * pix2

    data_mrsimulator = simulator_process(sim, data_object)
    return data_mrsimulator, data_source
