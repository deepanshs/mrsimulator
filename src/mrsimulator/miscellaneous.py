# -*- coding: utf-8 -*-
import json
import os
import numpy as np

__author__ = "Deepansh J. Srivastava"
__email__ = ["srivastava.89@osu.edu", "deepansh2012@gmail.com"]


def import_json(filename):
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
    data = {}
    for line_ in f:
        if "NP" in line_:
            data["count"] = int(line_.split("=")[1])
        if "X0" in line_:
            data["zero_index_coordinate"] = float(line_.split("=")[1])
        if "SW" in line_:
            data["increment"] = float(line_.split("=")[1]) / data["count"]
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
    data["skip_header"] = skip_header
    data["skip_footer"] = skip_footer
    return data


def read_dmfit_files(filename):
    """Load a dmfit output file"""
    # source data
    # data_object = import_json(filename)

    # test_data_object = data_object["test_data"]

    # source_file = test_data_object["filename"]
    path, extension = os.path.split(filename)
    # path_ = path.split(path_)[0]
    # source_file = path.join(path_, filename)

    data = _get_header_and_footer(filename)
    print(data)
    x, y = np.genfromtxt(
        filename,
        skip_header=data["skip_header"],
        # delimiter='tab',
        skip_footer=data["skip_footer"],
        unpack=True,
    )
    print(x)
    data_source = x + 1j * y
    data_source = data_source[::-1]
    data_source = np.roll(data_source, 1)
    freq = np.arange(data["count"]) * data["increment"] - data["zero_index_coordinate"]
    return freq, data_source
