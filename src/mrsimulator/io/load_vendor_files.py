# -*- coding: utf-8 -*-
import csdmpy as cp
import nmrglue as ng
import numpy as np

__author__ = "Alexis McCarthy"
__email__ = "mccarthy.677@osu.edu"


def load_bruker(filename):
    dic, data = ng.bruker.read(filename)
    udic = ng.bruker.guess_udic(dic, data)
    dims = [
        convert_to_dimension(value)
        for key, value in list(udic.items())
        if type(key) == int
    ]
    return cp.CSDM(
        dimensions=dims, dependent_variables=[cp.as_dependent_variable(data)]
    )


def load_tecmag(filename):
    dic, data = ng.tecmag.read(filename)
    udic = ng.tecmag.guess_udic(dic, data)
    car_list = list(udic[0]["car"])
    obs_list = list(udic[0]["obs"])
    for i in range(len(udic) - 1):
        if type(udic[i]) == dict:
            udic[i]["car"] = car_list[i]
            udic[i]["obs"] = obs_list[i]
    dims = [
        convert_to_dimension(value)
        for key, value in list(udic.items())
        if type(key) == int and value["size"] != 1
    ]
    return cp.CSDM(
        dimensions=dims,
        dependent_variables=[cp.as_dependent_variable((np.squeeze(data)).copy())],
    )


def convert_to_dimension(udic):
    return cp.LinearDimension(
        count=udic["size"],
        increment=(str(1 / udic["sw"]) + " s"),
        reciprocal={
            "coordinates_offset": f'{udic["car"]} Hz',
            "origin_offset": f'{udic["obs"]} MHz',
        },
        label=udic["label"],
    )
