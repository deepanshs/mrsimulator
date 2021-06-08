# -*- coding: utf-8 -*-
import csdmpy as cp
import nmrglue as ng
import numpy as np

__author__ = "Alexis McCarthy"
__email__ = "mccarthy.677@osu.edu"


def load(filename, spec="bruker"):
    """Load a vendor serialized file.

    Args:
        str filename: A file or a folder name where the vendor data is stored.
        str vendor: An enumeration containing the name of the vendor that
            serialized the file. The valid enumerations are
            ['bruker', 'tecmag', 'agilent'].
    """
    enumerations = ["bruker", "tecmag", "agilent"]
    if spec not in enumerations:
        raise ValueError(f"Invalid `spec`. Valid enumerations are {enumerations}.")

    mod = getattr(ng, spec)
    dic, data = mod.read(filename)
    udic = mod.guess_udic(dic, data)
    dims = to_dimension_objects(udic)
    return cp.CSDM(
        dimensions=dims, dependent_variables=[cp.as_dependent_variable(data)]
    )


# def load_bruker(filename):
#     dic, data = ng.bruker.read(filename)
#     udic = ng.bruker.guess_udic(dic, data)
#     dims = to_dimension_objects(udic)
#     # dims = [
#     #     convert_to_dimension(value)
#     #     for key, value in list(udic.items())
#     #     if type(key) == int
#     # ]
#     return cp.CSDM(
#         dimensions=dims, dependent_variables=[cp.as_dependent_variable(data)]
#     )


def load_tecmag(filename):
    dic, data = ng.tecmag.read(filename)
    udic = ng.tecmag.guess_udic(dic, data)
    car_list = list(udic[0]["car"])
    obs_list = list(udic[0]["obs"])
    for i in range(len(udic) - 1):
        if type(udic[i]) == dict:
            udic[i]["car"] = car_list[i]
            udic[i]["obs"] = obs_list[i]
    dims = to_dimension_objects(udic)
    # dims = [
    #     convert_to_dimension(value)
    #     for key, value in list(udic.items())
    #     if type(key) == int and value["size"] != 1
    # ]
    return cp.CSDM(
        dimensions=dims,
        dependent_variables=[cp.as_dependent_variable((np.squeeze(data)).copy())],
    )


# def load_agilent(filename):
#     dic, data = ng.agilent.read(filename)
#     udic = ng.agilent.guess_udic(dic, data)
#     dims = to_dimension_objects(udic)
#     # dims = [
#     #     convert_to_dimension(value)
#     #     for key, value in list(udic.items())
#     #     if type(key) == int
#     # ]
#     return cp.CSDM(
#         dimensions=dims, dependent_variables=[cp.as_dependent_variable(data)]
#     )


def convert_to_dimension(udic):
    return cp.LinearDimension(
        count=udic["size"],
        increment=(f'{(1 / udic["sw"])} s'),
        reciprocal={
            "coordinates_offset": f'{udic["car"]} Hz',
            "origin_offset": f'{udic["obs"]} MHz',
        },
        label=udic["label"],
    )


def to_dimension_objects(udic):
    obj = [udic[i] for i in range(udic["ndim"])]
    return [convert_to_dimension(value) for value in obj if value["size"] != 1]


# Things to consider.
def auto_fix(csdm_object):
    """Auto correct data, if possible.

    Operations include
    1: origin correction
    2: zeroth order phase correction.  # To do
    3: First order phase correction.  # To do

    Args:
        csdm_object: CSDM object
    """
    correct_origin(csdm_object)


def correct_origin(csdm_object):
    """Correct the origin of the reference frame bu setting origin as the coordinates
    corresponding to the maximum amplitude.

    Args:
        csdm_object: CSDM object
    """
    dat_ = csdm_object.y[0].components[0]
    index = np.where(dat_ == dat_.max())[0]
    _ = [
        setattr(dim, "coordinates_offset", -dim.coordinates[i])
        for dim, i in zip(csdm_object.x, index)
    ]


def zeroth_order_correction(csdm_object):
    pass
