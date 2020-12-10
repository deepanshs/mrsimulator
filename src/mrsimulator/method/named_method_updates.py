# -*- coding: utf-8 -*-
__author__ = "Deepansh J. Srivastava"
__email__ = "srivastava.89@osu.edu"

named_methods = [
    "BlochDecaySpectrum",
    "BlochDecayCentralTransitionSpectrum",
    "ThreeQ_VAS",
    "FiveQ_VAS",
    "SevenQ_VAS",
    "ST1_VAS",
    "ST2_VAS",
    "SSB2D",
]

shear_factor_MQ_MAS = {
    3: {1.5: 21 / 27, 2.5: 114 / 72, 3.5: 303 / 135, 4.5: 546 / 216},
    5: {2.5: 150 / 72, 3.5: 165 / 135, 4.5: 570 / 216},
    7: {3.5: 483 / 135, 4.5: 84 / 216},
    9: {4.5: 1116 / 216},
}
MQ_p_symmetry = {
    "ThreeQ_VAS": {"mq": 1.5},
    "FiveQ_VAS": {"mq": 2.5},
    "SevenQ_VAS": {"mq": 3.5},
}
shear_factor_ST_MAS = {
    3: {1.5: 24 / 27, 2.5: 21 / 72, 3.5: 84 / 135, 4.5: 165 / 216},
    5: {2.5: 132 / 72, 3.5: 69 / 135, 4.5: 12 / 216},
    7: {3.5: 324 / 135, 4.5: 243 / 216},
    9: {4.5: 600 / 216},
}
ST_p_symmetry = {"ST1_VAS": {"st": 1.5}, "ST2_VAS": {"st": 2.5}}


def Bloch_decay_update(method):
    method.spectral_dimensions[0].events[0].transition_query.P = {"channel-1": [[-1]]}
    return method


def Bloch_decay_CT_update(method):
    method.spectral_dimensions[0].events[0].transition_query.P = {"channel-1": [[-1]]}
    method.spectral_dimensions[0].events[0].transition_query.D = {"channel-1": [[0]]}

    return method


def SSB2D_update(method):
    for dim in method.spectral_dimensions:
        dim.events[0].transition_query.P = {"channel-1": [[-1]]}
        dim.events[0].transition_query.D = {"channel-1": [[0]]}

    # Add the affine matrix
    if method.affine_matrix is None:
        method.affine_matrix = [1, -1, 0, 1]
    method.description = "Simulate a 2D sideband separation method."

    return method


def MQ_VAS_update(method):
    mq = MQ_p_symmetry[method.name]["mq"]
    spin = method.channels[0].spin

    # select the coherence for the first event
    P = int(2 * mq)
    nQ = P
    P = -P if mq == spin else P

    # setting transition symmetry elements
    method.spectral_dimensions[0].events[0].transition_query.P = {"channel-1": [[P]]}
    method.spectral_dimensions[0].events[0].transition_query.D = {"channel-1": [[0]]}

    method.spectral_dimensions[1].events[0].transition_query.P = {"channel-1": [[-1]]}
    method.spectral_dimensions[1].events[0].transition_query.D = {"channel-1": [[0]]}

    # Add the affine matrix
    if method.affine_matrix is None:
        k = shear_factor_MQ_MAS[nQ][spin]
        method.affine_matrix = [1 / (1 + k), k / (1 + k), 0, 1]

    method.description = f"Simulate a {nQ}Q variable-angle spinning spectrum."
    return method


def ST_VAS_update(method):
    st = ST_p_symmetry[method.name]["st"]
    spin = method.channels[0].spin

    # select the coherence for the first event
    d = st ** 2 - (st - 1) ** 2
    D = [[d], [-d]]

    # setting transition symmetry elements
    method.spectral_dimensions[0].events[0].transition_query.P = {"channel-1": [[-1]]}
    method.spectral_dimensions[0].events[0].transition_query.D = {"channel-1": D}

    method.spectral_dimensions[1].events[0].transition_query.P = {"channel-1": [[-1]]}
    method.spectral_dimensions[1].events[0].transition_query.D = {"channel-1": [[0]]}

    # Add the affine matrix
    if method.affine_matrix is None:
        k = shear_factor_ST_MAS[int(2 * st)][spin]
        method.affine_matrix = [1 / (1 + k), k / (1 + k), 0, 1]

    method.description = (
        f"Simulate a {st} -> {st-1} and {-st+1} -> {-st} satellite-transition "
        "variable-angle spinning spectrum."
    )
    return method


# Generic method update
def update_method(method):
    if method.name == "BlochDecaySpectrum":
        return Bloch_decay_update(method)
    if method.name == "BlochDecayCentralTransitionSpectrum":
        return Bloch_decay_CT_update(method)
    if method.name in MQ_p_symmetry.keys():
        return MQ_VAS_update(method)
    if method.name in ST_p_symmetry.keys():
        return ST_VAS_update(method)
    if method.name == "SSB2D":
        return SSB2D_update(method)
    return method
