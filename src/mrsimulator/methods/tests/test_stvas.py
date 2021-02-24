# -*- coding: utf-8 -*-
import numpy as np
import pytest
from mrsimulator.method import Method
from mrsimulator.method.transition_query import TransitionQuery
from mrsimulator.methods import ST1_VAS
from mrsimulator.methods import ST2_VAS

__author__ = "Deepansh J. Srivastava"
__email__ = "srivastava.89@osu.edu"

methods = [ST1_VAS, ST2_VAS]
names = ["ST1_VAS", "ST2_VAS"]


def test_ST_VAS_rotor_freq():
    def error(name):
        return f"`rotor_frequency` value cannot be modified for {name} method."

    for name, method in zip(names, methods):
        e = error(name)
        with pytest.raises(ValueError, match=f".*{e}.*"):
            method(rotor_frequency=10, spectral_dimensions=[{}, {}])


def test_ST_VAS_spectral_dimension_count():
    e = "Method requires exactly 2 spectral dimensions, given 1."
    for name, method in zip(names, methods):
        with pytest.raises(ValueError, match=f".*{e}.*"):
            method(spectral_dimensions=[{}])


def test_ST_VAS_setting_transition_query():
    def error(name):
        return f"`transition_query` value cannot be modified for {name} method."

    for name, method in zip(names, methods):
        e = error(name)
        with pytest.raises(ValueError, match=f".*{e}.*"):
            method(
                spectral_dimensions=[
                    {"events": [{"transition_query": {"P": [-1]}}]},
                    {},
                ],
            )


def test_ST_VAS_affine():
    sites = ["87Rb", "27Al"]
    spins = [1.5, 2.5]
    k_ST_MAS = {
        3: {1.5: 24 / 27, 2.5: 21 / 72, 3.5: 84 / 135, 4.5: 165 / 216},
        5: {2.5: 132 / 72, 3.5: 69 / 135, 4.5: 12 / 216},
        7: {3.5: 324 / 135, 4.5: 243 / 216},
        9: {4.5: 600 / 216},
    }
    for j, method in enumerate(methods):
        for i, isotope in zip(spins[j:], sites[j:]):
            meth = method(channels=[isotope])
            k = k_ST_MAS[3 + 2 * j][i]
            assert meth.spectral_dimensions[0].events[0].fraction == 1
            assert meth.spectral_dimensions[1].events[0].fraction == 1
            assert np.allclose(meth.affine_matrix, [1 / (1 + k), k / (1 + k), 0, 1])


def test_ST1_VAS_general():
    """Inner satellite-transition variable-angle spinning method"""
    mth = ST1_VAS(
        channels=["87Rb"],
        magnetic_flux_density=9.4,  # in T
        spectral_dimensions=[
            {
                "count": 1024,
                "spectral_width": 5e4,  # in Hz
                "reference_offset": 0,  # in Hz
            },
            {
                "count": 1024,
                "spectral_width": 5e4,  # in Hz
                "reference_offset": 0,  # in Hz
            },
        ],
    )
    assert mth.name == "ST1_VAS"

    des = (
        "Simulate a 1.5 -> 0.5 and -0.5 -> -1.5 satellite-transition variable-angle "
        "spinning spectrum."
    )
    assert mth.description == des
    assert mth.spectral_dimensions[0].events[0].transition_query == TransitionQuery(
        P={"channel-1": [[-1]]}, D={"channel-1": [[2], [-2]]}
    )
    assert mth.spectral_dimensions[1].events[0].transition_query == TransitionQuery(
        P={"channel-1": [[-1]]}, D={"channel-1": [[0]]}
    )
    assert Method.parse_dict_with_units(mth.json()) == mth


def test_ST2_VAS_general():
    """Second to inner satellite-transition variable-angle spinning method"""
    mth = ST2_VAS(
        channels=["17O"],
        magnetic_flux_density=9.4,  # in T
        spectral_dimensions=[
            {
                "count": 1024,
                "spectral_width": 5e4,  # in Hz
                "reference_offset": 0,  # in Hz
            },
            {
                "count": 1024,
                "spectral_width": 5e4,  # in Hz
                "reference_offset": 0,  # in Hz
            },
        ],
    )
    assert mth.name == "ST2_VAS"

    des = (
        "Simulate a 2.5 -> 1.5 and -1.5 -> -2.5 satellite-transition variable-angle "
        "spinning spectrum."
    )
    assert mth.description == des
    assert mth.spectral_dimensions[0].events[0].transition_query == TransitionQuery(
        P={"channel-1": [[-1]]}, D={"channel-1": [[4], [-4]]}
    )
    assert mth.spectral_dimensions[1].events[0].transition_query == TransitionQuery(
        P={"channel-1": [[-1]]}, D={"channel-1": [[0]]}
    )
    assert Method.parse_dict_with_units(mth.json()) == mth
