import numpy as np
import pytest
from mrsimulator.method.lib import ST1_VAS
from mrsimulator.method.lib import ST2_VAS
from mrsimulator.method.query import TransitionQuery

__author__ = "Deepansh J. Srivastava"
__email__ = "srivastava.89@osu.edu"

methods = [ST1_VAS, ST2_VAS]
names = ["ST1_VAS", "ST2_VAS"]


def sample_test_output(n):
    return {
        "magnetic_flux_density": "9.4 T",
        "rotor_angle": "0.9553166181245 rad",
        "rotor_frequency": "1000000000000.0 Hz",
        "spectral_dimensions": [
            {
                "count": 1024,
                "spectral_width": "50000.0 Hz",
                "events": [
                    {
                        "transition_queries": [
                            {"ch1": {"P": [-1], "D": [i]}} for i in [n, -n]
                        ]
                    }
                ],
            },
            {
                "count": 1024,
                "spectral_width": "50000.0 Hz",
                "events": [{"transition_queries": [{"ch1": {"P": [-1], "D": [0]}}]}],
            },
        ],
    }


def test_ST_VAS_rotor_freq():
    e = "`rotor_frequency=1e12 Hz` is fixed for all 2D named Methods,"
    isotope = ["87Rb", "27Al"]
    for iso, method in zip(isotope, methods):
        with pytest.raises(ValueError, match=f".*{e}.*"):
            method(channels=[iso], rotor_frequency=10, spectral_dimensions=[{}, {}])


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
        magnetic_flux_density=11.7,  # in T
        spectral_dimensions=[
            {"count": 1024, "spectral_width": 3e4},
            {"count": 1024, "spectral_width": 2e4},
        ],
    )
    assert mth.name == "ST1_VAS"

    des = (
        "Simulate a P=-1 (ST) to P=-1 (CT) ST-CT variable-angle spinning "
        "correlation spectrum."
    )
    assert mth.description == des
    assert mth.spectral_dimensions[0].events[0].transition_queries == [
        TransitionQuery(ch1={"P": [-1], "D": [2]}),
        TransitionQuery(ch1={"P": [-1], "D": [-2]}),
    ]
    assert mth.spectral_dimensions[1].events[0].transition_queries == [
        TransitionQuery(ch1={"P": [-1], "D": [0]})
    ]
    assert ST1_VAS.parse_dict_with_units(mth.json()) == mth

    assert np.allclose(mth.affine_matrix, [0.52941176, 0.47058824, 0.0, 1.0])

    serialize = mth.json()
    _ = serialize.pop("affine_matrix")

    assert serialize == {
        "channels": [
            {
                "natural_abundance": 27.83,
                "gyromagnetic_ratio": 13.983992836343669,
                "quadrupole_moment": 0.127,
                "atomic_number": 37,
                "spin_multiplicity": 4,
                "symbol": "87Rb",
            }
        ],
        "description": des,
        "magnetic_flux_density": "11.7 T",
        "name": "ST1_VAS",
        "rotor_angle": "0.9553166181245 rad",
        "rotor_frequency": "1000000000000.0 Hz",
        "spectral_dimensions": [
            {
                "count": 1024,
                "spectral_width": "30000.0 Hz",
                "events": [
                    {
                        "transition_queries": [
                            {"ch1": {"P": [-1], "D": [i]}} for i in [2, -2]
                        ]
                    }
                ],
            },
            {
                "count": 1024,
                "spectral_width": "20000.0 Hz",
                "events": [{"transition_queries": [{"ch1": {"P": [-1], "D": [0]}}]}],
            },
        ],
    }


def test_ST2_VAS_general():
    """Second to inner satellite-transition variable-angle spinning method"""
    mth = ST2_VAS(
        channels=["17O"],
        magnetic_flux_density=9.4,  # in T
        spectral_dimensions=[
            {"count": 1024, "spectral_width": 5e4},
            {"count": 1024, "spectral_width": 5e4},
        ],
    )
    assert mth.name == "ST2_VAS"

    des = (
        "Simulate a P=-1 (ST) to P=-1 (CT) ST-CT variable-angle spinning "
        "correlation spectrum."
    )
    assert mth.description == des
    assert mth.spectral_dimensions[0].events[0].transition_queries == [
        TransitionQuery(ch1={"P": [-1], "D": [4]}),
        TransitionQuery(ch1={"P": [-1], "D": [-4]}),
    ]
    assert mth.spectral_dimensions[1].events[0].transition_queries == [
        TransitionQuery(ch1={"P": [-1], "D": [0]})
    ]
    assert ST2_VAS.parse_dict_with_units(mth.json()) == mth

    assert np.allclose(mth.affine_matrix, [0.35294118, 0.64705882, 0.0, 1.0])

    serialize = mth.json()
    _ = serialize.pop("affine_matrix")

    assert serialize == {
        "channels": [
            {
                "natural_abundance": 0.038,
                "gyromagnetic_ratio": -5.774236332534915,
                "quadrupole_moment": -0.02578,
                "atomic_number": 8,
                "spin_multiplicity": 6,
                "symbol": "17O",
            }
        ],
        "description": des,
        "name": "ST2_VAS",
        **sample_test_output(4),
    }
