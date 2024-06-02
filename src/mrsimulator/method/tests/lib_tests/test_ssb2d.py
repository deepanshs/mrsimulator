import numpy as np
import pytest
from mrsimulator.method.lib import SSB2D
from mrsimulator.method.query import TransitionQuery

__author__ = "Deepansh J. Srivastava"
__email__ = "srivastava.89@osu.edu"


def test_SSB_rotor_freq():
    e = "`rotor_frequency` cannot be zero for SSB2D method."
    with pytest.raises(ValueError, match=f".*{e}.*"):
        SSB2D(channels=["1H"], spectral_dimensions=[{}, {}])

    with pytest.raises(ValueError, match=f".*{e}.*"):
        SSB2D(channels=["1H"], rotor_frequency=0, spectral_dimensions=[{}, {}])


def test_SSB_affine():
    mth = SSB2D(channels=["13C"], rotor_frequency=1200)
    assert np.allclose(mth.affine_matrix, [1, -1, 0, 1])
    assert SSB2D.parse_dict_with_units(mth.json()) == mth


def test_SSB_general():
    """Inner satellite-transition variable-angle spinning method"""
    mth = SSB2D(
        channels=["87Rb"],
        magnetic_flux_density=9.4,  # in T
        rotor_frequency=1200,
        spectral_dimensions=[
            {"count": 1024, "spectral_width": 5e4},
            {"count": 1024, "spectral_width": 5e4},
        ],
    )
    assert mth.name == "SSB2D"

    assert mth.description == "Simulate a 2D sideband separation method."

    # test transition query
    tq = TransitionQuery(ch1={"P": [-1], "D": [0]})
    assert mth.spectral_dimensions[0].events[0].transition_queries[0] == tq
    assert mth.spectral_dimensions[1].events[0].transition_queries[0] == tq

    # test rotor_frequency
    assert mth.spectral_dimensions[0].events[0].rotor_frequency == 1200
    assert mth.spectral_dimensions[1].events[0].rotor_frequency == 1e12

    # check serialization
    assert SSB2D.parse_dict_with_units(mth.json()) == mth

    assert np.allclose(mth.affine_matrix, [1, -1, 0.0, 1.0])

    serialize = mth.json()
    _ = serialize.pop("affine_matrix")

    assert serialize == {
        "channels": ["87Rb"],
        "description": mth.description,
        "magnetic_flux_density": "9.4 T",
        "name": "SSB2D",
        "rotor_angle": "0.9553166181245 rad",
        "rotor_frequency": "1200.0 Hz",
        "spectral_dimensions": [
            {
                "count": 1024,
                "spectral_width": "50000.0 Hz",
                "events": [{"transition_queries": [{"ch1": {"P": [-1], "D": [0]}}]}],
            },
            {
                "count": 1024,
                "events": [
                    {
                        "rotor_frequency": "1000000000000.0 Hz",
                        "transition_queries": [{"ch1": {"P": [-1], "D": [0]}}],
                    }
                ],
                "spectral_width": "50000.0 Hz",
            },
        ],
    }
