# -*- coding: utf-8 -*-
"""Test for the base ConfigSimulator class."""
import pytest
from mrsimulator import Simulator

__author__ = "Deepansh Srivastava"
__email__ = "srivastava.89@osu.edu"


def test_config():

    # set config
    b = Simulator()

    error = "value is not a valid dict"
    with pytest.raises(ValueError, match=f".*{error}.*"):
        b.config = 5

    a = Simulator()

    assert b == a

    # number of sidebands
    assert a.config.number_of_sidebands == 64

    a.config.number_of_sidebands = 1.4
    assert a.config.number_of_sidebands == 1

    a.config.number_of_sidebands = 10
    assert a.config.number_of_sidebands == 10

    error = "ensure this value is greater than 0"
    with pytest.raises(ValueError, match=f".*{error}.*"):
        a.config.number_of_sidebands = 0

    # integration density
    assert a.config.integration_density == 70
    a.config.integration_density = 20
    assert a.config.integration_density == 20

    error = "ensure this value is greater than 0"
    with pytest.raises(ValueError, match=f".*{error}.*"):
        a.config.integration_density = -12

    error = "value is not a valid integer"
    with pytest.raises(ValueError, match=f".*{error}.*"):
        a.config.integration_density = {}

    # integration volume
    assert a.config.integration_volume == "octant"
    a.config.integration_volume = "hemisphere"
    assert a.config.integration_volume == "hemisphere"

    error = "unexpected value; permitted: 'octant', 'hemisphere'"
    with pytest.raises(ValueError, match=f".*{error}.*"):
        a.config.integration_volume = "sphere"

    # decompose spectrum
    assert a.config.decompose_spectrum == "none"
    a.config.decompose_spectrum = "spin_system"
    assert a.config.decompose_spectrum == "spin_system"

    error = "unexpected value; permitted: 'none', 'spin_system'"
    with pytest.raises(ValueError, match=f".*{error}.*"):
        a.config.decompose_spectrum = "haha"

    # overall
    assert a.config.dict() == {
        "decompose_spectrum": "spin_system",
        "number_of_sidebands": 10,
        "integration_volume": "hemisphere",
        "integration_density": 20,
    }

    assert a.config.get_int_dict() == {
        "decompose_spectrum": 1,
        "number_of_sidebands": 10,
        "integration_volume": 1,
        "integration_density": 20,
    }

    assert b != a

    # get orientation count
    assert a.config.get_orientations_count() == 4 * 21 * 22 / 2
