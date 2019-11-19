# -*- coding: utf-8 -*-
"""Test for the base ConfigSimulator class."""
import pytest
from mrsimulator import Simulator


def test_config():

    # set config
    b = Simulator()
    error = "instance of ConfigSimulator expected"
    with pytest.raises(ValueError, match=".*{0}.*".format(error)):
        b.config = {}

    a = Simulator()

    assert b == a

    # number of sidebands
    assert a.config.number_of_sidebands == 64
    a.config.number_of_sidebands = 10
    assert a.config.number_of_sidebands == 10

    error = "Expecting a positive integer, found"
    with pytest.raises(ValueError, match=".*{0}.*".format(error)):
        a.config.number_of_sidebands = "2"
    with pytest.raises(ValueError, match=".*{0}.*".format(error)):
        a.config.number_of_sidebands = 1.4

    # integration density
    assert a.config.integration_density == 70
    a.config.integration_density = 20
    assert a.config.integration_density == 20

    error = "Expecting a positive integer, found"
    with pytest.raises(ValueError, match=".*{0}.*".format(error)):
        a.config.integration_density = -12
    with pytest.raises(ValueError, match=".*{0}.*".format(error)):
        a.config.integration_density = {}

    # integration volume
    assert a.config.integration_volume == "octant"
    a.config.integration_volume = "hemisphere"
    assert a.config.integration_volume == "hemisphere"

    error = (
        "value is not a valid enumeration literal; permitted: 'octant', 'hemisphere'"
    )
    with pytest.raises(ValueError, match=".*{0}.*".format(error)):
        a.config.integration_volume = "sphere"

    # decompose
    assert a.config.decompose is False
    a.config.decompose = True
    assert a.config.decompose is True

    error = "Expecting a boolean."
    with pytest.raises(ValueError, match=".*{0}.*".format(error)):
        a.config.decompose = [5, 23]

    assert a.config.dict() == {
        "decompose": True,
        "number_of_sidebands": 10,
        "integration_volume": "hemisphere",
        "integration_density": 20,
    }

    assert a.config._dict == {
        "decompose": True,
        "number_of_sidebands": 10,
        "integration_volume": 1,
        "integration_density": 20,
    }

    assert str(a.config) == str(
        {
            "number_of_sidebands": 10,
            "integration_volume": "hemisphere",
            "integration_density": 20,
            "decompose": True,
        }
    )
    assert b != a
