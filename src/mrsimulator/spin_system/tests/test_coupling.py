"""Test for the base Site class."""
import pytest
from mrsimulator import Coupling
from pydantic.v1 import ValidationError

__author__ = "Deepansh Srivastava"
__email__ = "srivastava.89@osu.edu"


# Test Site ===========================================================================
def test_direct_init_coupling1():
    # test 1 --------------------------------------------------------------------------
    the_coupling = Coupling(site_index=[0, 1], isotropic_j=10)
    assert the_coupling.site_index == [0, 1]
    assert the_coupling.isotropic_j == 10
    assert the_coupling.property_units["isotropic_j"] == "Hz"

    assert the_coupling.j_symmetric is None
    assert the_coupling.j_antisymmetric is None
    assert the_coupling.dipolar is None

    # test 2 --------------------------------------------------------------------------
    the_coupling = Coupling(site_index=[0, 1], isotropic_j=10, j_symmetric=None)
    assert the_coupling.j_symmetric is None
    assert the_coupling.j_antisymmetric is None
    assert the_coupling.dipolar is None

    # test 3 --------------------------------------------------------------------------
    the_coupling = Coupling(
        site_index=[2, 3],
        isotropic_j=10,
        j_symmetric={"zeta": 12.1, "eta": 0.1},
    )
    assert the_coupling.site_index == [2, 3]
    assert the_coupling.isotropic_j == 10
    assert the_coupling.property_units["isotropic_j"] == "Hz"

    assert the_coupling.j_antisymmetric is None

    assert the_coupling.dipolar is None

    assert the_coupling.j_symmetric.zeta == 12.1
    assert the_coupling.j_symmetric.property_units["zeta"] == "Hz"
    assert the_coupling.j_symmetric.eta == 0.1

    # test 4 --------------------------------------------------------------------------
    error = "ensure this value is less than or equal to 1"
    with pytest.raises(ValidationError, match=f".*{error}.*"):
        Coupling(
            site_index=[3, 1],
            isotropic_j=10,
            j_symmetric={"zeta": 12.1, "eta": 1.5},
        )

    error = [
        "Site index must a list of two integers",
        "The two site indexes must be unique integers",
    ]
    with pytest.raises(ValidationError, match=".*{}.*".format(*error)):
        Coupling(site_index=[])

    with pytest.raises(ValidationError, match=".*{}.*".format(*error)):
        Coupling(site_index=[2, 3, 4])

    with pytest.raises(ValidationError, match=".*{}.*".format(*error)):
        Coupling(site_index=[2])

    with pytest.raises(ValidationError, match=".*{1}.*".format(*error)):
        Coupling(site_index=[1, 1])

    ax = Coupling.parse_dict_with_units({"site_index": [0, 1], "isotropic_j": "5 cHz"})
    assert ax.json() == {
        "site_index": [0, 1],
        "isotropic_j": "0.05 Hz",
    }


def test_parse_json_site():
    # test 1 --------------------------------------------------------------------------
    good_json_coupling = {
        "site_index": [0, 1],
        "isotropic_j": "5 Hz",
        "j_symmetric": {"zeta": "13.89 Hz", "eta": 0.25},
    }

    # testing method parse_dict_with_units()
    the_coupling = Coupling.parse_dict_with_units(good_json_coupling)
    assert the_coupling.site_index == [0, 1]
    assert the_coupling.isotropic_j == 5
    assert the_coupling.property_units["isotropic_j"] == "Hz"

    assert the_coupling.j_antisymmetric is None
    assert the_coupling.dipolar is None

    assert the_coupling.j_symmetric.zeta == 13.89
    assert the_coupling.j_symmetric.property_units["zeta"] == "Hz"
    assert the_coupling.j_symmetric.eta == 0.25
    assert the_coupling.j_symmetric.alpha is None
    assert the_coupling.j_symmetric.beta is None
    assert the_coupling.j_symmetric.gamma is None

    # test 2 --------------------------------------------------------------------------
    good_json_2 = {
        "site_index": [5, 3],
        "isotropic_j": "-10 Hz",
        "dipolar": {"D": "5.12 kHz", "eta": 0.5},
    }

    the_coupling = Coupling.parse_dict_with_units(good_json_2)
    assert the_coupling.site_index == [5, 3]
    assert the_coupling.isotropic_j == -10.0
    assert the_coupling.property_units["isotropic_j"] == "Hz"

    assert the_coupling.j_antisymmetric is None
    assert the_coupling.j_symmetric is None

    assert the_coupling.dipolar.D == 5120.0
    assert the_coupling.dipolar.eta == 0.5

    # test 3 bad input ----------------------------------------------------------------
    bad_json = {"zite_index": [0, 1], "isotropic_j": "5 ppm"}

    error = "Error enforcing units for isotropic_j: 5 ppm"
    with pytest.raises(Exception, match=f".*{error}.*"):
        Coupling.parse_dict_with_units(bad_json)

    error = "Error enforcing units for j_symmetric.zeta: ppm."
    with pytest.raises(ValueError, match=f".*{error}.*"):
        Coupling.parse_dict_with_units(
            {
                "site_index": [0, 1],
                "isotropic_j": "10 kHz",
                "j_symmetric": {"zeta": "12.1 ppm", "eta": 0.5},
            }
        )


def test_site_object_methods():
    good_json_2 = {"site_index": [0, 1], "isotropic_j": "-10 Hz"}
    the_coupling = Coupling.parse_dict_with_units(good_json_2)

    # testing method dict()
    result = {
        "site_index": [0, 1],
        "isotropic_j": -10.0,
        "property_units": {"isotropic_j": "Hz"},
        "name": None,
        "label": None,
        "description": None,
        "dipolar": None,
        "j_symmetric": None,
        "j_antisymmetric": None,
    }
    assert the_coupling.dict() == result, "Failed Coupling.dict()"

    # testing method json()
    result = {"site_index": [0, 1], "isotropic_j": "-10.0 Hz"}
    assert the_coupling.json() == result, "Failed Coupling.json()"

    result = {
        "site_index": [1, 2],
        "isotropic_j": "10.0 Hz",
        "j_symmetric": {"zeta": "12.1 Hz", "eta": 0.1, "alpha": "2.1 rad"},
        "j_antisymmetric": {
            "zeta": "-1.1 Hz",
            "alpha": "0.1 rad",
            "beta": "2.5 rad",
        },
        "dipolar": {"D": "10000.0 Hz", "eta": 0.6},
    }
    the_coupling = Coupling(
        site_index=[1, 2],
        isotropic_j=10,
        j_symmetric={"zeta": 12.1, "eta": 0.1, "alpha": 2.1},
        j_antisymmetric={"zeta": -1.1, "alpha": 0.1, "beta": 2.5},
        dipolar={"D": 10e3, "eta": 0.6},
    )
    assert the_coupling.json() == result, "Failed Coupling.json()"


def test_site_dipolar_set_to_None():
    a = Coupling(site_index=[0, 1], dipolar=None)
    assert a.site_index == [0, 1]
    assert a.isotropic_j == 0
    assert a.dipolar is None


def test_equality():
    a = Coupling(site_index=[0, 1])
    b = Coupling(site_index=[0, 1])
    assert a == b

    c = Coupling(site_index=[0, 1], isotropic_j=16)
    assert a != c
