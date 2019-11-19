# -*- coding: utf-8 -*-
"""Test for the base Dimension class."""
import numpy as np
import pytest
from mrsimulator import Dimension
from pydantic import ValidationError


def test_direct_init_dimension():
    the_dimension = Dimension(number_of_points=1024, spectral_width=100)

    assert the_dimension.number_of_points == 1024
    error = "ensure this value is greater than 0"
    with pytest.raises(ValidationError, match=".*{0}.*".format(error)):
        the_dimension.number_of_points = -1024

    # spectral width test
    assert the_dimension.spectral_width == 100
    with pytest.raises(ValidationError, match=".*{0}.*".format(error)):
        the_dimension.spectral_width = 0
    # ensure the default value is Hz
    assert the_dimension.property_units["spectral_width"] == "Hz"

    assert the_dimension.spin is None
    assert the_dimension.natural_abundance is None
    assert the_dimension.gyromagnetic_ratio is None
    assert the_dimension.quadrupole_moment is None
    assert the_dimension.atomic_number is None
    assert the_dimension.larmor_frequency is None

    # rotor angle test
    assert the_dimension.rotor_angle == 0.9553166
    # ensure the unit of rotor angle is radian
    assert the_dimension.property_units["rotor_angle"] == "rad"
    the_dimension.rotor_angle = 1.5707963268
    assert the_dimension.rotor_angle == 1.5707963268
    # rotor angle cannot ve greater than pi/2 and less that 0
    error = "ensure this value is less than or equal to 1.5707963268"
    with pytest.raises(ValidationError, match=".*{0}.*".format(error)):
        the_dimension.rotor_angle = 2.57079
    error = "ensure this value is greater than or equal to 0"
    with pytest.raises(ValidationError, match=".*{0}.*".format(error)):
        the_dimension.rotor_angle = -0.1243

    # rotor frequency test
    assert the_dimension.rotor_frequency == 0
    the_dimension.rotor_frequency = 1000
    assert the_dimension.rotor_frequency == 1000
    assert the_dimension.property_units["rotor_frequency"] == "Hz"

    # magnetic flux density
    assert the_dimension.magnetic_flux_density == 9.4
    assert the_dimension.property_units["magnetic_flux_density"] == "T"
    error = "ensure this value is greater than or equal to 0"
    with pytest.raises(ValidationError, match=".*{0}.*".format(error)):
        the_dimension.magnetic_flux_density = -1

    # isotope
    assert the_dimension.isotope is None
    the_dimension.isotope = "13C"
    print(the_dimension)
    assert the_dimension.isotope == "13C"
    assert the_dimension.gyromagnetic_ratio == 10.7084
    assert the_dimension.spin == 0.5
    assert the_dimension.larmor_frequency == -10.7084 * 9.4 * 1e6

    coordinate_Hz = (np.arange(1024) - 512) * (100.0 / 1024)
    assert np.allclose(the_dimension.coordinates_Hz, coordinate_Hz)

    coordinate_ppm = coordinate_Hz / (10.7084 * 9.4)
    assert np.allclose(the_dimension.coordinates_ppm, coordinate_ppm)

    error = "field required"
    with pytest.raises(ValidationError, match=".*{0}.*".format(error)):
        Dimension()

    # when isotope is None
    the_dimension.isotope = None
    assert the_dimension.isotope is None
    assert the_dimension.spin is None
    assert the_dimension.natural_abundance is None
    assert the_dimension.gyromagnetic_ratio is None
    assert the_dimension.quadrupole_moment is None
    assert the_dimension.atomic_number is None
    assert the_dimension.larmor_frequency is None

    assert np.allclose(the_dimension.coordinates_Hz, coordinate_Hz)
    assert the_dimension.coordinates_ppm is None

    error = "Could not parse isotope string"
    with pytest.raises(Exception, match=".*{0}.*".format(error)):
        Dimension(isotope="monkey", spectral_width=10000)

    with pytest.raises(Exception, match=".*{0}.*".format(error)):
        Dimension(isotope="15H", spectral_width=10000)

    result = {
        "number_of_points": 1024,
        "spectral_width": "100.0 Hz",
        "reference_offset": "0 Hz",
        "magnetic_flux_density": "9.4 T",
        "rotor_frequency": "1000.0 Hz",
        "rotor_angle": "1.5707963268 rad",
        "label": "",
    }
    assert the_dimension.to_dict_with_units() == result

    the_dimension.isotope = "27Al"
    result = {
        "isotope": "27Al",
        "number_of_points": 1024,
        "spectral_width": "100.0 Hz",
        "reference_offset": "0 Hz",
        "magnetic_flux_density": "9.4 T",
        "rotor_frequency": "1000.0 Hz",
        "rotor_angle": "1.5707963268 rad",
        "label": "",
    }
    assert the_dimension.to_dict_with_units() == result


def test_parse_json_spectrum():
    good_json = {
        "number_of_points": 1024,
        "spectral_width": "100 Hz",
        "reference_offset": "0 Hz",
        "magnetic_flux_density": "9.4 T",
        "rotor_frequency": "0 Hz",
        "rotor_angle": "54.935 degree",  # 54.935 degrees in radians
        "isotope": "1H",
    }

    spec = Dimension.parse_dict_with_units(good_json)
    assert spec.spin == 0.5
    assert spec.isotope == "1H"
    assert spec.gyromagnetic_ratio == 42.57748
    assert spec.natural_abundance == 99.985
    assert spec.quadrupole_moment == 0.0
    assert spec.atomic_number == 1
    assert np.allclose(spec.rotor_angle, 0.95879662)


def test_parse_json_spectrum2():
    good_json = {
        "number_of_points": 1024,
        "spectral_width": "100 Hz",
        "reference_offset": "0 Hz",
        "magnetic_flux_density": "9.4 T",
        "rotor_frequency": "0 Hz",
        "rotor_angle": "54.935 degree",  # 54.935 degrees in radians
        "rotor_phase": "0 rad",
        "isotope": "29Si",
    }

    spec = Dimension.parse_dict_with_units(good_json)
    assert spec.spin == 0.5
    assert spec.isotope == "29Si"
    assert spec.gyromagnetic_ratio == -8.465499
    assert spec.natural_abundance == 4.683
    assert spec.quadrupole_moment == 0.0
    assert spec.atomic_number == 14
    assert np.allclose(spec.rotor_angle, 0.95879662)
