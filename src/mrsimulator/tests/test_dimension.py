import numpy as np
from mrsimulator import Dimension


def test_direct_init_spectrum():
    the_spectrum = Dimension(number_of_points=1024, spectral_width=100)
    assert the_spectrum.number_of_points == 1024
    assert the_spectrum.spectral_width == 100
    # assert the_spectrum.property_units["spectral_width"] == 'Hz'

    Dimension(
        number_of_points=1024,
        spectral_width=100,
        reference_offset=0,
        magnetic_flux_density=9.4,
        rotor_frequency=0,
        rotor_angle=0.9553,  # 54.935 degrees in radians
        rotor_phase=0,
        isotope="1H",
    )


def test_parse_json_spectrum():
    good_json = {
        "number_of_points": 1024,
        "spectral_width": "100 Hz",
        "reference_offset": "0 Hz",
        "magnetic_flux_density": "9.4 T",
        "rotor_frequency": "0 Hz",
        "rotor_angle": "54.935 degree",  # 54.935 degrees in radians
        "rotor_phase": "0 rad",
        "isotope": "1H",
    }

    spec = Dimension.parse_dict_with_units(good_json)
    assert spec.spin == 1
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
    assert spec.spin == 1
    assert spec.isotope == "29Si"
    assert spec.gyromagnetic_ratio == -8.465499
    assert spec.natural_abundance == 4.683
    assert spec.quadrupole_moment == 0.0
    assert spec.atomic_number == 14
    assert np.allclose(spec.rotor_angle, 0.95879662)
