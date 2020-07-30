# -*- coding: utf-8 -*-
import numpy as np
from mrsimulator.methods import BlochDecaySpectrum


def test_BlochDecaySpectrum():
    # test-1
    m1 = BlochDecaySpectrum()

    event_dictionary_ = {
        "fraction": 1.0,
        "magnetic_flux_density": "9.4 T",
        "rotor_frequency": "0.0 Hz",
        "rotor_angle": "0.9553166 rad",
        "transition_query": {"P": {"channel-1": [[-1.0]]}},
        "user_variables": ["magnetic_flux_density", "rotor_frequency", "rotor_angle"],
    }
    dimension_dictionary_ = {
        "count": 1024,
        "spectral_width": "25000.0 Hz",
        "reference_offset": "0.0 Hz",
        "events": [event_dictionary_],
    }

    should_be = {
        "name": "BlochDecaySpectrum",
        "description": "Simulate a 1D Bloch decay spectrum.",
        "channels": ["1H"],
        "spectral_dimensions": [dimension_dictionary_],
    }
    assert m1.to_dict_with_units() == should_be

    # test-2
    m2_dict = {
        "channels": ["29Si"],
        "magnetic_flux_density": "11.7 T",
        "rotor_angle": "90 deg",
        "spectral_dimensions": [{}],
    }
    m2 = BlochDecaySpectrum.parse_dict_with_units(m2_dict)

    angle = 90 * np.pi / 180
    event_dictionary_ = {
        "fraction": 1.0,
        "magnetic_flux_density": "11.7 T",
        "rotor_frequency": "0.0 Hz",
        "rotor_angle": f"{angle} rad",
        "transition_query": {"P": {"channel-1": [[-1.0]]}},
        "user_variables": ["magnetic_flux_density", "rotor_frequency", "rotor_angle"],
    }
    dimension_dictionary_ = {
        "count": 1024,
        "spectral_width": "25000.0 Hz",
        "reference_offset": "0.0 Hz",
        "events": [event_dictionary_],
    }

    should_be = {
        "name": "BlochDecaySpectrum",
        "description": "Simulate a 1D Bloch decay spectrum.",
        "channels": ["29Si"],
        "spectral_dimensions": [dimension_dictionary_],
    }

    assert m2.to_dict_with_units() == should_be
