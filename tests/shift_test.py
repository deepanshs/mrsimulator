# -*- coding: utf-8 -*-
"""Test for shift and reference offset."""
import numpy as np
from mrsimulator import Isotopomer
from mrsimulator import Method
from mrsimulator import Simulator
from mrsimulator import Site


def pre_setup(isotope, shift, reference_offset):
    isotopomer = Isotopomer(
        sites=[Site(isotope=isotope, isotropic_chemical_shift=shift)]
    )
    method = Method.parse_dict_with_units(
        dict(
            channels=[isotope],
            spectral_dimensions=[
                {
                    "count": 2046,
                    "spectral_width": "25 kHz",
                    "reference_offset": f"{reference_offset} Hz",
                    "events": [{"magnetic_flux_density": "9.4 T"}],
                }
            ],
        )
    )
    sim = Simulator()
    sim.isotopomers.append(isotopomer)
    sim.methods += [method]
    sim.run()
    sim.methods[0].simulation.dimensions[0].to("ppm", "nmr_frequency_ratio")
    x, y = sim.methods[0].simulation.to_list()
    return x[np.argmax(y)], abs(x[1] - x[0])


# shift test for negative gyromagnetic ratio
def test_reference_offset_for_negative_gamma():
    shifts = [10, 0, -10, 21.34, 54.2]
    reference_offsets = [-1000, 0, 4231, -2301, 2352.124]
    for shift in shifts:
        for offset in reference_offsets:
            max_val, dx = pre_setup(
                isotope="29Si", shift=shift, reference_offset=offset
            )
            assert (
                abs(max_val.value - shift) <= dx.value
            ), f"failed at shift={shift}, offset={offset}"


# shift test for positive gyromagnetic ratio
def test_reference_offset_for_positive_gamma():
    shifts = [10, 0, -10, -8.34, 24.2]  # ppm
    reference_offsets = [-10, 0, 4231, -2301, 2352.124]  # kHz
    for shift in shifts:
        for offset in reference_offsets:
            max_val, dx = pre_setup(isotope="1H", shift=shift, reference_offset=offset)
            assert (
                abs(max_val.value - shift) <= dx.value
            ), f"failed at shift={shift}, offset={offset}"
