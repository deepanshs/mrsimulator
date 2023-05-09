"""Test for shift and reference offset."""
import numpy as np
from mrsimulator import Method
from mrsimulator import Simulator
from mrsimulator import Site
from mrsimulator import SpinSystem
from mrsimulator.method.lib import BlochDecaySpectrum


def pre_setup(isotope, shift, reference_offset):
    spin_system = SpinSystem(
        sites=[Site(isotope=isotope, isotropic_chemical_shift=shift)]
    )
    method = Method.parse_dict_with_units(
        dict(
            channels=[isotope],
            spectral_dimensions=[
                {
                    "count": 2048,
                    "spectral_width": "25 kHz",
                    "reference_offset": f"{reference_offset} Hz",
                    "events": [
                        {
                            "magnetic_flux_density": "9.4 T",
                            "transition_queries": [{"ch1": {"P": [-1]}}],
                        }
                    ],
                }
            ],
        )
    )
    sim = Simulator()
    sim.spin_systems.append(spin_system)
    sim.methods = [method]
    sim.run()
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


def test_large_shifts():
    shifts = [-3480.5, -3461.1, -3440.5]  # in ppm

    mth = BlochDecaySpectrum(
        channels=["207Pb"],
        magnetic_flux_density=11.744,  # in T
        rotor_frequency=20000,  # in Hz
        spectral_dimensions=[
            {
                "count": 2048,
                "spectral_width": 25000.0,
                "reference_offset": -366357.20126875,
                "label": "207Pb - Frequency",
                "origin_offset": 104692954,
            }
        ],
    )

    for shift in shifts:
        sys = SpinSystem(sites=[Site(isotope="207Pb", isotropic_chemical_shift=shift)])
        sim = Simulator(spin_systems=[sys], methods=[mth])
        sim.run()

        simulation = sim.methods[0].simulation
        simulation.x[0].to("ppm", "nmr_frequency_ratio")
        max_amp = simulation.x[0].coordinates[int(np.argmax(simulation))]

        np.testing.assert_almost_equal(max_amp.value, shift, decimal=1)
