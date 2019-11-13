from mrsimulator import Dimension
from mrsimulator import Isotopomer, Site
from mrsimulator import Simulator
from mrsimulator.methods import one_d_spectrum
import numpy as np


def pre_setup(isotope, shift, reference_offset):
    isotopomer = Isotopomer(
        sites=[Site(isotope=isotope, isotropic_chemical_shift=shift)]
    )
    dimension = Dimension.parse_dict_with_units(
        dict(
            number_of_points=2046,
            isotope=isotope,
            magnetic_flux_density="9.4 T",
            spectral_width="25 kHz",
            reference_offset=f"{reference_offset} Hz",
        )
    )
    sim = Simulator()
    sim.isotopomers.append(isotopomer)
    sim.dimensions.append(dimension)
    x, y = sim.run(method=one_d_spectrum)
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
