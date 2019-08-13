# -*- coding: utf-8 -*-
from numpy.fft import fft, fftfreq
import numpy as np
from numpy import dot, exp

from mrsimulator.python.angular_momentum import (
    wigner_d_matrix_cosines as wigner_matrices,
)
from mrsimulator.python.angular_momentum import wigner_dm0_vector as rotation_lab
from mrsimulator.python.angular_momentum import wigner_rotation as rotation


from mrsimulator.python.Hamiltonian import nuclear_shielding as NS
from mrsimulator.python.orientation import (
    cosine_of_polar_angles_and_amplitudes as polar_coordinates,
)

from mrsimulator.python.orientation import average_over_octant as averager
from mrsimulator.python.utils import pre_phase_components
import mrsimulator.python.transition_function as tf

# from timeit import default_timer
# import matplotlib.pyplot as plt

__author__ = "Deepansh J. Srivastava"
__email__ = ["srivastava.89@osu.edu", "deepansh2012@gmail.com"]


def simulator(
    spectrum, isotopomers, transitions=[[-0.5, 0.5]], nt=90, number_of_sidebands=128
):

    B0 = spectrum["magnetic_flux_density"]
    spin_frequency = spectrum["rotor_frequency"]
    rotor_angle = spectrum["rotor_angle"]

    frequency_scaling_factor = spectrum["gyromagnetic_ratio"] * B0

    number_of_points = spectrum["number_of_points"]
    spectral_width = spectrum["spectral_width"]
    reference_offset = spectrum["reference_offset"]
    frequency = (np.arange(number_of_points) / number_of_points) - 0.5
    frequency *= spectral_width
    frequency += reference_offset
    increment = frequency[1] - frequency[0]

    if spin_frequency < 1.0e-3:
        spin_frequency = 1.0e9
        rotor_angle = 0.0
        number_of_sidebands = 1

    shift_half_bin = 0.5

    # orientations
    cos_alpha, cos_beta, orientation_amp = polar_coordinates(nt)
    orientation_amp /= np.float64(number_of_sidebands)
    n_orientation = cos_beta.size

    # sideband freq
    vr_freq = fftfreq(number_of_sidebands, d=1.0 / number_of_sidebands)
    vr_freq *= spin_frequency / increment

    # wigner matrix
    wigner_2 = wigner_matrices(2, cos_beta)

    # rotor to lab frame transformation
    lab_vector_2 = rotation_lab(2, rotor_angle)

    # pre phase
    pre_phase = pre_phase_components(number_of_sidebands, spin_frequency)

    # allocate memory for calculations
    R2_out = np.empty((n_orientation, 5), dtype=np.complex128)
    spectrum = np.zeros(number_of_points)

    shape = (n_orientation, number_of_sidebands)
    temp = np.empty(shape, dtype=np.complex128)
    sideband_amplitude = np.empty(shape, dtype=np.float64)
    local_frequency = np.empty(n_orientation, dtype=np.float64)
    freq_offset = np.empty(n_orientation, dtype=np.float64)
    offset = np.empty(number_of_sidebands, dtype=np.float64)

    # start calculating the spectrum for every site in every isotopomer.
    for isotopomer in isotopomers:

        sites = isotopomer["sites"]
        spec = np.zeros(number_of_points)

        for transition in transitions:
            for site in sites:
                iso = site["isotropic_chemical_shift"]
                if iso.unit.physical_type == "dimensionless":
                    iso = iso.value * frequency_scaling_factor
                else:
                    iso = iso.value

                zeta = site["shielding_symmetric"]["anisotropy"]
                if zeta.unit.physical_type == "dimensionless":
                    zeta = zeta.value * frequency_scaling_factor
                else:
                    zeta = zeta.value

                eta = site["shielding_symmetric"]["asymmetry"]

                # Hailtonian
                # nuclear shielding
                R0, R2 = NS(iso, zeta, eta)
                scale = tf.p(transition[1], transition[0])
                R0 *= scale
                R2 *= scale

                local_frequency_offset = (
                    shift_half_bin + (R0 - frequency[0]) / increment
                )

                # rotation from PAS to Rotor frame over all orientations
                R2_out = rotation(2, R2, wigner_matrix=wigner_2, cos_alpha=cos_alpha)
                # rotation from rotor to lab frame over all orientations
                R2_out *= lab_vector_2

                # calculating side-band amplitudes
                temp[:] = fft(exp(dot(R2_out, pre_phase[2:7])), axis=-1)
                sideband_amplitude[:] = temp.real ** 2 + temp.imag ** 2
                sideband_amplitude *= orientation_amp[:, np.newaxis]

                # calculating local frequencies
                local_frequency[:] = R2_out[:, 2].real / increment

                # interpolate in-between the frequencies to generate a smooth spectrum.
                offset[:] = vr_freq + local_frequency_offset

                # print("before", default_timer() - start0)
                for j, shift in enumerate(offset):
                    if int(shift) >= 0 and int(shift) <= number_of_points:
                        freq_offset[:] = shift + local_frequency
                        # This is the slowest part of the code.
                        averager(spec, freq_offset, nt, sideband_amplitude[:, j])
                        # np.vectorize(averager(spec, freq_offset,
                        #                       nt, sideband_amplitude[:, j]))

                # print("time for computing site", default_timer() - start0)
        # average over all spins
        spectrum += spec * isotopomer["abundance"]

    return frequency, spectrum
