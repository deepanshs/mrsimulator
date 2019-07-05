from numpy.fft import fft, fftfreq
import numpy as np

from mrsimulator.python.angular_momentum import (
    wigner_d_matrix_cosines as wigner_matrices,
)
from mrsimulator.python.angular_momentum import wigner_dm0_vector as rotation_lab
from mrsimulator.python.angular_momentum import wigner_rotation as rotation


from mrsimulator.python.Hamiltonian import nuclear_shielding as NS
from mrsimulator.python.orientation import (
    trig_of_polar_angles_and_amplitudes as polar_coordinates,
)

# from mrsimulator.python.orientation import average_over_octant as averager
from mrsimulator.python.utils import pre_phase_components
from timeit import default_timer


def simulator(
    spectrum, isotopomers, transitions=[[-0.5, 0.5]], nt=90, number_of_sidebands=128
):

    start = default_timer()

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

    shift = 0.0
    if number_of_points % 2 == 0:
        shift = 0.5

    # powder setup
    # orientations
    cos_alpha, cos_beta, orientation_amp = polar_coordinates(nt)
    orientation_amp /= np.float64(number_of_sidebands)
    n_orientation = cos_beta.size

    # wigner matrix
    wigner_2 = wigner_matrices(2, cos_beta)

    # pre phase
    pre_phase = pre_phase_components(number_of_sidebands, spin_frequency)

    # sideband freq
    fft_freq = fftfreq(number_of_sidebands, d=1.0 / number_of_sidebands)
    print(fft_freq)
    vr_freq = np.tile(fft_freq * spin_frequency / increment, n_orientation)
    vr_freq.shape = (n_orientation, number_of_sidebands)

    # rotor to lab frame transformation
    lab_vector_2 = rotation_lab(2, rotor_angle)

    # allocate memory for calculations
    R2_out = np.empty((n_orientation, 5), dtype=np.complex128)
    spectrum = np.zeros(number_of_points)

    shape = (n_orientation, number_of_sidebands)
    sideband_amplitude = np.empty(shape, dtype=np.float64)
    local_frequency = np.empty(shape, dtype=np.float64)

    # start calculating the spectrum for every site in every isotopomer.
    for isotopomer in isotopomers:

        sites = isotopomer["sites"]
        # n_sites = len(sites)
        amp = np.zeros(number_of_points)

        for transition in transitions:
            for i, site in enumerate(sites):
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
                R0, R2 = NS(iso, zeta, eta, transition)

                local_frequency_offset = shift + (R0 - frequency[0]) / increment

                start0 = default_timer()
                # rotation from PAS to Rotor frame over all orientations
                R2_out = rotation(2, R2, wigner_matrix=wigner_2, cos_alpha=cos_alpha)
                # rotation from rotor to lab frame over all orientations
                R2_out *= lab_vector_2

                # calculating side-band amplitudes
                sideband_amplitude = np.exp(np.dot(R2_out, pre_phase[2:7]))
                sideband_amplitude = np.abs(fft(sideband_amplitude, axis=-1)) ** 2
                sideband_amplitude *= orientation_amp[:, np.newaxis]

                # calculating local frequencies
                local_frequency[:] = (
                    local_frequency_offset
                    + vr_freq
                    + (R2_out[:, 2].real[:, np.newaxis]) / increment
                )
                print("time for calculating the frequencies", default_timer() - start0)

                # interpolate in-between the frequencies to generate a smooth spectrum.
                # This is the slowest part of the code.
                # for j in range(number_of_sidebands):
                #     averager(amp, local_frequency[:, j],
                #              nt, sideband_amplitude[:, j], 0)

                amp, freq = np.histogram(
                    local_frequency,
                    weights=sideband_amplitude,
                    bins=number_of_points,
                    range=[0, frequency.size],
                )
                # average over all spins
                spectrum += amp * isotopomer["abundance"]

    print("time", default_timer() - start)
    return frequency, spectrum
