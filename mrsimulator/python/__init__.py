import numpy as np
from numpy.fft import fft, fftfreq
from numpy import dot, exp, zeros
from mrsimulator.python import angular_momentum
from mrsimulator.python import orientation
from mrsimulator.python import utils


from mrsimulator.python.Hamiltonian import nuclear_shielding as NS

coordinate_distribution_ = orientation.cosine_of_polar_angles_and_amplitudes


class MRSplan:
    def __init__(
        self,
        geodesic_polyhedron_frequency,
        number_of_sidebands,
        sample_rotation_frequency_in_Hz,
        rotor_angle_in_rad,
        increment,
        allow_fourth_rank,
    ):

        nt = geodesic_polyhedron_frequency
        # orientations
        self.cos_alpha, cos_beta, self.orientation_amp = coordinate_distribution_(nt)
        self.orientation_amp /= np.float64(number_of_sidebands)
        n_orientation = cos_beta.size

        # sideband freq
        self.vr_freq = fftfreq(number_of_sidebands, d=1.0 / number_of_sidebands)
        self.vr_freq *= sample_rotation_frequency_in_Hz / increment

        # wigner matrix
        self.wigner_2 = angular_momentum.wigner_d_matrix_cosines(2, cos_beta)

        # rotor to lab frame transformation
        self.lab_vector_2 = angular_momentum.wigner_dm0_vector(2, rotor_angle_in_rad)

        # pre phase
        self.pre_phase = utils.pre_phase_components(
            number_of_sidebands, sample_rotation_frequency_in_Hz
        )

        # # allocate memory for calculations
        # R2_out = np.empty((n_orientation, 5), dtype=np.complex128)
        # self.spectrum = np.zeros(number_of_points)

        shape = (n_orientation, number_of_sidebands)
        self.vector = np.empty(shape, dtype=np.complex128)
        self.sideband_amplitude = np.empty(shape, dtype=np.float64)
        self.local_frequency = np.empty(n_orientation, dtype=np.float64)
        self.freq_offset = np.empty(n_orientation, dtype=np.float64)
        self.offset = np.empty(number_of_sidebands, dtype=np.float64)

    def process(self, R0, R2, R4):
        # rotation from PAS to Rotor frame over all orientations
        R2_out = angular_momentum.wigner_rotation(
            2, R2, wigner_matrix=self.wigner_2, cos_alpha=self.cos_alpha
        )
        # rotation from rotor to lab frame over all orientations
        R2_out *= self.lab_vector_2

        # calculating side-band amplitudes
        self.vector[:] = fft(exp(dot(R2_out, self.pre_phase[2:7])), axis=-1)
        self.sideband_amplitude[:] = self.vector.real ** 2 + self.vector.imag ** 2
        self.sideband_amplitude *= self.orientation_amp[:, np.newaxis]
