cimport test as clib

cimport numpy as np
import numpy as np
import cython

from libcpp cimport bool as bool_t


__author__ = "Deepansh J. Srivastava"
__email__ = "srivastava.89@osu.edu"


clib.generate_table()

## wigner matrices

@cython.boundscheck(False)
@cython.wraparound(False)
def wigner_d_matrices(int l, np.ndarray[double] angle):
    cdef int n = angle.size
    cdef int n1 = (2*l+1)**2
    cdef np.ndarray[double] wigner = np.empty(n*n1, dtype=np.float64)
    clib.wigner_d_matrices(l, n, &angle[0], &wigner[0])
    return wigner


@cython.boundscheck(False)
@cython.wraparound(False)
def wigner_d_matrices_from_exp_I_beta(int l, bool_t half, np.ndarray[double complex] exp_I_beta):
    r"""
    Returns a :math:`(2l+1) \times (2l+1)` wigner-d(beta) matrix of rank $l$ at
    a given angle `beta` in the form of `exp(i\beta)`. Currently only rank l=2 and
    l=4 is supported.

    If `exp_I_beta` is a 1D-numpy array of size n, a
    `n x (2l+1) x (2l+1)` matrix is returned instead.

    :ivar l: The angular momentum quantum number.
    :ivar half: Compute only half of wigner matrix
    :ivar exp_I_beta: An 1D numpy array or a scalar representing $\exp\beta$.
    """
    n1 = (2 * l + 1)
    cdef int n = exp_I_beta.size
    size = n * n1*(l+1) if half else n * n1**2

    cdef np.ndarray[double, ndim=1] wigner = np.empty(size)

    clib.wigner_d_matrices_from_exp_I_beta(l, n, half, &exp_I_beta[0], &wigner[0])

    if half:
        return wigner.reshape(n, (l+1), n1)

    return wigner.reshape(n, n1, n1)

@cython.boundscheck(False)
@cython.wraparound(False)
def wigner_dm0_vector(int l, double beta):
    r"""

    """
    cdef int n1 = (2 * l + 1)
    cdef np.ndarray[double] R_out = np.zeros(n1, dtype=np.float64)
    clib.wigner_dm0_vector(l, beta, &R_out[0])
    return R_out


## wigner rotations

@cython.boundscheck(False)
@cython.wraparound(False)
def single_wigner_rotation(int l, np.ndarray[double] euler_angles, np.ndarray[double complex] R_in):
    cdef int n1 = (2 * l + 1)
    cdef np.ndarray[double complex] R_out = np.zeros(n1, dtype=np.complex128)
    clib.single_wigner_rotation(l, &euler_angles[0], &R_in[0], &R_out[0])
    return R_out


@cython.boundscheck(False)
@cython.wraparound(False)
def __wigner_rotation_2(int l, np.ndarray[double] cos_alpha,
                        np.ndarray[double] cos_beta,
                        np.ndarray[double complex] R_in):

    cdef int n1 = 2 * l + 1
    cdef int n = cos_alpha.size
    cdef np.ndarray[double, ndim=1] wigner
    cdef np.ndarray[double complex, ndim=1] exp_I_beta
    wigner = np.empty(n1 * (l+1) * n, dtype=np.float64)
    sin_beta = np.sqrt(1 - cos_beta**2)
    exp_I_beta = np.asarray(cos_beta + 1j*sin_beta, dtype=np.complex128)
    clib.wigner_d_matrices_from_exp_I_beta(l, n, True, &exp_I_beta[0], &wigner[0])

    cdef np.ndarray[double complex] exp_im_alpha
    exp_im_alpha = np.empty(4 * n, dtype=np.complex128)
    exp_im_alpha[3*n:] = cos_alpha + 1j*np.sqrt(1.0 - cos_alpha**2)
    clib.get_exp_Im_alpha(n, 1, &exp_im_alpha[0])

    cdef np.ndarray[complex] R_out = np.zeros((l + 1)*n, dtype=np.complex128)


    clib.__wigner_rotation_2(l, n, &wigner[0], &exp_im_alpha[0], &R_in[0], &R_out[0])
    return R_out.reshape(n, (l + 1))


@cython.boundscheck(False)
@cython.wraparound(False)
def get_exp_Im_alpha(int n, np.ndarray[double] cos_alpha, bool_t allow_fourth_rank):
    cdef unsigned int n_ = n
    cdef np.ndarray[double complex] exp_Im_alpha = np.empty(4*n, dtype=np.complex128)
    exp_Im_alpha[3*n:] = cos_alpha + 1j*np.sqrt(1.0 - cos_alpha**2)
    clib.get_exp_Im_alpha(n_, allow_fourth_rank, &exp_Im_alpha[0])
    return exp_Im_alpha


@cython.boundscheck(False)
@cython.wraparound(False)
def pre_phase_components(unsigned int number_of_sidebands, double rotor_frequency_in_Hz):
    cdef int n1 = 4 * number_of_sidebands
    cdef np.ndarray[double] pre_phase = np.zeros(2*n1, dtype=np.float64)
    clib.get_sideband_phase_components(number_of_sidebands, rotor_frequency_in_Hz, &pre_phase[0])
    return pre_phase.view(dtype=np.complex128).reshape(4, number_of_sidebands)


@cython.boundscheck(False)
@cython.wraparound(False)
def cosine_of_polar_angles_and_amplitudes(int integration_density=72):
    r"""
    Calculate the direction cosines and the related amplitudes for
    the positive quadrant of the sphere. The direction cosines corresponds to
    angle $\alpha$ and $\beta$, where $\alpha$ is the azimuthal angle and
    $\beta$ is the polar angle. The amplitudes are evaluated as

        `amp = 1/r**3`

    where `r` is the distance from the origin to the face of the unit
    octahedron in the positive quadrant along the line given by the values of
    $\alpha$ and $\beta$.

    :ivar integration_density:
        The value is an integer which represents the frequency of class I
        geodesic polyhedra. These polyhedra are used in calculating the
        spherical average. Presently we only use octahedral as the frequency1
        polyhedra. As the frequency of the geodesic polyhedron increases, the
        polyhedra approach a sphere geometry. A higher frequency will result in a
        better powder averaging. The default value is 72.
        Read more on the `Geodesic polyhedron
        <https://en.wikipedia.org/wiki/Geodesic_polyhedron>`_.

    :return cos_alpha: The cosine of the azimuthal angle.
    :return cos_beta: The cosine of the polar angle.
    :return amp: The amplitude at the given $\alpha$ and $\beta$.
    """
    nt = integration_density
    cdef unsigned int octant_orientations = int((nt+1) * (nt+2)/2)

    cdef np.ndarray[double complex] exp_I_alpha = np.empty(octant_orientations, dtype=np.complex128)
    cdef np.ndarray[double complex] exp_I_beta = np.empty(octant_orientations, dtype=np.complex128)
    cdef np.ndarray[double] amp = np.empty(octant_orientations, dtype=np.float64)

    clib.averaging_setup(nt, &exp_I_alpha[0], &exp_I_beta[0], &amp[0])

    return exp_I_alpha, exp_I_beta, amp


@cython.boundscheck(False)
@cython.wraparound(False)
def octahedronInterpolation(np.ndarray[double] spec, np.ndarray[double, ndim=2] freq, int nt, np.ndarray[double, ndim=2] amp, int stride=1):
    cdef int i
    cdef int number_of_sidebands = amp.shape[0]
    for i in range(number_of_sidebands):
        clib.octahedronInterpolation(&spec[0], &freq[i,0], nt, &amp[i,0], stride, spec.size)


@cython.boundscheck(False)
@cython.wraparound(False)
def triangle_interpolation1D(vector, np.ndarray[double, ndim=1] spectrum_amp,
                           double amp=1):
    r"""
    Given a vector of three points, this method interpolates the
    between the points to form a triangle. The height of the triangle is given
    as `2.0/(f[2]-f[1])` where `f` is the array `vector` sorted in an ascending
    order.

    :ivar vector: 1-D array of three points.
    :ivar spectrum_amp: A numpy array of amplitudes. This array is output.
    :ivar offset: A float specifying the offset. The points from array `vector`
                  are incremented or decremented based in this values. The
                  default value is 0.
    :ivar amp: A float specifying the offset. The points from array `vector`
               are incremented or decremented based in this values. The
               default value is 0.
    """
    cdef np.ndarray[int, ndim=1] points = np.asarray([spectrum_amp.size], dtype=np.int32)
    cdef np.ndarray[double, ndim=1] f_vector = np.asarray(vector, dtype=np.float64)

    cdef double *f1 = &f_vector[0]
    cdef double *f2 = &f_vector[1]
    cdef double *f3 = &f_vector[2]

    cdef np.ndarray[double, ndim=1] amp_ = np.asarray([amp])

    clib.triangle_interpolation1D(f1, f2, f3, &amp_[0], &spectrum_amp[0], &points[0])


@cython.boundscheck(False)
@cython.wraparound(False)
def triangle_interpolation2D(vector1, vector2, np.ndarray[double, ndim=2] spectrum_amp,
                            double amp=1):
    r"""
    Given a vector of three points, this method interpolates the
    between the points to form a triangle. The height of the triangle is given
    as `2.0/(f[2]-f[1])` where `f` is the array `vector` sorted in an ascending
    order.

    :ivar vector1: 1-D array of three points.
    :ivar vector2: 1-D array of three points.
    :ivar spectrum_amp: A numpy array of amplitudes. This array is the output.
    """
    shape = np.asarray([spectrum_amp.shape[0], spectrum_amp.shape[1]], dtype=np.int32)
    # cdef np.ndarray[int, ndim=1] points = shape
    cdef np.ndarray[double, ndim=1] f1_vector = np.asarray(vector1, dtype=np.float64)
    cdef np.ndarray[double, ndim=1] f2_vector = np.asarray(vector2, dtype=np.float64)

    cdef double *f11 = &f1_vector[0]
    cdef double *f12 = &f1_vector[1]
    cdef double *f13 = &f1_vector[2]

    cdef double *f21 = &f2_vector[0]
    cdef double *f22 = &f2_vector[1]
    cdef double *f23 = &f2_vector[2]

    cdef np.ndarray[double, ndim=1] amp_ = np.asarray([amp])

    clib.triangle_interpolation2D(f11, f12, f13, f21, f22, f23, &amp_[0],
                &spectrum_amp[0, 0], shape[0], shape[1])

@cython.boundscheck(False)
@cython.wraparound(False)
def __batch_wigner_rotation(unsigned int octant_orientations,
                            unsigned int n_octants,
                            np.ndarray[double] wigner_2j_matrices,
                            np.ndarray[double complex] R2,
                            np.ndarray[double] wigner_4j_matrices,
                            np.ndarray[double complex] R4,
                            np.ndarray[double complex] exp_Im_alpha):

    cdef np.ndarray[double complex] w2 = np.empty(3*octant_orientations*n_octants, dtype=np.complex128)
    cdef np.ndarray[double complex] w4 = np.empty(5*octant_orientations*n_octants, dtype=np.complex128)
    clib.__batch_wigner_rotation(octant_orientations, n_octants,
                            &wigner_2j_matrices[0], &R2[0], &wigner_4j_matrices[0],
                            &R4[0], &exp_Im_alpha[0], &w2[0], &w4[0])
    return w2, w4


# @cython.boundscheck(False)
# @cython.wraparound(False)
# def _one_d_simulator(
#         # spectrum information
#         double reference_offset,
#         double increment,
#         int number_of_points,

#         float spin_quantum_number = 0.5,
#         float larmor_frequency = 0.0,

#         # CSA tensor information
#         isotropic_chemical_shift = None,
#         shielding_anisotropy = None,
#         shielding_asymmetry = None,
#         shielding_orientations = None,

#         # quad tensor information
#         quadrupolar_coupling_constant = None,
#         quadrupolar_eta = None,
#         quadrupole_orientations = None,

#         second_order_quad = 1,
#         remove_2nd_order_quad_isotropic = 0,

#         # dipolar coupling
#         D = None,

#         # spin rate, spin angle and number spinning sidebands
#         int number_of_sidebands = 128,
#         double rotor_frequency_in_Hz = 0.0,
#         rotor_angle_in_rad = None,

#         m_final = 0.5,
#         m_initial = -0.5,

#         # Euler angle -> principal to molecular frame
#         # omega_PM=None,

#         # Euler angles for powder averaging scheme
#         int integration_density=90,
#         int integration_volume=0):



#     nt = integration_density
#     if isotropic_chemical_shift is None:
#         isotropic_chemical_shift = 0
#     isotropic_chemical_shift = np.asarray([isotropic_chemical_shift], dtype=np.float64).ravel()
#     cdef number_of_sites = isotropic_chemical_shift.size
#     cdef np.ndarray[double, ndim=1] isotropic_chemical_shift_c = isotropic_chemical_shift

#     cdef np.ndarray[float] spin = np.ones(number_of_sites, dtype=np.float32)*spin_quantum_number
#     cdef np.ndarray[double] gyromagnetic_ratio = np.ones(number_of_sites, dtype=np.float64)*larmor_frequency

#     if spin_quantum_number > 0.5 and larmor_frequency == 0.0:
#         raise Exception("'larmor_frequency' is required for quadrupole spins.")

#     # Shielding anisotropic values
#     if shielding_anisotropy is None:
#         shielding_anisotropy = np.ones(number_of_sites, dtype=np.float64).ravel() #*1e-4*increment
#     else:
#         shielding_anisotropy = np.asarray([shielding_anisotropy], dtype=np.float64).ravel()
#     if shielding_anisotropy.size != number_of_sites:
#         raise Exception("Number of shielding anisotropies are not consistent with the number of spins.")
#     cdef np.ndarray[double, ndim=1] shielding_anisotropy_c = shielding_anisotropy

#     # Shielding asymmetry values
#     if shielding_asymmetry is None:
#         shielding_asymmetry = np.zeros(number_of_sites, dtype=np.float64).ravel()
#     else:
#         shielding_asymmetry = np.asarray([shielding_asymmetry], dtype=np.float64).ravel()
#     if shielding_asymmetry.size != number_of_sites:
#         raise Exception("Number of shielding asymmetry are not consistent with the number of spins.")
#     cdef np.ndarray[double, ndim=1] shielding_asymmetry_c = shielding_asymmetry

#     # Shielding orientations
#     if shielding_orientations is None:
#         shielding_orientations = np.zeros(3*number_of_sites, dtype=np.float64).ravel()
#     else:
#         shielding_orientations = np.asarray([shielding_orientations], dtype=np.float64).ravel()
#     if shielding_orientations.size != 3*number_of_sites:
#         raise Exception("Number of euler angles are not consistent with the number of shielding tensors.")
#     cdef np.ndarray[double, ndim=1] shielding_orientations_c = shielding_orientations*np.pi/180.0

#     # Quad coupling constant
#     if quadrupolar_coupling_constant is None:
#         quadrupolar_coupling_constant = np.zeros(number_of_sites, dtype=np.float64).ravel()
#     else:
#         quadrupolar_coupling_constant = np.asarray([quadrupolar_coupling_constant], dtype=np.float64).ravel()
#     if quadrupolar_coupling_constant.size != number_of_sites:
#         raise Exception("Number of quad coupling constants are not consistent with the number of spins.")
#     cdef np.ndarray[double, ndim=1] quadrupolar_coupling_constant_c = quadrupolar_coupling_constant

#     # Quad asymmetry value
#     if quadrupolar_eta is None:
#         quadrupolar_eta = np.zeros(number_of_sites, dtype=np.float64).ravel()
#     else:
#         quadrupolar_eta = np.asarray([quadrupolar_eta], dtype=np.float64).ravel()
#     if quadrupolar_eta.size != number_of_sites:
#         raise Exception("Number of quad asymmetry are not consistent with the number of spins.")
#     cdef np.ndarray[double, ndim=1] quadrupole_asymmetry_c = quadrupolar_eta

#     # Quadrupolar orientations
#     if quadrupole_orientations is None:
#         quadrupole_orientations = np.zeros(3*number_of_sites, dtype=np.float64).ravel()
#     else:
#         quadrupole_orientations = np.asarray([quadrupole_orientations], dtype=np.float64).ravel()
#     if quadrupole_orientations.size != 3*number_of_sites:
#         raise Exception("Number of euler angles are not consistent with the number of quad tensors.")
#     cdef np.ndarray[double, ndim=1] quadrupole_orientations_c = quadrupole_orientations*np.pi/180.0

#     # Dipolar coupling constant
#     if D is None:
#         D = np.zeros(number_of_sites, dtype=np.float64).ravel()
#     else:
#         D = np.asarray([D], dtype=np.float64).ravel()
#     if D.size != number_of_sites:
#         raise Exception("Number of dipolar coupling are not consistent with the number of spins.")
#     cdef np.ndarray[double, ndim=1] D_c = D

#     # if rotor_angle is None:
#     #     rotor_angle = 54.735
#     cdef double rotor_angle_in_rad_c = rotor_angle_in_rad
#     cdef second_order_quad_c = second_order_quad

#     cdef np.ndarray[double, ndim=1] transition_c = np.asarray([m_initial, m_final], dtype=np.float64)

#     cdef np.ndarray[double, ndim=1] amp = np.zeros(number_of_points * number_of_sites)

#     cdef clib.site_struct sites_c

#     sites_c.number_of_sites = number_of_sites
#     sites_c.spin = &spin[0]
#     sites_c.gyromagnetic_ratio = &gyromagnetic_ratio[0]

#     sites_c.isotropic_chemical_shift_in_ppm = &isotropic_chemical_shift_c[0]
#     sites_c.shielding_symmetric_zeta_in_ppm = &shielding_anisotropy_c[0]
#     sites_c.shielding_asymmetry = &shielding_asymmetry_c[0]
#     sites_c.shielding_orientation = &shielding_orientations_c[0]

#     sites_c.quadrupolar_Cq_in_Hz = &quadrupolar_coupling_constant_c[0]
#     sites_c.quadrupolar_eta = &quadrupole_asymmetry_c[0]
#     sites_c.quadrupolar_orientation = &quadrupole_orientations_c[0]


#     cdef bool_t remove_second_order_quad_isotropic_c = remove_2nd_order_quad_isotropic

#     cdef clib.MRS_dimension *dimension[1]
#     clib.mrsimulator_core(
#             # spectrum information and related amplitude
#             &amp[0],
#             reference_offset,
#             increment,
#             number_of_points,

#             &sites_c,
#             dimension,

#             second_order_quad_c,
#             remove_second_order_quad_isotropic_c,

#             # spin rate, spin angle and number spinning sidebands
#             number_of_sidebands,
#             rotor_frequency_in_Hz,
#             rotor_angle_in_rad_c,

#             &transition_c[0],
#             integration_density,
#             integration_volume,           # 0-octant, 1-hemisphere, 2-sphere.
#             1
#             )


#     freq = np.arange(number_of_points)*increment + reference_offset

#     return freq, amp
