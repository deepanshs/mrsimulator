cimport sandbox as clib

cimport numpy as np
import numpy as np
import cython

__author__ = "Deepansh J. Srivastava"
__email__ = ["srivastava.89@osu.edu", "deepansh2012@gmail.com"]

# @cython.boundscheck(False)
# @cython.wraparound(False)
# def MRS_plan(int geodesic_polyhedron_frequency,
#         int number_of_sidebands)

@cython.boundscheck(False)
@cython.wraparound(False)
def pre_phase_components(int number_of_sidebands, double sample_rotation_frequency_in_Hz):
    r"""

    """
    cdef int n1 = 9 * number_of_sidebands
    cdef np.ndarray[complex] pre_phase = np.zeros(n1, dtype=np.complex128)
    clib.__get_components(number_of_sidebands, sample_rotation_frequency_in_Hz, &pre_phase[0])
    return pre_phase.reshape(9, number_of_sidebands)


@cython.boundscheck(False)
@cython.wraparound(False)
def wigner_dm0_vector(int l, double beta):
    r"""

    """
    cdef int n1 = (2 * l + 1)
    cdef np.ndarray[double] R_out = np.zeros(n1, dtype=np.float64)
    clib.wigner_dm0_vector(l, beta, &R_out[0])
    return R_out


@cython.boundscheck(False)
@cython.wraparound(False)
def wigner_rotation(int l, np.ndarray[complex] R_in,
                    cos_alpha = None, cos_beta = None,
                    wigner_matrix=None, phase_alpha=None):
    r"""

    """
    cdef int n1 = 2 * l + 1
    cdef np.ndarray[double, ndim=1] wigner, cos_alpha_c
    cdef np.ndarray[double complex] exp_I_beta_c
    cos_alpha_c = np.asarray(cos_alpha, dtype=np.float64)

    if wigner_matrix is None:
        n = cos_beta.size
        sin_beta = np.sqrt(1 - cos_beta**2)
        wigner = np.empty(n1**2 * n)
        exp_I_beta_c = np.asarray(cos_beta+1j*sin_beta, dtype=np.complex128)
        clib.wigner_d_matrices_from_exp_I_beta(l, n, &exp_I_beta_c[0], &wigner[0])
    else:
        n = wigner_matrix.shape[0]
        wigner = np.asarray(wigner_matrix.ravel(), dtype=np.float64)

    # if phase_alpha is None:
    #     alpha = np.arccos(cos_alpha)
    #     temp = np.tile(np.arange(9)-4.0, n).reshape(n,9)
    #     phase_alpha = np.exp(-1j*temp*alpha[:, np.newaxis])

    # cdef np.ndarray[complex] phase_alpha_c = np.asarray(phase_alpha.ravel(), dtype=np.complex128)

    cdef np.ndarray[complex] R_out = np.zeros(n1*n, dtype=np.complex128)

    clib.__wigner_rotation(l, n, &wigner[0],
                           &cos_alpha_c[0], &R_in[0], &R_out[0])
    return R_out.reshape(n, n1)


@cython.boundscheck(False)
@cython.wraparound(False)
def wigner_d_matrix_cosines(int l, np.ndarray[double] cos_beta):
    r"""
    Returns a :math:`(2l+1) \times (2l+1)` wigner-d(cos_beta) matrix for rank $l$ at
    a given `cos_beta`. Currently only rank l=2 and l=4 is supported.

    If `cos_beta` is a 1D-numpy array of size n, a
    `n x (2l+1) x (2l+1)` matrix is returned instead.

    :ivar l: The angular momentum quantum number.
    :ivar cos_beta: An 1D numpy array or a scalar representing the cosine of $\beta$ angles.
    """
    n1 = (2 * l + 1)
    cdef int n = cos_beta.size
    cdef np.ndarray[double complex] exp_I_beta_c
    sin_beta = np.sqrt(1 - cos_beta**2)
    exp_I_beta_c = np.asarray(cos_beta+1j*sin_beta, dtype=np.complex128)
    cdef np.ndarray[double, ndim=1] wigner = np.empty(n * n1**2)
    clib.wigner_d_matrices_from_exp_I_beta(l, n, &exp_I_beta_c[0], &wigner[0])
    return wigner.reshape(n, n1, n1)


@cython.boundscheck(False)
@cython.wraparound(False)
def cosine_of_polar_angles_and_amplitudes(int geodesic_polyhedron_frequency=72):
    r"""
    Calculate the direction cosines and the related amplitudes for
    the positive quadrant of the sphere. The direction cosines corresponds to
    angle $\alpha$ and $\beta$, where $\alpha$ is the azimuthal angle and
    $\beta$ is the polar angle. The amplitudes are evaluated as

        `amp = 1/r**3`

    where `r` is the distance from the origin to the face of the unit
    octahedron in the positive quadrant along the line given by the values of
    $\alpha$ and $\beta$.

    :ivar geodesic_polyhedron_frequency:
        The value is an integer which represents the frequency of class I
        geodesic polyhedra. These polyhedra are used in calculating the
        spherical average. Presently we only use octahedral as the frequency1
        polyhedra. As the frequency of the geodesic polyhedron increases, the
        polyhedra approach a sphere geometry. For line-shape simulation, a higher
        frequency will result in a better powder averaging.
        The default value is 72.
        Read more on the `Geodesic polyhedron
        <https://en.wikipedia.org/wiki/Geodesic_polyhedron>`_.

    :return cos_alpha: The cosine of the azimuthal angle.
    :return cos_beta: The cosine of the polar angle.
    :return amp: The amplitude at the given $\alpha$ and $\beta$.
    """
    nt = geodesic_polyhedron_frequency
    cdef unsigned int n_orientations = int((nt+1) * (nt+2)/2)

    cdef np.ndarray[double complex] exp_I_alpha = np.empty(n_orientations, dtype=np.complex128)
    cdef np.ndarray[double complex] exp_I_beta = np.empty(n_orientations, dtype=np.complex128)
    cdef np.ndarray[double] amp = np.empty(n_orientations, dtype=np.float64)


    clib.__powder_averaging_setup(nt, &exp_I_alpha[0], &exp_I_beta[0], &amp[0], 1)

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
def triangle_interpolation(vector, np.ndarray[double, ndim=1] spectrum_amp,
                           double amp=1):
    r"""
    Given a vector of three points, this method interpolates the
    between the points to form a triangle. The height of the triangle is given
    as `2.0/(f[2]-f[1])` where `f` is the array `vector` sorted in an ascending
    order.

    :ivar vector: 1-D array of three points.
    :ivar spectrum_amp: A numpy array of amplitudes. This array is updated.
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

    clib.triangle_interpolation(f1, f2, f3, &amp_[0], &spectrum_amp[0], &points[0])



@cython.boundscheck(False)
@cython.wraparound(False)
def _one_d_simulator(
        # spectrum information
        double reference_offset,
        double increment,
        int number_of_points,

        float spin_quantum_number = 0.5,
        float larmor_frequency = 0.0,

        # CSA tensor information
        isotropic_chemical_shift = None,
        chemical_shift_anisotropy = None,
        chemical_shift_asymmetry = None,

        # quad tensor information
        quadrupolar_coupling_constant = None,
        quadrupolar_asymmetry = None,
        second_order_quad = 1,
        remove_second_order_quad_isotropic = 0,

        # dipolar coupling
        D = None,

        # spin rate, spin angle and number spinning sidebands
        int number_of_sidebands = 128,
        double sample_rotation_frequency_in_Hz = 0.0,
        rotor_angle = None,

        m_final = 0.5,
        m_initial = -0.5,

        # Euler angle -> principal to molecular frame
        # omega_PM=None,

        # Euler angles for powder averaging scheme
        int geodesic_polyhedron_frequency=90):

    nt = geodesic_polyhedron_frequency
    if isotropic_chemical_shift is None:
        isotropic_chemical_shift = 0
    isotropic_chemical_shift = np.asarray([isotropic_chemical_shift], dtype=np.float64).ravel()
    cdef number_of_sites = isotropic_chemical_shift.size
    cdef np.ndarray[double, ndim=1] isotropic_chemical_shift_c = isotropic_chemical_shift

    if spin_quantum_number > 0.5 and larmor_frequency == 0.0:
        raise Exception("'larmor_frequency' is required for quadrupole spins.")

    # Shielding anisotropic values
    if chemical_shift_anisotropy is None:
        chemical_shift_anisotropy = np.ones(number_of_sites, dtype=np.float64).ravel() #*1e-4*increment
    else:
        chemical_shift_anisotropy = np.asarray([chemical_shift_anisotropy], dtype=np.float64).ravel()
        # chemical_shift_anisotropy[np.where(chemical_shift_anisotropy==0.)] = 1e-4*increment
    if chemical_shift_anisotropy.size != number_of_sites:
        raise Exception("Number of shielding anisotropies are not consistent with the number of spins.")
    cdef np.ndarray[double, ndim=1] chemical_shift_anisotropy_c = chemical_shift_anisotropy

    # Shielding asymmetry values
    if chemical_shift_asymmetry is None:
        chemical_shift_asymmetry = np.zeros(number_of_sites, dtype=np.float64).ravel()
    else:
        chemical_shift_asymmetry = np.asarray([chemical_shift_asymmetry], dtype=np.float64).ravel()
    if chemical_shift_asymmetry.size != number_of_sites:
        raise Exception("Number of shielding asymmetry are not consistent with the number of spins.")
    cdef np.ndarray[double, ndim=1] chemical_shift_asymmetry_c = chemical_shift_asymmetry

    # Quad coupling constant
    if quadrupolar_coupling_constant is None:
        quadrupolar_coupling_constant = np.zeros(number_of_sites, dtype=np.float64).ravel()
    else:
        quadrupolar_coupling_constant = np.asarray([quadrupolar_coupling_constant], dtype=np.float64).ravel()
    if quadrupolar_coupling_constant.size != number_of_sites:
        raise Exception("Number of quad coupling constants are not consistent with the number of spins.")
    cdef np.ndarray[double, ndim=1] quadrupolar_coupling_constant_c = quadrupolar_coupling_constant

    # Quad asymmetry value
    if quadrupolar_asymmetry is None:
        quadrupolar_asymmetry = np.zeros(number_of_sites, dtype=np.float64).ravel()
    else:
        quadrupolar_asymmetry = np.asarray([quadrupolar_asymmetry], dtype=np.float64).ravel()
    if quadrupolar_asymmetry.size != number_of_sites:
        raise Exception("Number of quad asymmetry are not consistent with the number of spins.")
    cdef np.ndarray[double, ndim=1] quadrupolar_asymmetry_c = quadrupolar_asymmetry


    # Dipolar coupling constant
    if D is None:
        D = np.zeros(number_of_sites, dtype=np.float64).ravel()
    else:
        D = np.asarray([D], dtype=np.float64).ravel()
    if D.size != number_of_sites:
        raise Exception("Number of dipolar coupling are not consistent with the number of spins.")
    cdef np.ndarray[double, ndim=1] D_c = D

    if rotor_angle is None:
        rotor_angle = 54.735
    cdef double rotor_angle_in_rad_c = np.pi*rotor_angle/180.

    cdef second_order_quad_c = second_order_quad

    cdef np.ndarray[double, ndim=1] transition_c = np.asarray([m_initial, m_final], dtype=np.float64)

    cdef np.ndarray[double, ndim=1] amp = np.zeros(number_of_points * number_of_sites)

    cdef np.ndarray[double] ori_n = np.zeros(3*number_of_sites, dtype=np.float64)
    cdef np.ndarray[double] ori_e = np.zeros(3*number_of_sites, dtype=np.float64)

    cdef clib.isotopomer_ravel isotopomer_struct

    isotopomer_struct.number_of_sites = number_of_sites
    isotopomer_struct.spin = spin_quantum_number
    isotopomer_struct.larmor_frequency = larmor_frequency

    isotopomer_struct.isotropic_chemical_shift_in_Hz = &isotropic_chemical_shift_c[0]
    isotopomer_struct.shielding_anisotropy_in_Hz = &chemical_shift_anisotropy_c[0]
    isotopomer_struct.shielding_asymmetry = &chemical_shift_asymmetry_c[0]
    isotopomer_struct.shielding_orientation = &ori_n[0]

    isotopomer_struct.quadrupolar_constant_in_Hz = &quadrupolar_coupling_constant_c[0]
    isotopomer_struct.quadrupolar_asymmetry = &quadrupolar_asymmetry_c[0]
    isotopomer_struct.quadrupolar_orientation = &ori_e[0]

    isotopomer_struct.dipolar_couplings = &D_c[0]

    cdef int remove_second_order_quad_isotropic_c = remove_second_order_quad_isotropic
    clib.spinning_sideband_core(
            # spectrum information and related amplitude
            &amp[0],
            reference_offset,
            increment,
            number_of_points,

            &isotopomer_struct,

            second_order_quad_c,
            remove_second_order_quad_isotropic_c,

            # spin rate, spin angle and number spinning sidebands
            number_of_sidebands,
            sample_rotation_frequency_in_Hz,
            rotor_angle_in_rad_c,

            &transition_c[0],
            geodesic_polyhedron_frequency)


    freq = np.arange(number_of_points)*increment + reference_offset

    return freq, amp
