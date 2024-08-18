cimport test as clib

cimport numpy as np
import numpy as np
import cython

from libcpp cimport bool as bool_t


__author__ = "Deepansh J. Srivastava"
__email__ = "srivastava.89@osu.edu"


clib.generate_tables()

## wigner d elements
@cython.boundscheck(False)
@cython.wraparound(False)
def wigner_d_element(float l, float m1, float m2, double beta):
    return clib.wigner_d_element(l, m1, m2, beta)

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
    clib.get_exp_Im_angle(n, 1, &exp_im_alpha[0], 0.0)

    cdef np.ndarray[complex] R_out = np.zeros((l + 1)*n, dtype=np.complex128)


    clib.__wigner_rotation_2(l, n, &wigner[0], &exp_im_alpha[0], &R_in[0], &R_out[0])
    return R_out.reshape(n, (l + 1))


@cython.boundscheck(False)
@cython.wraparound(False)
def get_exp_Im_angle(int n, np.ndarray[double] cos_alpha, bool_t allow_4th_rank):
    cdef unsigned int n_ = n
    cdef np.ndarray[double complex] exp_Im_angle = np.empty(4*n, dtype=np.complex128)
    exp_Im_angle[3*n:] = cos_alpha + 1j*np.sqrt(1.0 - cos_alpha**2)
    clib.get_exp_Im_angle(n_, allow_4th_rank, &exp_Im_angle[0], 0.0)
    return exp_Im_angle


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
    cdef bool_t interploation = True

    cdef np.ndarray[double complex] exp_I_alpha = np.empty(octant_orientations, dtype=np.complex128)
    cdef np.ndarray[double complex] exp_I_beta = np.empty(octant_orientations, dtype=np.complex128)
    cdef np.ndarray[double] amp = np.empty(octant_orientations, dtype=np.float64)

    clib.averaging_setup(nt, &exp_I_alpha[0], &exp_I_beta[0], &amp[0], interploation)

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
                           double amp=1, str type="linear"):
    r"""Given a vector of three points, this method interpolates the
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
    :ivar type: Linear or Gaussian interpolation for delta functions.
    """
    cdef np.ndarray[int, ndim=1] points = np.asarray([spectrum_amp.size/2], dtype=np.int32)
    cdef np.ndarray[double, ndim=1] f_vector = np.asarray(vector, dtype=np.float64)

    cdef double *f1 = &f_vector[0]
    cdef double *f2 = &f_vector[1]
    cdef double *f3 = &f_vector[2]

    cdef np.ndarray[double, ndim=1] amp_ = np.asarray([amp])

    iso_intrp = 0 if type == "linear" else 1
    clib.triangle_interpolation1D(f1, f2, f3, &amp_[0], &spectrum_amp[0],
                &points[0], iso_intrp)


@cython.boundscheck(False)
@cython.wraparound(False)
def triangle_interpolation2D(vector1, vector2, np.ndarray[double, ndim=2] spectrum_amp,
                            double amp=1, str type="linear"):
    r"""Given a vector of three points, this method interpolates the
    between the points to form a triangle. The height of the triangle is given
    as `2.0/(f[2]-f[1])` where `f` is the array `vector` sorted in an ascending
    order.

    :ivar vector1: 1-D array of three points.
    :ivar vector2: 1-D array of three points.
    :ivar spectrum_amp: A numpy array of amplitudes. This array is the output.
    """
    shape = np.asarray([spectrum_amp.shape[0], spectrum_amp.shape[1]/2], dtype=np.int32)
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

    iso_intrp = 0 if type == "linear" else 1
    clib.triangle_interpolation2D(f11, f12, f13, f21, f22, f23, &amp_[0],
                &spectrum_amp[0, 0], shape[0], shape[1], iso_intrp)

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

@cython.boundscheck(False)
@cython.wraparound(False)
def rank_2_tensor_products(np.ndarray[double complex] tensor_a, np.ndarray[double complex] tensor_b):

    cdef np.ndarray[double] R_a = tensor_a.view(dtype=float)
    cdef np.ndarray[double] R_b = tensor_b.view(dtype=float)
    cdef np.ndarray[double] Delta_0 = np.zeros(2, dtype=float)
    cdef np.ndarray[double] Delta_2 = np.zeros(10, dtype=float)
    cdef np.ndarray[double] Delta_4 = np.zeros(18, dtype=float)
    clib.rank_2_tensor_products(&R_a[0], &R_b[0], &Delta_0[0], &Delta_2[0], &Delta_4[0])

    return Delta_0.view(dtype=complex), Delta_2.view(dtype=complex), Delta_4.view(dtype=complex)


@cython.boundscheck(False)
@cython.wraparound(False)
def vm_absd(double a):
    return clib.test_vm_absd(a)

@cython.boundscheck(False)
@cython.wraparound(False)
def vm_add(np.ndarray[double] A, np.ndarray[double] B):
    cdef np.ndarray[double] res = np.zeros(A.size, dtype=float)
    clib.test_vm_double_add(A.size, &A[0], &B[0], &res[0])
    return res

@cython.boundscheck(False)
@cython.wraparound(False)
def vm_add_inplace(np.ndarray[double] A, np.ndarray[double] B):
    clib.test_vm_double_add_inplace(A.size, &A[0], &B[0])

@cython.boundscheck(False)
@cython.wraparound(False)
def vm_sub(np.ndarray[double] A, np.ndarray[double] B):
    cdef np.ndarray[double] res = np.zeros(A.size, dtype=float)
    clib.test_vm_double_subtract(A.size, &A[0], &B[0], &res[0])
    return res

@cython.boundscheck(False)
@cython.wraparound(False)
def vm_sub_inplace(np.ndarray[double] A, np.ndarray[double] B):
    clib.test_vm_double_subtract_inplace(A.size, &A[0], &B[0])

@cython.boundscheck(False)
@cython.wraparound(False)
def vm_mult(np.ndarray[double] A, np.ndarray[double] B):
    cdef np.ndarray[double] res = np.zeros(A.size, dtype=float)
    clib.test_vm_double_multiply(A.size, &A[0], 1, &B[0], &res[0])
    return res

@cython.boundscheck(False)
@cython.wraparound(False)
def vm_mult_inplace(np.ndarray[double] A, np.ndarray[double] B):
    clib.test_vm_double_multiply_inplace(A.size, &A[0], 1, &B[0], 1)

@cython.boundscheck(False)
@cython.wraparound(False)
def vm_div(np.ndarray[double] A, np.ndarray[double] B):
    cdef np.ndarray[double] res = np.zeros(A.size, dtype=float)
    clib.test_vm_double_divide(A.size, &A[0], &B[0], &res[0])
    return res

@cython.boundscheck(False)
@cython.wraparound(False)
def vm_div_inplace(np.ndarray[double] A, np.ndarray[double] B):
    clib.test_vm_double_divide_inplace(A.size, &A[0], &B[0])

@cython.boundscheck(False)
@cython.wraparound(False)
def vm_cmult(np.ndarray[complex] A, np.ndarray[complex] B):
    cdef np.ndarray[complex] res = np.zeros(A.size, dtype=complex)
    clib.test_vm_double_complex_multiply(A.size, &A[0], &B[0], &res[0])
    return res

@cython.boundscheck(False)
@cython.wraparound(False)
def vm_cmult_conj(np.ndarray[complex] A, np.ndarray[complex] B):
    cdef np.ndarray[complex] res = np.zeros(A.size, dtype=complex)
    clib.test_vm_double_complex_conj_multiply(A.size, &A[0], &B[0], &res[0])
    return res

@cython.boundscheck(False)
@cython.wraparound(False)
def vm_sq(np.ndarray[double] A):
    cdef np.ndarray[double] res = np.zeros(A.size, dtype=float)
    clib.test_vm_double_square(A.size, &A[0], &res[0])
    return res

@cython.boundscheck(False)
@cython.wraparound(False)
def vm_sq_inplace(np.ndarray[double] A):
    cdef np.ndarray[double] res = np.zeros(A.size, dtype=float)
    clib.test_vm_double_square_inplace(A.size, &A[0])

@cython.boundscheck(False)
@cython.wraparound(False)
def vm_sqrt(np.ndarray[double] A):
    cdef np.ndarray[double] res = np.zeros(A.size, dtype=float)
    clib.test_vm_double_square_root(A.size, &A[0], &res[0])
    return res

@cython.boundscheck(False)
@cython.wraparound(False)
def vm_sqrt_inplace(np.ndarray[double] A):
    cdef np.ndarray[double] res = np.zeros(A.size, dtype=float)
    clib.test_vm_double_square_root_inplace(A.size, &A[0])


@cython.boundscheck(False)
@cython.wraparound(False)
def vm_cos(np.ndarray[double] A):
    cdef np.ndarray[double] res = np.zeros(A.size, dtype=float)
    clib.test_vm_double_cosine(A.size, &A[0], &res[0])
    return res

@cython.boundscheck(False)
@cython.wraparound(False)
def vm_sin(np.ndarray[double] A):
    cdef np.ndarray[double] res = np.zeros(A.size, dtype=float)
    clib.test_vm_double_sine(A.size, &A[0], &res[0])
    return res

@cython.boundscheck(False)
@cython.wraparound(False)
def vm_cos_I_sin(np.ndarray[double] A):
    cdef np.ndarray[complex] res = np.zeros(A.size, dtype=complex)
    clib.test_vm_cosine_I_sine(A.size, &A[0], &res[0])
    return res

@cython.boundscheck(False)
@cython.wraparound(False)
def vm_exp(np.ndarray[double] A):
    cdef np.ndarray[double] res = np.zeros(A.size, dtype=float)
    clib.test_vm_double_exp(A.size, &A[0], &res[0], 1)
    return res

@cython.boundscheck(False)
@cython.wraparound(False)
def vm_I_exp(np.ndarray[complex] A):
    cdef np.ndarray[complex] res = np.zeros(A.size, dtype=complex)
    clib.test_vm_double_complex_exp(A.size, &A[0], &res[0])
    return res
