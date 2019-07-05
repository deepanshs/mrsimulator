cimport sandbox as clib

cimport numpy as np
import numpy as np
import cython
from .utils import __get_spin_attribute__

__author__ = "Deepansh J. Srivastava"
__email__ = ["srivastava.89@osu.edu", "deepansh2012@gmail.com"]


@cython.boundscheck(False)
@cython.wraparound(False)
def pre_phase_components(int number_of_sidebands, double sample_rotation_frequency):
    r"""

    """
    cdef int n1 = 9 * number_of_sidebands
    cdef np.ndarray[complex] pre_phase = np.zeros(n1, dtype=np.complex128)
    clib.__get_pre_phase_components(number_of_sidebands, sample_rotation_frequency, &pre_phase[0])
    return pre_phase.reshape(9, number_of_sidebands)


@cython.boundscheck(False)
@cython.wraparound(False)
def wigner_dm0_vector(int l, double beta):
    r"""

    """
    cdef int n1 = (2 * l + 1)
    cdef np.ndarray[double] R_out = np.zeros(n1, dtype=np.float64)
    clib.__wigner_dm0_vector(l, beta, &R_out[0])
    return R_out


@cython.boundscheck(False)
@cython.wraparound(False)
def wigner_rotation(int l, np.ndarray[complex] R_in,
                    cos_alpha = None, cos_beta = None,
                    wigner_matrix=None, phase_alpha=None):
    r"""

    """
    cdef int n1 = 2 * l + 1
    cdef np.ndarray[double, ndim=1] wigner, cos_alpha_c, cos_beta_c
    cos_alpha_c = np.asarray(cos_alpha, dtype=np.float64)

    if wigner_matrix is None:
        n = cos_beta.size
        wigner = np.empty(n1**2 * n)
        cos_beta_c = np.asarray(cos_beta, dtype=np.float64)
        clib.__wigner_d_matrix_cosine(l, n, &cos_beta_c[0], &wigner[0])
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
    Returns a $(2l+1) \times (2l+1)$ wigner-d(cos_beta) matrix for rank $l$ at
    a given `cos_beta`. Currently only rank l=2 and l=4 is supported.

    If `cos_beta` is a 1D-numpy array of size n, a
    `n x (2l+1) x (2l+1)` matrix is returned instead.

    :ivar l: The angular momentum quantum number.
    :ivar cos_beta: An 1D numpy array or a scalar representing the cosine of $\beta$ angles.
    """
    n1 = (2 * l + 1)
    cdef int n = cos_beta.size
    cdef np.ndarray[double, ndim=1] wigner = np.empty(n * n1**2)
    clib.__wigner_d_matrix_cosine(l, n, &cos_beta[0], &wigner[0])
    return wigner.reshape(n, n1, n1)


@cython.boundscheck(False)
@cython.wraparound(False)
def trig_of_polar_angles_and_amplitudes(int geodesic_polyhedron_frequency=72):
    r"""
    Calculates and return the direction cosines and the related amplitudes for
    the positive quadrant of the sphere. The direction cosines corresponds to
    angle $\alpha$ and $\beta$, where $\alpha$ is the azimuthal angle and
    $\beta$ is the polar angle. The amplitudes are evaluated as

        >>> amp = 1/r**3

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

    cdef np.ndarray[double] cos_alpha = np.empty(n_orientations, dtype=np.float64)
    cdef np.ndarray[double] cos_beta = np.empty(n_orientations, dtype=np.float64)
    cdef np.ndarray[double] amp = np.empty(n_orientations, dtype=np.float64)

    clib.__powder_averaging_setup(nt, &cos_alpha[0], &cos_beta[0], &amp[0], 1)

    return cos_alpha, cos_beta, amp


@cython.boundscheck(False)
@cython.wraparound(False)
def triangle_interpolation(vector, np.ndarray[double, ndim=1] spectrum_amp,
                           double offset=0, double amp=1):
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
    cdef np.ndarray[double, ndim=1] f_vector = np.asarray(vector + offset, dtype=np.float64)
    f_vector+=0.5
    cdef double *f1 = &f_vector[0]
    cdef double *f2 = &f_vector[1]
    cdef double *f3 = &f_vector[2]

    cdef np.ndarray[double, ndim=1] amp_ = np.asarray([amp])

    clib.triangle_interpolation(f1, f2, f3, &amp_[0], &spectrum_amp[0], &points[0])
