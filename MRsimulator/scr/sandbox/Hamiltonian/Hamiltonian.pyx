

import numpy as np

class Hamiltonian:
    def __init__(self):
        self.__first_order_nuclear_shielding__ = 0
        self.__first_order_electric_quadrupole_coupling__ = 0
        self.__second_order_electric_quadrupole_coupling__ = 0

class FirstOrderNuclearShielding(Hamiltonian):
    def __init__(self):
        self.__first_order_nuclear_shielding__ = 1

class FirstOrderElectricQuadrupoleCoupling(Hamiltonian):
    def __init__(self):
        self.__first_order_electric_quadrupole_coupling__ = 1

class SecondOrderElectricQuadrupoleCoupling(Hamiltonian):
    def __init__(self):
        self.__second_order_electric_quadrupole_coupling__ = 1


@cython.boundscheck(False)
@cython.wraparound(False)
def rasterize(
            int number_of_points_x,
            double increment_x,
            double reference_offset_x,
            int number_of_points_y,
            double increment_y,
            double reference_offset_y,
            a,
            b,
            c
        ):

    a = np.asarray(a, dtype=np.float64).ravel()
    a[::2] -= reference_offset_x
    a[::2] /= increment_x

    a[1::2] -= reference_offset_y
    a[1::2] /= increment_y

    b = np.asarray(b, dtype=np.float64).ravel()
    b[::2] -= reference_offset_x
    b[::2] /= increment_x

    b[1::2] -= reference_offset_y
    b[1::2] /= increment_y

    c = np.asarray(c, dtype=np.float64).ravel()
    c[::2] -= reference_offset_x
    c[::2] /= increment_x

    c[1::2] -= reference_offset_y
    c[1::2] /= increment_y

    # print (a, b, c)

    cdef np.ndarray[double, ndim=1, mode='c'] amp = np.zeros((number_of_points_x*number_of_points_y), dtype=np.float64)
    cdef np.ndarray[double, ndim=1, mode='c'] v0, v1, v2

    cdef int i, j= int(a.size/2), i_

    for i in range(j):
        i_ = 2*i
        v0 = a[i_ : i_+1]
        v1 = b[i_ : i_+1]
        v2 = c[i_ : i_+1]
        # print (v0, v1, v2)

        clib.rasterization(&amp[0],
                        &v0[0],
                        &v1[0],
                        &v2[0],
                        number_of_points_x,
                        number_of_points_y)

    x = np.arange(number_of_points_x)*increment_x + reference_offset_x
    y = np.arange(number_of_points_y)*increment_y + reference_offset_y
    
    return x, y, amp.reshape(number_of_points_y, number_of_points_x)


@cython.boundscheck(False)
@cython.wraparound(False)
def CSA_static_lineshape(
        np.ndarray[double, ndim=1, mode='c'] haeberlen_values,
        int number_of_points,
        double start_frequency, 
        double frequency_bandwidth,
        int octa=1,
        int nt=64):
    """
    The method computes a static chemical shielding anisotropy (CSA) NMR
    lineshape spectrum by applying the powder averaging scheme to the CSA
    tensor in the principal axis system (PAS). Note, the CSA tensor is diagonal
    in the PAS with three principal components. The following code uses the
    Haeberlen convention for the principal components.

    The amplitude of the spectrum is evaluated at frequencies which are given by
    ``freq = np.arange(number_of_points)/number_of_points * frequency_bandwidth + start_frequency``

    The code implements the powder averaging scheme by
    Alderman, Solum and Grant, J. Chem. Phys, 84, 1985. DOI: 10.1063/1.450211
    
    :attr:haeberlen_values: namedtuple: A namedTuple HaeberlenNotation from
                                        pymatgen.analysis.nmr
    :attr:number_of_points: int: The number of points in the frequency dimension.
    :attr:start_frequency: float:The starting frequency.
    :attr:frequency_bandwidth float: The spectral width of the frequency spectrum.
    :attr:nt: int: The number of equilateral triangle segments along the edge of
                   an octahedron face. A higher number results in better averaging.
                   The default value is 48.

    :returns:freq: A ``Numpy array`` of frequencies.
    :returns:amp: A ``Numpy array`` of amplitudes corresponding to the frequencies.
    """

    cdef np.ndarray[double, ndim=1, mode='c'] freq = np.arange(number_of_points)/number_of_points * \
                frequency_bandwidth + start_frequency

    cdef np.ndarray[double, ndim=1, mode='c'] amp = np.zeros(number_of_points)

    cdef double cpu_time_
    cdef double iso = haeberlen_values[0]
    cdef double zeta = haeberlen_values[1]
    cdef double eta = haeberlen_values[2]

    clib.lineshape_csa_static(
                &amp[0],
                &cpu_time_,
                number_of_points,
                nt,
                start_frequency,
                frequency_bandwidth,
                iso,
                zeta,
                eta,
                octa,
                1)

    return freq, amp, cpu_time_

