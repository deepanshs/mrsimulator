
cimport nmr_methods as clib
cimport numpy as np
import numpy as np
# from nmr.lib import EulerAnglesInRadians
import cython

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


@cython.boundscheck(False)
@cython.wraparound(False)
def CSA_spinning_sideband(
        # spectrum information
        double spectral_start,
        double spectral_increment,
        int number_of_points,

        qunatum_number,
        wo,

        # CSA tensor information
        iso,
        aniso,
        eta,

        # quad tensor information
        Cq,
        eta_e,
        quad_second_order,

        # spin rate, spin angle and number spinning sidebands
        int ph_step,
        double spin_frequency,
        double rotor_angle,

        m_final,
        m_initial,

        # Euler angle -> principal to molecular frame
        omega_PM=None,

        # Euler angles for powder averaging scheme
        int averaging_scheme=2,
        int averaging_size=31):

    """
    The firs
    
    """
    iso = np.asarray(iso)
    cdef number_of_sites = iso.size

    cdef np.ndarray[double, ndim=1, mode="c"] amp = np.zeros(number_of_points * number_of_sites) 
    cdef double cpu_time_
    

    cdef np.ndarray[double, ndim=1, mode="c"] transition_c = np.asarray([m_initial, m_final], dtype=np.float64)


    cdef np.ndarray[double, ndim=1, mode="c"] qunatum_number_c = np.asarray([qunatum_number], dtype=np.float64).ravel()
    cdef np.ndarray[double, ndim=1, mode="c"] wo_c = np.asarray([wo], dtype=np.float64).ravel()

    cdef np.ndarray[double, ndim=1, mode="c"] iso_n_c = np.asarray([iso], dtype=np.float64).ravel()

    aniso = np.asarray([aniso], dtype=np.float64).ravel()
    aniso[np.where(aniso==0.)] = 1e-4*spectral_increment
    cdef np.ndarray[double, ndim=1, mode="c"] aniso_n_c = aniso
    cdef np.ndarray[double, ndim=1, mode="c"] eta_n_c = np.asarray([eta], dtype=np.float64).ravel()

    cdef np.ndarray[double, ndim=1, mode="c"] Cq_e_c = np.asarray([Cq], dtype=np.float64).ravel()
    cdef np.ndarray[double, ndim=1, mode="c"] eta_e_c = np.asarray([eta_e], dtype=np.float64).ravel()

    if omega_PM is None:
        omega_PM = np.zeros(3)

    cdef np.ndarray[double, ndim=1, mode="c"] omega_PM_c = np.asarray(omega_PM, dtype=np.float64).ravel()
    
    if averaging_scheme == 0:
        clib.lineshape_cas_spinning_sideband_core(
                # spectrum information and related amplitude
                &amp[0],
                &cpu_time_,
                spectral_start,
                spectral_increment,
                number_of_points,

                &qunatum_number_c[0],
                &wo_c[0],

                # CSA tensor information
                &iso_n_c[0],
                &aniso_n_c[0],
                &eta_n_c[0],

                # quad tensor information
                &Cq_e_c[0],
                &eta_e_c[0],
                quad_second_order,

                # spin rate, spin angle and number spinning sidebands
                ph_step,
                spin_frequency,
                rotor_angle,

                &transition_c[0],
                # Euler angle -> principal to molecular frame
                # &omega_PM_c[0],

                # Euler angles for powder averaging scheme
                averaging_size,
                number_of_sites
                )
    
    else:
        clib.lineshape_cas_spinning_sideband_cython_angles(
                # spectrum information and related amplitude
                &amp[0],
                &cpu_time_,
                spectral_start,
                spectral_increment,
                number_of_points,

                # CSA tensor information
                iso,
                aniso,
                eta,

                # spin rate, spin angle and number spinning sidebands
                ph_step,
                spin_frequency,
                rotor_angle,

                # Euler angle -> principal to molecular frame
                &omega_PM_c[0],

                # Euler angles for powder averaging scheme
                averaging_scheme,
                averaging_size
                )

    freq = np.arange(number_of_points)*spectral_increment + spectral_start   

    return freq, amp, cpu_time_