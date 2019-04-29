cimport nmr_methods as clib
cimport numpy as np
import numpy as np
# from nmr.lib import EulerAnglesInRadians
import cython

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


@cython.boundscheck(False)
@cython.wraparound(False)
def one_d_simulator(
        # spectrum information
        double spectral_start,
        double spectral_increment,
        int number_of_points,

        qunatum_number = None,
        wo = None,

        # CSA tensor information
        iso = None,
        aniso = None,
        eta = None,

        # quad tensor information
        Cq = None,
        eta_e = None,
        quad_second_order = False,

        # dipolar coupling
        D = None,

        # spin rate, spin angle and number spinning sidebands
        int ph_step = 16,
        double spin_frequency = 0,
        rotor_angle = None,

        m_final = 0.5,
        m_initial = -0.5,

        # Euler angle -> principal to molecular frame
        omega_PM=None,

        # Euler angles for powder averaging scheme
        int averaging_scheme=0,
        int averaging_size=32):

    """
    The firs
    
    """

    if iso is None:
        iso = 0
    iso = np.asarray([iso], dtype=np.float64).ravel()
    cdef number_of_sites = iso.size
    cdef np.ndarray[double, ndim=1] iso_n_c = iso
    
    # quantum numbers
    if qunatum_number is None:
        qunatum_number = np.ones(number_of_sites, dtype=np.float64)*0.5
    else:
        qunatum_number = np.asarray([qunatum_number], dtype=np.float64).ravel()
    if qunatum_number.size != number_of_sites:
        raise Exception("Number of spin quantum number are not consistent with the number of spins.")
    cdef np.ndarray[double, ndim=1] qunatum_number_c = qunatum_number

    # Larmor frequency
    if wo is None:
        wo = np.ones(number_of_sites, dtype=np.float64)
    else:
        wo = np.asarray([wo], dtype=np.float64).ravel()
    if wo.size != number_of_sites:
        raise Exception("Number of Larmor frequencies are not consistent with the number of spins.")
    cdef np.ndarray[double, ndim=1] wo_c = wo

    # Shielding anisotropic values
    if aniso is None:
        aniso = np.ones(number_of_sites, dtype=np.float64).ravel()*1e-4*spectral_increment
    else:
        aniso = np.asarray([aniso], dtype=np.float64).ravel()
        aniso[np.where(aniso==0.)] = 1e-4*spectral_increment
    if aniso.size != number_of_sites:
        raise Exception("Number of shielding anisotropies are not consistent with the number of spins.")
    cdef np.ndarray[double, ndim=1] aniso_n_c = aniso

    # Shielding asymmetry values
    if eta is None:
        eta = np.zeros(number_of_sites, dtype=np.float64).ravel()
    else:
        eta = np.asarray([eta], dtype=np.float64).ravel()
    if eta.size != number_of_sites:
        raise Exception("Number of shielding asymmetry are not consistent with the number of spins.")
    cdef np.ndarray[double, ndim=1] eta_n_c = eta

    # Quad coupling constant
    if Cq is None:
        Cq = np.zeros(number_of_sites, dtype=np.float64).ravel()
    else:
        Cq = np.asarray([Cq], dtype=np.float64).ravel()
    if Cq.size != number_of_sites:
        raise Exception("Number of quad coupling constants are not consistent with the number of spins.")
    cdef np.ndarray[double, ndim=1] Cq_e_c = Cq

    # Quad asymmetry value
    if eta_e is None:
        eta_e = np.zeros(number_of_sites, dtype=np.float64).ravel()
    else:
        eta_e = np.asarray([eta_e], dtype=np.float64).ravel()
    if eta_e.size != number_of_sites:
        raise Exception("Number of quad asymmetry are not consistent with the number of spins.")
    cdef np.ndarray[double, ndim=1] eta_e_c = eta_e

    
    # Dipolar coupling constant
    if D is None:
        D = np.zeros(number_of_sites, dtype=np.float64).ravel()
    else:
        D = np.asarray([D], dtype=np.float64).ravel()
    if D.size != number_of_sites:
        raise Exception("Number of dipolar coupling are not consistent with the number of spins.")
    cdef np.ndarray[double, ndim=1] D_c = D


    if omega_PM is None:
        omega_PM = np.zeros(3)
    cdef np.ndarray[double, ndim=1] omega_PM_c = np.asarray(omega_PM, dtype=np.float64).ravel()
    
    if rotor_angle is None:
        if spin_frequency == 0:
            rotor_angle = 0
        else:
            rotor_angle = 54.735
    cdef double rotor_angle_c = np.pi*rotor_angle/180.

    quad_second_order_c = 0
    if quad_second_order:
        quad_second_order_c = 1

    cdef np.ndarray[double, ndim=1] transition_c = np.asarray([m_initial, m_final], dtype=np.float64)
    
    cdef np.ndarray[double, ndim=1] amp = np.zeros(number_of_points * number_of_sites) 
    cdef double cpu_time_
    # print('rotor_angle in rad', rotor_angle_in_rad)
    # print('transition', transition_c)

    if averaging_scheme == 0:
        clib.spinning_sideband_core(
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
                quad_second_order_c,

                # Dipolar couplings 
                &D_c[0],

                # spin rate, spin angle and number spinning sidebands
                ph_step,
                spin_frequency,
                rotor_angle_c,

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
                rotor_angle_c,

                # Euler angle -> principal to molecular frame
                &omega_PM_c[0],

                # Euler angles for powder averaging scheme
                averaging_scheme,
                averaging_size
                )

    freq = np.arange(number_of_points)*spectral_increment + spectral_start   

    return freq, amp, cpu_time_


    # # print('spectral_start', spectral_start)
    # # print('spectral_increment', spectral_increment)
    # # print('number_of_points', number_of_points)

    # # print('qunatum_number', qunatum_number)
    # # print('wo', wo)
    # # print('iso', iso)
    # # print('aniso', aniso)
    # # print('eta', eta)


    # # print('Cq', Cq)
    # # print('eta_e', eta_e)
    # # print('quad_second_order', quad_second_order)

    # # print('ph_step', ph_step)
    # # print('spin_frequency', spin_frequency)
    # # print('rotor_angle', rotor_angle)

    # # print('m_final', m_final)
    # # print('m_initial', m_initial)

    # return _one_d_simulator(
    #         # spectrum information
    #         spectral_start,
    #         spectral_increment,
    #         number_of_points,

    #         qunatum_number_c,
    #         wo_c,

    #         # CSA tensor information
    #         iso_n_c,
    #         aniso_n_c,
    #         eta_n_c,

    #         # quad tensor information
    #         Cq_e_c,
    #         eta_e_c,
    #         quad_second_order_c,

    #         # spin rate, spin angle and number spinning sidebands
    #         ph_step,
    #         spin_frequency,
    #         rotor_angle_c,  ## in radians

    #         transition_c,

    #         # Euler angle -> principal to molecular frame
    #         omega_PM_c,

    #         # Euler angles for powder averaging scheme
    #         averaging_scheme,
    #         averaging_size)







@cython.boundscheck(False)
@cython.wraparound(False)
def _one_d_simulator(
            # spectrum information
            double spectral_start,
            double spectral_increment,
            int number_of_points,

            np.ndarray[double, ndim=1] qunatum_number,
            np.ndarray[double, ndim=1] wo,

            # CSA tensor information
            np.ndarray[double, ndim=1] iso,
            np.ndarray[double, ndim=1] aniso,
            np.ndarray[double, ndim=1] eta,

            # quad tensor information
            np.ndarray[double, ndim=1] Cq,
            np.ndarray[double, ndim=1] eta_e,
            int quad_second_order,

            # Dipolar coupling
            np.ndarray[double, ndim=1] D,

            # spin rate, spin angle and number spinning sidebands
            int ph_step,
            double spin_frequency,
            double rotor_angle_in_rad,

            np.ndarray[double, ndim=1] transition_c,

            # Euler angle -> principal to molecular frame
            np.ndarray[double, ndim=1] omega_PM,

            # Euler angles for powder averaging scheme
            int averaging_scheme,
            int averaging_size):

    # print('spectral_start', spectral_start)
    # print('spectral_increment', spectral_increment)
    # print('number_of_points', number_of_points)

    # print('qunatum_number', qunatum_number)
    # print('wo', wo)
    # print('iso', iso)
    # print('aniso', aniso)
    # print('eta', eta)


    # print('Cq', Cq)
    # print('eta_e', eta_e)
    # print('quad_second_order', quad_second_order)

    # print('ph_step', ph_step)
    # print('spin_frequency', spin_frequency)
    # print('rotor_angle in rad', rotor_angle_in_rad)
    # print('transition', transition_c)
    
    
    cdef number_of_sites = iso.size
    
    cdef np.ndarray[double, ndim=1] amp = np.zeros(number_of_points * number_of_sites) 
    cdef double cpu_time_
    # print('rotor_angle in rad', rotor_angle_in_rad)
    # print('transition', transition_c)

    if averaging_scheme == 0:
        clib.spinning_sideband_core(
                # spectrum information and related amplitude
                &amp[0],
                &cpu_time_,
                spectral_start,
                spectral_increment,
                number_of_points,

                &qunatum_number[0],
                &wo[0],

                # CSA tensor information
                &iso[0],
                &aniso[0],
                &eta[0],

                # quad tensor information
                &Cq[0],
                &eta_e[0],
                quad_second_order,

                # dipolar coupling
                &D[0],

                # spin rate, spin angle and number spinning sidebands
                ph_step,
                spin_frequency,
                rotor_angle_in_rad,

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
                rotor_angle_in_rad,

                # Euler angle -> principal to molecular frame
                &omega_PM[0],

                # Euler angles for powder averaging scheme
                averaging_scheme,
                averaging_size
                )

    freq = np.arange(number_of_points)*spectral_increment + spectral_start   

    return freq, amp, cpu_time_
