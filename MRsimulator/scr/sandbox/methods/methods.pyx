
@cython.boundscheck(False)
@cython.wraparound(False)
def one_d_simulator(
        # spectrum information
        double reference_offset,
        double increment,
        int number_of_points,

        # spin info
        double spin_quantum_number = 0.5,
        double larmor_frequency = 1.0,

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
        int number_of_sidebands = 64,
        double spin_frequency = 0.0,
        rotor_angle = None,

        double m_final = 0.5,
        double m_initial = -0.5,

        # Euler angle -> principal to molecular frame
        omega_PM=None,

        # Euler angles for powder averaging scheme
        int averaging_scheme=0,
        int averaging_size=128):

    """
    A one-dimensional NMR spectrum simulator. 
    
    The NMR spectrum is simulated over a range of frequencies evaluated as
        freq = np.arange(number_of_points)*increment + reference_offset

    The input can either be a number or a one-dimensional array of numbers.
    When the input is provided as a a one-dimensional array of size N,  

    :params: number_of_points: The number of frequency points over which the lineshape is evaluated. 
    :params: increment: The frequency increment.
    :params: reference_offset: The start frequency.
    

    :params: spin_quantum_number: The input can be a number or an array of `N` nummber
                                  representing the spin quantum numbers of N sites.
    :params: larmor_frequency: The larmor frequency or an array of `N` larmor frequencies
                                for N sites.
    :params: isotropic_chemical_shift: 
    :params: increment

    
    """

    if iso is None:
        iso = 0
    iso = np.asarray([iso], dtype=np.float64).ravel()
    cdef number_of_sites = iso.size
    cdef np.ndarray[double, ndim=1] iso_n_c = iso
    
    # quantum numbers
    # if spin_quantum_number is None:
    #     spin_quantum_number = np.ones(number_of_sites, dtype=np.float64)*0.5
    # else:
    #     spin_quantum_number = np.asarray([spin_quantum_number], dtype=np.float64).ravel()
    # if spin_quantum_number.size != number_of_sites:
    #     raise Exception("Number of spin quantum number are not consistent with the number of spins.")
    # cdef np.ndarray[double, ndim=1] qunatum_number_c = spin_quantum_number

    # Larmor frequency
    # if larmor_frequency is None:
    #     larmor_frequency = np.ones(number_of_sites, dtype=np.float64)
    # else:
    #     larmor_frequency = np.asarray([larmor_frequency], dtype=np.float64).ravel()
    # if larmor_frequency.size != number_of_sites:
    #     raise Exception("Number of Larmor frequencies are not consistent with the number of spins.")
    # cdef np.ndarray[double, ndim=1] wo_c = larmor_frequency

    # Shielding anisotropic values
    if aniso is None:
        aniso = np.ones(number_of_sites, dtype=np.float64).ravel()*1e-4*increment
    else:
        aniso = np.asarray([aniso], dtype=np.float64).ravel()
        aniso[np.where(aniso==0.)] = 1e-4*increment
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
            rotor_angle = 0.
            spin_frequency = 1e-5*increment
        else:
            rotor_angle = 54.7356
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
                reference_offset,
                increment,
                number_of_points,

                spin_quantum_number,
                larmor_frequency,

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
                number_of_sidebands,
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
                reference_offset,
                increment,
                number_of_points,

                # CSA tensor information
                iso,
                aniso,
                eta,

                # spin rate, spin angle and number spinning sidebands
                number_of_sidebands,
                spin_frequency,
                rotor_angle_c,

                # Euler angle -> principal to molecular frame
                &omega_PM_c[0],

                # Euler angles for powder averaging scheme
                averaging_scheme,
                averaging_size
                )

    freq = np.arange(number_of_points)*increment + reference_offset   

    return freq, amp, cpu_time_


    # # print('reference_offset', reference_offset)
    # # print('increment', increment)
    # # print('number_of_points', number_of_points)

    # # print('spin_quantum_number', spin_quantum_number)
    # # print('larmor_frequency', larmor_frequency)
    # # print('iso', iso)
    # # print('aniso', aniso)
    # # print('eta', eta)


    # # print('Cq', Cq)
    # # print('eta_e', eta_e)
    # # print('quad_second_order', quad_second_order)

    # # print('number_of_sidebands', number_of_sidebands)
    # # print('spin_frequency', spin_frequency)
    # # print('rotor_angle', rotor_angle)

    # # print('m_final', m_final)
    # # print('m_initial', m_initial)

    # return _one_d_simulator(
    #         # spectrum information
    #         reference_offset,
    #         increment,
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
    #         number_of_sidebands,
    #         spin_frequency,
    #         rotor_angle_c,  ## in radians

    #         transition_c,

    #         # Euler angle -> principal to molecular frame
    #         omega_PM_c,

    #         # Euler angles for powder averaging scheme
    #         averaging_scheme,
    #         averaging_size)






# # @cython.language_level(3)
# @cython.boundscheck(False)
# @cython.wraparound(False)
# def _one_d_simulator(
#             # spectrum information
#             double reference_offset,
#             double increment,
#             int number_of_points,

#             np.ndarray[double, ndim=1] spin_quantum_number,
#             np.ndarray[double, ndim=1] larmor_frequency,

#             # CSA tensor information
#             np.ndarray[double, ndim=1] iso,
#             np.ndarray[double, ndim=1] aniso,
#             np.ndarray[double, ndim=1] eta,

#             # quad tensor information
#             np.ndarray[double, ndim=1] Cq,
#             np.ndarray[double, ndim=1] eta_e,
#             int quad_second_order,

#             # Dipolar coupling
#             np.ndarray[double, ndim=1] D,

#             # spin rate, spin angle and number spinning sidebands
#             int number_of_sidebands,
#             double spin_frequency,
#             double rotor_angle_in_rad,

#             np.ndarray[double, ndim=1] transition_c,

#             # Euler angle -> principal to molecular frame
#             np.ndarray[double, ndim=1] omega_PM,

#             # Euler angles for powder averaging scheme
#             int averaging_scheme,
#             int averaging_size):

#     # print('reference_offset', reference_offset)
#     # print('increment', increment)
#     # print('number_of_points', number_of_points)

#     # print('spin_quantum_number', spin_quantum_number)
#     # print('larmor_frequency', larmor_frequency)
#     # print('iso', iso)
#     # print('aniso', aniso)
#     # print('eta', eta)


#     # print('Cq', Cq)
#     # print('eta_e', eta_e)
#     # print('quad_second_order', quad_second_order)

#     # print('number_of_sidebands', number_of_sidebands)
#     # print('spin_frequency', spin_frequency)
#     # print('rotor_angle in rad', rotor_angle_in_rad)
#     # print('transition', transition_c)
    
    
#     cdef number_of_sites = iso.size
    
#     cdef np.ndarray[double, ndim=1] amp = np.zeros(number_of_points * number_of_sites) 
#     cdef double cpu_time_
#     # print('rotor_angle in rad', rotor_angle_in_rad)
#     # print('transition', transition_c)

#     if averaging_scheme == 0:
#         clib.spinning_sideband_core(
#                 # spectrum information and related amplitude
#                 &amp[0],
#                 &cpu_time_,
#                 reference_offset,
#                 increment,
#                 number_of_points,

#                 &spin_quantum_number[0],
#                 &larmor_frequency[0],

#                 # CSA tensor information
#                 &iso[0],
#                 &aniso[0],
#                 &eta[0],

#                 # quad tensor information
#                 &Cq[0],
#                 &eta_e[0],
#                 quad_second_order,

#                 # dipolar coupling
#                 &D[0],

#                 # spin rate, spin angle and number spinning sidebands
#                 number_of_sidebands,
#                 spin_frequency,
#                 rotor_angle_in_rad,

#                 &transition_c[0],
#                 # Euler angle -> principal to molecular frame
#                 # &omega_PM_c[0],

#                 # Euler angles for powder averaging scheme
#                 averaging_size,
#                 number_of_sites
#                 )
    
#     else:
#         clib.lineshape_cas_spinning_sideband_cython_angles(
#                 # spectrum information and related amplitude
#                 &amp[0],
#                 &cpu_time_,
#                 reference_offset,
#                 increment,
#                 number_of_points,

#                 # CSA tensor information
#                 iso,
#                 aniso,
#                 eta,

#                 # spin rate, spin angle and number spinning sidebands
#                 number_of_sidebands,
#                 spin_frequency,
#                 rotor_angle_in_rad,

#                 # Euler angle -> principal to molecular frame
#                 &omega_PM[0],

#                 # Euler angles for powder averaging scheme
#                 averaging_scheme,
#                 averaging_size
#                 )

#     freq = np.arange(number_of_points)*increment + reference_offset   

#     return freq, amp, cpu_time_
