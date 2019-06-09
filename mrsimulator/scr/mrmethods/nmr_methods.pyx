cimport nmr_methods as clib
from nmr_methods cimport (
    __powder_averaging_setup
)
cimport numpy as np
import numpy as np
import cython
from .utils import __get_spin_attribute__


__author__ = "Deepansh J. Srivastava"
__email__ = ["srivastava.89@osu.edu", "deepansh2012@gmail.com"]


@cython.boundscheck(False)
@cython.wraparound(False)
def one_d_spectrum(dict spectrum,
       list isotopomers,
       int verbose=0,
       int number_of_sidebands=90,
       int geodesic_polyhedron_frequency=72):
    """

    :ivar verbose:
        The allowed values are 0, 1, and 11. When the value is 1, the output is
        printed on the screen. When the value is 11, in addition to the output
        from 1, execution time is also printed on the screen.
        The default value is 0.
    :ivar number_of_sidebands:
        The value is an integer which corresponds to the number of sidebands
        simulated in the spectrum. The default value is 90. Note, when the
        sample spinning frequency is low, computation of more sidebands may be
        required for an acceptable result. The user is advised to ensure that
        enough sidebands are requested for computation.
    :ivar geodesic_polyhedron_frequency:
        The value is an integer which represents the frequency of class I
        geodesic polyhedra. These polyhedra are used in calculating the
        spherical average. Presently we only use octahedral as the frequency1
        polyhedra. As the frequency of the geodesic polyhedron increases, the
        polyhedra approach a sphere geometry. In lineshape simulation, a higher
        frequency will result in a better powder averaging.
        The default value is 72.
        Read more on the `Geodesic polyhedron <https://en.wikipedia.org/wiki/Geodesic_polyhedron>`_.
    """
# ---------------------------------------------------------------------
# spectrum ________________________________________________________
    cdef double sample_rotation_frequency = spectrum['rotor_frequency']
    cdef double rotor_angle = spectrum['rotor_angle']
    B0 = spectrum['magnetic_flux_density']

    if verbose in [1, 11]:
        print ('Setting up the virtual NMR spectrometer')
        print ('---------------------------------------')
        print (f'Adjusting the magnetic flux density to {B0} T.')
        _angle = spectrum['rotor_angle']
        print (f'Setting rotation angle to {_angle} rad.')
        print (f'Setting rotation frequency to {sample_rotation_frequency} Hz.')

# ---------------------------------------------------------------------
# spin observed _______________________________________________________
    # obs_dict = __get_spin_attribute__[detect]
    isotope = spectrum['isotope']

    # spin quantum number of the observed spin
    cdef double spin_quantum_number = spectrum['spin']

    # transitions of the observed spin
    number_of_energy_levels = int(2*spin_quantum_number+1)
    cdef int number_of_transitions = number_of_energy_levels-1
    energy_states = np.arange(number_of_energy_levels) - spin_quantum_number
    transitions = [ [energy_states[i], energy_states[i+1]] for i in range(number_of_transitions)]
    # print(transitions)
    cdef np.ndarray[double, ndim=1] transition_array = np.asarray(transitions).ravel()

    # gyromagnetic ratio
    cdef double larmor_frequency = spectrum['gyromagnetic_ratio']*B0

# ---------------------------------------------------------------------
# dimension axis ______________________________________________________
    cdef int number_of_points = spectrum['number_of_points']
    cdef double spectral_width = spectrum['spectral_width']
    cdef double increment = spectral_width/number_of_points
    cdef double reference_offset = spectrum['reference_offset']

    if verbose in [1, 11]:
        print ((f'Detecting {isotope}(I={spin_quantum_number}, '
                f'precession frequency = {larmor_frequency} MHz) isotope.'))
        print ((f'Recording {isotope} spectrum with {number_of_points} '
                f'points over a {spectral_width} Hz bandwidth and a '
                f'reference offset of {reference_offset} Hz.'))

    reference_offset -= spectral_width/2.0

    # CSA
    cdef int number_of_sites
    cdef np.ndarray[double] iso_n
    cdef np.ndarray[double] aniso_n
    cdef np.ndarray[double] eta_n

    # quad
    cdef np.ndarray[double] Cq_e
    cdef np.ndarray[double] eta_e

    cdef np.ndarray[double] D_c

    cdef int i, trans__
    cdef np.ndarray[double, ndim=1] amp
    cdef double cpu_time_ = 0.0
    amp1 = np.zeros(number_of_points, dtype=np.float64)

    # print (amp1)
    # octahedral power orientation averaging
    cdef unsigned int n_orientations = int((geodesic_polyhedron_frequency+1) * (geodesic_polyhedron_frequency+2)/2)
    cdef np.ndarray[double] cosAlpha = np.zeros(n_orientations, dtype=np.float64)
    cdef np.ndarray[double] cosBeta = np.zeros(n_orientations, dtype=np.float64)
    cdef np.ndarray[double] amp_orientation = np.zeros(n_orientations, dtype=np.float64)

    # print (cosAlpha)

    __powder_averaging_setup(geodesic_polyhedron_frequency, &cosAlpha[0],
                             &cosBeta[0], &amp_orientation[0], 1)
    amp_orientation*=(1./number_of_sidebands)

    list_index_isotopomer = []
    # print (cosAlpha)
    # print(cosBeta)
    # print(amp_orientation)
    # ---------------------------------------------------------------------
    # sample _______________________________________________________________
    for index_isotopomer, isotopomer in enumerate(isotopomers):
        abundance = isotopomer['abundance']
        sub_sites = [site for site in isotopomer['sites'] if site['isotope_symbol'] == isotope]

        number_of_sites= len(sub_sites)

        # site specification
        # CSA
        iso_n = np.empty(number_of_sites, dtype=np.float64)
        aniso_n = np.empty(number_of_sites, dtype=np.float64)
        eta_n = np.empty(number_of_sites, dtype=np.float64)

        # quad
        Cq_e = np.zeros(number_of_sites, dtype=np.float64)
        eta_e = np.zeros(number_of_sites, dtype=np.float64)

        D_c = np.zeros(number_of_sites, dtype=np.float64)

        # cdef int i, size = len(sample)

        amp = np.zeros(number_of_points*number_of_sites)

        if number_of_sites > 0: list_index_isotopomer.append(index_isotopomer)
        for i in range(number_of_sites):
            site = sub_sites[i]

            # CSA tensor
            iso = site['isotropic_chemical_shift']
            aniso = site['shielding_symmetric']['anisotropy']
            eta = site['shielding_symmetric']['asymmetry']

            if verbose in [1, 11]:
                text = ((
                    f"\n{isotope} site {i} from isotopomer {index_isotopomer} "
                    f"@ {abundance*100}% abundance"
                ))
                len_ = len(text)
                print(text)
                print(f"{'-'*(len_-1)}")
                print(f'isotropic chemical shift = {str(iso)}')
                print(f'chemical shift anisotropy = {str(aniso)}')
                print(f'chemical shift asymmetry = {eta}')

            # CSA tensor to Hz
            if iso.unit.physical_type == 'frequency':
                iso_n[i] = iso.to('Hz').value
            else:
                iso_n[i] = iso.value*larmor_frequency
            if aniso.unit.physical_type == 'frequency':
                aniso_n[i] = aniso.to('Hz').value
            else:
                aniso_n[i] = aniso.value*larmor_frequency
            eta_n[i] = eta

            # quad tensor
            # if spin.electric_quadrupole_tensor is not ():
            #     Cq_e[i] = spin.electric_quadrupole_tensor.Cq.to('Hz').value
            #     eta_e[i] = spin.electric_quadrupole_tensor.eta

        # print('\n')
        # print(transition_array)
        # print(number_of_transitions)



        for trans__ in range(number_of_transitions):
            # spin transitions
            # m_initial = transition_array[i][0]
            # m_final   = transition_array[i][1]
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
                    &iso_n[0],
                    &aniso_n[0],
                    &eta_n[0],

                    # quad tensor information
                    &Cq_e[0],
                    &eta_e[0],
                    1, # quad_second_order_c,

                    # Dipolar couplings
                    &D_c[0],

                    # spin rate, spin angle and number spinning sidebands
                    number_of_sidebands,
                    sample_rotation_frequency,
                    rotor_angle,

                    &transition_array[2*trans__],
                    # Euler angle -> principal to molecular frame
                    # &omega_PM_c[0],

                    # Euler angles for powder averaging scheme

                    # powder orientation averager
                    n_orientations,
                    &cosAlpha[0],
                    &cosBeta[0],
                    &amp_orientation[0],
                    geodesic_polyhedron_frequency,

                    number_of_sites
                    )


        amp1 += amp.reshape(number_of_sites, number_of_points).sum(axis=0)*abundance

    if verbose in [11]:
        print(f'\nExecution time {cpu_time_} s')
    freq = (np.arange(number_of_points)*increment + reference_offset)
    return freq, amp1, larmor_frequency, list_index_isotopomer




@cython.boundscheck(False)
@cython.wraparound(False)
def _one_d_simulator(
        # spectrum information
        double reference_offset,
        double increment,
        int number_of_points,

        float quantum_number = 0.5,
        float larmor_frequency = 0.0,

        # CSA tensor information
        isotropic_chemical_shift = None,
        chemical_shift_anisotropy = None,
        chemical_shift_asymmetry = None,

        # quad tensor information
        quadrupolar_coupling_constant = None,
        quadrupolar_asymmetry = None,
        second_order_quad = 1,

        # dipolar coupling
        D = None,

        # spin rate, spin angle and number spinning sidebands
        int number_of_sidebands = 96,
        double sample_rotation_frequency = 0.0,
        rotor_angle = None,

        m_final = 0.5,
        m_initial = -0.5,

        # Euler angle -> principal to molecular frame
        # omega_PM=None,

        # Euler angles for powder averaging scheme
        int averaging_scheme=0,
        int averaging_size=64):


    if isotropic_chemical_shift is None:
        isotropic_chemical_shift = 0
    isotropic_chemical_shift = np.asarray([isotropic_chemical_shift], dtype=np.float64).ravel()
    cdef number_of_sites = isotropic_chemical_shift.size
    cdef np.ndarray[double, ndim=1] isotropic_chemical_shift_c = isotropic_chemical_shift

    if quantum_number > 0.5 and larmor_frequency == 0.0:
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


    # if omega_PM is None:
    #     omega_PM = np.zeros(3)
    # cdef np.ndarray[double, ndim=1] omega_PM_c = np.asarray(omega_PM, dtype=np.float64).ravel()

    if rotor_angle is None:
        rotor_angle = 54.735
    cdef double rotor_angle_c = np.pi*rotor_angle/180.

    # if sample_rotation_frequency == 0:
    #     rotor_angle_c = 0
    #     sample_rotation_frequency = 1000000
    print(rotor_angle_c, sample_rotation_frequency)

    cdef second_order_quad_c = second_order_quad
    # if quad_second_order:
    #     quad_second_order_c = 1

    cdef np.ndarray[double, ndim=1] transition_c = np.asarray([m_initial, m_final], dtype=np.float64)

    cdef np.ndarray[double, ndim=1] amp = np.zeros(number_of_points * number_of_sites)
    cdef double cpu_time_
    # # print('rotor_angle in rad', rotor_angle_in_rad)
    # # print('transition', transition_c)

    # octahedral power orientation averaging
    cdef unsigned int n_orientations = int((averaging_size+1) * (averaging_size+2)/2)
    cdef np.ndarray[double, ndim=1] cosAlpha = np.empty(n_orientations, dtype=np.float64)
    cdef np.ndarray[double, ndim=1] cosBeta = np.empty(n_orientations, dtype=np.float64)
    cdef np.ndarray[double, ndim=1] amp_orientation = np.empty(n_orientations, dtype=np.float64)
    __powder_averaging_setup(averaging_size, &cosAlpha[0],
                             &cosBeta[0], &amp_orientation[0], 1)


    if averaging_scheme == 0:
        clib.spinning_sideband_core(
                # spectrum information and related amplitude
                &amp[0],
                &cpu_time_,
                reference_offset,
                increment,
                number_of_points,

                quantum_number,
                larmor_frequency,

                # CSA tensor information
                &isotropic_chemical_shift_c[0],
                &chemical_shift_anisotropy_c[0],
                &chemical_shift_asymmetry_c[0],

                # quad tensor information
                &quadrupolar_coupling_constant_c[0],
                &quadrupolar_asymmetry_c[0],
                second_order_quad_c,

                # Dipolar couplings
                &D_c[0],

                # spin rate, spin angle and number spinning sidebands
                number_of_sidebands,
                sample_rotation_frequency,
                rotor_angle_c,

                &transition_c[0],
                # Euler angle -> principal to molecular frame
                # &omega_PM_c[0],

                # Euler angles for powder averaging scheme

                # powder orientation averager
                n_orientations,
                &cosAlpha[0],
                &cosBeta[0],
                &amp_orientation[0],
                averaging_size,

                number_of_sites
                )

    # else:
    #     clib.lineshape_cas_spinning_sideband_cython_angles(
    #             # spectrum information and related amplitude
    #             &amp[0],
    #             &cpu_time_,
    #             spectral_start,
    #             spectral_increment,
    #             number_of_points,

    #             # CSA tensor information
    #             iso,
    #             aniso,
    #             eta,

    #             # spin rate, spin angle and number spinning sidebands
    #             ph_step,
    #             spin_frequency,
    #             rotor_angle_c,

    #             # Euler angle -> principal to molecular frame
    #             &omega_PM_c[0],

    #             # Euler angles for powder averaging scheme
    #             averaging_scheme,
    #             averaging_size
    #             )

    freq = np.arange(number_of_points)*increment + reference_offset

    return freq, amp, cpu_time_
