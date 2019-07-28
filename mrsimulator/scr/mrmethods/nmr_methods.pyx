cimport nmr_methods as clib
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
        polyhedra approach a sphere geometry. For line-shape simulation, a higher
        frequency will result in a better powder averaging.
        The default value is 72.
        Read more on the `Geodesic polyhedron <https://en.wikipedia.org/wiki/Geodesic_polyhedron>`_.
    """
# ---------------------------------------------------------------------
# spectrum ________________________________________________________
    cdef double sample_rotation_frequency_in_Hz = spectrum['rotor_frequency']
    cdef double rotor_angle_in_rad = spectrum['rotor_angle']
    B0 = spectrum['magnetic_flux_density']

    if verbose in [1, 11]:
        print ('Setting up the virtual NMR spectrometer')
        print ('---------------------------------------')
        print (f'Adjusting the magnetic flux density to {B0} T.')
        _angle = spectrum['rotor_angle']
        print (f'Setting rotation angle to {_angle} rad.')
        print (f'Setting rotation frequency to {sample_rotation_frequency_in_Hz} Hz.')

# ---------------------------------------------------------------------
# spin observed _______________________________________________________
    # obs_dict = __get_spin_attribute__[detect]
    isotope = spectrum['isotope']

    # spin quantum number of the observed spin
    cdef double spin_quantum_number = spectrum['spin']

    # transitions of the observed spin
    cdef np.ndarray[double, ndim=1] transition_array
    cdef int number_of_transitions
    transition_array = np.asarray(spectrum['transitions']).ravel()
    number_of_transitions = int(transition_array.size/2)
    # else:
    #     energy_level_count = int(2*spin_quantum_number+1)
    #     number_of_transitions = energy_level_count-1
    #     energy_states = np.arange(energy_level_count) - spin_quantum_number
    #     transitions = [ [energy_states[i], energy_states[i+1]] for i in range(number_of_transitions)]
    #     transition_array = np.asarray(transitions).ravel()

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
    cdef np.ndarray[double] zeta_n
    cdef np.ndarray[double] eta_n
    cdef np.ndarray[double] ori_n

    # quad
    cdef np.ndarray[double] Cq_e
    cdef np.ndarray[double] eta_e
    cdef np.ndarray[double] ori_e

    cdef np.ndarray[double] D_c

    cdef int i, trans__
    cdef np.ndarray[double, ndim=1] amp
    cdef double cpu_time_ = 0.0
    amp1 = np.zeros(number_of_points, dtype=np.float64)

    cdef clib.isotopomer_ravel isotopomer_struct
    # cdef clib.isotopomers_list *isotopomers_list_c

    list_index_isotopomer = []
    # ---------------------------------------------------------------------
    # sample _______________________________________________________________
    for index_isotopomer, isotopomer in enumerate(isotopomers):
        abundance = isotopomer['abundance']
        sub_sites = [site for site in isotopomer['sites'] if site['isotope_symbol'] == isotope]

        number_of_sites= len(sub_sites)

        # site specification
        # CSA
        iso_n = np.empty(number_of_sites, dtype=np.float64)
        zeta_n = np.empty(number_of_sites, dtype=np.float64)
        eta_n = np.empty(number_of_sites, dtype=np.float64)
        ori_n = np.zeros(3*number_of_sites, dtype=np.float64)

        # quad
        Cq_e = np.zeros(number_of_sites, dtype=np.float64)
        eta_e = np.zeros(number_of_sites, dtype=np.float64)
        ori_e = np.zeros(3*number_of_sites, dtype=np.float64)

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
                print(f'Isotropic chemical shift = {str(iso)}')
                print(f'Shielding anisotropy = {str(aniso)}')
                print(f'Shielding asymmetry = {eta}')

            # CSA tensor in Hz
            if iso.unit.physical_type == 'frequency':
                iso_n[i] = iso.to('Hz').value
            else:
                iso_n[i] = iso.value*larmor_frequency
            if aniso.unit.physical_type == 'frequency':
                zeta_n[i] = aniso.to('Hz').value
            else:
                zeta_n[i] = aniso.value*larmor_frequency
            eta_n[i] = eta

            # quad tensor
            if spin_quantum_number > 0.5:
                Cq_e[i] = site['quadrupolar']['anisotropy']
                eta_e[i] = site['quadrupolar']['asymmetry']

                if verbose in [1, 11]:
                    print(f'Quadrupolar coupling constant = {Cq_e[i]/1e6} MHz')
                    print(f'Quadrupolar asymmetry = {eta}')

        # print(transition_array)
        # print(number_of_transitions)

        isotopomer_struct.number_of_sites = number_of_sites
        isotopomer_struct.spin = spin_quantum_number
        isotopomer_struct.larmor_frequency = larmor_frequency

        isotopomer_struct.isotropic_chemical_shift_in_Hz = &iso_n[0]
        isotopomer_struct.shielding_anisotropy_in_Hz = &zeta_n[0]
        isotopomer_struct.shielding_asymmetry = &eta_n[0]
        isotopomer_struct.shielding_orientation = &ori_n[0]

        isotopomer_struct.quadrupolar_constant_in_Hz = &Cq_e[0]
        isotopomer_struct.quadrupolar_asymmetry = &eta_e[0]
        isotopomer_struct.quadrupolar_orientation = &ori_e[0]

        isotopomer_struct.dipolar_couplings = &D_c[0]

        # isotopomers_list_c[i] = isotopomer_struct
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

                    &isotopomer_struct,

                    1, # quad_second_order_c,
                    0, # turn off quad second order isotropic contribution

                    # spin rate, spin angle and number spinning sidebands
                    number_of_sidebands,
                    sample_rotation_frequency_in_Hz,
                    rotor_angle_in_rad,

                    &transition_array[2*trans__],
                    geodesic_polyhedron_frequency)


        amp1 += amp.reshape(number_of_sites, number_of_points).sum(axis=0)*abundance

    if verbose in [11]:
        print(f'\nExecution time {cpu_time_} s')
    freq = (np.arange(number_of_points)*increment + reference_offset)
    return freq, amp1, larmor_frequency, list_index_isotopomer
