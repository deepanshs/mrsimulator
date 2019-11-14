cimport nmr_methods as clib
from libcpp cimport bool as bool_t
from numpy cimport ndarray
import numpy as np
import cython

__author__ = "Deepansh J. Srivastava"
__email__ = ["srivastava.89@osu.edu", "deepansh2012@gmail.com"]


@cython.boundscheck(False)
@cython.wraparound(False)
def one_d_spectrum(dimension,
       list isotopomers,
       int verbose=0,
       int number_of_sidebands=90,
       int integration_density=72,
       bool_t individual_spectrum=False,
       int integration_volume=1):
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
    :ivar integration_density:
        The value is an integer which represents the frequency of class I
        geodesic polyhedra. These polyhedra are used in calculating the
        spherical average. Presently we only use octahedral as the frequency1
        polyhedra. As the frequency of the geodesic polyhedron increases, the
        polyhedra approach a sphere geometry. For line-shape simulation, a higher
        frequency will result in a better powder averaging.
        The default value is 72.
        Read more on the `Geodesic polyhedron <https://en.wikipedia.org/wiki/Geodesic_polyhedron>`_.
    :ivar individual_spectrum:
        A boolean. If true, returns an ordered list of spectrum corresponding
        to the ordered list of isotopomers.
    """

# ---------------------------------------------------------------------
# observed spin _______________________________________________________
    dimension = dimension[0]
    isotope = dimension.isotope
    # spin quantum number of the observed spin
    cdef double spin_quantum_number = dimension.spin

    # gyromagnetic ratio
    cdef double larmor_frequency = dimension.larmor_frequency
    cdef double factor = 1.0
    if larmor_frequency < 0.0:
        factor = -1.0

    # if verbose in [1, 11]:
    #     print(f'Simulating {isotope} (I={spin_quantum_number})')
    #     print(f'Larmor frequency (ω0 = - γ B0) = {larmor_frequency/1.0e6} MHz')
    #     print((f'Recording {isotope} spectrum with {number_of_points} '
    #             f'points over {spectral_width} Hz bandwidth'))
    #     print((f"and a reference offset of {dimension['reference_offset']} Hz."))

    # transitions of the observed spin
    cdef ndarray[double, ndim=1] transition_array
    cdef int number_of_transitions
    transition_array = np.asarray([-0.5, 0.5])
    number_of_transitions = int(transition_array.size/2)
    # else:
    #     energy_level_count = int(2*spin_quantum_number+1)
    #     number_of_transitions = energy_level_count-1
    #     energy_states = np.arange(energy_level_count) - spin_quantum_number
    #     transitions = [ [energy_states[i], energy_states[i+1]] for i in range(number_of_transitions)]
    #     transition_array = np.asarray(transitions).ravel()


# Generate the dimension coordinates __________________________________________________
    cdef int number_of_points = dimension.number_of_points
    cdef double spectral_width = dimension.spectral_width
    cdef double increment = spectral_width/number_of_points
    cdef double reference_offset = dimension.reference_offset *factor

    offset = increment*int(number_of_points/2.0)
    reference_offset -= offset
    # freq = dimension.coordinates_ppm

# create MRS_plan _____________________________________________________________________
    cdef double sample_rotation_frequency_in_Hz = dimension.rotor_frequency
    cdef double rotor_angle_in_rad = dimension.rotor_angle

    # B0 = dimension.magnetic_flux_density

    # if verbose in [1, 11]:
    #     text = "`one_d_spectrum` method simulation parameters."
    #     len_ = len(text)
    #     print(text)
    #     print(f"{'-'*(len_-1)}")
    #     print (f'Macroscopic magnetic flux density (B0) = {B0} T')
    #     print (f'Sample rotation angle is (θ) = {rotor_angle_in_rad} rad')
    #     print (f'Sample rotation frequency (𝜈r) = {sample_rotation_frequency_in_Hz} Hz')

# sites _______________________________________________________________________________
    # CSA
    cdef int number_of_sites
    cdef ndarray[double] iso_n
    cdef ndarray[double] zeta_n
    cdef ndarray[double] eta_n
    cdef ndarray[double] ori_n

    # quad
    cdef ndarray[double] Cq_e
    cdef ndarray[double] eta_e
    cdef ndarray[double] ori_e

    cdef ndarray[double] D_c

    cdef int i, trans__
    cdef ndarray[double, ndim=1] amp
    amp1 = np.zeros(number_of_points, dtype=np.float64)
    amp_individual = []

    cdef clib.isotopomer_ravel isotopomer_struct
    # cdef clib.isotopomers_list *isotopomers_list_c

    # ---------------------------------------------------------------------
    # sample _______________________________________________________________
    for index_isotopomer, isotopomer in enumerate(isotopomers):
        abundance = isotopomer['abundance']
        sub_sites = [site for site in isotopomer['sites'] if site['isotope'] == isotope]

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

        for i in range(number_of_sites):
            site = sub_sites[i]


            # CSA tensor
            iso = site['isotropic_chemical_shift']
            if iso is None: iso = 0.0

            shielding = site['shielding_symmetric']
            if shielding is not None:
                zeta = shielding['zeta']
                if zeta is None: zeta = 0.0
                eta = shielding['eta']
                if eta is None: eta = 0.0
                alpha, beta, gamma = shielding['alpha'], shielding['beta'], shielding['gamma']
                if alpha is None: alpha = 0.0
                if beta is None: beta = 0.0
                if gamma is None: gamma = 0.0
            else:
                zeta, eta, alpha, beta, gamma = 0.0, 0.0, 0.0, 0.0, 0.0

            # if verbose in [1, 11]:
            #     text = ((
            #         f"\n{isotope} site {i} from isotopomer {index_isotopomer} "
            #         f"@ {abundance}% abundance"
            #     ))
            #     len_ = len(text)
            #     print(text)
            #     print(f"{'-'*(len_-1)}")
            #     print(f'Isotropic chemical shift (δ) = {str(1e6*iso/larmor_frequency)} ppm')
            #     print(f'Shielding anisotropy (ζ) = {str(1e6*zeta/larmor_frequency)} ppm')
            #     print(f'Shielding asymmetry (η) = {eta}')
            #     print(f'Shielding orientation = [alpha = {alpha}, beta = {beta}, gamma = {gamma}]')

            # CSA tensor in Hz
            iso_n[i] = iso #* factor
            zeta_n[i] = zeta #* factor
            eta_n[i] = eta
            ori_n[3*i:3*i+3] = [alpha, beta, gamma]

            # quad tensor
            if spin_quantum_number > 0.5:
                quad = site['quadrupolar']
                if quad is not None:
                    Cq = quad['Cq']
                    if Cq is None: Cq = 0.0
                    eta = quad['eta']
                    if eta is None: eta = 0.0
                    alpha, beta, gamma = quad['alpha'], quad['beta'], quad['gamma']
                    if alpha is None: alpha = 0.0
                    if beta is None: beta = 0.0
                    if gamma is None: gamma = 0.0
                else:
                    Cq, eta, alpha, beta, gamma = 0.0, 0.0, 0.0, 0.0, 0.0

                Cq_e[i] = Cq #* factor
                eta_e[i] = eta
                ori_e[3*i:3*i+3] = [alpha, beta, gamma]

                # if verbose in [1, 11]:
                #     print(f'Quadrupolar coupling constant (Cq) = {Cq_e[i]/1e6} MHz')
                #     print(f'Quadrupolar asymmetry (η) = {eta}')
                #     print(f'Quadrupolar orientation = [alpha = {alpha}, beta = {beta}, gamma = {gamma}]')


        if number_of_sites != 0:
            isotopomer_struct.number_of_sites = number_of_sites
            isotopomer_struct.spin = spin_quantum_number
            isotopomer_struct.larmor_frequency = larmor_frequency

            isotopomer_struct.isotropic_chemical_shift_in_Hz = &iso_n[0]
            isotopomer_struct.shielding_anisotropy_in_Hz = &zeta_n[0]
            isotopomer_struct.shielding_asymmetry = &eta_n[0]
            isotopomer_struct.shielding_orientation = &ori_n[0]

            isotopomer_struct.quadrupole_coupling_constant_in_Hz = &Cq_e[0]
            isotopomer_struct.quadrupole_asymmetry = &eta_e[0]
            isotopomer_struct.quadrupole_orientation = &ori_e[0]

            isotopomer_struct.dipolar_couplings = &D_c[0]

            # isotopomers_list_c[i] = isotopomer_struct
            for trans__ in range(number_of_transitions):
                # spin transitions
                # m_initial = transition_array[i][0]
                # m_final   = transition_array[i][1]
                clib.spinning_sideband_core(
                    # spectrum information and related amplitude
                    &amp[0],
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
                    integration_density,
                    integration_volume,          # 0-octant, 1-hemisphere, 2-sphere.
                    )

            temp = amp.reshape(number_of_sites, number_of_points).sum(axis=0)*abundance

            ## reverse the spectrum if larmor_frequency is negative
            if larmor_frequency < 0:
                if number_of_points % 2 == 0:
                    temp[1:] = temp[1:][::-1]
                else:
                    temp = temp[::-1]

            if individual_spectrum:
                amp_individual.append(temp)
            else:
                amp1 += temp
        else:
            if individual_spectrum:
                amp_individual.append([])

    if individual_spectrum:
        amp1 = amp_individual


    return amp1
