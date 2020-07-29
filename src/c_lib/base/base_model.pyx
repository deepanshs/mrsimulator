cimport base_model as clib
from libcpp cimport bool as bool_t
from numpy cimport ndarray
import numpy as np
import cython
from mrsimulator import sandbox as sb

__author__ = "Deepansh J. Srivastava"
__email__ = "deepansh2012@gmail.com"


@cython.boundscheck(False)
@cython.wraparound(False)
def one_d_spectrum(method,
       list spin_systems,
       int verbose=0,
       int number_of_sidebands=90,
       unsigned int integration_density=72,
       unsigned int decompose_spectrum=0,
       unsigned int integration_volume=1,
       bool_t interpolation=True):
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
        polyhedra approach a sphere geometry. A higher frequency will result in a
        better powder averaging. The default value is 72.
        Read more on the `Geodesic polyhedron <https://en.wikipedia.org/wiki/Geodesic_polyhedron>`_.
    :ivar decompose_spectrum:
        An unsigned integer. When value is 0, the spectum is a sum of spectrum from all
        spin systems. If value is 1, spectrum from individual spin systems is stored
        separately.
    """

    cdef int transition_increment


# ---------------------------------------------------------------------
# observed spin _______________________________________________________
    # dimension = dimension[0]
    channel = method.channels[0].symbol
    # spin quantum number of the observed spin
    cdef double spin_quantum_number = method.channels[0].spin

    # gyromagnetic ratio
    # cdef double larmor_frequency =
    cdef gyromagnetic_ratio = method.channels[0].gyromagnetic_ratio
    cdef double factor = 1.0
    if gyromagnetic_ratio > 0.0:
        factor = -1.0

    # if verbose in [1, 11]:
    #     print(f'Simulating {isotope} (I={spin_quantum_number})')
    #     print(f'Larmor frequency (ω0 = - γ B0) = {larmor_frequency/1.0e6} MHz')
    #     print((f'Recording {isotope} spectrum with {number_of_points} '
    #             f'points over {spectral_width} Hz bandwidth'))
    #     print((f"and a reference offset of {dimension['reference_offset']} Hz."))

    # transitions of the observed spin
    cdef ndarray[float, ndim=1] transition_array
    cdef int number_of_transitions
    # transition_array = np.asarray([-0.5, 0.5]).ravel()
    # number_of_transitions = int(transition_array.size/2)
    # else:
    #     energy_level_count = int(2*spin_quantum_number+1)
    #     number_of_transitions = energy_level_count-1
    #     energy_states = np.arange(energy_level_count) - spin_quantum_number
    #     transitions = [ [energy_states[i], energy_states[i+1]] for i in range(number_of_transitions)]
    #     transition_array = np.asarray(transitions).ravel()


# -------------------------------------------------------------------------------------
# schemes

    cdef bool_t allow_fourth_rank = 0
    if spin_quantum_number > 0.5:
        allow_fourth_rank = 1

    # if sample_rotation_frequency_in_Hz < 1.0e-3:
    #     sample_rotation_frequency_in_Hz = 1.0e9
    #     rotor_angle_in_rad = 0.0
    #     number_of_sidebands = 1

# create averaging scheme _____________________________________________________
    cdef clib.MRS_averaging_scheme *the_averaging_scheme
    the_averaging_scheme = clib.MRS_create_averaging_scheme(
        integration_density=integration_density, allow_fourth_rank=allow_fourth_rank, integration_volume=integration_volume
    )

    max_n_sidebands = number_of_sidebands

# create sequences____________________________________________________________

    cdef int n_sequence = len(method.spectral_dimensions)
    total_n_points = 1
    cdef ndarray[int] n_event
    cdef ndarray[double] magnetic_flux_density_in_T
    cdef ndarray[double] srfiH
    cdef ndarray[double] rair
    cdef ndarray[int] cnt
    cdef ndarray[double] coord_off
    cdef ndarray[double] incre

    Bo = []
    vr = []
    th = []
    event_i = []
    count = []
    increment = []
    coordinates_offset = []

    prev_n_sidebands = 0
    for i, seq in enumerate(method.spectral_dimensions):
        for event in seq.events:
            if event.rotor_frequency < 1.0e-3:
                sample_rotation_frequency_in_Hz = 1.0e9
                rotor_angle_in_rad = 0.0
                number_of_sidebands = 1
                if prev_n_sidebands == 0: prev_n_sidebands = 1
            else:
                sample_rotation_frequency_in_Hz = event.rotor_frequency
                rotor_angle_in_rad = event.rotor_angle
                if prev_n_sidebands == 0: prev_n_sidebands = number_of_sidebands

            if prev_n_sidebands != number_of_sidebands:
                raise ValueError(
                    (
                        'The library does not support spectral dimensions containing '
                        'both zero and non-zero rotor frequencies. Consider using a '
                        'smaller value instead of zero.'
                    )
                )

            Bo.append(event.magnetic_flux_density)  # in T
            vr.append(sample_rotation_frequency_in_Hz) # in Hz
            th.append(rotor_angle_in_rad) # in rad

        total_n_points *= seq.count

        count.append(seq.count)
        offset = seq.spectral_width / 2.0
        coordinates_offset.append(-seq.reference_offset * factor - offset)
        increment.append(seq.spectral_width / seq.count)
        event_i.append(len(seq.events))

        seq.origin_offset = np.abs(Bo[0] * gyromagnetic_ratio * 1e6)

    magnetic_flux_density_in_T = np.asarray(Bo, dtype=np.float64)
    srfiH = np.asarray(vr, dtype=np.float64)
    rair = np.asarray(th, dtype=np.float64)
    cnt = np.asarray(count, dtype=np.int32)
    incre = np.asarray(increment, dtype=np.float64)
    coord_off = np.asarray(coordinates_offset, dtype=np.float64)
    n_event = np.asarray(event_i, dtype=np.int32)
    # magnetic_flux_density_in_T = np.asarray([dimension.magnetic_flux_density], dtype=np.float64)
    # srfiH = np.asarray([sample_rotation_frequency_in_Hz], dtype=np.float64)
    # rair = np.asarray([rotor_angle_in_rad], dtype=np.float64)
    # create spectral_dimensions



    the_sequence = clib.MRS_create_plans_for_sequence(
        the_averaging_scheme, &cnt[0], &coord_off[0],
        &incre[0], &magnetic_flux_density_in_T[0], &srfiH[0],
        &rair[0], &n_event[0], n_sequence, number_of_sidebands)

# normalization factor for the spectrum
    norm = np.prod(incre)

# create fftw scheme __________________________________________________________

    cdef clib.MRS_fftw_scheme *the_fftw_scheme
    the_fftw_scheme = clib.create_fftw_scheme(the_averaging_scheme.total_orientations, number_of_sidebands)
# # _____________________________________________________________________________


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
    cdef ndarray[float] spin_i
    cdef ndarray[double] gyromagnetic_ratio_i

    cdef ndarray[double] iso_n
    cdef ndarray[double] zeta_n
    cdef ndarray[double] eta_n
    cdef ndarray[double] ori_n

    # quad
    cdef ndarray[double] Cq_e
    cdef ndarray[double] eta_e
    cdef ndarray[double] ori_e

    cdef ndarray[double] D_c

    cdef int trans__, pathway_increment, pathway_count, transition_count_per_pathway
    cdef ndarray[double, ndim=1] amp
    amp1 = np.zeros(total_n_points, dtype=np.float64)
    amp_individual = []

    cdef clib.isotopomer_ravel isotopomer_struct

    index_ = []
    # cdef clib.isotopomers_list *isotopomers_list_c

    # ---------------------------------------------------------------------
    # sample _______________________________________________________________
    for index, spin_sys in enumerate(spin_systems):
        abundance = spin_sys.abundance
        isotopes = [site.isotope.symbol for site in spin_sys.sites]
        if channel not in isotopes:
            continue

        # sub_sites = [site for site in spin_sys.sites if site.isotope.symbol == isotope]
        index_.append(index)
        number_of_sites= len(spin_sys.sites)

        if number_of_sites > 2:
            continue
        # site specification
        # CSA

        spin_i = np.empty(number_of_sites, dtype=np.float32)
        gyromagnetic_ratio_i = np.empty(number_of_sites, dtype=np.float64)

        iso_n = np.empty(number_of_sites, dtype=np.float64)
        zeta_n = np.empty(number_of_sites, dtype=np.float64)
        eta_n = np.empty(number_of_sites, dtype=np.float64)
        ori_n = np.zeros(3*number_of_sites, dtype=np.float64)

        # quad
        Cq_e = np.zeros(number_of_sites, dtype=np.float64)
        eta_e = np.zeros(number_of_sites, dtype=np.float64)
        ori_e = np.zeros(3*number_of_sites, dtype=np.float64)

        # for n sites, coupling grows as sum_{i=1}^{n-1}(i)
        D_c = np.zeros(number_of_sites, dtype=np.float64)

        # cdef int i, size = len(sample)

        amp = np.zeros(total_n_points)

        for i in range(number_of_sites):
            site = spin_sys.sites[i]
            spin_i[i] = site.isotope.spin
            gyromagnetic_ratio_i[i] = site.isotope.gyromagnetic_ratio


            # CSA tensor
            iso = site.isotropic_chemical_shift
            if iso is None: iso = 0.0

            shielding = site.shielding_symmetric
            if shielding is not None:
                zeta = shielding.zeta
                if zeta is None: zeta = 0.0
                eta = shielding.eta
                if eta is None: eta = 0.0
                alpha, beta, gamma = shielding.alpha, shielding.beta, shielding.gamma
                if alpha is None: alpha = 0.0
                if beta is None: beta = 0.0
                if gamma is None: gamma = 0.0
            else:
                zeta, eta, alpha, beta, gamma = 0.0, 0.0, 0.0, 0.0, 0.0

            # if verbose in [1, 11]:
            #     text = ((
            #         f"\n{isotope} site {i} from spin system {index_isotopomer} "
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
            iso_n[i] = iso
            zeta_n[i] = zeta
            eta_n[i] = eta
            ori_n[3*i:3*i+3] = [alpha, beta, gamma]

            # quad tensor
            if spin_quantum_number > 0.5:
                quad = site.quadrupolar
                if quad is not None:
                    Cq = quad.Cq
                    if Cq is None: Cq = 0.0
                    eta = quad.eta
                    if eta is None: eta = 0.0
                    alpha, beta, gamma = quad.alpha, quad.beta, quad.gamma
                    if alpha is None: alpha = 0.0
                    if beta is None: beta = 0.0
                    if gamma is None: gamma = 0.0
                else:
                    Cq, eta, alpha, beta, gamma = 0.0, 0.0, 0.0, 0.0, 0.0

                Cq_e[i] = Cq
                eta_e[i] = eta
                ori_e[3*i:3*i+3] = [alpha, beta, gamma]

                # if verbose in [1, 11]:
                #     print(f'Quadrupolar coupling constant (Cq) = {Cq_e[i]/1e6} MHz')
                #     print(f'Quadrupolar asymmetry (η) = {eta}')
                #     print(f'Quadrupolar orientation = [alpha = {alpha}, beta = {beta}, gamma = {gamma}]')

        if number_of_sites != 0:
            transition_pathway = spin_sys.transition_pathways
            if transition_pathway is None:
                transition_pathway = method.get_transition_pathways(spin_sys)
            transition_pathway = np.asarray(transition_pathway)
            pathway_count, transition_count_per_pathway = transition_pathway.shape

            # convert transition objects to list
            lst = []
            for item in transition_pathway.ravel():
                lst += item.tolist()

            transition_array = np.asarray(lst, dtype=np.float32).ravel()
            pathway_increment = 2*number_of_sites*transition_count_per_pathway

            # if spin_sys.transitions is not None:
            #     transition_array = np.asarray(
            #         spin_sys.transitions, dtype=np.float32
            #     ).ravel()
            # else:
            #     transition_array = np.asarray([0.5, -0.5], dtype=np.float32)

            # the number 2 is because of single site transition [mi, mf]
            # it dose not work for coupled sites.
            # transition_increment = 2*number_of_sites
            # number_of_transitions = int((transition_array.size)/transition_increment)

            isotopomer_struct.number_of_sites = number_of_sites
            isotopomer_struct.spin = &spin_i[0]
            isotopomer_struct.gyromagnetic_ratio = &gyromagnetic_ratio_i[0]

            isotopomer_struct.isotropic_chemical_shift_in_ppm = &iso_n[0]
            isotopomer_struct.shielding_symmetric_zeta_in_ppm = &zeta_n[0]
            isotopomer_struct.shielding_symmetric_eta = &eta_n[0]
            isotopomer_struct.shielding_orientation = &ori_n[0]

            isotopomer_struct.quadrupolar_Cq_in_Hz = &Cq_e[0]
            isotopomer_struct.quadrupolar_eta = &eta_e[0]
            isotopomer_struct.quadrupolar_orientation = &ori_e[0]

            isotopomer_struct.dipolar_couplings = &D_c[0]

            # isotopomers_list_c[i] = isotopomer_struct
            for trans__ in range(pathway_count):
                clib.__mrsimulator_core(
                    # spectrum information and related amplitude
                    &amp[0],
                    &isotopomer_struct,
                    0, # turn off quad second order isotropic contribution
                    &transition_array[pathway_increment*trans__],
                    the_sequence,
                    n_sequence,
                    the_fftw_scheme,
                    the_averaging_scheme,
                    interpolation,
                    )

            temp = amp*abundance/norm

            ## reverse the spectrum if gyromagnetic ratio is positive.
            if gyromagnetic_ratio < 0:
                if total_n_points % 2 == 0:
                    temp[1:] = temp[1:][::-1]
                else:
                    temp = temp[::-1]

            if decompose_spectrum == 1:
                amp_individual.append(temp)
            else:
                amp1 += temp
        else:
            if decompose_spectrum == 1:
                amp_individual.append([])

    if decompose_spectrum == 1 and len(amp_individual) != 0:
        amp1 = amp_individual

    return amp1, index_
