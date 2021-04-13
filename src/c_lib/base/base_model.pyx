cimport base_model as clib
from libcpp cimport bool as bool_t
from numpy cimport ndarray
import numpy as np
import cython

__author__ = "Deepansh J. Srivastava"
__email__ = "srivastava.89@osu.edu"

clib.generate_table()

@cython.profile(False)
@cython.boundscheck(False)
@cython.wraparound(False)
def one_d_spectrum(method,
       list spin_systems,
       int verbose=0,
       unsigned int number_of_sidebands=90,
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
        An unsigned integer. When value is 0, the spectrum is a sum of spectrum from all
        spin systems. If value is 1, spectrum from individual spin systems is stored
        separately.
    """

# observed spin _______________________________________________________
    channel = method.channels[0].symbol
    # spin quantum number of the observed spin
    cdef double spin_quantum_number = method.channels[0].spin

    # gyromagnetic ratio
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
    cdef int transition_increment
    cdef ndarray[float, ndim=1] transition_pathway_c
    cdef int number_of_transitions
    # transition_pathway_c = np.asarray([-0.5, 0.5]).ravel()
    # number_of_transitions = int(transition_pathway_c.size/2)
    # else:
    #     energy_level_count = int(2*spin_quantum_number+1)
    #     number_of_transitions = energy_level_count-1
    #     energy_states = np.arange(energy_level_count) - spin_quantum_number
    #     transitions = [ [energy_states[i], energy_states[i+1]] for i in range(number_of_transitions)]
    #     transition_pathway_c = np.asarray(transitions).ravel()

    cdef bool_t allow_fourth_rank = 0
    if spin_quantum_number > 0.5:
        allow_fourth_rank = 1

# create averaging scheme _____________________________________________________
    cdef clib.MRS_averaging_scheme *the_averaging_scheme
    the_averaging_scheme = clib.MRS_create_averaging_scheme(
        integration_density=integration_density, allow_fourth_rank=allow_fourth_rank,
        integration_volume=integration_volume
    )

# create spectral dimensions _______________________________________________
    cdef int n_dimension = len(method.spectral_dimensions)
    # if n_sequence > 1:
    #     number_of_sidebands = 1

    max_n_sidebands = number_of_sidebands

    total_n_points = 1
    cdef ndarray[int] n_event
    cdef ndarray[double] magnetic_flux_density_in_T, frac
    cdef ndarray[double] srfiH
    cdef ndarray[double] rair
    cdef ndarray[int] cnt
    cdef ndarray[double] coord_off
    cdef ndarray[double] incre
    freq_contrib = np.asarray([])

    fr = []
    Bo = []
    vr = []
    th = []
    event_i = []
    count = []
    increment = []
    coordinates_offset = []

    prev_n_sidebands = 0
    for i, dim in enumerate(method.spectral_dimensions):
        for event in dim.events:
            freq_contrib = np.append(freq_contrib, event._freq_contrib_flags())
            if event.rotor_frequency < 1.0e-3:
                rotor_frequency_in_Hz = 1.0e9
                rotor_angle_in_rad = 0.0
                number_of_sidebands = 1
                if prev_n_sidebands == 0: prev_n_sidebands = 1
            else:
                rotor_frequency_in_Hz = event.rotor_frequency
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

            fr.append(event.fraction) # fraction
            Bo.append(event.magnetic_flux_density)  # in T
            vr.append(rotor_frequency_in_Hz) # in Hz
            th.append(rotor_angle_in_rad) # in rad

        total_n_points *= dim.count

        count.append(dim.count)
        offset = dim.spectral_width / 2.0
        coordinates_offset.append(-dim.reference_offset * factor - offset)
        increment.append(dim.spectral_width / dim.count)
        event_i.append(len(dim.events))

        dim.origin_offset = np.abs(Bo[0] * gyromagnetic_ratio * 1e6)

    frac = np.asarray(fr, dtype=np.float64)
    magnetic_flux_density_in_T = np.asarray(Bo, dtype=np.float64)
    srfiH = np.asarray(vr, dtype=np.float64)
    rair = np.asarray(th, dtype=np.float64)
    cnt = np.asarray(count, dtype=np.int32)
    incre = np.asarray(increment, dtype=np.float64)
    coord_off = np.asarray(coordinates_offset, dtype=np.float64)
    n_event = np.asarray(event_i, dtype=np.int32)

    # create spectral_dimensions
    dimensions = clib.MRS_create_dimensions(the_averaging_scheme, &cnt[0],
        &coord_off[0], &incre[0], &frac[0], &magnetic_flux_density_in_T[0],
        &srfiH[0], &rair[0], &n_event[0], n_dimension, number_of_sidebands)

# normalization factor for the spectrum
    norm = np.prod(incre)

# create fftw scheme __________________________________________________________
    cdef clib.MRS_fftw_scheme *the_fftw_scheme
    the_fftw_scheme = clib.create_fftw_scheme(the_averaging_scheme.total_orientations, number_of_sidebands)

# _____________________________________________________________________________

# frequency contrib
    cdef ndarray[bool_t] freq_contrib_c = np.asarray(freq_contrib, dtype=bool)

# affine transformation
    cdef ndarray[double] affine_matrix_c
    if method.affine_matrix is None:
        affine_matrix_c = np.asarray([1, 0, 0, 1], dtype=np.float64)
    else:
        increment_fraction = [incre/item for item in incre]
        matrix = method.affine_matrix.ravel() * np.asarray(increment_fraction).ravel()
        affine_matrix_c = np.asarray(matrix, dtype=np.float64)
        if affine_matrix_c[2] != 0:
            affine_matrix_c[2] /= affine_matrix_c[0]
            affine_matrix_c[3] -=  affine_matrix_c[1]*affine_matrix_c[2]

    # B0 = dimension.magnetic_flux_density

    # if verbose in [1, 11]:
    #     text = "`one_d_spectrum` method simulation parameters."
    #     len_ = len(text)
    #     print(text)
    #     print(f"{'-'*(len_-1)}")
    #     print (f'Macroscopic magnetic flux density (B0) = {B0} T')
    #     print (f'Sample rotation angle is (θ) = {rotor_angle_in_rad} rad')
    #     print (f'Sample rotation frequency (𝜈r) = {rotor_frequency_in_Hz} Hz')

# sites _______________________________________________________________________________
    p_isotopes = None

    cdef int number_of_sites, number_of_couplings, p_number_of_sites=0
    cdef ndarray[int] spin_index_ij
    cdef ndarray[float] spin_i
    cdef ndarray[double] gyromagnetic_ratio_i

    # CSA
    cdef ndarray[double] iso_n
    cdef ndarray[double] zeta_n
    cdef ndarray[double] eta_n
    cdef ndarray[double] ori_n

    # quad
    cdef ndarray[double] Cq_e
    cdef ndarray[double] eta_e
    cdef ndarray[double] ori_e

    # J-coupling
    cdef ndarray[double] iso_j
    cdef ndarray[double] zeta_j
    cdef ndarray[double] eta_j
    cdef ndarray[double] ori_j

    # quad
    cdef ndarray[double] D_d
    cdef ndarray[double] eta_d
    cdef ndarray[double] ori_d



    cdef int trans__, pathway_increment, pathway_count, transition_count_per_pathway
    cdef ndarray[double, ndim=1] amp
    amp1 = np.zeros(total_n_points, dtype=np.float64)
    amp_individual = []

    cdef clib.site_struct sites_c
    cdef clib.coupling_struct couplings_c

    index_ = []

    # -------------------------------------------------------------------------
    # sample __________________________________________________________________
    for index, spin_sys in enumerate(spin_systems):
        abundance = spin_sys.abundance
        isotopes = [site.isotope.symbol for site in spin_sys.sites]
        if channel not in isotopes:
            if decompose_spectrum == 1:
                amp_individual.append(np.zeros(method.shape()))
            continue

        # sub_sites = [site for site in spin_sys.sites if site.isotope.symbol == isotope]
        index_.append(index)
        number_of_sites = len(spin_sys.sites)

        # ------------------------------------------------------------------------
        #                          Site specification
        # ------------------------------------------------------------------------
        # CSA
        spin_i = np.empty(number_of_sites, dtype=np.float32)
        gyromagnetic_ratio_i = np.empty(number_of_sites, dtype=np.float64)

        iso_n = np.zeros(number_of_sites, dtype=np.float64)
        zeta_n = np.zeros(number_of_sites, dtype=np.float64)
        eta_n = np.zeros(number_of_sites, dtype=np.float64)
        ori_n = np.zeros(3*number_of_sites, dtype=np.float64)

        # Quad
        Cq_e = np.zeros(number_of_sites, dtype=np.float64)
        eta_e = np.zeros(number_of_sites, dtype=np.float64)
        ori_e = np.zeros(3*number_of_sites, dtype=np.float64)

        # Extract and assign site information from Site objects to C structure
        # ---------------------------------------------------------------------
        for i in range(number_of_sites):
            site = spin_sys.sites[i]
            spin_i[i] = site.isotope.spin
            gyromagnetic_ratio_i[i] = site.isotope.gyromagnetic_ratio
            i3 = 3*i

            # CSA tensor
            if site.isotropic_chemical_shift is not None:
                iso_n[i] = site.isotropic_chemical_shift

            shielding = site.shielding_symmetric
            if shielding is not None:
                if shielding.zeta is not None:
                    zeta_n[i] = shielding.zeta
                if shielding.eta is not None:
                    eta_n[i] = shielding.eta
                if shielding.alpha is not None:
                    ori_n[i3] = shielding.alpha
                if shielding.beta is not None:
                    ori_n[i3+1] = shielding.beta
                if shielding.gamma is not None:
                    ori_n[i3+2] = shielding.gamma

            # if verbose in [1, 11]:
            #     text = ((
            #         f"\n{isotope} site {i} from spin system {index} "
            #         f"@ {abundance}% abundance"
            #     ))
            #     len_ = len(text)
            #     print(text)
            #     print(f"{'-'*(len_-1)}")
            #     print(f'Isotropic chemical shift (δ) = {str(1e6*iso/larmor_frequency)} ppm')
            #     print(f'Shielding anisotropy (ζ) = {str(1e6*zeta/larmor_frequency)} ppm')
            #     print(f'Shielding asymmetry (η) = {eta}')
            #     print(f'Shielding orientation = [alpha = {alpha}, beta = {beta}, gamma = {gamma}]')

            # quad tensor
            if spin_quantum_number > 0.5:
                quad = site.quadrupolar
                if quad is not None:
                    if quad.Cq is not None:
                        Cq_e[i] = quad.Cq
                    if quad.eta is not None:
                        eta_e[i] = quad.eta
                    if quad.alpha is not None:
                        ori_e[i3] = quad.alpha
                    if quad.beta is not None:
                        ori_e[i3+1] = quad.beta
                    if quad.gamma is not None:
                        ori_e[i3+2] = quad.gamma

                # if verbose in [1, 11]:
                #     print(f'Quadrupolar coupling constant (Cq) = {Cq_e[i]/1e6} MHz')
                #     print(f'Quadrupolar asymmetry (η) = {eta}')
                #     print(f'Quadrupolar orientation = [alpha = {alpha}, beta = {beta}, gamma = {gamma}]')

        # sites packed as c struct
        sites_c.number_of_sites = number_of_sites
        sites_c.spin = &spin_i[0]
        sites_c.gyromagnetic_ratio = &gyromagnetic_ratio_i[0]

        sites_c.isotropic_chemical_shift_in_ppm = &iso_n[0]
        sites_c.shielding_symmetric_zeta_in_ppm = &zeta_n[0]
        sites_c.shielding_symmetric_eta = &eta_n[0]
        sites_c.shielding_orientation = &ori_n[0]

        sites_c.quadrupolar_Cq_in_Hz = &Cq_e[0]
        sites_c.quadrupolar_eta = &eta_e[0]
        sites_c.quadrupolar_orientation = &ori_e[0]
        # ------------------------------------------------------------------------
        #                           Coupling specification
        # ------------------------------------------------------------------------
        # J-coupling
        couplings_c.number_of_couplings = 0
        if spin_sys.couplings is not None:
            number_of_couplings = len(spin_sys.couplings)
            spin_index_ij = np.zeros(2*number_of_couplings, dtype=np.int32)

            iso_j = np.zeros(number_of_couplings, dtype=np.float64)
            zeta_j = np.zeros(number_of_couplings, dtype=np.float64)
            eta_j = np.zeros(number_of_couplings, dtype=np.float64)
            ori_j = np.zeros(3*number_of_couplings, dtype=np.float64)

            # Dipolar
            D_d = np.zeros(number_of_couplings, dtype=np.float64)
            eta_d = np.zeros(number_of_couplings, dtype=np.float64)
            ori_d = np.zeros(3*number_of_couplings, dtype=np.float64)

            # Extract and assign coupling information from Site objects to C structure
            for i in range(number_of_couplings):
                coupling = spin_sys.couplings[i]
                spin_index_ij[2*i: 2*i+2] = coupling.site_index
                i3 = 3*i

                # J tensor
                if coupling.isotropic_j is not None:
                    iso_j[i] = coupling.isotropic_j

                J_sym = coupling.j_symmetric
                if J_sym is not None:
                    if J_sym.zeta is not None:
                        zeta_j[i] = J_sym.zeta
                    if J_sym.eta is not None:
                        eta_j[i] = J_sym.eta
                    if J_sym.alpha is not None:
                        ori_j[i3] = J_sym.alpha
                    if J_sym.beta is not None:
                        ori_j[i3+1] = J_sym.beta
                    if J_sym.gamma is not None:
                        ori_j[i3+2] = J_sym.gamma

                # dipolar tensor
                dipolar = coupling.dipolar
                if dipolar is not None:
                    if dipolar.D is not None:
                        D_d[i] = dipolar.D
                    if dipolar.eta is not None:
                        eta_d[i] = dipolar.eta
                    if dipolar.alpha is not None:
                        ori_d[i3] = dipolar.alpha
                    if dipolar.beta is not None:
                        ori_d[i3+1] = dipolar.beta
                    if dipolar.gamma is not None:
                        ori_d[i3+2] = dipolar.gamma

            # if verbose in [1, 11]:
            #     print(f'N couplings = {number_of_couplings}')
            #     print(f'site index J = {spin_index_ij}')
            #     print(f'Isotropic J = {iso_j} Hz')
            #     print(f'J anisotropy = {zeta_j} Hz')
            #     print(f'J asymmetry = {eta_j}')
            #     print(f'J orientation = {ori_j}')

            #     print(f'Dipolar coupling constant = {D_d} Hz')
            #     print(f'Dipolar asymmetry = {eta_d}')
            #     print(f'Dipolar orientation = {ori_d}')

            # couplings packed as c struct
            couplings_c.number_of_couplings = number_of_couplings
            couplings_c.site_index = &spin_index_ij[0]

            couplings_c.isotropic_j_in_Hz = &iso_j[0]
            couplings_c.j_symmetric_zeta_in_Hz = &zeta_j[0]
            couplings_c.j_symmetric_eta = &eta_j[0]
            couplings_c.j_orientation = &ori_j[0]

            couplings_c.dipolar_coupling_in_Hz = &D_d[0]
            couplings_c.dipolar_eta = &eta_d[0]
            couplings_c.dipolar_orientation = &ori_d[0]


        # Spectrum amplitude vector -------------------------------------------
        amp = np.zeros(total_n_points)

        # if number_of_sites == 0:
        #     if decompose_spectrum == 1:
        #         amp_individual.append([])
        #     continue

        if number_of_sites != p_number_of_sites and isotopes != p_isotopes:
            transition_pathway = spin_sys.transition_pathways
            if transition_pathway is None:
                transition_pathway = np.asarray(method._get_transition_pathways_np(spin_sys))
                transition_pathway_c = np.asarray(transition_pathway, dtype=np.float32).ravel()
            else:
                transition_pathway = np.asarray(transition_pathway)
                # convert transition objects to list
                lst = [item.tolist() for item in transition_pathway.ravel()]
                transition_pathway_c = np.asarray(lst, dtype=np.float32).ravel()

            pathway_count, transition_count_per_pathway = transition_pathway.shape[:2]
            pathway_increment = 2*number_of_sites*transition_count_per_pathway

            p_number_of_sites = number_of_sites
            p_isotopes = isotopes

        # if spin_sys.transitions is not None:
        #     transition_pathway_c = np.asarray(
        #         spin_sys.transitions, dtype=np.float32
        #     ).ravel()
        # else:
        #     transition_pathway_c = np.asarray([0.5, -0.5], dtype=np.float32)

        # the number 2 is because of single site transition [mi, mf]
        # it dose not work for coupled sites.
        # transition_increment = 2*number_of_sites
        # number_of_transitions = int((transition_pathway_c.size)/transition_increment)

        for trans__ in range(pathway_count):
            clib.__mrsimulator_core(
                # spectrum information and related amplitude
                &amp[0],
                &sites_c,
                &couplings_c,
                &transition_pathway_c[pathway_increment*trans__],
                n_dimension,          # The total number of spectroscopic dimensions.
                dimensions,           # Pointer to MRS_dimension structure
                the_fftw_scheme,      # Pointer to the fftw scheme.
                the_averaging_scheme, # Pointer to the powder averaging scheme.
                interpolation,
                &freq_contrib_c[0],
                &affine_matrix_c[0],
                )

        temp = amp*abundance/norm

        if decompose_spectrum == 1:
            amp_individual.append(temp.reshape(method.shape()))
        else:
            amp1 += temp
        # else:
        #     if decompose_spectrum == 1:
        #         amp_individual.append([])

    # reverse the spectrum if gyromagnetic ratio is positive.
    if decompose_spectrum == 1 and len(amp_individual) != 0:
        if gyromagnetic_ratio < 0:
            amp1 = [np.fft.fftn(np.fft.ifftn(item).conj()).real for item in amp_individual]
        else:
            amp1 = amp_individual
    else:
        amp1.shape = method.shape()
        if gyromagnetic_ratio < 0:
            amp1 = np.fft.fftn(np.fft.ifftn(amp1).conj()).real

    clib.MRS_free_dimension(dimensions, n_dimension)
    clib.MRS_free_averaging_scheme(the_averaging_scheme)
    clib.MRS_free_fftw_scheme(the_fftw_scheme)
    return amp1


@cython.profile(False)
@cython.boundscheck(False)
@cython.wraparound(False)
def get_zeeman_states(sys):
    cdef int i, j, n_site = len(sys.sites)

    two_Ip1 = [int(2 * site.isotope.spin + 1) for site in sys.sites]
    spin_quantum_numbers = [
        np.arange(two_Ip1[i]) - site.isotope.spin for i, site in enumerate(sys.sites)
    ]

    lst = []
    for j in range(n_site):
        k = 1
        for i in range(n_site):
            if i == j:
                k = np.kron(k, spin_quantum_numbers[i])
            else:
                k = np.kron(k, np.ones(two_Ip1[i]))
        lst.append(k)
    return np.asarray(lst).T
