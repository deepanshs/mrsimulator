cimport base_model as clib
cimport numpy as cnp
from libcpp cimport bool as bool_t
from libc.stdint cimport int32_t
from numpy cimport ndarray
import numpy as np
import cython

__author__ = "Deepansh J. Srivastava"
__email__ = "srivastava.89@osu.edu"

clib.generate_tables()

@cython.profile(False)
@cython.boundscheck(False)
@cython.wraparound(False)
def core_simulator(method,
       list spin_systems,
       list transition_pathways,          # Same length as spin_systems list and
       list transition_weights,  # elements coorespond to each system
       int verbose=0,  # for debug purpose only.
       unsigned int number_of_sidebands=90,
       unsigned int integration_density=72,
       unsigned int decompose_spectrum=0,
       unsigned int integration_volume=1,
       unsigned int isotropic_interpolation=0,
       unsigned int number_of_gamma_angles=1,
       bool_t is_complex=True,
       bool_t interpolation=True,
       bool_t auto_switch=True,
       debug=False,
       bool_t user_defined=False,
       ndarray[double, ndim=1] alpha=np.zeros(1, dtype=float),
       ndarray[double, ndim=1] beta=np.zeros(1, dtype=float),
       ndarray[double, ndim=1] weight=np.zeros(1, dtype=float),
       ndarray[int32_t, ndim=1] positions=np.ones(1, dtype=np.intc)):
    """core simulator init"""

# initialization and config
    # observed spin is always channel at index 0_______________________________________
    channel = method.channels[0].symbol
    cdef double spin_quantum_number = method.channels[0].spin

    # gyromagnetic ratio and reverse axis factor
    cdef gyromagnetic_ratio = method.channels[0].gyromagnetic_ratio
    cdef double gyro_factor = 1.0
    if gyromagnetic_ratio > 0.0:
        gyro_factor = -1.0

    # config for spin I=0.5
    cdef bool_t allow_4th_rank = 1

    # transitions of the observed spin
    cdef int transition_increment, number_of_transitions, i
    cdef ndarray[float, ndim=1] transition_pathway_c
    cdef ndarray[double, ndim=1] transition_pathway_weight_c

# create averaging scheme _____________________________________________________
    cdef clib.MRS_averaging_scheme *averaging_scheme
    cdef unsigned int position_size = 0

    if user_defined:
        if interpolation:
            position_size = np.uint32(positions.size / 3) if positions is not None else 0
        else:
            positions = None
        averaging_scheme = clib.MRS_create_averaging_scheme_from_alpha_beta(
            alpha=&alpha[0], beta=&beta[0], weight=&weight[0], n_angles=alpha.size,
            allow_4th_rank=allow_4th_rank, n_gamma=number_of_gamma_angles,
            position_size=position_size, positions=&positions[0], interpolation=interpolation,
            is_complex=is_complex
        )
    else:
        averaging_scheme = clib.MRS_create_averaging_scheme(
            integration_density=integration_density, allow_4th_rank=allow_4th_rank,
            n_gamma=number_of_gamma_angles, integration_volume=integration_volume,
            interpolation=interpolation, is_complex=is_complex
        )

# create C spectral dimensions ________________________________________________
    cdef int n_dimension = len(method.spectral_dimensions)

    # Define and allocate C numpy arrays
    n_points = 1
    cdef int n_ev
    cdef ndarray[int] n_event
    cdef ndarray[double] magnetic_flux_density_in_T
    cdef ndarray[double] frac   # for SpectralEvent
    cdef ndarray[double] dur    # for DelayEvent
    cdef ndarray[unsigned char] is_spectral
    cdef ndarray[double] srfiH
    cdef ndarray[double] rair
    cdef ndarray[int] cnt
    cdef ndarray[double] coord_off
    cdef ndarray[double] incre
    cdef ndarray[unsigned int] n_dim_sidebands

    # Loop through dimensions and grab attributes/values in python
    freq_contrib = np.asarray([])

    fr = []  # fraction attribute
    du = []  # duration attribute
    is_spec = []  # if event is a SpectralEvent
    Bo = []
    vr = []
    th = []
    event_i = []
    count = []
    increment = []
    coordinates_offset = []
    dim_sidebands = []

    for i, dim in enumerate(method.spectral_dimensions):
        n_ev = 0
        track = []
        for event in dim.events:
            if event.__class__.__name__ not in ["MixingEvent", "RotationEvent"]:
                freq_contrib = np.append(freq_contrib, event._freq_contrib_flags())

                if event.rotor_frequency < 1.0e-3:
                    rotor_frequency_in_Hz = 1.0e-6
                else:
                    rotor_frequency_in_Hz = event.rotor_frequency

                rotor_angle_in_rad = event.rotor_angle
                track.append(event.rotor_frequency < 1e12 and event.rotor_frequency != 0)

                # Update event attribute depending on event type
                if event.__class__.__name__ == "SpectralEvent":
                    fr.append(event.fraction)
                    du.append(0)
                    is_spec.append(1)
                elif event.__class__.__name__ == "DelayEvent":
                    fr.append(0)
                    du.append(event.duration)
                    is_spec.append(0)

                Bo.append(event.magnetic_flux_density)  # in T
                vr.append(rotor_frequency_in_Hz) # in Hz
                th.append(rotor_angle_in_rad) # in rad
                n_ev +=1

        n_points *= dim.count

        count.append(dim.count)
        offset = dim.spectral_width / 2.0
        coordinates_offset.append(-dim.reference_offset * gyro_factor - offset)
        increment.append(dim.spectral_width / dim.count)
        event_i.append(n_ev)

        dim_sidebands.append(number_of_sidebands if np.any(track) else 1)

    # Assing values to previously defined C numpy arrays
    frac = np.asarray(fr, dtype=np.float64)
    dur = np.asarray(du, dtype=np.float64)
    is_spectral = np.asarray(is_spec, dtype=np.uint8)
    magnetic_flux_density_in_T = np.asarray(Bo, dtype=np.float64)
    srfiH = np.asarray(vr, dtype=np.float64)
    rair = np.asarray(th, dtype=np.float64)
    cnt = np.asarray(count, dtype=np.int32)
    incre = np.asarray(increment, dtype=np.float64)
    coord_off = np.asarray(coordinates_offset, dtype=np.float64)
    n_event = np.asarray(event_i, dtype=np.int32)
    n_dim_sidebands = np.asarray(dim_sidebands, dtype=np.uint32)

    # # special 1D case with 1 event.
    # if np.all(srfiH == 1e-3) and np.all(rair - rair[0] == 0):
    #     # rair[:] = 0.0
    #     n_dim_sidebands[0] = 1
    # if np.all(srfiH == 1e12):
    #     n_dim_sidebands[0] = 1

    if srfiH.size == 1 and srfiH[0] == 1e-6 and auto_switch:
        rair[0] = 0.0
        n_dim_sidebands[0] = 1

    # create spectral_dimensions
    dimensions = clib.MRS_create_dimensions(averaging_scheme, &cnt[0],
        &coord_off[0], &incre[0], &frac[0], &dur[0], &is_spectral[0],
        &magnetic_flux_density_in_T[0], &srfiH[0], &rair[0],
        &n_event[0], n_dimension, &n_dim_sidebands[0])

# normalization factor for the spectrum
    norm = np.abs(np.prod(incre))

# create fftw scheme __________________________________________________________
    cdef unsigned int max_sidebands = n_dim_sidebands.max()
    cdef clib.MRS_fftw_scheme *fftw_scheme
    fftw_scheme = clib.create_fftw_scheme(averaging_scheme.total_orientations, max_sidebands)
# _____________________________________________________________________________

# frequency contrib
    # numpy uint8 corresponds to an unsigned char C type
    # More types found here: https://numpy.org/doc/stable/user/basics.types.html#array-types-and-conversions-between-types
    cdef ndarray[cnp.uint8_t] f_contrib = np.asarray(freq_contrib, dtype=np.uint8)

# affine transformation
    cdef ndarray[double] affine_matrix_c
    if method.affine_matrix is None:
        affine_matrix_c = np.asarray([1, 0, 0, 1], dtype=np.float64)
    else:
        increment_fraction = [incre/item for item in incre]
        matrix = np.asarray(method.affine_matrix).ravel() * np.asarray(increment_fraction).ravel()
        affine_matrix_c = np.asarray(matrix, dtype=np.float64)
        # affine matrix = [[a, b], [c, d]] represented as [[a, b], [c/a, d - bc/a]]
        if affine_matrix_c[2] != 0:
            affine_matrix_c[2] /= affine_matrix_c[0]
            affine_matrix_c[3] -=  affine_matrix_c[1]*affine_matrix_c[2]

# sites _______________________________________________________________________________
    cdef int number_of_sites, number_of_couplings
    cdef ndarray[int] spin_index_ij
    cdef ndarray[float] spin_i
    cdef ndarray[double] gyromagnetic_ratio_i
    cdef ndarray[double]  one_minus_sigma_iso_ref_i

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
    cdef ndarray[double, ndim=1] amp = np.zeros(2 * n_points, dtype=np.float64)
    amp1 = np.zeros(n_points, dtype=np.complex128)
    amp_individual = []

    cdef clib.site_struct sites_c
    cdef clib.coupling_struct couplings_c

    # -------------------------------------------------------------------------
    # sample __________________________________________________________________
    for idx in range(len(spin_systems)):
        spin_sys = spin_systems[idx]
        segments = transition_pathways[idx]
        weights = transition_weights[idx]

        abundance = spin_sys.abundance
        isotopes = [site.isotope.symbol for site in spin_sys.sites]
        if channel not in isotopes:
            if decompose_spectrum == 1:
                amp_individual.append(np.zeros(method.shape()))
            continue

        # sub_sites = [site for site in spin_sys.sites if site.isotope.symbol == isotope]
        # index_.append(index)
        number_of_sites = len(spin_sys.sites)

        # ------------------------------------------------------------------------
        #                          Site specification
        # ------------------------------------------------------------------------
        # CSA
        spin_i = np.empty(number_of_sites, dtype=np.float32)
        gyromagnetic_ratio_i = np.empty(number_of_sites, dtype=np.float64)
        one_minus_sigma_iso_ref_i = np.empty(number_of_sites, dtype=np.float64)

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
            one_minus_sigma_iso_ref_i[i] = site.isotope.ref_larmor_ratio
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

            if debug:
                print(f'Isotropic chemical shift (δ) = {str(iso_n)} ppm')
                print(f'Shielding anisotropy (ζ) = {str(zeta_n)} ppm')
                print(f'Shielding asymmetry (η) = {eta_n}')
                print(f'Shielding orientation (alpha/beta/gamma) = {ori_n}')

            # quad tensor
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

            if debug:
                print(f'Quadrupolar coupling constant (Cq) = {Cq_e[i]/1e6} MHz')
                print(f'Quadrupolar asymmetry (η) = {eta_e}')
                print(f'Quadrupolar orientation (alpha/beta/gamma) = {ori_e}]')

        # sites packed as c struct
        sites_c.number_of_sites = number_of_sites
        sites_c.spin = &spin_i[0]
        sites_c.gyromagnetic_ratio = &gyromagnetic_ratio_i[0]
        sites_c.one_minus_sigma_iso_ref = &one_minus_sigma_iso_ref_i[0]

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

            # if debug:
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


        # if number_of_sites == 0:
        #     if decompose_spectrum == 1:
        #         amp_individual.append([])
        #     continue

        transition_pathway_c = np.asarray(segments, dtype=np.float32).ravel()  # NOTE: Why 32 then 64 bits?
        transition_pathway_weight_c = weights.view(dtype=np.float64)

        pathway_count, transition_count_per_pathway = segments.shape[:2]
        pathway_increment = 2*number_of_sites*transition_count_per_pathway

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

        # print('pathway', transition_pathway_c)
        # print('weight', transition_pathway_weight_c)
        # print('pathway_count, inc', pathway_count, pathway_increment)
        for trans__ in range(pathway_count):
            clib.__mrsimulator_core(
                &amp[0],  # as complex array
                &sites_c,
                &couplings_c,
                &transition_pathway_c[pathway_increment*trans__],
                &transition_pathway_weight_c[2*trans__],
                n_dimension,      # The total number of spectroscopic dimensions.
                dimensions,       # Pointer to MRS_dimension structure
                fftw_scheme,      # Pointer to the fftw scheme.
                averaging_scheme, # Pointer to the powder averaging scheme.
                isotropic_interpolation,
                &f_contrib[0],
                &affine_matrix_c[0],
            )

        temp = amp.view(dtype=np.complex128)
        temp *= abundance/norm

        if decompose_spectrum == 1:
            amp_individual.append(temp.copy().reshape(method.shape()))
        else:
            amp1 += temp
        amp[:] = 0

    # reverse the spectrum if gyromagnetic ratio is positive.
    if decompose_spectrum == 1 and len(amp_individual) != 0:
        if gyromagnetic_ratio < 0:
            amp1 = [np.fft.fftn(np.fft.ifftn(item).conj()) for item in amp_individual]
        else:
            amp1 = amp_individual
    else:
        amp1.shape = method.shape()
        if gyromagnetic_ratio < 0:
            amp1 = np.fft.fftn(np.fft.ifftn(amp1).conj())

    clib.MRS_free_dimension(dimensions, n_dimension)
    clib.MRS_free_averaging_scheme(averaging_scheme)
    clib.MRS_free_fftw_scheme(fftw_scheme)
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


# @cython.profile(False)
# @cython.boundscheck(False)
# @cython.wraparound(False)
# def transition_connect_factor(float l, float m1_f, float m1_i, float m2_f,
#                         float m2_i, double theta, double phi):
#     """Evaluate the probability of connecting two transitions driven by an external rf
#     pulse of phase phi and angle theta. The connected transitions are
#     | m1_f >< m1_i | --> | m2_f > < m2_i |.

#     Args:
#         float l: The angular momentum quantum number of the spin involved in the transition.
#         float m1_f Final quantum number of the starting transition.
#         float m1_i Initial quantum number of the starting transition.
#         float m2_f Final quantum number of the connecting transition.
#         float m2_i Initial quantum number of the connecting transition.
#         float theta The tip-angle of the rf pulse.
#         float phi The phase of the rf pulse.

#     Return: A complex amplitude.
#     """
#     cdef ndarray[double] factor = np.asarray([1, 0], dtype=np.float64)
#     clib.transition_connect_factor(l, m1_f, m1_i, m2_f, m2_i, theta, phi, &factor[0])
#     factor = np.around(factor, decimals=12)
#     return complex(factor[0], factor[1])


@cython.profile(False)
@cython.boundscheck(False)
@cython.wraparound(False)
def transition_connect_factor(float l, float m1_f, float m1_i, float m2_f,
                        float m2_i, double alpha, double beta, double gamma):
    """Evaluate the probability of connecting two transitions driven by a rotation
    defined by the euler angles alpha, beta, and gamma in the ZYZ convention.
    The connected transitions are | m1_f >< m1_i | --> | m2_f > < m2_i |.

    Args:
        float l: The angular momentum quantum number of the spin involved in the transition.
        float m1_f Final quantum number of the starting transition.
        float m1_i Initial quantum number of the starting transition.
        float m2_f Final quantum number of the connecting transition.
        float m2_i Initial quantum number of the connecting transition.
        float alpha The first angle of rotation about the Z axis
        float beta The second angle of rotation about the transformed Y axis
        float gamma the third angle of rotation about the transformed Z axis

    Return: A complex amplitude.
    """
    cdef ndarray[double] factor = np.asarray([1, 0], dtype=np.float64)
    clib.transition_connect_factor(l, m1_f, m1_i, m2_f, m2_i, alpha, beta, gamma, &factor[0])
    factor = np.around(factor, decimals=12)
    return complex(factor[0], factor[1])


@cython.profile(False)
@cython.boundscheck(False)
@cython.wraparound(False)
def calculate_transition_connect_weight(
        ndarray[float, ndim=2] trans1,
        ndarray[float, ndim=2] trans2,
        ndarray[float, ndim=1] spin,
        ndarray[double, ndim=1] alpha,
        ndarray[double, ndim=1] beta,
        ndarray[double, ndim=1] gamma
    ):
    """Evaluate the probability of connecting two transitions driven by a rotation described
    by the Euler angles alpha, beta, gamma in the ZYZ convention. The connected transitions are
    | m1_f >< m1_i | --> | m2_f > < m2_i |.

    Args:
        float l: The angular momentum quantum number of the spin involved in the transition.
        float m1_f Final quantum number of the starting transition.
        float m1_i Initial quantum number of the starting transition.
        float m2_f Final quantum number of the connecting transition.
        float m2_i Initial quantum number of the connecting transition.
        float alpha First euler angle.
        float beta Second euler angle.
        float gamma Third euler angle.

    Return: A complex amplitude.
    """
    cdef ndarray[double] factor = np.asarray([1, 0], dtype=np.float64)
    cdef int i, n_sites = spin.size
    cdef float m1_f, m1_i, m2_f, m2_i
    for i in range(n_sites):
        m1_f = trans1[1][i]  # starting transition final state
        m1_i = trans1[0][i]  # starting transition initial state
        m2_f = trans2[1][i]  # landing transition final state
        m2_i = trans2[0][i]  # landing transition initial state

        clib.transition_connect_factor(
            spin[i], m1_f, m1_i, m2_f, m2_i, alpha[i], beta[i], gamma[i], &factor[0]
        )
    return complex(factor[0], factor[1])


@cython.profile(False)
@cython.boundscheck(False)
@cython.wraparound(False)
def get_Haeberlen_components(
        ndarray[double, ndim=2] expr_base_p,
        ndarray[double, ndim=2] expr_base_q,
        double zeta,
        double eta,
        double rho):
    """Return random extended czjzek tensors in Haeberlen convention"""
    cdef int n = expr_base_p.shape[1]
    cdef ndarray[double, ndim=2] param = np.empty((n, 2), dtype=float)
    clib.vm_haeberlen_components(
        n, &expr_base_p[0, 0], &expr_base_q[0, 0], zeta, eta, rho, &param[0, 0]
    )
    return param[:, 0], param[:, 1]


# @cython.profile(False)
# @cython.boundscheck(False)
# @cython.wraparound(False)
# def pathway_rotation_factor(float l, float *pathway, float m2_a, float m1_b,
#                         float m2_b, double theta, double phi):
#     cdef ndarray[double] factor = np.zeros(2, dtype=np.float64)
#     clib.transition_connect_factor(l, m1_a, m2_a, m1_b, m2_b, theta, phi, &factor[0])
#     return complex(factor[0], factor[1])
