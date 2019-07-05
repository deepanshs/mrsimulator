cdef extern from "angular_momentum.h":
    void __wigner_d_matrix_cosine(int l, int n, double *cos_angle,
                                  double *wigner)

cdef extern from "powder_setup.h":
    void __powder_averaging_setup(int nt, double *cosAlpha, double *cosBeta,
                                  double *amp, int space)

cdef extern from "isotopomer_ravel.h":
    ctypedef struct isotopomer_ravel:
        int number_of_sites;                    # Number of sites
        float spin;                             # The spin quantum number
        double larmor_frequency;                # Larmor frequency (MHz)
        double *isotropic_chemical_shift_in_Hz; # Isotropic chemical shift (Hz)
        double *shielding_anisotropy_in_Hz;     # Nuclear shielding anisotropy (Hz)
        double *shielding_asymmetry;            # Nuclear shielding asymmetry parameter
        double *shielding_orientation;          # Nuclear shielding PAS to CRS euler angles (rad.)
        double *quadrupolar_constant_in_Hz;     # Quadrupolar coupling constant (Hz)
        double *quadrupolar_asymmetry;          # Quadrupolar asymmetry parameter
        double *quadrupolar_orientation;        # Quadrupolar PAS to CRS euler angles (rad.)
        double *dipolar_couplings;              # dipolar coupling stored as list of lists

    ctypedef struct isotopomers_list:
        isotopomer_ravel *isotopomers

cdef extern from "spinning_sidebands.h":
    void __get_pre_phase_components(
        int number_of_sidebands,
        double spin_frequency,
        double complex *pre_phase)

    void __powder_averaging_setup(
        int nt,
        double *cos_alpha,
        double *cos_beta,
        double *amp,
        int space)   # 1 for octant, 2 for hemisphere and 4 for sphere

    void spinning_sideband_core(
        # spectrum information and related amplitude
        double * spec,
        double * cpu_time_,
        double spectral_start,
        double spectral_increment,
        int number_of_points,

        isotopomer_ravel *ravel_isotopomer,

        int quadSecondOrder,                # Quad theory for second order,
        int remove_second_order_quad_iso,   # remove the isotropic contribution from the
                                            # second order quad Hamiltonian.

        # spin rate, spin angle and number spinning sidebands
        int number_of_sidebands,
        double sample_rotation_frequency,
        double rotor_angle,

        # The transition as transition[0] = mi and transition[1] = mf
        double *transition,
        int geodesic_polyhedron_frequency)
        # # Euler angle -> principal to molecular frame
        # # double *omega_PM,

        # # Euler angles for powder averaging scheme
        # # powder orientation average
        # unsigned int n_orientations,
        # double *cos_alpha,
        # double *amp,
        # int nt,

        # unsigned int number_of_site,
        # double *wigner_2j_matrices,
        # double *wigner_4j_matrices,
        # double complex *pre_phase,
        # double *rotor_lab_2,
        # double *rotor_lab_4)
