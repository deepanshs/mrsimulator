from libcpp cimport bool as bool_t

cdef extern from "angular_momentum.h":
    void wigner_d_matrix(const int l, const int n, const double *angle, double *wigner)

    void __wigner_d_matrix_cosine(const int l, const int n, const double *cos_angle,
                                  double *wigner)

    void __wigner_rotation(const int l, const int n, const double *wigner, const double *cos_alpha,
                           const double complex *R_in, double complex *R_out)

    void __wigner_rotation_2(const int l, const int n, const double *wigner, const double complex *exp_Im_alpha,
                           const double complex *R_in, double complex *R_out)

    void single_wigner_rotation(const int l, const double *euler_angles,
                            const double complex *R_in, double complex *R_out)

    void __wigner_dm0_vector(const int l, const double beta, double *R_out)

    void get_exp_Im_alpha(const unsigned int n_orientations,
                             const double *cos_alpha,
                             const bool_t allow_fourth_rank,
                             double complex *exp_Im_alpha)

    void __batch_wigner_rotation(const unsigned int octant_orientations,
                             const unsigned int n_octants,
                             const double *wigner_2j_matrices,
                             const double complex *R2,
                             const double *wigner_4j_matrices,
                             const double complex *R4,
                             double complex *exp_Im_alpha,
                             double complex *w2, double complex *w4)


cdef extern from "powder_setup.h":
    void __powder_averaging_setup(
        int nt,
        double *cosAlpha,
        double *cosBeta,
        double *amp)

cdef extern from "interpolation.h":
    void triangle_interpolation(
        double *freq1,
        double *freq2,
        double *freq3,
        double *amp,
        double *spec,
        int *points)

cdef extern from "octahedron.h":
    void octahedronInterpolation(
        double *spec,
        double *freq,
        int nt,
        double *amp,
        int stride,
        int m)

cdef extern from "mrsimulator.h":
    void __get_components(
        int number_of_sidebands,
        double spin_frequency,
        double complex *pre_phase)

    ctypedef struct MRS_plan

    MRS_plan *MRS_create_plan(
        unsigned int geodesic_polyhedron_frequency,
        int number_of_sidebands,
        double sample_rotation_frequency_in_Hz,
        double rotor_angle_in_rad, double increment,
        bool_t allow_fourth_rank)


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
    void spinning_sideband_core(
        # spectrum information and related amplitude
        double * spec,
        double spectral_start,
        double spectral_increment,
        int number_of_points,

        isotopomer_ravel *ravel_isotopomer,

        int quadSecondOrder,                # Quad theory for second order,
        int remove_second_order_quad_iso,   # remove the isotropic contribution from the
                                            # second order quad Hamiltonian.

        # spin rate, spin angle and number spinning sidebands
        int number_of_sidebands,
        double sample_rotation_frequency_in_Hz,
        double rotor_angle_in_rad,

        # The transition as transition[0] = mi and transition[1] = mf
        double *transition,
        int geodesic_polyhedron_frequency)
