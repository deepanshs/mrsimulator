from libcpp cimport bool as bool_t

cdef extern from "angular_momentum.h":
    void __wigner_d_matrix(int l, int n, double *angle, double *wigner)

    void __wigner_d_matrix_cosine(int l, int n, double *cos_angle,
                                  double *wigner)

    void __wigner_rotation(int l, int n, double *wigner, double *cos_alpha,
                           double complex *R_in, double complex *R_out)

    void __wigner_dm0_vector(int l, double beta, double *R_out)


cdef extern from "powder_setup.h":
    void __powder_averaging_setup(
        int nt,
        double *cosAlpha,
        double *cosBeta,
        double *amp,
        int space)   # 1 for octant, 2 for hemisphere and 4 for sphere

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

        int quad_second_order,                # Quad theory for second order,
        int remove_second_order_quad_isotropic,   # remove the isotropic contribution from the
                                            # second order quad Hamiltonian.

        # spin rate, spin angle and number spinning sidebands
        int number_of_sidebands,
        double sample_rotation_frequency_in_Hz,
        double rotor_angle_in_rad,

        # The transition as transition[0] = mi and transition[1] = mf
        double *transition,
        int geodesic_polyhedron_frequency)
