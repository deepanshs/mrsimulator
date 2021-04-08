# -*- coding: utf-8 -*-
#
#  test.pxd
#
#  @copyright Deepansh J. Srivastava, 2019-2021.
#  Created by Deepansh J. Srivastava.
#  Contact email = srivastava.89@osu.edu
#

from libcpp cimport bool as bool_t

cdef extern from "tables/trig.h":
    void generate_table()

cdef extern from "angular_momentum.h":
    void wigner_d_matrices(const int l, const int n, const double *angle, double *wigner)

    void wigner_d_matrices_from_exp_I_beta(const int l, const int n, const bool_t half,
                            const double complex *exp_I_beta, double *wigner)

    # void __wigner_rotation(const int l, const int n, const double *wigner, const double *cos_alpha,
    #                        const double complex *R_in, double complex *R_out)

    void __wigner_rotation_2(const int l, const int n, const double *wigner,
                            const void *exp_Im_alpha, const void *R_in, void *R_out)

    void single_wigner_rotation(const int l, const double *euler_angles, const void *R_in,
                            void *R_out)

    void wigner_dm0_vector(const int l, const double beta, double *R_out)

    void get_exp_Im_alpha(const unsigned int octant_orientations,
                            const bool_t allow_fourth_rank, void *exp_Im_alpha)

    void __batch_wigner_rotation(const unsigned int octant_orientations,
                            const unsigned int n_octants, double *wigner_2j_matrices, void *R2,
                            double *wigner_4j_matrices, void *R4, void *exp_Im_alpha, void *w2,
                            void *w4)


cdef extern from "octahedron.h":
    void averaging_setup(
        int nt,
        double complex *exp_I_alpha,
        double complex *exp_I_beta,
        double *amp)

cdef extern from "interpolation.h":
    void triangle_interpolation1D(
        double *freq1,
        double *freq2,
        double *freq3,
        double *amp,
        double *spec,
        int *points)

    void triangle_interpolation2D(
        double *freq11,
        double *freq12,
        double *freq13,
        double *freq21,
        double *freq22,
        double *freq23,
        double *amp,
        double *spec,
        int m0,
        int m1)

    void octahedronInterpolation(
        double *spec,
        double *freq,
        int nt,
        double *amp,
        int stride,
        int m)

cdef extern from "mrsimulator.h":
    void get_sideband_phase_components(
        unsigned int number_of_sidebands,
        double spin_frequency,
        double *pre_phase)

#     ctypedef struct MRS_plan

#     MRS_plan *MRS_create_plan(
#         unsigned int integration_density,
#         int number_of_sidebands,
#         double rotor_frequency_in_Hz,
#         double rotor_angle_in_rad, double increment,
#         bool_t allow_fourth_rank)


cdef extern from "object_struct.h":
    ctypedef struct site_struct:
        int number_of_sites;                     # Number of sites
        float *spin;                             # The spin quantum number
        double *gyromagnetic_ratio;              # Larmor frequency (MHz)
        double *isotropic_chemical_shift_in_ppm; # Isotropic chemical shift (Hz)
        double *shielding_symmetric_zeta_in_ppm; # Nuclear shielding anisotropy (Hz)
        double *shielding_symmetric_eta;         # Nuclear shielding asymmetry parameter
        double *shielding_orientation;           # Nuclear shielding PAS to CRS euler angles (rad.)
        double *quadrupolar_Cq_in_Hz;            # Quadrupolar coupling constant (Hz)
        double *quadrupolar_eta;                 # Quadrupolar asymmetry parameter
        double *quadrupolar_orientation;         # Quadrupolar PAS to CRS euler angles (rad.)

cdef extern from "method.h":
    ctypedef struct MRS_event:
        double fraction                    # The weighted frequency contribution from the event.
        double magnetic_flux_density_in_T  #  he magnetic flux density in T.
        double rotor_angle_in_rad          # The rotor angle in radians.
        double rotor_frequency_in_Hz       # The sample rotation frequency in Hz.

    ctypedef struct MRS_dimension:
        int count                       #  The number of coordinates along the dimension.
        double increment                # Increment of coordinates along the dimension.
        double coordinates_offset       #  Start coordinate of the dimension.
        MRS_event *events               # Holds a list of events.
        unsigned int n_events           # The number of events.

    # MRS_dimension *MRS_create_dimensions(
    #     MRS_averaging_scheme *scheme,
    #     int count,
    #     double coordinates_offset,
    #     double increment,
    #     double *magnetic_flux_density_in_T,
    #     double *rotor_frequency_in_Hz,
    #     double *rotor_angle_in_rad,
    #     unsigned int n_events,
    #     int number_of_sidebands)

cdef extern from "simulation.h":
    void mrsimulator_core(
        # spectrum information and related amplitude
        double * spec,
        double spectral_start,
        double spectral_increment,
        int number_of_points,

        site_struct *sites,
        MRS_dimension *dimensions[],              # the dimensions in the method.

        int quad_second_order,                    # Quad theory for second order,
        bool_t remove_2nd_order_quad_isotropic,   # remove the isotropic contribution from the
                                                  # second order quad Hamiltonian.

        # spin rate, spin angle and number spinning sidebands
        unsigned int number_of_sidebands,
        double rotor_frequency_in_Hz,
        double rotor_angle_in_rad,

        # The transition as transition[0] = mi and transition[1] = mf
        double *transition,
        int integration_density,
        unsigned int integration_volume,      # 0-octant, 1-hemisphere, 2-sphere.
        bool_t interpolation
        )
