# -*- coding: utf-8 -*-
#
#  nmr_method.pxd
#
#  @copyright Deepansh J. Srivastava, 2019-2021.
#  Created by Deepansh J. Srivastava.
#  Contact email = srivastava.89@osu.edu
#
from libcpp cimport bool as bool_t
from libc.stdint cimport int32_t


cdef extern from "angular_momentum/wigner_element.h":
    void transition_connect_factor(const float l, const float m1_f,
                            const float m1_i, const float m2_f, const float m2_i,
                            const double alpha, const double beta, const double gamma,
                            double *factor)

cdef extern from "tables.h":
    void generate_tables()

cdef extern from "angular_momentum/wigner_matrix.h":
    void wigner_d_matrices_from_exp_I_beta(int l, int n, bool_t half,
                                void *exp_I_beta, double *wigner)


cdef extern from "schemes.h":
    ctypedef struct MRS_averaging_scheme:
        unsigned int total_orientations

    ctypedef struct MRS_fftw_scheme:
        pass

    MRS_averaging_scheme *MRS_create_averaging_scheme(
                            unsigned int integration_density,
                            bool_t allow_4th_rank,
                            unsigned int n_gamma,
                            unsigned int integration_volume,
                            bool_t interpolation,
                            bool_t is_complex)

    MRS_averaging_scheme *MRS_create_averaging_scheme_from_alpha_beta(
                            double *alpha,
                            double *beta,
                            double *weight,
                            unsigned int n_angles,
                            bool_t allow_4th_rank,
                            unsigned int n_gamma,
                            unsigned int position_size,
                            int32_t *positions,
                            bool_t interpolation,
                            bool_t is_complex)

    void MRS_free_averaging_scheme(MRS_averaging_scheme *scheme)
    MRS_fftw_scheme *create_fftw_scheme(unsigned int total_orientations,
                                    unsigned int number_of_sidebands)
    void MRS_free_fftw_scheme(MRS_fftw_scheme *fftw_scheme)


cdef extern from "mrsimulator.h":
    ctypedef struct MRS_plan:
        unsigned int number_of_sidebands
        double rotor_frequency_in_Hz
        double rotor_angle_in_rad

    MRS_plan *MRS_create_plan(MRS_averaging_scheme *scheme, unsigned int number_of_sidebands,
                          double rotor_frequency_in_Hz,
                          double rotor_angle_in_rad, double increment,
                          bool_t allow_4th_rank)
    void MRS_free_plan(MRS_plan *plan)
    void MRS_get_amplitudes_from_plan(MRS_plan *plan, bool_t refresh)
    void MRS_get_frequencies_from_plan(MRS_plan *plan, double R0, double complex *R2,
                                  double complex *R4, bool_t refresh)

cdef extern from "object_struct.h":
    ctypedef struct site_struct:
        int number_of_sites                     # Number of sites
        float *spin                             # The spin quantum number
        double *gyromagnetic_ratio              # gyromagnetic ratio in (MHz/T)
        double *one_minus_sigma_iso_ref;        # Ref freq / Larmor freq ratio
        double *isotropic_chemical_shift_in_ppm # Isotropic chemical shift (Hz)
        double *shielding_symmetric_zeta_in_ppm # Nuclear shielding anisotropy (Hz)
        double *shielding_symmetric_eta         # Nuclear shielding asymmetry
        double *shielding_orientation           # Nuclear shielding PAS to CRS euler angles (rad.)
        double *quadrupolar_Cq_in_Hz            # Quadrupolar coupling constant (Hz)
        double *quadrupolar_eta                 # Quadrupolar asymmetry parameter
        double *quadrupolar_orientation         # Quadrupolar PAS to CRS euler angles (rad.)

    ctypedef struct coupling_struct:
        int number_of_couplings            # Number of couplings
        int *site_index                    # The site indexes of the coupled sites.
        double *isotropic_j_in_Hz          # isotropic J-coupling (Hz).
        double *j_symmetric_zeta_in_Hz     # J-coupling anisotropy (Hz).
        double *j_symmetric_eta            # J-coupling asymmetry.
        double *j_orientation              # J tensor PAS to CRS euler angles (rad.)
        double *dipolar_coupling_in_Hz     # Dipolar coupling constant (Hz)
        double *dipolar_eta                # Dipolar asymmetry parameter
        double *dipolar_orientation        # Dipolar tensor PAS to CRS euler angles (rad.)


cdef extern from "method.h":
    ctypedef struct MRS_event:
        double fraction                    # The weighted frequency contribution from the event.
        double duration                    # The duration of a delay event in µs.
        unsigned char is_spectral                   # True if the event is a SpectralEvent, False if a DelayEvent
        double magnetic_flux_density_in_T  # The magnetic flux density in T.
        double rotor_angle_in_rad          # The rotor angle in radians.
        double rotor_frequency_in_Hz       # The sample rotation frequency in Hz.

    ctypedef struct MRS_dimension:
        int count                       # The number of coordinates along the dimension.
        double increment                # Increment of coordinates along the dimension.
        double coordinates_offset       # Start coordinate of the dimension.
        MRS_event *events               # Holds a list of events.
        unsigned int n_events           # The number of events.

    MRS_dimension *MRS_create_dimensions(
        MRS_averaging_scheme *scheme,
        int *count,
        double *coordinates_offset,
        double *increment,
        double *fraction,
        double *duration,
        unsigned char *is_spectral,
        double *magnetic_flux_density_in_T,
        double *rotor_frequency_in_Hz,
        double *rotor_angle_in_rad,
        int *n_events,
        unsigned int n_dim,
        unsigned int *number_of_sidebands)

    void MRS_free_dimension(MRS_dimension *dimensions, int n)


cdef extern from "simulation.h":
    void mrsimulator_core(
        # spectrum information and related amplitude
        double *spec,
        double spectral_start,
        double spectral_increment,
        int number_of_points,

        site_struct *sites,
        coupling_struct *couplings,

        MRS_dimension *dimensions[],
        int n_dimension,

        int quad_second_order,                    # Quad theory for second order,

        # spin rate, spin angle and number spinning sidebands
        unsigned int number_of_sidebands,
        double rotor_frequency_in_Hz,
        double rotor_angle_in_rad,

        float *transition_pathway, # Pointer to a list of transitions.
        int integration_density,
        unsigned int integration_volume,  # 0-octant, 1-hemisphere, 2-sphere
        unsigned int interpolate_type,
        unsigned char *freq_contrib,
        double *affine_matrix,
        )

    void __mrsimulator_core(
        # spectrum information and related amplitude
        double *spec,
        site_struct *sites,
        coupling_struct *couplings,
        float *transition_pathway,    # Pointer to a list of transitions.
        double *transition_pathway_weight,  # The complex weight of transition pathway.
        int n_dimension,              # the number of dimensions.
        MRS_dimension *dimensions,    # the dimensions within method.
        MRS_fftw_scheme *fftw_scheme, # the fftw scheme
        MRS_averaging_scheme *scheme, # the powder averaging scheme
        unsigned int interpolate_type,
        unsigned char *freq_contrib,
        double *affine_matrix,
        )


cdef extern from "vm_linalg.h":
    void vm_haeberlen_components(
        int n,
        double *expr_base_p,
        double *expr_base_q,
        double zeta,
        double eta,
        double rho,
        double *param)
