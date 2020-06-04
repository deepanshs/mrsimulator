# -*- coding: utf-8 -*-
#
#  nmr_method.pxd
#
#  @copyright Deepansh J. Srivastava, 2019-2020.
#  Created by Deepansh J. Srivastava.
#  Contact email = deepansh2012@gmail.com
#
from libcpp cimport bool as bool_t

cdef extern from "angular_momentum.h":
    void wigner_d_matrices_from_exp_I_beta(int l, int n, void *exp_I_beta,
                                  double *wigner)

cdef extern from "schemes.h":
    ctypedef struct MRS_averaging_scheme:
        unsigned int total_orientations

    MRS_averaging_scheme *MRS_create_averaging_scheme(
                            unsigned int integration_density,
                            bool_t allow_fourth_rank,
                            unsigned int integration_volume)

    void MRS_free_averaging_scheme(MRS_averaging_scheme *scheme)

    ctypedef struct MRS_fftw_scheme:
        pass

    MRS_fftw_scheme *create_fftw_scheme(unsigned int total_orientations,
                                    int number_of_sidebands)

cdef extern from "mrsimulator.h":

    ctypedef struct MRS_plan:
        MRS_averaging_scheme *averaging_scheme
        int number_of_sidebands
        double sample_rotation_frequency_in_Hz
        double rotor_angle_in_rad
        # double complex *vector

    MRS_plan *MRS_create_plan(MRS_averaging_scheme *scheme, int number_of_sidebands,
                          double sample_rotation_frequency_in_Hz,
                          double rotor_angle_in_rad, double increment,
                          bool_t allow_fourth_rank)
    void MRS_free_plan(MRS_plan *plan)
    void MRS_get_amplitudes_from_plan(MRS_plan *plan, bool_t refresh)
    void MRS_get_frequencies_from_plan(MRS_plan *plan, double R0, double complex *R2,
                                  double complex *R4, bool_t refresh)

cdef extern from "isotopomer_ravel.h":
    ctypedef struct isotopomer_ravel:
        int number_of_sites                    # Number of sites
        float *spin                            # The spin quantum number
        double *gyromagnetic_ratio             # gyromagnetic ratio in (MHz/T)
        double *isotropic_chemical_shift_in_ppm # Isotropic chemical shift (Hz)
        double *shielding_symmetric_zeta_in_ppm     # Nuclear shielding anisotropy (Hz)
        double *shielding_symmetric_eta            # Nuclear shielding asymmetry
        double *shielding_orientation          # Nuclear shielding PAS to CRS euler angles (rad.)
        double *quadrupolar_Cq_in_Hz     # Quadrupolar coupling constant (Hz)
        double *quadrupolar_eta          # Quadrupolar asymmetry parameter
        double *quadrupolar_orientation        # Quadrupolar PAS to CRS euler angles (rad.)
        double *dipolar_couplings             # dipolar coupling stored as list of lists

cdef extern from "method.h":
    ctypedef struct MRS_event:
        double fraction                    # The weighted frequency contribution from the event.
        double magnetic_flux_density_in_T  #  he magnetic flux density in T.
        double rotor_angle_in_rad          # The rotor angle in radians.
        double sample_rotation_frequency_in_Hz # The sample rotation frequency in Hz.

    ctypedef struct MRS_sequence:
        int count                       #  The number of coordinates along the dimension.
        double increment                # Increment of coordinates along the dimension.
        double coordinates_offset       #  Start coordinate of the dimension.
        MRS_event *events               # Holds a list of events.
        unsigned int n_events           # The number of events.

    MRS_sequence *MRS_sequence_malloc(int n)

    MRS_sequence *MRS_create_plans_for_sequence(
        MRS_averaging_scheme *scheme,
        int *count,
        double *coordinates_offset,
        double *increment,
        double *magnetic_flux_density_in_T,
        double *sample_rotation_frequency_in_Hz,
        double *rotor_angle_in_rad,
        int *n_events,
        unsigned int n_seq,
        int number_of_sidebands)

cdef extern from "simulation.h":
    void mrsimulator_core(
        # spectrum information and related amplitude
        double * spec,
        double spectral_start,
        double spectral_increment,
        int number_of_points,

        isotopomer_ravel *ravel_isotopomer,
        MRS_sequence *the_sequence[],
        int n_sequence,

        int quad_second_order,                    # Quad theory for second order,
        bool_t remove_2nd_order_quad_isotropic,   # remove the isotropic contribution from the
                                                  # second order quad Hamiltonian.

        # spin rate, spin angle and number spinning sidebands
        int number_of_sidebands,
        double sample_rotation_frequency_in_Hz,
        double rotor_angle_in_rad,

        # The transition as transition[0] = mi and transition[1] = mf
        float *transition,
        int integration_density,
        unsigned int integration_volume,  # 0-octant, 1-hemisphere, 2-sphere
        bool_t interpolation
        )

    void __mrsimulator_core(
        # spectrum information and related amplitude
        double * spec,

        isotopomer_ravel *ravel_isotopomer,
        bool_t remove_2nd_order_quad_isotropic,   # remove the isotropic contribution from the
                                                  # second order quad Hamiltonian.

        # The transition as transition[0] = mi and transition[1] = mf
        float *transition,
        MRS_sequence *the_sequence, # the sequences within method.
        int n_sequence, # the number of sequences.
        MRS_fftw_scheme *fftw_scheme, # the fftw scheme
        MRS_averaging_scheme *scheme, # the powder averaging scheme
        bool_t interpolation
        )
