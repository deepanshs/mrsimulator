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

cdef extern from "averaging_scheme.h":
    ctypedef struct MRS_averaging_scheme:
        unsigned int total_orientations

    MRS_averaging_scheme * MRS_create_averaging_scheme(
                            unsigned int integration_density,
                            bool_t allow_fourth_rank,
                            unsigned int integration_volume)

    void MRS_free_averaging_scheme(MRS_averaging_scheme *scheme)

cdef extern from "mrsimulator.h":

    ctypedef struct MRS_plan:
        MRS_averaging_scheme *averaging_scheme
        int number_of_sidebands
        double sample_rotation_frequency_in_Hz
        double rotor_angle_in_rad
        double complex *vector

    MRS_plan *MRS_create_plan(MRS_averaging_scheme *scheme, int number_of_sidebands,
                          double sample_rotation_frequency_in_Hz,
                          double rotor_angle_in_rad, double increment,
                          bool_t allow_fourth_rank)
    void MRS_free_plan(MRS_plan *plan)
    void MRS_get_amplitudes_from_plan(MRS_plan *plan, double complex *R2,
                                  double complex *R4)
    void MRS_get_frequencies_from_plan(MRS_plan *plan, double R0)

    ctypedef struct MRS_dimension:
        pass

    MRS_dimension *MRS_create_dimension(int count, double coordinates_offset,
                                    double increment)

# cdef extern from "site.h":
#     ctypedef struct site:
#         unsigned int index
#         double isotropic_chemical_shift_in_Hz # Isotropic chemical shift (Hz)
#         double shielding_symmetric_anisotropy_in_Hz     # Nuclear shielding anisotropy (Hz)
#         double shielding_symmetric_asymmetry            # Nuclear shielding asymmetry
#         double shielding_symmetric_orientation[3]          # Nuclear shielding PAS to CRS euler angles (rad.)
#         double quadrupole_coupling_constant_in_Hz     # Quadrupolar coupling constant (Hz)
#         double quadrupole_asymmetry          # Quadrupolar asymmetry parameter
#         double quadrupole_orientation[3]        # Quadrupolar PAS to CRS euler angles (rad.)
#         double dipolar_coupling              # dipolar coupling sof the site


cdef extern from "isotopomer_ravel.h":
    ctypedef struct isotopomer_ravel:
        int number_of_sites                    # Number of sites
        float spin                             # The spin quantum number
        double larmor_frequency                # Larmor frequency (MHz)
        double *isotropic_chemical_shift_in_Hz # Isotropic chemical shift (Hz)
        double *shielding_anisotropy_in_Hz     # Nuclear shielding anisotropy (Hz)
        double *shielding_asymmetry            # Nuclear shielding asymmetry
        double *shielding_orientation          # Nuclear shielding PAS to CRS euler angles (rad.)
        double *quadrupole_coupling_constant_in_Hz     # Quadrupolar coupling constant (Hz)
        double *quadrupole_asymmetry          # Quadrupolar asymmetry parameter
        double *quadrupole_orientation        # Quadrupolar PAS to CRS euler angles (rad.)
        double *dipolar_couplings              # dipolar coupling stored as list of lists

cdef extern from "simulation.h":
    void mrsimulator_core(
        # spectrum information and related amplitude
        double * spec,
        double spectral_start,
        double spectral_increment,
        int number_of_points,

        isotopomer_ravel *ravel_isotopomer,

        int quad_second_order,                    # Quad theory for second order,
        int remove_second_order_quad_isotropic,   # remove the isotropic contribution from the
                                                  # second order quad Hamiltonian.

        # spin rate, spin angle and number spinning sidebands
        int number_of_sidebands,
        double sample_rotation_frequency_in_Hz,
        double rotor_angle_in_rad,

        # The transition as transition[0] = mi and transition[1] = mf
        double *transition,
        int integration_density,
        unsigned int integration_volume,  # 0-octant, 1-hemisphere, 2-sphere
        bool_t interpolation
        )

    void __mrsimulator_core(
        # spectrum information and related amplitude
        double * spec,
        # double spectral_start,
        # double spectral_increment,
        # int number_of_points,

        isotopomer_ravel *ravel_isotopomer,

        # int quad_second_order,                    # Quad theory for second order,
        int remove_second_order_quad_isotropic,   # remove the isotropic contribution from the
                                                  # second order quad Hamiltonian.

        # spin rate, spin angle and number spinning sidebands
        # int number_of_sidebands,
        # double sample_rotation_frequency_in_Hz,
        # double rotor_angle_in_rad,

        # The transition as transition[0] = mi and transition[1] = mf
        double *transition,

        # int integration_density,
        # unsigned int integration_volume             # 0-octant, 1-hemisphere, 2-sphere
        MRS_plan *plan,
        MRS_dimension *dimension,
        bool_t interpolation
        )
