# -*- coding: utf-8 -*-
#
#  sandbox.pxd
#
#  @copyright Deepansh J. Srivastava, 2019-2021.
#  Created by Deepansh J. Srivastava
#  Contact email = srivastava.89@osu.edu
#

from libcpp cimport bool as bool_t

cdef extern from "schemes.h":
    ctypedef struct MRS_averaging_scheme:
        unsigned int total_orientations
        unsigned int integration_density
        unsigned int integration_volume

    MRS_averaging_scheme * MRS_create_averaging_scheme(
                            unsigned int integration_density,
                            bool_t allow_fourth_rank,
                            unsigned int integration_volume)

    MRS_averaging_scheme *MRS_create_averaging_scheme_from_alpha_beta(
                            double *alpha, double *beta,
                            double *weight, unsigned int n_angles,
                            bool_t allow_fourth_rank)

    void MRS_free_averaging_scheme(MRS_averaging_scheme *scheme)

cdef extern from "mrsimulator.h":

    ctypedef struct MRS_plan:
        MRS_averaging_scheme *averaging_scheme
        unsigned int number_of_sidebands
        double rotor_frequency_in_Hz
        double rotor_angle_in_rad
        # double complex *vector

    MRS_plan *MRS_create_plan(MRS_averaging_scheme *scheme, unsigned int number_of_sidebands,
                          double rotor_frequency_in_Hz,
                          double rotor_angle_in_rad, double increment,
                          bool_t allow_fourth_rank)
    void MRS_free_plan(MRS_plan *plan)
    void MRS_get_amplitudes_from_plan(MRS_plan *plan, double complex *R2,
                                  double complex *R4)
    void MRS_get_frequencies_from_plan(MRS_plan *plan, double R0)
