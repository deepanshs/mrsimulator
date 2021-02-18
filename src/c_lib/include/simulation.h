// -*- coding: utf-8 -*-
//
//  simulation.h
//
//  @copyright Deepansh J. Srivastava, 2019-2020.
//  Created by Deepansh J. Srivastava, Apr 11, 2019.
//  Contact email = srivastava.89@osu.edu
//

#include "method.h"
#include "mrsimulator.h"
#include "octahedron.h"
// headerdoc

extern void mrsimulator_core(
    // spectrum information and related amplitude
    double *spec,               // The amplitude of the spectrum.
    double spectral_start,      // The start of the frequency spectrum.
    double spectral_increment,  // The increment of the frequency spectrum.
    int number_of_points,       // Number of points on the frequency spectrum.
    site_struct *sites,         // Pointer to sites within a spin system.
    coupling_struct *couplings, // Pointer to couplings within spin system.
    MRS_dimension *dimensions,  // Pointer to the spectral dimension.
    int n_dimension,            // Number of dimensions.
    int quad_second_order,      // Quad theory for second order,

    // spin rate, spin angle and number spinning sidebands
    unsigned int number_of_sidebands,       // The number of sidebands
    double sample_rotation_frequency_in_Hz, // The rotor spin frequency
    double rotor_angle_in_rad,              // The rotor angle relative to lab-frame z-axis

    // Pointer to the transition.
    float *transition,

    // powder orientation average
    // The number of triangle along the edge of octahedron.
    int integration_density,
    unsigned int integration_volume, // 0-octant, 1-hemisphere, 2-sphere.
    bool interpolation, bool *freq_contrib, double *affine_matrix);

extern void __mrsimulator_core(
    // spectrum information and related amplitude
    double *spec,               // The pointer to the vector representing the spectrum.
    site_struct *sites,         // A pointer to a list of sites within a spin system.
    coupling_struct *couplings, // A pointer to a list of couplings within a spin system.

    // Pointer to the transitions. transition[0] = mi and transition[1] = mf
    float *transition, MRS_dimension *dimensions, int n_dimension, MRS_fftw_scheme *fftw_scheme,
    MRS_averaging_scheme *scheme, bool interpolation, bool *freq_contrib, double *affine_matrix);
