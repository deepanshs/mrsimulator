// -*- coding: utf-8 -*-
//
//  simulation.h
//
//  @copyright Deepansh J. Srivastava, 2019-2020.
//  Created by Deepansh J. Srivastava, Apr 11, 2019.
//  Contact email = deepansh2012@gmail.com
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

    isotopomer_ravel *ravel_isotopomer,  // isotopomer structure
    MRS_sequence *the_sequence,          // the transition sequence.
    int n_sequence,                      // number of sequences.

    int quad_second_order,                 // Quad theory for second order,
    bool remove_2nd_order_quad_isotropic,  // remove the isotropic
                                           // contribution from the second
                                           // order quad Hamiltonian.

    // spin rate, spin angle and number spinning sidebands
    int number_of_sidebands,                 // The number of sidebands
    double sample_rotation_frequency_in_Hz,  // The rotor spin frequency
    double rotor_angle_in_rad,  // The rotor angle relative to lab-frame z-axis

    // Pointer to the transition.
    float *transition,

    // powder orientation average
    // The number of triangle along the edge of octahedron.
    int integration_density,
    unsigned int integration_volume,  // 0-octant, 1-hemisphere, 2-sphere.
    bool interpolation);

extern void __mrsimulator_core(
    // spectrum information and related amplitude
    double *
        spec,  // The pointer to the amplitude vector representing the spectrum.
    isotopomer_ravel *ravel_isotopomer,    // isotopomer structure
    bool remove_2nd_order_quad_isotropic,  // remove the isotropic
                                           // contribution from the second
                                           // order quad Hamiltonian.

    // Pointer to the transitions. transition[0] = mi and transition[1] = mf
    float *transition, MRS_sequence *the_sequence, int n_sequence,
    MRS_fftw_scheme *fftw_scheme, MRS_averaging_scheme *scheme,
    bool interpolation);
