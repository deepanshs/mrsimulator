// -*- coding: utf-8 -*-
//
//  simulation.h
//
//  @copyright Deepansh J. Srivastava, 2019-2021.
//  Created by Deepansh J. Srivastava, Apr 11, 2019.
//  Contact email = srivastava.89@osu.edu
//

#include "method.h"
#include "mrsimulator.h"
#include "octahedron.h"
// headerdoc

extern void mrsimulator_core(
    // spectrum information and related amplitude
    double *spec,                // The amplitude of the spectrum.
    double spectral_start,       // The start of the frequency spectrum.
    double spectral_increment,   // The increment of the frequency spectrum.
    int number_of_points,        // Number of points on the frequency spectrum.
    site_struct *sites,          // Pointer to sites within a spin system.
    coupling_struct *couplings,  // Pointer to couplings within spin system.
    MRS_dimension *dimensions,   // Pointer to the spectral dimension.
    int n_dimension,             // Number of dimensions.
    int quad_second_order,       // Quad theory for second order,

    unsigned int number_of_sidebands,  // The number of sidebands.
    double rotor_frequency_in_Hz,      // The rotor spin frequency.
    double rotor_angle_in_rad,         // The rotor angle relative to lab-frame z-axis.

    float *transition_pathway,  // Pointer to a list of transitions.

    // powder orientation average
    int integration_density,  // The number of triangle along the edge of octahedron.
    unsigned int integration_volume,  // 0-octant, 1-hemisphere, 2-sphere.
    bool interpolation, bool *freq_contrib, double *affine_matrix);

extern void __mrsimulator_core(
    // spectrum information and related amplitude
    double *spec,                // Pointer to the spectrum array.
    site_struct *sites,          // Pointer to a list of sites within a spin system.
    coupling_struct *couplings,  // Pointer to a list of couplings within a spin system.

    // A pointer to a spin transition pathway packed as a series of transitions. Each
    // transition is a list of quantum numbers packed as quantum numbers from the
    // initial energy state followed by the quantum numbers from the final energy state.
    // The energy states are given in Zeeman basis.
    float *transition_pathway,
    int n_dimension,               // The total number of spectroscopic dimensions.
    MRS_dimension *dimensions,     // Pointer to MRS_dimension structure.
    MRS_fftw_scheme *fftw_scheme,  // Pointer to the fftw scheme.
    MRS_averaging_scheme *scheme,  // Pointer to the powder averaging scheme.
    bool interpolation,            // If true, perform a 1D interpolation.

    /**
     * Each event consists of the following freq contrib ordered as
     * 1. Shielding 1st order 0th rank
     * 2. Shielding 1st order 2th rank
     * 3. Quad 1st order 2th rank
     * 4. Quad 2st order 0th rank
     * 5. Quad 2st order 2th rank
     * 6. Quad 2st order 4th rank
     *
     * The freq contrib from each event is a list of boolean, where 1 mean allow
     * frequency contribution and 0 means remove. The `freq_contrib` variable is
     * a stack of boolean list, where the stack is ordered according to the
     * events.
     */
    bool *freq_contrib,
    double *affine_matrix  // Affine transformation matrix.
);
