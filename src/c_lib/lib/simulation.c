// -*- coding: utf-8 -*-
//
//  sideband_simulator.c
//
//  @copyright Deepansh J. Srivastava, 2019-2021.
//  Created by Deepansh J. Srivastava, Apr 11, 2019.
//  Contact email = srivastava.89@osu.edu
//

#include "simulation.h"

#include "frequency_averaging.h"

/**
 * Each event consists of the following freq contrib ordered as
 * 1. Shielding 1st order 0th rank
 * 2. Shielding 1st order 2th rank
 * 3. Quad 1st order 2th rank
 * 4. Quad 2st order 0th rank
 * 5. Quad 2st order 2th rank
 * 6. Quad 2st order 4th rank
 *
 * The freq contrib from each event is a list of boolean, where 1 mean include frequency
 * contribution and 0 means exclude. The `freq_contrib` variable is a stack of boolean
 * list, where the stack is ordered according to the events. The variable
 * `FREQ_CONTRIB_INCREMENT` is the length of the freq contribs.
 */
int FREQ_CONTRIB_INCREMENT = 18;

static inline void __zero_components(double *R0, complex128 *R2, complex128 *R4) {
  *R0 = 0.0;
  vm_double_zeros(10, (double *)R2);
  vm_double_zeros(18, (double *)R4);
}

// Calculate spectrum from the spin systems for a single transition.
void __mrsimulator_core(
    // spectrum information and related amplitude
    double *spec,                // Pointer to the spectrum array (complex).
    site_struct *sites,          // Pointer to a list of sites within a spin system.
    coupling_struct *couplings,  // Pointer to a list of couplings within a spin system.

    // A pointer to a spin transition pathway packed as a series of transitions. Each
    // transition is a list of quantum numbers packed as quantum numbers from the
    // initial energy state followed by the quantum numbers from the final energy state.
    // The energy states are given in Zeeman basis.
    float *transition_pathway,          // Pointer to the transition pathway,
    double *transition_pathway_weight,  // The complex weight of transition pathway.
    int n_dimension,                    // The total number of spectroscopic dimensions.
    MRS_dimension *dimensions,          // Pointer to MRS_dimension structure.
    MRS_fftw_scheme *fftw_scheme,       // Pointer to the fftw scheme.
    MRS_averaging_scheme *scheme,       // Pointer to the powder averaging scheme.
    unsigned int iso_intrp,       // Isotropic interpolation scheme (linear | Gaussian)
    unsigned char *freq_contrib,  // A list of freq_contrib booleans.
    double *affine_matrix         // Affine transformation matrix.
) {
  /*
  The sideband computation is based on the method described by Eden and Levitt
  et al. `Computation of Orientational Averages in Solid-State NMR by Gaussian
  Spherical Quadrature` JMR, 132, 1998. https://doi.org/10.1006/jmre.1998.1427
  */
  unsigned char is_spectral;
  unsigned int evt;
  int dim, total_pts = scheme->n_gamma * scheme->total_orientations;
  double B0_in_T, fraction, duration;
  vm_double_zeros(total_pts, scheme->phase);

  // Allocate memory for zeroth, second, and fourth-rank tensor components.
  // variable with _temp allocate temporary memory for tensor components
  double R0 = 0.0, R0_temp = 0.0;
  complex128 R2[5], R4[9], R2_temp[5], R4_temp[9];

  // `transition_increment` is the step size to the next transition within the pathway.
  int transition_increment = 2 * sites->number_of_sites;
  float *transition = transition_pathway;

  MRS_plan *plan = NULL;
  MRS_event *event;

  // openblas_set_num_threads(1);

  // Loop over the dimension.
  for (dim = 0; dim < n_dimension; dim++) {
    // Reset the freqs to zero at the start of each spectral dimension.
    vm_double_zeros(total_pts, dimensions[dim].local_frequency);
    dimensions[dim].R0_offset = 0.0;

    plan = dimensions[dim].events->plan;
    vm_double_ones(plan->size, dimensions[dim].freq_amplitude);
    // Loop over the events per dimension.
    for (evt = 0; evt < dimensions[dim].n_events; evt++) {
      event = &dimensions[dim].events[evt];
      plan = event->plan;
      B0_in_T = event->magnetic_flux_density_in_T;
      fraction = event->fraction;
      duration = event->duration;
      is_spectral = event->is_spectral;

      /* Initialize with zeroing all spatial components */
      __zero_components(&R0, R2, R4);

      /* Rotate all frequency components from PAS to a common frame */
      MRS_rotate_components_from_PAS_to_common_frame(
          sites,                 // Pointer to a list of sites within a spin system.
          couplings,             // Pointer to a list of couplings within a spin system.
          transition,            // Pointer to a single transition.
          plan->allow_4th_rank,  // If 1, prepare for 4th rank computation.
          &R0,                   // The R0 components.
          R2,                    // The R2 components.
          R4,                    // The R4 components.
          &R0_temp,              // The temporary R0 components.
          R2_temp,               // The temporary R2 components.
          R4_temp,               // The temporary R4 components.
          B0_in_T,               // Magnetic flux density in T.
          freq_contrib           // The pointer to freq contribs boolean.
      );

      // The number 6 comes from the six types of pre-listed freq contributions.
      freq_contrib += FREQ_CONTRIB_INCREMENT;

      /* Get frequencies and amplitudes per octant .................................. */
      /* IMPORTANT: Always evaluate the frequencies before the amplitudes. */
      // NOTE: How to incorporate both "fraction" and "duration" into this function?
      // Possibly calculate normalized frequencies first, then decide if frac or dur
      MRS_get_normalized_frequencies_from_plan(
          scheme, plan, R0, R2, R4, &dimensions[dim], fraction, is_spectral, duration);
      MRS_get_amplitudes_from_plan(scheme, plan, fftw_scheme,
                                   event->event_freq_amplitude, 1);

      if (plan->number_of_sidebands != 1) {
        vm_double_multiply_inplace(plan->size, (double *)fftw_scheme->vector, 2,
                                   dimensions[dim].freq_amplitude, 1);
      }
      transition += transition_increment;  // increment to next transition
    }  // end events
  }  // end dimensions

  // calculate phase exponent of delay events
  vm_cosine_I_sine(total_pts, scheme->phase, scheme->exp_I_phase);
  cblas_zscal(total_pts, transition_pathway_weight, (double *)scheme->exp_I_phase, 1);

  /* ---------------------------------------------------------------------
   *              Delta and triangle tenting interpolation
   */

  switch (n_dimension) {
  case 1:
    one_dimensional_averaging(dimensions, scheme, spec, iso_intrp, scheme->exp_I_phase);
    break;
  case 2:
    two_dimensional_averaging(dimensions, scheme, spec, affine_matrix, iso_intrp,
                              scheme->exp_I_phase);
    break;
  }
}

void mrsimulator_core(
    // spectrum information and related amplitude
    double *spec,                // Pointer to the spectrum array (complex).
    double coordinates_offset,   // The start of the frequency spectrum.
    double increment,            // The increment of the frequency spectrum.
    int count,                   // Number of points on the frequency spectrum.
    site_struct *sites,          // Pointer to a list of sites within a spin system.
    coupling_struct *couplings,  // Pointer to a list of couplings within a spin system.
    MRS_dimension *dimensions,   // the dimensions in the method.
    int n_dimension,             // The number of dimension.
    int quad_second_order,       // Quad theory for second order,

    // spin rate, spin angle and number spinning sidebands
    unsigned int number_of_sidebands,  // The number of sidebands
    double rotor_frequency_in_Hz,      // The rotor spin frequency
    double rotor_angle_in_rad,         // The rotor angle relative to lab-frame z-axis

    // Pointer to the a list of transitions.
    float *transition_pathway, double *transition_pathway_weight,
    // powder orientation average
    int integration_density,  // The number of triangle along the edge of octahedron
    unsigned int integration_volume,  // 0-octant, 1-hemisphere, 2-sphere.
    unsigned int interpolate_type, unsigned char *freq_contrib, double *affine_matrix) {
  // int num_process = openblas_get_num_procs();
  // int num_threads = openblas_get_num_threads();
  // openblas_set_num_threads(1);
  // printf("%d processors", num_process);
  // printf("%d threads", num_threads);
  // int parallel = openblas_get_parallel();
  // printf("%d parallel", parallel);

  bool allow_4th_rank = false;
  bool interpolation = true;
  bool is_complex = true;

  if (sites[0].spin[0] > 0.5 && quad_second_order == 1) {
    allow_4th_rank = true;
  }

  // check for spinning speed
  if (rotor_frequency_in_Hz < 1.0e-3) {
    rotor_frequency_in_Hz = 1.0e9;
    rotor_angle_in_rad = 0.0;
    number_of_sidebands = 1;
  }

  MRS_averaging_scheme *scheme =
      MRS_create_averaging_scheme(integration_density, allow_4th_rank, 9,
                                  integration_volume, interpolation, is_complex);

  MRS_fftw_scheme *fftw_scheme =
      create_fftw_scheme(scheme->total_orientations, number_of_sidebands);

  // gettimeofday(&all_site_time, NULL);
  __mrsimulator_core(
      // spectrum information and related amplitude
      spec,                // Pointer to the spectrum array (complex).
      sites,               // Pointer to a list of sites within the spin system.
      couplings,           // Pointer to a list of couplings within a spin system.
      transition_pathway,  // Pointer to a list of transition.
      transition_pathway_weight, n_dimension, dimensions, fftw_scheme, scheme,
      interpolate_type, freq_contrib, affine_matrix);

  // gettimeofday(&end, NULL);
  // clock_time = (double)(end.tv_usec - begin.tv_usec) / 1000000. +
  //              (double)(end.tv_sec - begin.tv_sec);
  // printf("time %f s\n", clock_time);
  // cpu_time_[0] += clock_time;

  /* clean up */
  MRS_free_fftw_scheme(fftw_scheme);
  MRS_free_averaging_scheme(scheme);
  // MRS_free_plan(plan);
}
