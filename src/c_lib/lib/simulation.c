// -*- coding: utf-8 -*-
//
//  sideband_simulator.c
//
//  @copyright Deepansh J. Srivastava, 2019-2020.
//  Created by Deepansh J. Srivastava, Apr 11, 2019
//  Contact email = deepansh2012@gmail.com
//

#include "simulation.h"

static inline void __zero_components(double *R0, complex128 *R2,
                                     complex128 *R4) {
  R0[0] = 0.0;
  vm_double_zeros(10, (double *)R2);
  vm_double_zeros(18, (double *)R4);
}

void __mrsimulator_core(
    // spectrum information and related amplitude
    double *spec,  // amplitude vector representing the spectrum.

    // A pointer to the isotopomer_ravel structure containing information about
    // the sites within an isotopomer.
    isotopomer_ravel *ravel_isotopomer,

    // remove the isotropic contribution from the second order quad Hamiltonian.
    bool remove_2nd_order_quad_isotropic,

    // A pointer to a spin transition packed as quantum numbers from the initial
    // energy state followed by the quantum numbers from the final energy state.
    // The energy states are given in Zeeman basis.
    float *transition,

    // A pointer to MRS_sequence structure containing information on the events
    // per spectroscopic dimension.
    MRS_sequence *the_sequence,

    // The total number of spectroscopic dimensions.
    int n_sequence,

    // A pointer to the fftw scheme.
    MRS_fftw_scheme *fftw_scheme,

    // A pointer to the powder averaging scheme.
    MRS_averaging_scheme *scheme,

    // if true, perform a 1D interpolation
    bool interpolation) {
  /*
  The sideband computation is based on the method described by Eden and Levitt
  et. al. `Computation of Orientational Averages in Solid-State NMR by Gaussian
  Spherical Quadrature` JMR, 132, 1998. https://doi.org/10.1006/jmre.1998.1427
  */

  unsigned int j, evt, step_vector = 0, address;
  int i, seq;
  double offset, B0_in_T;

  double R0 = 0.0;
  complex128 *R2 = malloc_complex128(5);
  complex128 *R4 = malloc_complex128(9);

  double R0_temp = 0.0;
  complex128 *R2_temp = malloc_complex128(5);
  complex128 *R4_temp = malloc_complex128(9);

  double *spec_site_ptr;
  int transition_increment = 2 * ravel_isotopomer->number_of_sites;
  MRS_plan *plan;

  // spec_site = site * the_sequence[0].count;
  spec_site_ptr = &spec[0];

  // Loop over the sequence.
  for (seq = 0; seq < n_sequence; seq++) {
    // Loop over the events per sequence.
    for (evt = 0; evt < the_sequence[seq].n_events; evt++) {
      plan = the_sequence[seq].events[evt].plan;
      B0_in_T = the_sequence[seq].events[evt].magnetic_flux_density_in_T;

      /* Initialize with zeroing all spatial components */
      __zero_components(&R0, R2, R4);

      /* Rotate all frequency components from PAS to a common frame */
      MRS_rotate_components_from_PAS_to_common_frame(
          ravel_isotopomer,         // isotopomer structure
          transition,               // the transition
          plan->allow_fourth_rank,  // if 1, prepare for 4th rank computation
          &R0,                      // the R0 components
          R2,                       // the R2 components
          R4,                       // the R4 components
          &R0_temp,                 // the temporary R0 components
          R2_temp,                  // the temporary R2 components
          R4_temp,                  // the temporary R4 components
          remove_2nd_order_quad_isotropic,  // if true, remove second order quad
                                            // isotropic shift
          B0_in_T                           // magnetic flux density in T.
      );

      // Add a loop over all couplings.. here

      /* Get frequencies and amplitudes per octant ......................... */
      /* Always evalute the frequencies before the amplitudes. */
      MRS_get_normalized_frequencies_from_plan(
          scheme, plan, R0, R2, R4, 1, the_sequence[seq].normalize_offset,
          the_sequence[seq].inverse_increment);
      MRS_get_amplitudes_from_plan(scheme, plan, fftw_scheme, 1);

      transition += transition_increment;
    }  // end events

    /* ---------------------------------------------------------------------
     *              Calculating the tent for every sideband
     * Allowing only sidebands that are within the spectral bandwidth
     *
     * for (i = 0; i < plan->number_of_sidebands; i++) {
     *   offset = plan->vr_freq[i] + plan->isotropic_offset;
     *   if ((int)offset >= 0 && (int)offset <= dimension->count) {

     *     vm_double_ramp(plan->octant_orientations,
     plan->local_frequency, 1.0,
     *                    offset, plan->freq_offset);
     *     octahedronInterpolation(
     *         spec_site_ptr, plan->freq_offset,
     *         plan->integration_density,
     *         (double *)&plan->vector[i * plan->octant_orientations], 2,
     *         dimension->count);
     *   }
     * }
     */
    if (interpolation) {
      for (i = 0; i < plan->number_of_sidebands; i++) {
        offset = plan->vr_freq[i] + plan->isotropic_offset;
        if ((int)offset >= 0 && (int)offset <= the_sequence[0].count) {
          step_vector = i * scheme->total_orientations;
          for (j = 0; j < plan->n_octants; j++) {
            address = j * scheme->octant_orientations;

            // Add offset(isotropic + sideband_order) to the local frequency
            // from [n to n+octant_orientation]
            vm_double_ramp(scheme->octant_orientations,
                           &scheme->local_frequency[address], 1.0, offset,
                           scheme->freq_offset);
            // Perform tenting on every sideband order over all orientations
            octahedronInterpolation(spec_site_ptr, scheme->freq_offset,
                                    scheme->integration_density,
                                    (double *)&fftw_scheme->vector[step_vector],
                                    2, the_sequence[0].count);
            step_vector += scheme->octant_orientations;
          }
        }
      }
    }

    // gettimeofday(&end_site_time, NULL);
    // clock_time =
    //     (double)(end_site_time.tv_usec - start_site_time.tv_usec) /
    //     1000000.
    //     + (double)(end_site_time.tv_sec - start_site_time.tv_sec);
    // printf("Total time per site %f \n", clock_time);
  }
}

void mrsimulator_core(
    // spectrum information and related amplitude
    double *spec,               // The amplitude of the spectrum.
    double coordinates_offset,  // The start of the frequency spectrum.
    double increment,           // The increment of the frequency spectrum.
    int count,                  // Number of points on the frequency spectrum.
    isotopomer_ravel *ravel_isotopomer,    // SpinSystem structure
    MRS_sequence *the_sequence,            // the sequences in the method.
    int n_sequence,                        // The number of sequence.
    int quad_second_order,                 // Quad theory for second order,
    bool remove_2nd_order_quad_isotropic,  // remove the isotropic
                                           // contribution from the second
                                           // order quad interaction.

    // spin rate, spin angle and number spinning sidebands
    int number_of_sidebands,                 // The number of sidebands
    double sample_rotation_frequency_in_Hz,  // The rotor spin frequency
    double rotor_angle_in_rad,  // The rotor angle relative to lab-frame z-axis

    // Pointer to the transitions. transition[0] = mi and transition[1] = mf
    float *transition,

    // powder orientation average
    int integration_density,          // The number of triangle along the edge
                                      // of octahedron
    unsigned int integration_volume,  // 0-octant, 1-hemisphere, 2-sphere.
    bool interpolation) {
  // int num_process = openblas_get_num_procs();
  // int num_threads = openblas_get_num_threads();
  // // openblas_set_num_threads(1);
  // printf("%d processors", num_process);
  // printf("%d threads", num_threads);
  // int parallel = openblas_get_parallel();
  // printf("%d parallel", parallel);

  bool allow_fourth_rank = false;
  if (ravel_isotopomer[0].spin[0] > 0.5 && quad_second_order == 1) {
    allow_fourth_rank = true;
  }

  // check for spinning speed
  if (sample_rotation_frequency_in_Hz < 1.0e-3) {
    sample_rotation_frequency_in_Hz = 1.0e9;
    rotor_angle_in_rad = 0.0;
    number_of_sidebands = 1;
  }

  MRS_averaging_scheme *scheme = MRS_create_averaging_scheme(
      integration_density, allow_fourth_rank, integration_volume);

  MRS_fftw_scheme *fftw_scheme =
      create_fftw_scheme(scheme->total_orientations, number_of_sidebands);

  // gettimeofday(&all_site_time, NULL);
  __mrsimulator_core(
      // spectrum information and related amplitude
      spec,  // The amplitude of the spectrum.

      ravel_isotopomer,  // isotopomer structure

      remove_2nd_order_quad_isotropic,  // remove the isotropic contribution
                                        // from the second order quad
                                        // Hamiltonian.

      // Pointer to the transitions. transition[0] = mi and transition[1] = mf
      transition,

      the_sequence, n_sequence, fftw_scheme, scheme, interpolation);

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
