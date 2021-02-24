// -*- coding: utf-8 -*-
//
//  sideband_simulator.c
//
//  @copyright Deepansh J. Srivastava, 2019-2021.
//  Created by Deepansh J. Srivastava, Apr 11, 2019
//  Contact email = srivastava.89@osu.edu
//

#include "simulation.h"

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
int FREQ_CONTRIB_INCREMENT = 6;

static inline void __zero_components(double *R0, complex128 *R2, complex128 *R4) {
  R0[0] = 0.0;
  vm_double_zeros(10, (double *)R2);
  vm_double_zeros(18, (double *)R4);
}

static inline void one_dimensional_averaging(MRS_dimension *dimensions,
                                             MRS_averaging_scheme *scheme,
                                             MRS_fftw_scheme *fftw_scheme, double *spec,
                                             unsigned int number_of_sidebands) {
  unsigned int i, j, evt, step_vector = 0, address;
  MRS_plan *plan;
  MRS_event *event;
  int size = scheme->total_orientations * number_of_sidebands;
  double *freq_amp = malloc_double(size);
  double offset, offset1;

  vm_double_ones(size, freq_amp);

  // offset = plan->vr_freq[i] + plan->isotropic_offset +
  //          dimensions[dim].normalize_offset;
  offset = dimensions[0].normalize_offset + dimensions[0].R0_offset;
  for (evt = 0; evt < dimensions[0].n_events; evt++) {
    event = &dimensions[0].events[evt];
    plan = event->plan;
    // offset += plan->R0_offset;
    vm_double_multiply_inplace(size, event->freq_amplitude, 1, freq_amp, 1);
  }

  for (j = 0; j < scheme->octant_orientations; j++) {
    cblas_dscal(plan->n_octants * number_of_sidebands, plan->norm_amplitudes[j],
                &freq_amp[j], scheme->octant_orientations);
  }

  for (i = 0; i < number_of_sidebands; i++) {
    offset1 = offset + plan->vr_freq[i] * dimensions[0].inverse_increment;
    if ((int)offset1 >= 0 && (int)offset1 <= dimensions[0].count) {
      step_vector = i * scheme->total_orientations;
      for (j = 0; j < plan->n_octants; j++) {
        address = j * scheme->octant_orientations;
        // Add offset(isotropic + sideband_order) to the local frequency
        // from [n to n+octant_orientation]
        vm_double_add_offset(scheme->octant_orientations,
                             &dimensions[0].local_frequency[address], offset1,
                             dimensions[0].freq_offset);
        // Perform tenting on every sideband order over all orientations
        octahedronInterpolation(spec, dimensions[0].freq_offset,
                                scheme->integration_density, &freq_amp[step_vector], 1,
                                dimensions[0].count);
        step_vector += scheme->octant_orientations;
      }
    }
  }
  free(freq_amp);
}

static inline void two_dimensional_averaging(MRS_dimension *dimensions,
                                             MRS_averaging_scheme *scheme,
                                             MRS_fftw_scheme *fftw_scheme, double *spec,
                                             unsigned int number_of_sidebands,
                                             double *affine_matrix) {
  unsigned int i, k, j, evt;
  unsigned int step_vector_i = 0, step_vector_k = 0, address;
  MRS_plan *planA, *planB;
  MRS_event *event;
  int size = scheme->total_orientations * number_of_sidebands;
  double *freq_ampA = malloc_double(size);
  double *freq_ampB = malloc_double(size);
  double *freq_amp = malloc_double(scheme->total_orientations);
  double offset0, offset1, offsetA, offsetB;
  double *dim0, *dim1;
  double norm0, norm1;

  vm_double_ones(size, freq_ampA);
  vm_double_ones(size, freq_ampB);

  dim0 = dimensions[0].local_frequency;
  dim1 = dimensions[1].local_frequency;

  // scale and shear the first dimension.
  if (affine_matrix[0] != 1) {
    cblas_dscal(scheme->total_orientations, affine_matrix[0], dim0, 1);
    // dimensions[0].R0_offset *= affine_matrix[0];
  }
  if (affine_matrix[1] != 0) {
    cblas_daxpy(scheme->total_orientations, affine_matrix[1], dim1, 1, dim0, 1);
    // dimensions[0].R0_offset += affine_matrix[1] *
    // dimensions[1].R0_offset;
  }

  // scale and shear the second dimension.
  if (affine_matrix[3] != 1) {
    cblas_dscal(scheme->total_orientations, affine_matrix[3], dim1, 1);
    // dimensions[1].R0_offset *= affine_matrix[3];
  }
  if (affine_matrix[2] != 0) {
    cblas_daxpy(scheme->total_orientations, affine_matrix[2], dim0, 1, dim1, 1);
    // dimensions[1].R0_offset += affine_matrix[2] *
    // dimensions[0].R0_offset;
  }

  // offset = plan->vr_freq[i] + plan->isotropic_offset +
  //          dimensions[dim].normalize_offset;

  offset0 = dimensions[0].R0_offset;
  for (evt = 0; evt < dimensions[0].n_events; evt++) {
    event = &dimensions[0].events[evt];
    planA = event->plan;
    // offset0 += plan->R0_offset;
    vm_double_multiply_inplace(size, event->freq_amplitude, 1, freq_ampA, 1);
  }

  offset1 = dimensions[1].R0_offset;
  for (evt = 0; evt < dimensions[1].n_events; evt++) {
    event = &dimensions[1].events[evt];
    planB = event->plan;
    // offset1 += plan->R0_offset;
    vm_double_multiply_inplace(size, event->freq_amplitude, 1, freq_ampB, 1);
  }

  for (j = 0; j < scheme->octant_orientations; j++) {
    cblas_dscal(planA->n_octants * number_of_sidebands, planA->norm_amplitudes[j],
                &freq_ampB[j], scheme->octant_orientations);
  }

  for (i = 0; i < number_of_sidebands; i++) {
    offsetA = offset0 + planA->vr_freq[i] * dimensions[0].inverse_increment;
    for (k = 0; k < number_of_sidebands; k++) {
      offsetB = offset1 + planB->vr_freq[k] * dimensions[1].inverse_increment;

      norm0 = offsetA;
      norm1 = offsetB;

      // scale and shear the offsets
      norm0 *= affine_matrix[0];
      norm0 += affine_matrix[1] * offsetB;

      norm1 *= affine_matrix[3];
      norm1 += affine_matrix[2] * norm0;

      norm0 += dimensions[0].normalize_offset;
      norm1 += dimensions[1].normalize_offset;

      if ((int)norm0 >= 0 && (int)norm0 <= dimensions[0].count) {
        step_vector_i = i * scheme->total_orientations;
        // for (k = 0; k < number_of_sidebands; k++) {
        //   offsetB =
        //       offset1 + plan->vr_freq[k] * dimensions[1].inverse_increment;
        //   norm1 = offsetB + dimensions[1].normalize_offset;
        if ((int)norm1 >= 0 && (int)norm1 <= dimensions[1].count) {
          step_vector_k = k * scheme->total_orientations;

          // step_vector = 0;
          for (j = 0; j < planA->n_octants; j++) {
            address = j * scheme->octant_orientations;
            // Add offset(isotropic + sideband_order) to the local frequency
            // from [n to n+octant_orientation]
            vm_double_add_offset(scheme->octant_orientations, &dim0[address], norm0,
                                 dimensions[0].freq_offset);
            vm_double_add_offset(scheme->octant_orientations, &dim1[address], norm1,
                                 dimensions[1].freq_offset);

            vm_double_multiply(scheme->total_orientations,
                               &freq_ampA[step_vector_i + address],
                               &freq_ampB[step_vector_k + address], freq_amp);
            // Perform tenting on every sideband order over all orientations
            octahedronInterpolation2D(spec, dimensions[0].freq_offset,
                                      dimensions[1].freq_offset,
                                      scheme->integration_density, freq_amp, 1,
                                      dimensions[0].count, dimensions[1].count);
            // step_vector += scheme->octant_orientations;
          }
        }
      }
    }
  }
  free(freq_amp);
  free(freq_ampA);
  free(freq_ampB);
}

// Calculate spectrum from the spin systems for a single transition.
void __mrsimulator_core(
    // spectrum information and related amplitude
    double *spec,               // Pointer to the spectrum array.
    site_struct *sites,         // Pointer to a list of sites within a spin system.
    coupling_struct *couplings, // Pointer to a list of couplings within a spin system.

    // A pointer to a spin transition pathway packed as a series of transitions. Each
    // transition is a list of quantum numbers packed as quantum numbers from the
    // initial energy state followed by the quantum numbers from the final energy state.
    // The energy states are given in Zeeman basis.
    float *transition_pathway,
    int n_dimension,              // The total number of spectroscopic dimensions.
    MRS_dimension *dimensions,    // Pointer to MRS_dimension structure.
    MRS_fftw_scheme *fftw_scheme, // Pointer to the fftw scheme.
    MRS_averaging_scheme *scheme, // Pointer to the powder averaging scheme.
    bool interpolation,           // If true, perform a 1D interpolation.

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
    double *affine_matrix // Affine transformation matrix.
) {
  /*
  The sideband computation is based on the method described by Eden and Levitt
  et. al. `Computation of Orientational Averages in Solid-State NMR by Gaussian
  Spherical Quadrature` JMR, 132, 1998. https://doi.org/10.1006/jmre.1998.1427
  */
  bool refresh;
  unsigned int evt;
  int dim;
  double B0_in_T, fraction;

  // Allocate memory for zeroth, second, and fourth-rank tensor components.
  double R0 = 0.0;
  complex128 *R2 = malloc_complex128(5);
  complex128 *R4 = malloc_complex128(9);

  // Allocate memory for zeroth, second, and fourth-rank temporary tensor components.
  double R0_temp = 0.0;
  complex128 *R2_temp = malloc_complex128(5);
  complex128 *R4_temp = malloc_complex128(9);

  double *spec_site_ptr;
  // `transition_increment` is the step size to the next transition within the pathway.
  int transition_increment = 2 * sites->number_of_sites;

  MRS_plan *plan;
  MRS_event *event;

  openblas_set_num_threads(1);

  // spec_site = site * dimensions[0].count;
  spec_site_ptr = &spec[0];

  // Loop over the dimensionn.
  for (dim = 0; dim < n_dimension; dim++) {
    refresh = 1;
    // Loop over the events per dimension.
    for (evt = 0; evt < dimensions[dim].n_events; evt++) {
      event = &dimensions[dim].events[evt];
      plan = event->plan;
      B0_in_T = event->magnetic_flux_density_in_T;
      fraction = event->fraction;

      /* Initialize with zeroing all spatial components */
      __zero_components(&R0, R2, R4);

      /* Rotate all frequency components from PAS to a common frame */
      MRS_rotate_components_from_PAS_to_common_frame(
          sites,              // Pointer to a list of sites within a spin system.
          couplings,          // Pointer to a list of couplings within a spin system.
          transition_pathway, // Pointer to a list of transition.
          plan->allow_fourth_rank, // If 1, prepare for 4th rank computation.
          &R0,                     // The R0 components.
          R2,                      // The R2 components.
          R4,                      // The R4 components.
          &R0_temp,                // The temporary R0 components.
          R2_temp,                 // The temporary R2 components.
          R4_temp,                 // The temporary R4 components.
          B0_in_T,                 // Magnetic flux density in T.
          freq_contrib             // The pointer to freq contribs boolean.
      );

      // The number 6 comes from the six types of pre-listed freq contributions.
      freq_contrib += FREQ_CONTRIB_INCREMENT;

      /* Get frequencies and amplitudes per octant .................................. */
      /* IMPORTANT: Always evalute the frequencies before the amplitudes. */
      MRS_get_normalized_frequencies_from_plan(scheme, plan, R0, R2, R4, refresh,
                                               &dimensions[dim], fraction);
      MRS_get_amplitudes_from_plan(scheme, plan, fftw_scheme, 1);

      /* Copy the amplitudes from the `fftw_scheme->vector` to the
       * `event->freq_amplitude` for each event within the dimension.*/
      if (plan->number_of_sidebands != 1) {
        cblas_dcopy(plan->size, (double *)fftw_scheme->vector, 2, event->freq_amplitude,
                    1);
      }
      transition_pathway += transition_increment;
      refresh = 0;
    } // end events
  }   // end dimensions

  free(R2);
  free(R4);
  free(R2_temp);
  free(R4_temp);

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
  unsigned int i, j, step_vector, address;
  double offset, offset0;
  if (n_dimension == 1 && dimensions[0].n_events == 1) {
    /**
     * If the number of sidebands is 1, the sideband amplitude at every
     * sideband order is one. In this case, update the `fftw_scheme->vector` is
     * the same as the weights from the orientation averaging,
     */
    if (plan->number_of_sidebands == 1) {
      /* Copy the plan->norm_amplitudes to fftw_scheme->vector. */
      for (j = 0; j < plan->n_octants; j++) {
        cblas_dcopy(scheme->octant_orientations, plan->norm_amplitudes, 1,
                    (double *)&fftw_scheme->vector[j * scheme->octant_orientations], 2);
      }
    } else {
      /**
       * Scale the absolute value square with the powder scheme weights. Only
       * the real part is scaled and the imaginary part is left as is.
       */
      for (j = 0; j < scheme->octant_orientations; j++) {
        cblas_dscal(plan->n_octants * plan->number_of_sidebands,
                    plan->norm_amplitudes[j], (double *)&fftw_scheme->vector[j],
                    2 * scheme->octant_orientations);
      }
    }

    offset0 = dimensions[0].normalize_offset + dimensions[0].R0_offset;

    for (i = 0; i < plan->number_of_sidebands; i++) {
      offset = plan->vr_freq[i] * dimensions[0].inverse_increment + offset0;
      if ((int)offset >= 0 && (int)offset <= dimensions[0].count) {
        step_vector = i * scheme->total_orientations;
        for (j = 0; j < plan->n_octants; j++) {
          address = j * scheme->octant_orientations;

          // Add offset(isotropic + sideband_order) to the local frequency
          // from [n to n+octant_orientation]
          vm_double_add_offset(scheme->octant_orientations,
                               &dimensions[0].local_frequency[address], offset,
                               dimensions[0].freq_offset);
          // Perform tenting on every sideband order over all orientations.
          octahedronInterpolation(
              spec_site_ptr, dimensions[0].freq_offset, scheme->integration_density,
              (double *)&fftw_scheme->vector[step_vector], 2, dimensions[0].count);
          step_vector += scheme->octant_orientations;
        }
      }
    }
    return;
  }

  if (interpolation) {
    if (n_dimension == 1) {
      one_dimensional_averaging(dimensions, scheme, fftw_scheme, spec,
                                plan->number_of_sidebands);
      return;
    }

    if (n_dimension == 2) {
      two_dimensional_averaging(dimensions, scheme, fftw_scheme, spec,
                                plan->number_of_sidebands, affine_matrix);
      return;
    }
  }
}

void mrsimulator_core(
    // spectrum information and related amplitude
    double *spec,               // The amplitude of the spectrum.
    double coordinates_offset,  // The start of the frequency spectrum.
    double increment,           // The increment of the frequency spectrum.
    int count,                  // Number of points on the frequency spectrum.
    site_struct *sites,         // Pointer to a list of sites wiithin a spin system.
    coupling_struct *couplings, // Pointer to a list of couplings within a spin system.
    MRS_dimension *dimensions,  // the dimensions in the method.
    int n_dimension,            // The number of dimension.
    int quad_second_order,      // Quad theory for second order,

    // spin rate, spin angle and number spinning sidebands
    unsigned int number_of_sidebands,       // The number of sidebands
    double sample_rotation_frequency_in_Hz, // The rotor spin frequency
    double rotor_angle_in_rad, // The rotor angle relative to lab-frame z-axis

    // Pointer to the a list of transitions.
    float *transition_pathway,

    // powder orientation average
    int integration_density, // The number of triangle along the edge of octahedron
    unsigned int integration_volume, // 0-octant, 1-hemisphere, 2-sphere.
    bool interpolation, bool *freq_contrib, double *affine_matrix) {
  // int num_process = openblas_get_num_procs();
  // int num_threads = openblas_get_num_threads();
  openblas_set_num_threads(1);
  // printf("%d processors", num_process);
  // printf("%d threads", num_threads);
  // int parallel = openblas_get_parallel();
  // printf("%d parallel", parallel);

  bool allow_fourth_rank = false;
  if (sites[0].spin[0] > 0.5 && quad_second_order == 1) {
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
      spec, // The amplitude of the spectrum.

      sites,              // Pointer to a list of sites within the spin system.
      couplings,          // Pointer to a list of couplings within a spin system.
      transition_pathway, // Pointer to a list of transition.

      n_dimension, dimensions, fftw_scheme, scheme, interpolation, freq_contrib,
      affine_matrix);

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
