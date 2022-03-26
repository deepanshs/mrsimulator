// -*- coding: utf-8 -*-
//
//  dimensional_averaging.c
//
//  @copyright Deepansh J. Srivastava, 2019-2021.
//  Created by Deepansh J. Srivastava, Mar 2, 2021.
//  Contact email = srivastava.89@osu.edu
//

#include "frequency_averaging.h"

// Multiply the amplitudes from each event to the amplitudes from the first event.
// static inline void get_multi_event_amplitudes(int n_events, MRS_event *restrict
// event,
//                                               int size) {
//   double *amp = event->freq_amplitude;
//   n_events--;
//   while (n_events-- > 0) {
//     event++;
//     vm_double_multiply_inplace(size, event->freq_amplitude, 1, amp, 1);
//   }
// }

static inline void averaging_1D(MRS_dimension *dimensions, MRS_averaging_scheme *scheme,
                                MRS_fftw_scheme *fftw_scheme, double *spec,
                                double transition_pathway_weight,
                                unsigned int iso_intrp) {
  unsigned int i, j, k1, address;
  unsigned int nt = scheme->integration_density, npts = scheme->octant_orientations;

  double offset_0, offset;
  double *freq = dimensions->local_frequency, *amps = dimensions->freq_amplitude;

  bool delta_interpolation = false;
  MRS_plan *plan = dimensions->events->plan;

  /**
   * If the number of sidebands is 1, the sideband amplitude at every
   * sideband order is one. In this case, update the `fftw_scheme->vector` is
   * the same as the weights from the orientation averaging,
   */
  switch (plan->number_of_sidebands) {
  case 1:
    /* Copy the plan->norm_amplitudes to fftw_scheme->vector. */
    for (j = 0; j < plan->n_octants; j++) {
      cblas_dcopy(npts, plan->norm_amplitudes, 1, &amps[j * npts], 1);
    }
    break;
  default:
    /**
     * Scale the absolute value square with the powder scheme weights. Only
     * the real part is scaled and the imaginary part is left as is. */
    for (j = 0; j < npts; j++) {
      cblas_dscal(plan->n_octants * plan->number_of_sidebands, plan->norm_amplitudes[j],
                  &amps[j], npts);
    }
  }

  cblas_dscal(dimensions->events->plan->size, transition_pathway_weight, amps, 1);

  offset_0 = dimensions->normalize_offset + dimensions->R0_offset;
  if (fabs(*freq - freq[nt]) < TOL && fabs(*freq - freq[npts - 1]) < TOL)
    delta_interpolation = true;

  if (delta_interpolation) {
    offset_0 += *freq;
    for (i = 0; i < plan->number_of_sidebands; i++) {
      offset = offset_0 + plan->vr_freq[i];
      if ((int)offset >= 0 && (int)offset <= dimensions->count) {
        k1 = i * scheme->total_orientations;
        j = 0;
        while (j++ < plan->n_octants) {
          octahedronDeltaInterpolation(nt, &offset, &amps[k1], 1, dimensions->count,
                                       spec, iso_intrp);
          k1 += npts;
        }
      }
    }
    return;
  }

  for (i = 0; i < plan->number_of_sidebands; i++) {
    offset = offset_0 + plan->vr_freq[i];
    if ((int)offset >= 0 && (int)offset <= dimensions->count) {
      k1 = i * scheme->total_orientations;
      address = 0;
      for (j = 0; j < plan->n_octants; j++) {
        // Add offset(isotropic + sideband_order) to the local frequencies.
        vm_double_add_offset(npts, &freq[address], offset, dimensions->freq_offset);
        // Perform tenting on every sideband order over all orientations.
        octahedronInterpolation(spec, dimensions->freq_offset, nt, &amps[k1], 1,
                                dimensions->count);
        k1 += npts;
        address += npts;
      }
    }
  }
}

void one_dimensional_averaging(MRS_dimension *dimensions, MRS_averaging_scheme *scheme,
                               MRS_fftw_scheme *fftw_scheme, double *spec,
                               double transition_pathway_weight,
                               unsigned int iso_intrp) {
  averaging_1D(dimensions, scheme, fftw_scheme, spec, transition_pathway_weight,
               iso_intrp);
}

void two_dimensional_averaging(MRS_dimension *dimensions, MRS_averaging_scheme *scheme,
                               MRS_fftw_scheme *fftw_scheme, double *spec,
                               double transition_pathway_weight,
                               unsigned int number_of_sidebands, double *affine_matrix,
                               unsigned int iso_intrp) {
  unsigned int i, k, j;
  unsigned int step_vector_i = 0, step_vector_k = 0, address;
  MRS_plan *planA, *planB;
  int size = scheme->total_orientations * number_of_sidebands;
  double *freq_ampA, *freq_ampB, *freq_amp = scheme->scrach;
  double offset0, offset1, offsetA, offsetB;
  double *freq0, *freq1;
  double norm0, norm1;

  freq0 = dimensions[0].local_frequency;
  freq1 = dimensions[1].local_frequency;

  // scale and shear the first dimension.
  if (affine_matrix[0] != 1) {
    cblas_dscal(scheme->total_orientations, affine_matrix[0], freq0, 1);
  }
  if (affine_matrix[1] != 0) {
    cblas_daxpy(scheme->total_orientations, affine_matrix[1], freq1, 1, freq0, 1);
  }

  // scale and shear the second dimension.
  if (affine_matrix[3] != 1) {
    cblas_dscal(scheme->total_orientations, affine_matrix[3], freq1, 1);
  }
  if (affine_matrix[2] != 0) {
    cblas_daxpy(scheme->total_orientations, affine_matrix[2], freq0, 1, freq1, 1);
  }

  offset0 = dimensions[0].R0_offset;
  freq_ampA = dimensions[0].freq_amplitude;
  planA = dimensions[0].events[0].plan;

  offset1 = dimensions[1].R0_offset;
  freq_ampB = dimensions[1].freq_amplitude;
  planB = dimensions[1].events[0].plan;

  for (j = 0; j < scheme->octant_orientations; j++) {
    cblas_dscal(planA->n_octants * number_of_sidebands, planA->norm_amplitudes[j],
                &freq_ampB[j], scheme->octant_orientations);
  }
  cblas_dscal(size, transition_pathway_weight, freq_ampA, 1);

  for (i = 0; i < number_of_sidebands; i++) {
    offsetA = offset0 + planA->vr_freq[i];
    for (k = 0; k < number_of_sidebands; k++) {
      offsetB = offset1 + planB->vr_freq[k];

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
        //       offset1 + plan->vr_freq[k];
        //   norm1 = offsetB + dimensions[1].normalize_offset;
        if ((int)norm1 >= 0 && (int)norm1 <= dimensions[1].count) {
          step_vector_k = k * scheme->total_orientations;

          // step_vector = 0;
          for (j = 0; j < planA->n_octants; j++) {
            address = j * scheme->octant_orientations;
            // Add offset(isotropic + sideband_order) to the local frequency
            // from [n to n+octant_orientation]
            vm_double_add_offset(scheme->octant_orientations, &freq0[address], norm0,
                                 dimensions[0].freq_offset);
            vm_double_add_offset(scheme->octant_orientations, &freq1[address], norm1,
                                 dimensions[1].freq_offset);

            vm_double_multiply(scheme->octant_orientations,
                               &freq_ampA[step_vector_i + address],
                               &freq_ampB[step_vector_k + address], freq_amp);
            // Perform tenting on every sideband order over all orientations
            octahedronInterpolation2D(
                spec, dimensions[0].freq_offset, dimensions[1].freq_offset,
                scheme->integration_density, freq_amp, 1, dimensions[0].count,
                dimensions[1].count, iso_intrp);
            // step_vector += scheme->octant_orientations;
          }
        }
      }
    }
  }
}
