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

void one_dimensional_averaging(MRS_dimension *dimensions, MRS_averaging_scheme *scheme,
                               double *spec, unsigned int iso_intrp,
                               complex128 *exp_I_phase) {
  unsigned int i, j, k1, address, ptr, gamma_idx;
  unsigned int nt = scheme->integration_density, npts = scheme->octant_orientations;

  double offset_0, offset;
  double *freq, *phase_ptr, *amps = dimensions->freq_amplitude;
  double *amps_real = scheme->amps_real, *amps_imag = scheme->amps_imag;

  bool delta_interpolation = false;
  MRS_plan *planA = dimensions->events->plan;

  // get amplitudes for the interpolation
  switch (planA->number_of_sidebands) {
  case 1:
    for (j = 0; j < planA->n_octants; j++) {
      cblas_dcopy(npts, planA->norm_amplitudes, 1, &amps[j * npts], 1);
    }
    break;
  default:
    /**
     * Scale the absolute value square with the powder scheme weights. Only
     * the real part is scaled and the imaginary part is left as is. */
    for (j = 0; j < npts; j++) {
      cblas_dscal(planA->n_octants * planA->number_of_sidebands,
                  planA->norm_amplitudes[j], &amps[j], npts);
    }
  }

  offset_0 = dimensions->normalize_offset + dimensions->R0_offset;

  // gamma averaging
  for (gamma_idx = 0; gamma_idx < scheme->n_gamma; gamma_idx++) {
    ptr = scheme->total_orientations * gamma_idx;
    freq = &dimensions->local_frequency[ptr];
    phase_ptr = &(((double *)exp_I_phase)[2 * ptr]);

    if (fabs(*freq - freq[nt]) < TOL && fabs(*freq - freq[npts - 1]) < TOL)
      delta_interpolation = true;

    if (delta_interpolation) {
      offset_0 += *freq;
      for (i = 0; i < planA->number_of_sidebands; i++) {
        offset = offset_0 + planA->vr_freq[i];
        if ((int)offset >= 0 && (int)offset <= dimensions->count) {
          k1 = i * scheme->total_orientations;
          address = 0;
          j = 0;
          // multiply phase to the amplitudes.
          vm_double_multiply(scheme->total_orientations, phase_ptr, 2, &amps[k1],
                             amps_real);
          vm_double_multiply(scheme->total_orientations, phase_ptr + 1, 2, &amps[k1],
                             amps_imag);
          while (j++ < planA->n_octants) {
            octahedronDeltaInterpolation(nt, &offset, &amps_real[address], 1,
                                         dimensions->count, spec, iso_intrp);
            octahedronDeltaInterpolation(nt, &offset, &amps_imag[address], 1,
                                         dimensions->count, spec + 1, iso_intrp);
            address += npts;
          }
        }
      }
    }

    else {
      for (i = 0; i < planA->number_of_sidebands; i++) {
        offset = offset_0 + planA->vr_freq[i];
        if ((int)offset >= 0 && (int)offset <= dimensions->count) {
          k1 = i * scheme->total_orientations;
          address = 0;
          // multiply phase to the amplitudes.
          vm_double_multiply(scheme->total_orientations, phase_ptr, 2, &amps[k1],
                             amps_real);
          vm_double_multiply(scheme->total_orientations, phase_ptr + 1, 2, &amps[k1],
                             amps_imag);
          for (j = 0; j < planA->n_octants; j++) {
            // Add offset(isotropic + sideband_order) to the local frequencies.
            vm_double_add_offset(npts, &freq[address], offset, dimensions->freq_offset);
            // Perform tenting on every sideband order over all orientations.
            octahedronInterpolation(spec, dimensions->freq_offset, nt,
                                    &amps_real[address], 1, dimensions->count);
            octahedronInterpolation(spec + 1, dimensions->freq_offset, nt,
                                    &amps_imag[address], 1, dimensions->count);
            address += npts;
          }
        }
      }
    }
  }
}

void two_dimensional_averaging(MRS_dimension *dimensions, MRS_averaging_scheme *scheme,
                               double *spec, double *affine_matrix,
                               unsigned int iso_intrp, complex128 *exp_I_phase) {
  unsigned int i, k, j, index, step_vector_i, step_vector_k, address, gamma_idx;
  unsigned int npts = scheme->octant_orientations, ptr;

  MRS_plan *planA, *planB, *avg_plan;
  double *freq_ampA, *freq_ampB, *freq_amp = scheme->scrach, *avg_freq;
  double *ampsA_real = scheme->amps_real, *ampsA_imag = scheme->amps_imag;
  double offset0, offset1, offsetA, offsetB;
  double *freq0, *freq1, *phase_ptr;
  double norm0, norm1;

  offset0 = dimensions[0].R0_offset;
  freq_ampA = dimensions[0].freq_amplitude;
  planA = dimensions[0].events[0].plan;

  offset1 = dimensions[1].R0_offset;
  freq_ampB = dimensions[1].freq_amplitude;
  planB = dimensions[1].events[0].plan;

  // printf("%f %f\n", freq_ampA[0], freq_ampB[0]);
  index = planA->number_of_sidebands > 1;
  avg_plan = (index) ? planA : planB;
  avg_freq = (index) ? freq_ampA : freq_ampB;

  for (j = 0; j < npts; j++) {
    cblas_dscal(avg_plan->n_octants * avg_plan->number_of_sidebands,
                avg_plan->norm_amplitudes[j], &avg_freq[j], npts);
  }

  // gamma averaging
  for (gamma_idx = 0; gamma_idx < scheme->n_gamma; gamma_idx++) {
    ptr = scheme->total_orientations * gamma_idx;
    // printf("gma_idx %d ", gamma_idx);
    freq0 = &dimensions[0].local_frequency[ptr];
    freq1 = &dimensions[1].local_frequency[ptr];
    phase_ptr = &(((double *)exp_I_phase)[2 * ptr]);

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

    for (i = 0; i < planA->number_of_sidebands; i++) {
      offsetA = offset0 + planA->vr_freq[i];
      for (k = 0; k < planB->number_of_sidebands; k++) {
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

        // printf("f0 %f %f, (%f, %f), %f, %f\n", freq0[0] + norm0, freq1[0] + norm1,
        //  norm0, norm1, freq_ampA[0], freq_ampB[0]);
        if ((int)norm0 >= 0 && (int)norm0 <= dimensions[0].count) {
          step_vector_i = i * scheme->total_orientations;
          if ((int)norm1 >= 0 && (int)norm1 <= dimensions[1].count) {
            step_vector_k = k * scheme->total_orientations;

            // multiply phase to the amplitudes.
            vm_double_multiply(scheme->total_orientations, phase_ptr, 2,
                               &freq_ampA[step_vector_i], ampsA_real);
            vm_double_multiply(scheme->total_orientations, phase_ptr + 1, 2,
                               &freq_ampA[step_vector_i], ampsA_imag);

            // step_vector = 0;
            for (j = 0; j < planA->n_octants; j++) {
              address = j * npts;
              // Add offset(isotropic + sideband_order) to the local frequency
              // from [n to n+octant_orientation]
              vm_double_add_offset(npts, &freq0[address], norm0,
                                   dimensions[0].freq_offset);
              vm_double_add_offset(npts, &freq1[address], norm1,
                                   dimensions[1].freq_offset);

              vm_double_multiply(npts, &ampsA_real[address], 1,
                                 &freq_ampB[step_vector_k + address], freq_amp);
              // Perform tenting on every sideband order over all orientations
              octahedronInterpolation2D(
                  spec, dimensions[0].freq_offset, dimensions[1].freq_offset,
                  scheme->integration_density, freq_amp, 1, dimensions[0].count,
                  dimensions[1].count, iso_intrp);

              vm_double_multiply(npts, &ampsA_imag[address], 1,
                                 &freq_ampB[step_vector_k + address], freq_amp);
              // Perform tenting on every sideband order over all orientations
              octahedronInterpolation2D(
                  spec + 1, dimensions[0].freq_offset, dimensions[1].freq_offset,
                  scheme->integration_density, freq_amp, 1, dimensions[0].count,
                  dimensions[1].count, iso_intrp);
            }
          }
        }
      }
    }
  }
}
