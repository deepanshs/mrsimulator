// -*- coding: utf-8 -*-
//
//  dimensional_averaging.c
//
//  @copyright Deepansh J. Srivastava, 2019-2021.
//  Created by Deepansh J. Srivastava, Mar 2, 2021.
//  Contact email = srivastava.89@osu.edu
//

#include "frequency_averaging.h"

/**
 * npts Number of array points
 * n_octant Number of octant for powder averaging.
 * a11 complex array of size npts, sideband amplitude for event 11
 * a21 complex array of size npts, sideband amplitude for event 21
 * n1 array index along dimenion 1 where sideband amplitude is evaluated
 * n2 array index along dimenion 2 where sideband amplitude is evaluated
 * n1_sidebands Number of sidebands along dimension 1
 * n2_sidebands Number of sidebands along dimension 2
 * res complex array of size npts, stored result sideband amplitude for (n1, n2)
 * res_t temp vector to store intermediate complex array
 * fft1_index stores the fft sideband order index relative to array index for dim1
 * fft2_index stores the fft sideband order index relative to array index for dim2
 * n_min minimum sideband order for dim 1
 * n_max maximum sideband order for dim 1
 *
 * Evaluates: sum n2p a11(n1) * a11.conj()(n1') * a21(n2) * a21.conj()(n2')
 */
static inline void sideband_amplitude(int npts, int n_octant, complex128 *a11,
                                      complex128 *a21, int n1p, int n2,
                                      int n1_sidebands, int n2_sidebands,
                                      complex128 *res, complex128 *res_t,
                                      int *fft1_index, int *fft2_index, int n_min,
                                      int n_max) {
  // array a11 and a21 are packed as (sidebands, total_orientation) with total
  // orientation as the leading dimension.
  int n1, n2p, total_orientation = npts * n_octant;
  int n1p_ = n1p * total_orientation, n2_ = n2 * total_orientation, n1_idx, n12f;

  // zero the result array
  cblas_dscal(2 * npts, 0.0, (double *)res, 1);

  // calculate n1 = n1p + (n2p - n2) using fft_index (sideband order) for corresponding
  // dimension
  n12f = fft1_index[n1p] - fft2_index[n2];

  // compute sum_n2p (a11[n1] * conj(a21[n2p]))
  for (n2p = 0; n2p < n2_sidebands; n2p++) {
    n1 = n12f + fft2_index[n2p];
    if (n1 >= n_min && n1 <= n_max) {  // check if within dim1 sideband order
      // convert sideband order to array index
      n1_idx = (n1 >= 0) ? n1 : n1_sidebands + n1;
      vm_double_complex_conj_multiply_inplace(npts, &a11[n1_idx * total_orientation],
                                              &a21[n2p * total_orientation], res);
    }
  }
  // calculate conj(a11[n1p]) * a21[n2]
  vm_double_complex_conj_multiply(npts, &a21[n2_], &a11[n1p_], res_t);

  // final product
  vm_double_complex_multiply(npts, res_t, res, res);
  // printf("amps[0]=%.6e %.6e\n", ((double *)res)[0], ((double *)res)[1]);
}

void one_dimensional_averaging(MRS_dimension *dimensions, MRS_averaging_scheme *scheme,
                               double *spec, unsigned int iso_intrp,
                               complex128 *exp_I_phase) {
  unsigned int i, j, k1, address, ptr, gamma_idx;
  unsigned int nt = scheme->integration_density, npts = scheme->octant_orientations;

  double offset_0, offset;
  double *freq, *phase_ptr, *amps = dimensions->freq_amplitude;
  double *amps_real = scheme->amps_real, *amps_imag = scheme->amps_imag;

  bool delta_interpolation = false;
  bool user_defined = scheme->user_defined, interpolation = scheme->interpolation;
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

    if (absd(*freq - freq[nt]) < TOL && absd(*freq - freq[npts - 1]) < TOL)
      if (interpolation) delta_interpolation = true;

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
          if (scheme->is_complex) {
            vm_double_multiply(scheme->total_orientations, phase_ptr + 1, 2, &amps[k1],
                               amps_imag);
          }
          while (j++ < planA->n_octants) {
            octahedronDeltaInterpolation(nt, &offset, &amps_real[address], 1,
                                         dimensions->count, spec, iso_intrp);
            if (scheme->is_complex) {
              octahedronDeltaInterpolation(nt, &offset, &amps_imag[address], 1,
                                           dimensions->count, spec + 1, iso_intrp);
            }
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
          if (scheme->is_complex) {
            vm_double_multiply(scheme->total_orientations, phase_ptr + 1, 2, &amps[k1],
                               amps_imag);
          }
          for (j = 0; j < planA->n_octants; j++) {
            // Add offset(isotropic + sideband_order) to the local frequencies.
            vm_double_add_offset(npts, &freq[address], offset, dimensions->freq_offset);
            // Perform tenting on every sideband order over all orientations.
            one_d_averaging(spec, npts, dimensions->freq_offset, &amps_real[address],
                            &amps_imag[address], dimensions->count,
                            scheme->position_size, scheme->positions, nt, user_defined,
                            interpolation, scheme->is_complex);
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
  unsigned int i, k, j, address, gamma_idx;
  unsigned int npts = scheme->octant_orientations, ptr;
  int *fft1_index, *fft2_index, n_min, n_max;

  MRS_plan *planA, *planB;
  complex128 *freq_ampA, *freq_ampB, *phase_ptr;
  double *freq_amp = scheme->scratch;
  double offset0, offset1, offsetA, offsetB;
  double *freq0, *freq1;
  double norm0, norm1;
  complex128 *res_t = malloc_complex128(npts);
  bool user_defined = scheme->user_defined, interpolation = scheme->interpolation;

  offset0 = dimensions[0].R0_offset;
  freq_ampA = dimensions[0].events[0].event_freq_amplitude;
  planA = dimensions[0].events[0].plan;

  offset1 = dimensions[1].R0_offset;
  freq_ampB = dimensions[1].events[0].event_freq_amplitude;
  planB = dimensions[1].events[0].plan;

  fft1_index = get_FFT_index(planA->number_of_sidebands);
  fft2_index = get_FFT_index(planB->number_of_sidebands);

  n_min = -(int)(planA->number_of_sidebands / 2);
  n_max = (int)((planA->number_of_sidebands - 1) / 2);

  // gamma averaging
  for (gamma_idx = 0; gamma_idx < scheme->n_gamma; gamma_idx++) {
    ptr = scheme->total_orientations * gamma_idx;
    // printf("gma_idx %d ", gamma_idx);
    freq0 = &dimensions[0].local_frequency[ptr];
    freq1 = &dimensions[1].local_frequency[ptr];
    phase_ptr = &exp_I_phase[ptr];

    for (j = 0; j < planA->n_octants; j++) {
      vm_double_multiply_inplace(npts, planA->norm_amplitudes, 1,
                                 (double *)(&phase_ptr[j * npts]), 2);
      vm_double_multiply_inplace(npts, planA->norm_amplitudes, 1,
                                 (double *)(&phase_ptr[j * npts]) + 1, 2);
    }

    // scale and shear the first dimension.
    // matrix = [[a, b], [c/a, d - bc/a]]
    // f0' = a f0 + b * f1
    if (affine_matrix[0] != 1) {
      cblas_dscal(scheme->total_orientations, affine_matrix[0], freq0, 1);
    }
    if (affine_matrix[1] != 0) {
      cblas_daxpy(scheme->total_orientations, affine_matrix[1], freq1, 1, freq0, 1);
    }

    // scale and shear the second dimension.
    // f1' = (d - bc/a) f1 + c/a f0'
    //     = d f1 - bc/a f1 + c f0 + bc/a f1
    //     = c f0 + d f1
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
          if ((int)norm1 >= 0 && (int)norm1 <= dimensions[1].count) {
            for (j = 0; j < planA->n_octants; j++) {
              address = j * npts;
              // Add offset(isotropic + sideband_order) to the local frequency
              // from [n to n+octant_orientation]
              vm_double_add_offset(npts, &freq0[address], norm0,
                                   dimensions[0].freq_offset);
              vm_double_add_offset(npts, &freq1[address], norm1,
                                   dimensions[1].freq_offset);

              sideband_amplitude(npts, planA->n_octants, &freq_ampA[address],
                                 &freq_ampB[address], i, k, planA->number_of_sidebands,
                                 planB->number_of_sidebands, (complex128 *)freq_amp,
                                 res_t, fft1_index, fft2_index, n_min, n_max);

              vm_double_complex_multiply(npts, &phase_ptr[address],
                                         (complex128 *)freq_amp,
                                         (complex128 *)freq_amp);

              // real part
              two_d_averaging(spec, npts, dimensions[0].freq_offset,
                              dimensions[1].freq_offset, freq_amp, 2,
                              scheme->position_size, scheme->positions,
                              dimensions[0].count, dimensions[1].count, iso_intrp,
                              scheme->integration_density, user_defined, interpolation);

              if (scheme->is_complex) {
                // imaginary part
                two_d_averaging(spec + 1, npts, dimensions[0].freq_offset,
                                dimensions[1].freq_offset, freq_amp + 1, 2,
                                scheme->position_size, scheme->positions,
                                dimensions[0].count, dimensions[1].count, iso_intrp,
                                scheme->integration_density, user_defined,
                                interpolation);
              }
            }
          }
        }
      }
    }
  }
  free(res_t);
}
