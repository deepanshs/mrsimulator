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
 * a11 complex array of size npts, sideband amplitude event 11
 * a21 complex array of size npts, sideband amplitude event 21
 * res complex array of size npts, result  sideband amplitude.
 *
 * Evaluates: sum n2p a11(n1) * a11.conj()(n1') * a21(n2) * a21.conj()(n2')
 */
static inline void sideband_amplitude(int npts, int n_octant, complex128 *a11,
                                      complex128 *a21, int n1, int n2, int n1_sidebands,
                                      int n2_sidebands, complex128 *res,
                                      complex128 *res_t) {
  // array is packed as (sidebands, 1, orientations[npts])
  int i, n2p, n1p, d_npts = 2 * npts, n_min, n_max, size = npts * n_octant;
  int n1_ = n1 * size, n2_ = n2 * size, n1p_idx;
  double *fft1_index, *fft2_index;

  vm_double_zeros(2 * npts, (double *)res);
  fft1_index = get_FFT_order_freq(n1_sidebands, 1.0);
  fft2_index = get_FFT_order_freq(n2_sidebands, 1.0);

  n_min = -(int)(n1_sidebands / 2);
  n_max = (int)((n1_sidebands - 1) / 2);
  for (i = 0; i < n2_sidebands; i++) {
    n2p = (int)fft2_index[i];
    n1p = (int)fft1_index[n1] - (n2p - (int)fft2_index[n2]);
    if (n1p >= n_min && n1p <= n_max) {
      n1p_idx = (n1p >= 0) ? n1p : n1_sidebands + n1p;
      // printf(
      //     "sideband amps add n1=%i, n2=%i, n1p=%i, n2p=%i, n1d=%i, n2d=%i, n1pd=%i, "
      //     "n2pd=%i\n",
      //     (int)fft1_index[n1], (int)fft2_index[n2], n1p, n2p, n1, n2, n1p_idx, i);

      vm_double_complex_conj_multiply(npts, &a11[n1_], &a11[n1p_idx * size], res_t);
      vm_double_complex_conj_multiply(npts, res_t, &a21[i * size], res_t);
      vm_double_complex_multiply(npts, res_t, &a21[n2_], res_t);
      vm_double_add_inplace(d_npts, (double *)res_t, (double *)res);
    }
  }
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
            one_d_averaging(spec, npts, dimensions->freq_offset, &amps_real[address],
                            &amps_imag[address], dimensions->count,
                            scheme->position_size, scheme->positions, nt, user_defined,
                            interpolation);
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

  // gamma averaging
  for (gamma_idx = 0; gamma_idx < scheme->n_gamma; gamma_idx++) {
    ptr = scheme->total_orientations * gamma_idx;
    // printf("gma_idx %d ", gamma_idx);
    freq0 = &dimensions[0].local_frequency[ptr];
    freq1 = &dimensions[1].local_frequency[ptr];
    phase_ptr = &exp_I_phase[ptr];

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
                                 res_t);

              vm_double_complex_multiply(npts, &phase_ptr[address],
                                         (complex128 *)freq_amp,
                                         (complex128 *)freq_amp);

              vm_double_multiply_inplace(npts, planA->norm_amplitudes, 1, freq_amp, 2);
              // real part
              two_d_averaging(spec, npts, dimensions[0].freq_offset,
                              dimensions[1].freq_offset, freq_amp, 2,
                              scheme->position_size, scheme->positions,
                              dimensions[0].count, dimensions[1].count, iso_intrp,
                              scheme->integration_density, user_defined, interpolation);

              // imaginary part
              vm_double_multiply_inplace(npts, planA->norm_amplitudes, 1, freq_amp + 1,
                                         2);
              two_d_averaging(spec + 1, npts, dimensions[0].freq_offset,
                              dimensions[1].freq_offset, freq_amp + 1, 2,
                              scheme->position_size, scheme->positions,
                              dimensions[0].count, dimensions[1].count, iso_intrp,
                              scheme->integration_density, user_defined, interpolation);
            }
          }
        }
      }
    }
  }
  free(res_t);
}
