//
//  mrsimulator.h
//
//  Created by Deepansh J. Srivastava, Jun 9, 2019
//  Copyright © 2019 Deepansh J. Srivastava. All rights reserved.
//  Contact email = srivastava.89@osu.edu, deepansh2012@gmail.com
//

#include "mrsimulator.h"

/* free the buffer and pre-calculated tables from the mrsimulator plan.      */
void mrsimulator_free_plan(mrsimulator_plan *the_plan) {
  fftw_destroy_plan(the_plan->the_fftw_plan);
  fftw_free(the_plan->vector);
  //   DftiFreeDescriptor(&plan);
  free_double(the_plan->vr_freq);
  free_double(the_plan->rotor_lab_2);
  free_double(the_plan->rotor_lab_4);
  free_double(the_plan->freq_offset);
  free_double(the_plan->amplitudes);
  free_double(the_plan->local_frequency);
  free_double(the_plan->wigner_2j_matrices);
  free_double(the_plan->wigner_4j_matrices);
  free_double_complex(the_plan->w2);
  free_double_complex(the_plan->w4);
  free_double_complex(the_plan->pre_phase_2);
  free_double_complex(the_plan->pre_phase_4);
}

/* Create a new mrsimulator plan. */
mrsimulator_plan *
mrsimulator_create_plan(unsigned int geodesic_polyhedron_frequency,
                        int number_of_sidebands,
                        double sample_rotation_frequency, double rotor_angle,
                        double spectral_increment, bool ALLOW_FOURTH_RANK) {

  unsigned int nt = geodesic_polyhedron_frequency, j, i, size_2, size_4;
  unsigned int n_orientations = (nt + 1) * (nt + 2) / 2;

  mrsimulator_plan *plan = malloc(sizeof(mrsimulator_plan));
  plan->number_of_sidebands = number_of_sidebands;
  plan->n_orientations = n_orientations;
  plan->size = n_orientations * number_of_sidebands;
  plan->geodesic_polyhedron_frequency = nt;

  // setup the fftw routine
  plan->vector = (fftw_complex *)fftw_malloc(sizeof(fftw_complex) * plan->size);
  // gettimeofday(&fft_setup_time, NULL);
  plan->the_fftw_plan =
      fftw_plan_many_dft(1, &number_of_sidebands, n_orientations, plan->vector,
                         NULL, n_orientations, 1, plan->vector, NULL,
                         n_orientations, 1, FFTW_FORWARD, FFTW_ESTIMATE);

  // char *filename = "128_sidebands.wisdom";
  // int status = fftw_export_wisdom_to_filename(filename);
  // printf("file save status %i \n", status);

  // double complex *vector = malloc_double_complex(size);
  // DFTI_DESCRIPTOR_HANDLE plan;
  // MKL_LONG status;
  // MKL_LONG stride[2] = {0, n_orientations};
  // status = DftiCreateDescriptor(&plan, DFTI_DOUBLE, DFTI_COMPLEX, 1,
  //                               number_of_sidebands);
  // status = DftiSetValue(plan, DFTI_NUMBER_OF_TRANSFORMS, n_orientations);
  // status = DftiSetValue(plan, DFTI_INPUT_DISTANCE, 1);
  // status = DftiSetValue(plan, DFTI_OUTPUT_DISTANCE, 1);
  // status = DftiSetValue(plan, DFTI_INPUT_STRIDES, stride);
  // status = DftiSetValue(plan, DFTI_OUTPUT_STRIDES, stride);
  // DftiCommitDescriptor(plan);
  // gettimeofday(&fft_setup_time_end, NULL);
  // clock_time =
  //     (double)(fft_setup_time_end.tv_usec - fft_setup_time.tv_usec) /
  //     1000000. + (double)(fft_setup_time_end.tv_sec - fft_setup_time.tv_sec);
  // printf("fft time %f \n", clock_time);
  // fftw routine end

  // gettimeofday(&start1, NULL);

  // Sideband order frequency relative to fft output order.
  double spectral_increment_inverse = 1.0 / spectral_increment;
  plan->vr_freq = __get_frequency_in_FFT_order(number_of_sidebands,
                                               sample_rotation_frequency);
  cblas_dscal(number_of_sidebands, spectral_increment_inverse, plan->vr_freq,
              1);

  // orientation set-up
  double *cos_alpha = malloc_double(n_orientations);
  double *cos_beta = malloc_double(n_orientations);
  plan->amplitudes = malloc_double(n_orientations);
  __powder_averaging_setup(geodesic_polyhedron_frequency, cos_alpha, cos_beta,
                           plan->amplitudes, 1);
  // normalizing amplitudes
  double number_of_sideband_inverse = (1.0 / (double)number_of_sidebands);
  cblas_dscal(n_orientations, number_of_sideband_inverse, plan->amplitudes, 1);

  double *sin_alpha = malloc_double(n_orientations);
  for (i = 0; i < n_orientations; i++) {
    sin_alpha[i] = sqrt(1.0 - pow(cos_alpha[i], 2));
  }

  plan->exp_Im_alpha = malloc_double_complex(5 * n_orientations);
  int s_2 = 2 * n_orientations, s_3 = 3 * n_orientations;
  // s_5 = 5 * n_orientations, s_6 = 6 * n_orientations;
  int s_0 = 0 * n_orientations, s_1 = 1 * n_orientations;
  // s_7 = 7 * n_orientations, s_8 = 8 * n_orientations;
  // int s_4 = 4 * n_orientations;

  // cblas_zdscal(n_orientations, 0.0, &exp_Im_alpha[s_4], 1);

  cblas_dcopy(n_orientations, cos_alpha, 1, (double *)&plan->exp_Im_alpha[s_3],
              2);
  cblas_dcopy(n_orientations, sin_alpha, 1,
              (double *)&plan->exp_Im_alpha[s_3] + 1, 2);
  // vmzConj(n_orientations, &exp_Im_alpha[s_3], &exp_Im_alpha[s_5], VML_EP);

  vmzMul(n_orientations, &plan->exp_Im_alpha[s_3], &plan->exp_Im_alpha[s_3],
         &plan->exp_Im_alpha[s_2], VML_EP);
  // vmzConj(n_orientations, &exp_Im_alpha[s_2], &exp_Im_alpha[s_6], VML_EP);

  if (ALLOW_FOURTH_RANK) {
    vmzMul(n_orientations, &plan->exp_Im_alpha[s_2], &plan->exp_Im_alpha[s_3],
           &plan->exp_Im_alpha[s_1], VML_EP);
    // vmzConj(n_orientations, &exp_Im_alpha[s_1], &exp_Im_alpha[s_7], VML_EP);

    vmzMul(n_orientations, &plan->exp_Im_alpha[s_1], &plan->exp_Im_alpha[s_3],
           &plan->exp_Im_alpha[s_0], VML_EP);
    // vmzConj(n_orientations, &exp_Im_alpha[s_0], &exp_Im_alpha[s_8], VML_EP);
  }

  // Setup for second rank tensor.
  // buffer for storing calculated frequencies from wigner 2j.
  plan->w2 = malloc_double_complex(5 * n_orientations);

  // wigner 2j matrices for every orientation.
  plan->wigner_2j_matrices = malloc_double(25 * n_orientations);
  __wigner_d_matrix_cosine(2, n_orientations, cos_beta,
                           plan->wigner_2j_matrices);

  // wigner 2j dm0 vector for rotating rotor-frame to lab-frame.
  plan->rotor_lab_2 = malloc_double(5);
  __wigner_dm0_vector(2, rotor_angle, plan->rotor_lab_2);

  // phase for calculating side-band amplitudes from wigner 2j frequencies.
  size_4 = 9 * number_of_sidebands;
  double complex *pre_phase = malloc_double_complex(size_4);
  __get_pre_phase_components(number_of_sidebands, sample_rotation_frequency,
                             pre_phase);

  size_2 = 5 * number_of_sidebands;
  plan->pre_phase_2 = malloc_double_complex(size_2);

  cblas_zcopy(size_2, &pre_phase[2 * number_of_sidebands], 1, plan->pre_phase_2,
              1);
  j = 0;
  for (i = 0; i < 5; i++) {
    cblas_zdscal(number_of_sidebands, plan->rotor_lab_2[i],
                 &plan->pre_phase_2[j], 1);
    j += number_of_sidebands;
  }

  plan->w4 = NULL;
  plan->wigner_4j_matrices = NULL;
  plan->rotor_lab_4 = NULL;
  plan->pre_phase_4 = NULL;

  // Setup for fourth rank tensor
  if (ALLOW_FOURTH_RANK) {
    // buffer for storing calculated frequencies from wigner 4j.
    plan->w4 = malloc_double_complex(9 * n_orientations);

    // wigner 4j matrices for every orientation.
    plan->wigner_4j_matrices = malloc_double(81 * n_orientations);
    __wigner_d_matrix_cosine(4, n_orientations, cos_beta,
                             plan->wigner_4j_matrices);

    /* wigner 4j dm0 vector for rotating rotor-frame to lab-frame.           */
    plan->rotor_lab_4 = malloc_double(9);
    __wigner_dm0_vector(4, rotor_angle, plan->rotor_lab_4);

    /* phase for calculating side-band amplitudes from wigner 2j frequencies.*/
    plan->pre_phase_4 = pre_phase;
    j = 0;
    for (i = 0; i < 9; i++) {
      cblas_zdscal(number_of_sidebands, plan->rotor_lab_4[i],
                   &plan->pre_phase_4[j], 1);
      j += number_of_sidebands;
    }
  } else {
    free_double_complex(pre_phase);
  }
  free_double(cos_beta);
  free_double(cos_alpha);
  free_double(sin_alpha);

  /* buffer to hold the local frequencies and frequency offset. The buffer   *
   * is useful when the rotor angle is off magic angle (54.735 deg).         */
  plan->local_frequency = malloc_double(n_orientations);
  plan->freq_offset = malloc_double(n_orientations);

  return plan;
}

// void mrsimulator_simulate(mrsimulator_plan *plan, double *R0,
//                           double complex *R2, double complex *R4) {
//   rotate
// }

/* pre-calculating the phase step exponents. ------------------------------- *
 * ------------------------------------------------------------------------- *
 *   phi = exp(sum_m v_cs_m * I 2pi [(exp(I m wr t) - 1)/(I m wr)])          *
 *   pre_phase(m, t) = I 2 pi [(exp(I m wr t) - 1)/(I m wr)]                 *
 *                   = (2 pi / m wr) (exp(I m wr t) - 1)                     *
 *                     ----scale----                                         *
 *                   = scale * (exp(I m wr t) - 1)                           *
 * ------------------------------------------------------------------------- */
void __get_pre_phase_components(int number_of_sidebands,
                                double sample_rotation_frequency,
                                double complex *pre_phase) {
  double spin_angular_freq, tau, wrt, pht, scale;
  int step, i, m;

  // Calculate the spin angular frequency
  spin_angular_freq = sample_rotation_frequency * PI2;

  // Calculate tau increments, where tau = (rotor period / number of phase
  // steps)
  tau = 1.0 / ((double)number_of_sidebands * sample_rotation_frequency);

  // pre-calculate the m omega spinning frequencies
  double m_wr[9] = {-4., -3., -2., -1., 0., 1., 2., 3., 4.};
  cblas_dscal(9, spin_angular_freq, m_wr, 1);

  i = 0;
  for (m = 0; m <= 8; m++) {
    if (m != 4) {
      wrt = m_wr[m] * tau;
      pht = 0.0;
      scale = PI2 / m_wr[m];
      for (step = 0; step < number_of_sidebands; step++) {
        pre_phase[i++] = scale * (cexp(I * pht) - 1.0);
        pht += wrt;
      }
    } else {
      for (step = 0; step < number_of_sidebands; step++) {
        pre_phase[i++] = 0.0;
      }
    }
  }
}
