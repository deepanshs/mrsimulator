// -*- coding: utf-8 -*-
//
//  mrsimulator.c
//
//  @copyright Deepansh J. Srivastava, 2019-2020.
//  Created by Deepansh J. Srivastava, Jun 9, 2019
//  Contact email = deepansh2012@gmail.com
//

#include "mrsimulator.h"

/* free the buffer and pre-calculated tables from the mrsimulator plan. */
void MRS_free_plan(MRS_plan *the_plan) {
  fftw_destroy_plan(the_plan->the_fftw_plan);
  fftw_free(the_plan->vector);
  free(the_plan->vr_freq);
  free(the_plan->wigner_d2m0_vector);
  free(the_plan->wigner_d4m0_vector);
  // free(the_plan->freq_offset);
  free(the_plan->norm_amplitudes);
  // free(the_plan->local_frequency);
  // free(the_plan->wigner_2j_matrices);
  // free(the_plan->wigner_4j_matrices);
  // free(the_plan->w2);
  // free(the_plan->w4);
  free(the_plan->pre_phase_2);
  free(the_plan->pre_phase_4);
  MRS_free_averaging_scheme(the_plan->averaging_scheme);
}

/* Create a new mrsimulator plan.
 * A plan for mrsimulator contains buffers and tabulated values to produce
 * faster simulation. The plan includes,
 * 1) calculating an array of orientations over the surface of a sphere. Each
 * orientation is described by an azimuthal angle, (α), a polar angle, (β),
 * and a weighting factor describing the spherical average.
 * 2) calculating wigner-2j(β) and wigner-4j(β) matrices at every orientation
 * angle β,
 * 3) pre-calculating the exponent of the sideband order phase,
 * exp(-Imα), at every orientation angle α,
 * 4) creating the fftw plan,
 * 4) allocating buffer for storing the evaluated frequencies and their
 * respective amplitudes.
 */

MRS_plan *MRS_create_plan(MRS_averaging_scheme *scheme, int number_of_sidebands,
                          double sample_rotation_frequency_in_Hz,
                          double rotor_angle_in_rad, double increment,
                          bool allow_fourth_rank) {
  unsigned int j, i, size_2, size_4;

  MRS_plan *plan = malloc(sizeof(MRS_plan));
  plan->number_of_sidebands = number_of_sidebands;
  plan->sample_rotation_frequency_in_Hz = sample_rotation_frequency_in_Hz;

  plan->allow_fourth_rank = allow_fourth_rank;
  plan->one[0] = 1.0;
  plan->one[1] = 0.0;
  plan->zero[0] = 0.0;
  plan->zero[1] = 0.0;

  /**
   * Update the mrsimulator plan with the given spherical averaging scheme. We
   * create the coordinates on the surface of the unit sphere by projecting the
   * points on the face of the octahedron to a unit sphere.
   * Usually, before updating the averaging scheme, the memory allocated by the
   * previous scheme must be freed. Since, we are creating the scheme for this
   * plan for the very first time, there is no need to call
   * MRS_free_averaging_plan() method.
   */
  // plan->averaging_scheme = MRS_create_averaging_scheme(
  //     integration_density, allow_fourth_rank);
  plan->averaging_scheme = scheme;
  plan->n_octants = 1;
  if (scheme->integration_volume == 1) plan->n_octants = 4;
  if (scheme->integration_volume == 2) plan->n_octants = 8;

  /* Normalizing amplitudes from the spherical averaging scheme by the number
   * of sidebands */
  plan->norm_amplitudes = malloc_double(scheme->octant_orientations);
  cblas_dcopy(scheme->octant_orientations, scheme->amplitudes, 1,
              plan->norm_amplitudes, 1);
  double scale = (1.0 / (double)(plan->number_of_sidebands *
                                 plan->number_of_sidebands * plan->n_octants));
  cblas_dscal(scheme->octant_orientations, scale, plan->norm_amplitudes, 1);
  /* ----------------------------------------------------------------------- */
  /* fftw routine setup .................................................... */
  /* ....................................................................... */
  plan->size = scheme->total_orientations * plan->number_of_sidebands;
  plan->vector = (fftw_complex *)fftw_malloc(sizeof(fftw_complex) * plan->size);
  // malloc_complex128(plan->size);
  // gettimeofday(&fft_setup_time, NULL);

  int fftw_thread = fftw_init_threads();
  if (fftw_thread == 0) {
    printf("failed to initialize fftw threading");
  }
  fftw_plan_with_nthreads(4);

  plan->the_fftw_plan = fftw_plan_many_dft(
      1, &plan->number_of_sidebands, scheme->total_orientations, plan->vector,
      NULL, scheme->total_orientations, 1, plan->vector, NULL,
      scheme->total_orientations, 1, FFTW_FORWARD, FFTW_ESTIMATE);

  // char *filename = "128_sidebands.wisdom";
  // int status = fftw_export_wisdom_to_filename(filename);
  // printf("file save status %i \n", status);
  /* ----------------------------------------------------------------------- */

  /**
   * Update the mrsimulator plan with the given rotor angle in radian.
   * This method updates the wigner d^l_{m,0}(rotor_angle_in_rad) vectors used
   * in tranforming the l-rank tensors from the rotor frame to lab frame. Here l
   * is either 2 or 4.
   */
  MRS_plan_update_rotor_angle_in_rad(plan, rotor_angle_in_rad,
                                     allow_fourth_rank);

  /* sideband order frequency in fft output order. */
  double increment_inverse = 1.0 / increment;
  plan->vr_freq = __get_frequency_in_FFT_order(number_of_sidebands,
                                               sample_rotation_frequency_in_Hz);
  cblas_dscal(number_of_sidebands, increment_inverse, plan->vr_freq, 1);

  /**
   * calculating the sideband phase multiplier.
   *    pre_phase(m,t) =  I 2π [(exp(I m wr t) - 1)/(I m wr)].
   * @see __get_components()
   */
  size_4 = 9 * number_of_sidebands;
  complex128 *pre_phase = malloc_complex128(size_4);
  __get_components(number_of_sidebands, sample_rotation_frequency_in_Hz,
                   (double *)pre_phase);

  size_2 = 5 * number_of_sidebands;
  plan->pre_phase_2 = malloc_complex128(size_2);

  cblas_zcopy(size_2, (double *)(pre_phase[2 * number_of_sidebands]), 1,
              (double *)(plan->pre_phase_2), 1);
  /* multiplying the wigner 2j d^2_{m,0}(rotor_angle_in_rad) vector to the
   * sideband phase multiplier. This multiplication absorbs the rotation of the
   * second rank tensors from the rotor frame to the lab frame, thereby,
   * reducing the number of calculations involved per site. This step assumes
   * that the Euler angles invloved in the rotation of the 2nd rank tensors to
   * the lab frame is (0, rotor_angle_in_rad, 0).
   */

  j = 0;
  for (i = 0; i < 5; i++) {
    cblas_zdscal(number_of_sidebands, plan->wigner_d2m0_vector[i],
                 (double *)(plan->pre_phase_2[j]), 1);
    j += number_of_sidebands;
  }

  plan->pre_phase_4 = NULL;

  /* Setup for processing the fourth rank tensors. */
  if (allow_fourth_rank) {
    /* multiplying the wigner 4j d^4_{m,0} vector to the sideband phase
     * multiplier. This multiplication absorbs the rotation of the fourth rank
     * tensors from the rotor frame to the lab frame, thereby, reducing the
     * number of calculations involved per site. This step assumes that the
     * Euler angles involved in the rotation of the 4th rank tensors to the lab
     * frame is (0, rotor_angle_in_rad, 0).
     */
    plan->pre_phase_4 = pre_phase;
    j = 0;
    for (i = 0; i < 9; i++) {
      cblas_zdscal(number_of_sidebands, plan->wigner_d4m0_vector[i],
                   (double *)(plan->pre_phase_4[j]), 1);
      j += number_of_sidebands;
    }
  } else {
    free(pre_phase);
  }
  return plan;
}

/* Free the memory from the mrsimulator plan associated with the wigner
 * d^l_{m,0}(rotor_angle_in_rad) vectors. Here, l=2 or 4.
 *  */
void MRS_plan_free_rotor_angle_in_rad(MRS_plan *plan) {
  free(plan->wigner_d2m0_vector);
  free(plan->wigner_d4m0_vector);
  plan->wigner_d2m0_vector = NULL;
  plan->wigner_d4m0_vector = NULL;
}

/* Update the MRS plan for the given rotor angle in radians. */
void MRS_plan_update_rotor_angle_in_rad(MRS_plan *plan,
                                        double rotor_angle_in_rad,
                                        bool allow_fourth_rank) {
  plan->rotor_angle_in_rad = rotor_angle_in_rad;
  /**
   *  calculating wigner 2j d^2_{m,0} vector where m ∈ [-2, 2]. This vector is
   * used to rotate the second tank tensors from the rotor frame to the lab
   * frame.
   * @see wigner_dm0_vector()
   */
  plan->wigner_d2m0_vector = malloc_double(5);
  wigner_dm0_vector(2, rotor_angle_in_rad, plan->wigner_d2m0_vector);

  plan->wigner_d4m0_vector = NULL;
  if (allow_fourth_rank) {
    /* calculating wigner 4j d^4_{m,0} vector where m ∈ [-4, 4]. This vector is
     * used to rotate the fourth tank tensors from the rotor frame to the lab
     * frame.*/
    plan->wigner_d4m0_vector = malloc_double(9);
    wigner_dm0_vector(4, rotor_angle_in_rad, plan->wigner_d4m0_vector);
  }
}

/* Return a copy of thr mrsimulator plan. */
MRS_plan *MRS_copy_plan(MRS_plan *plan) {}

/**
 * @func MRS_get_amplitudes_from_plan
 *
 * The function evaluates the amplitudes at every orientation and at every
 * sideband per orientation. This is done in two steps.
 * 1) Rotate R2 and R4, given in the crystal or common frame to w2 and w4 in
 *    the lab frame using wigner 2j and 4j rotation matrices, respectively,
 *    at all orientations.
 * 2) Evalute the sideband amplitudes using equation [39] of the reference
 *    https://doi.org/10.1006/jmre.1998.1427.
 */
void MRS_get_amplitudes_from_plan(MRS_plan *plan, complex128 *R2,
                                  complex128 *R4) {
  unsigned int orientation;
  MRS_averaging_scheme *scheme = plan->averaging_scheme;
  __batch_wigner_rotation(scheme->octant_orientations, plan->n_octants,
                          scheme->wigner_2j_matrices, R2,
                          scheme->wigner_4j_matrices, R4, scheme->exp_Im_alpha,
                          scheme->w2, scheme->w4);

  /* Evaluating the exponent of the sideband phase w.r.t the second rank
   * tensors. The exponent is given as,
   * w2(Θ)*d^2_{m,0}(rotor_angle_in_rad) * 2πI [(exp(I m ωr t) - 1)/(I m ωr)]
   * |----lab frame 2nd rank tensors---|
   *       |------------------------- pre_phase_2 ---------------------------|
   * where `pre_phase_2` is pre-calculated and stored in the plan. The result
   * is stored in the plan as a complex double array under the variable
   * `vector`, which is interpreted as a row major matrix of shape
   * `number_of_sidebands` x `total_orientations` with `total_orientations`
   * as the leading dimension. */
  cblas_zgemm(CblasRowMajor, CblasTrans, CblasTrans, plan->number_of_sidebands,
              scheme->total_orientations, 5, (double *)(plan->one),
              (double *)(plan->pre_phase_2), plan->number_of_sidebands,
              (double *)(scheme->w2), 5, (double *)(plan->zero),
              (double *)(plan->vector), scheme->total_orientations);

  if (scheme->w4 != NULL) {
    /* Evaluating the exponent of the sideband phase w.r.t the fourth rank
     * tensors. The exponent is given as,
     * w4(Θ)*d^4_{m, 0}(rotor_angle_in_rad) * 2πI[(exp(I m ωr t) - 1)/(I m ωr)]
     * |----lab frame 4th rank tensors----|
     *       |-------------------------- pre_phase_4 -------------------------|
     * where `pre_phase_4` is pre-calculated and stored in the plan.
     * This operation will update the values in variable `vector`.
     */
    cblas_zgemm(CblasRowMajor, CblasTrans, CblasTrans,
                plan->number_of_sidebands, scheme->total_orientations, 9,
                (double *)(plan->one), (double *)(plan->pre_phase_4),
                plan->number_of_sidebands, (double *)(scheme->w4), 9,
                (double *)(plan->one), (double *)(plan->vector),
                scheme->total_orientations);
  }
  if (plan->number_of_sidebands == 1) {
    cblas_dscal(plan->size, 2.0, (double *)(plan->vector), 2);
    vm_double_exp(plan->size, (double *)(plan->vector),
                  (double *)(plan->vector), 2);
  } else {
    /* Evaluating the sideband phase exp(vector) */
    vm_double_complex_exp(plan->size, plan->vector, plan->vector);

    /* Evaluate the Fourier transform of vector, fft(vector). The fft operation
     * updates the value of the array, `vector` */
    fftw_execute(plan->the_fftw_plan);

    /* Taking the absolute value square of the vector array. The absolute value
     * square is stores as the real part of the `vector` array. The imaginary
     * part is garbage. This method avoids creating new arrays. */

    // cblas_dsbmv(CblasRowMajor, CblasUpper, 2 * plan->size, 0, 1.0,
    //             (double *)plan->vector, 1, (double *)plan->vector, 1, 0.0,
    //             (double *)plan->vector, 1);

    vm_double_square_inplace(2 * plan->size, (double *)plan->vector);
    cblas_daxpy(plan->size, 1.0, (double *)plan->vector + 1, 2,
                (double *)plan->vector, 2);
  }
  /* Scaling the absolute value square with the powder scheme weights. Only
   * the real part is scaled and the imaginary part is left as is.
   */
  for (orientation = 0; orientation < scheme->octant_orientations;
       orientation++) {
    cblas_dscal(plan->n_octants * plan->number_of_sidebands,
                plan->norm_amplitudes[orientation],
                (double *)&plan->vector[orientation],
                2 * scheme->octant_orientations);
  }
}

/**
 * @func MRS_get_frequencies_from_plan
 *
 * The complex128 w2 and w4 vectors are accessed from the plan as plan->w2
 * and plan->w4
 */
void MRS_get_frequencies_from_plan(MRS_plan *plan, double R0) {
  MRS_averaging_scheme *scheme = plan->averaging_scheme;
  /* Setting isotropic frequency contribution from the zeroth rank tensor. */
  plan->isotropic_offset = R0;

  /* Calculating the local anisotropic frequency contributions from the      *
   * second rank tensor. */
  plan->buffer = plan->wigner_d2m0_vector[2];
  cblas_daxpby(scheme->total_orientations, plan->buffer,
               (double *)&scheme->w2[2], 10, 0.0, scheme->local_frequency, 1);
  if (plan->allow_fourth_rank) {
    /* Calculating the local anisotropic frequency contributions from the    *
     * fourth rank tensor. */
    plan->buffer = plan->wigner_d4m0_vector[4];
    cblas_daxpy(scheme->total_orientations, plan->buffer,
                (double *)&scheme->w4[4], 18, scheme->local_frequency, 1);
  }
}

void MRS_get_normalized_frequencies_from_plan(MRS_plan *plan,
                                              MRS_dimension *dim, double R0) {
  MRS_averaging_scheme *scheme = plan->averaging_scheme;
  /* Calculating the normalized isotropic frequency contribution from the    *
   * zeroth rank tensor. */
  plan->isotropic_offset = dim->normalize_offset + R0 * dim->inverse_increment;

  /* Calculating the normalized local anisotropic frequency contributions    *
   * from the second rank tensor. */
  plan->buffer = dim->inverse_increment * plan->wigner_d2m0_vector[2];
  cblas_daxpby(scheme->total_orientations, plan->buffer,
               (double *)&scheme->w2[2], 10, 0.0, scheme->local_frequency, 1);
  if (plan->allow_fourth_rank) {
    /* Calculating the normalized local anisotropic frequency contributions  *
     * from the fourth rank tensor. */
    plan->buffer = dim->inverse_increment * plan->wigner_d4m0_vector[4];
    cblas_daxpy(scheme->total_orientations, plan->buffer,
                (double *)&scheme->w4[4], 18, scheme->local_frequency, 1);
  }
}

/**
 *  The function calculates the following.
 *   pre_phase(m, t) = I 2π [(exp(I m ωr t) - 1)/(I m ωr)]
 *                   = (2π / m ωr) (exp(I m ωr t) - 1)
 *                     |--scale--|
 *                   = scale (exp(I m ωr t) - 1)
 *                   = scale [[cos(m ωr t) -1] +Isin(m ωr t)]
 * where ωr is the sample spinning frequency in Hz, m goes from -4 to 4, and
 * t is a vector of length `number_of_sidebands` given as
 *    t = [0, 1, ... number_of_sidebands-1]/(ωr*number_of_sidebands)
 * `pre_phase` is a matrix of shape, `9 x number_of_sidebands`.
 *
 * Also,
 *   pre_phase(-m, t) = (-2π / m ωr) (exp(-I m ωr t) - 1)
 *                    = -scale [[cos(m ωr t) -1] -Isin(m ωr t)]
 *                    = scale [-[cos(m ωr t) -1] +Isin(m ωr t)]
 * That is pre_phase[-m] = -Re(pre_phase[m]) + Im(pre_phase[m])
 */
void __get_components_2(int number_of_sidebands,
                        double sample_rotation_frequency_in_Hz,
                        complex128 *pre_phase) {
  int m, i;
  double spin_angular_freq, tau, scale;

  // memset(pre_phase, 0, 2 * number_of_sidebands * sizeof(double));

  double *input = malloc_double(number_of_sidebands);
  double *ones = malloc_double(number_of_sidebands);
  double *phase = malloc_double(number_of_sidebands);

  vm_double_ones(number_of_sidebands, ones);
  // for (i = 0; i < number_of_sidebands; i++) {
  //   printf("%f", ones[i]);
  // }
  vm_double_arrange(number_of_sidebands, input);

  // Calculate the spin angular frequency
  spin_angular_freq = sample_rotation_frequency_in_Hz * PI2;

  // Calculate tau, where tau = (rotor period / number of phase steps)
  tau = 1.0 / ((double)number_of_sidebands * sample_rotation_frequency_in_Hz);

  // pre-calculate the m omega spinning frequencies
  double m_wr[9] = {-4., -3., -2., -1., 0., 1., 2., 3., 4.};
  cblas_dscal(9, spin_angular_freq, m_wr, 1);

  for (m = 0; m <= 3; m++) {
    /**
     * evaluate pre_phase = scale * (cexp(I * phase) - 1.0)
     * where phase = m_wr[m] * tau * [0 .. number_of_sidebands-1]
     * and scale = 2π/m_wr[m].
     */
    i = m * number_of_sidebands;
    scale = PI2 / m_wr[m];

    // step 1. calculate phase
    vm_double_ramp(number_of_sidebands, input, m_wr[m] * tau, 0.0, phase);

    // step 2. evaluate cexp(I * phase) = cos(phase) + I sin(phase)
    vm_cosine_I_sine(number_of_sidebands, phase, &pre_phase[i]);

    // step 3. subtract 1.0 from pre_phase
    cblas_daxpy(number_of_sidebands, -1.0, ones, 1, (double *)(pre_phase[i]),
                2);

    // step 4. scale pre_phase with factor `scale`
    cblas_zdscal(number_of_sidebands, scale, (double *)(pre_phase[i]), 1);

    /* The expression pre_phase[m] = scale * (cexp(I * phase) - 1.0) given
     * above from m is related to -m as pre_phase[-m] = -Re(pre_phase[m]) +
     * Im(pre_phase[m])
     */

    cblas_zcopy(number_of_sidebands, (double *)(pre_phase[i]), 1,
                (double *)(pre_phase[(8 - m) * number_of_sidebands]), 1);
    cblas_dscal(number_of_sidebands, -1.0, (double *)(pre_phase[i]), 2);
  }
  vm_double_zeros(2 * number_of_sidebands,
                  (double *)(pre_phase[4 * number_of_sidebands]));
  // memset((double *)&pre_phase[4 * number_of_sidebands], 0,
  //        2 * number_of_sidebands * sizeof(double));

  free(input);
  free(phase);
  free(ones);
}

/**
 *  The function calculates the following.
 *   pre_phase(m, t) = I 2π [(exp(I m ωr t) - 1)/(I m ωr)]
 *                   = (2π / m ωr) (exp(I m ωr t) - 1)
 *                     |--scale--|
 *                   = scale * (exp(I m ωr t) - 1)
 * where ωr is the sample spinning frequency in Hz, m goes from -4 to 4, and
 * t is a vector of length `number_of_sidebands` given as
 *    t = [0, 1, ... number_of_sidebands-1]/(ωr*number_of_sidebands)
 * `pre_phase` is a matrix of shape, `9 x number_of_sidebands`.
 */
void __get_components(int number_of_sidebands, double sample_rotation_frequency,
                      double *restrict pre_phase) {
  double spin_angular_freq, tau, wrt, pht, scale;
  int step, m;
  // double *pre_phase_ = (double *)pre_phase;

  // Calculate the spin angular frequency
  spin_angular_freq = sample_rotation_frequency * PI2;

  // Calculate tau increments, where tau = (rotor period / number of phase
  // steps)
  tau = 1.0 / ((double)number_of_sidebands * sample_rotation_frequency);

  // pre-calculate the m omega spinning frequencies
  double m_wr[9] = {-4., -3., -2., -1., 0., 1., 2., 3., 4.};
  cblas_dscal(9, spin_angular_freq, m_wr, 1);

  for (m = 0; m <= 8; m++) {
    if (m != 4) {
      wrt = m_wr[m] * tau;
      pht = 0.0;
      scale = PI2 / m_wr[m];
      for (step = 0; step < number_of_sidebands; step++) {
        *pre_phase++ = scale * (cos(pht) - 1.0);
        *pre_phase++ = scale * sin(pht);
        pht += wrt;
      }
    } else {
      vm_double_zeros(2 * number_of_sidebands, &pre_phase[0]);
      pre_phase += 2 * number_of_sidebands;
    }
  }
}

MRS_dimension *MRS_create_dimension(int count, double coordinates_offset,
                                    double increment) {
  MRS_dimension *dim = malloc(sizeof(MRS_dimension));
  dim->count = count;
  dim->coordinates_offset = coordinates_offset;
  dim->increment = increment;

  dim->inverse_increment = 1.0 / increment;
  dim->normalize_offset = 0.5 - (coordinates_offset * dim->inverse_increment);
  return dim;
}
