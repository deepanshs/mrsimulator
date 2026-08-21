// -*- coding: utf-8 -*-
//
//  mrsimulator.c
//
//  @copyright Deepansh J. Srivastava, 2019-2021.
//  Created by Deepansh J. Srivastava, Jun 9, 2019.
//  Contact email = srivastava.89@osu.edu
//

#include "mrsimulator.h"

double ONE[] = {1.0, 0.0};
double ZERO[] = {0.0, 0.0};

/**
 * Free the buffers and pre-calculated tables from the mrsimulator plan.
 */
void MRS_free_plan(MRS_plan *the_plan) {
  if (DEBUG)
    printf("plan_copy %d, plan copy_for_rotor_angle %d plan copy_for_rotor_freq %d\n",
           the_plan->copy, the_plan->copy_for_rotor_angle,
           the_plan->copy_for_rotor_freq);
  if (!the_plan->copy) {
    MRS_free_plan_for_rotor_angle_copy(the_plan);
    MRS_free_plan_for_rotor_freq_copy(the_plan);
    free(the_plan->norm_amplitudes);
    if (DEBUG) printf("freed norm_amplitudes\n");
    free(the_plan);
    return;
  }
  if (the_plan->copy_for_rotor_angle) {
    MRS_free_plan_for_rotor_angle_copy(the_plan);
  }
  if (the_plan->copy_for_rotor_freq) {
    MRS_free_plan_for_rotor_freq_copy(the_plan);
  }
  free(the_plan);
}

/**
 * Free the memory from the mrsimulator plan associated with the wigner
 * d^l_{m,0}(rotor_angle_in_rad) vectors at rotor_angle and all powder
 * orientations. Here, l=2 or 4.
 */
void MRS_free_plan_for_rotor_angle_copy(MRS_plan *the_plan) {
  free(the_plan->wigner_d2m0_vector);
  free(the_plan->wigner_d4m0_vector);
  free(the_plan->pre_phase_2);
  free(the_plan->pre_phase_4);
  if (DEBUG)
    printf("freed wigner_d2m0_vector, wigner_d4m0_vector, pre_phase_2, pre_phase_4\n");
}

/**
 * Free the memory from the mrsimulator plan associated with the sideband freq.
 */
void MRS_free_plan_for_rotor_freq_copy(MRS_plan *the_plan) {
  free(the_plan->vr_freq);
  if (DEBUG) printf("freed vr_freq\n");
}

/**
 * Release temporary MRS plan storage
 */
void MRS_plan_release_temp_storage(MRS_plan *the_plan) {
  free(the_plan->pre_phase);
  if (DEBUG) printf("freed pre_phase\n");
}

/**
 * Create a new mrsimulator plan.
 *
 * A plan for mrsimulator contains buffers and tabulated values to produce faster
 * simulation. The plan includes,
 * 1) calculating an array of orientations over the surface of a sphere. Each
 *    orientation is described by an azimuthal angle, (α), a polar angle, (β), and a
 *    weighting factor describing the spherical average.
 * 2) calculating wigner-2j(β) and wigner-4j(β) matrices at every orientation angle β,
 * 3) pre-calculating the exponent of the sideband order phase, exp(-Imα), at every
 *    orientation angle α,
 * 4) creating the fftw plan, 4) allocating buffer for storing the evaluated frequencies
 *    and their respective amplitudes.
 */
MRS_plan *MRS_create_plan(MRS_averaging_scheme *scheme,
                          unsigned int number_of_sidebands,
                          double rotor_frequency_in_Hz, double rotor_angle_in_rad,
                          double increment, bool allow_4th_rank) {
  MRS_plan *plan = malloc(sizeof(MRS_plan));
  plan->number_of_sidebands = number_of_sidebands;
  plan->rotor_frequency_in_Hz = rotor_frequency_in_Hz;
  plan->rotor_angle_in_rad = rotor_angle_in_rad;

  plan->allow_4th_rank = allow_4th_rank;

  plan->copy = false;
  plan->copy_for_rotor_angle = false;
  plan->copy_for_rotor_freq = false;
  plan->is_static = (rotor_frequency_in_Hz < 1.0e-3) ? true : false;
  /**
   * Update the mrsimulator plan with the given spherical averaging scheme. We create
   * the coordinates on the surface of the unit sphere by projecting the points on the
   * face of the octahedron to a unit sphere. Usually, before updating the averaging
   * scheme, the memory allocated by the previous scheme must be freed. Since, we are
   * creating the scheme for this plan for the very first time, there is no need to call
   * MRS_free_averaging_plan() method.
   */

  plan->n_octants = 1;
  if (scheme->integration_volume == 1) plan->n_octants = 4;
  if (scheme->integration_volume == 2) plan->n_octants = 8;

  /**
   * Normalizing amplitudes from the spherical averaging scheme by the number of
   * sidebands square times the number of octants.
   */
  plan->norm_amplitudes = malloc_double(scheme->octant_orientations);
  cblas_dcopy(scheme->octant_orientations, scheme->amplitudes, 1, plan->norm_amplitudes,
              1);
  double scale = (1.0 / (double)(plan->number_of_sidebands * plan->number_of_sidebands *
                                 plan->n_octants * scheme->n_gamma));
  cblas_dscal(scheme->octant_orientations, scale, plan->norm_amplitudes, 1);

  plan->size = scheme->total_orientations * plan->number_of_sidebands;

  /** Update the mrsimulator plan with the given rotor frequency in Hz. */
  MRS_plan_update_from_rotor_frequency_in_Hz(plan, rotor_frequency_in_Hz);

  /** Update the mrsimulator plan with the given rotor angle in radian. */
  MRS_plan_update_from_rotor_angle_in_rad(plan, rotor_angle_in_rad, allow_4th_rank);
  return plan;
}

/**
 * Update the MRS plan for the given sample rotation frequency in Hz.
 */
void MRS_plan_update_from_rotor_frequency_in_Hz(MRS_plan *plan,
                                                double rotor_frequency_in_Hz) {
  unsigned int size_4;
  // double increment_inverse = 1.0 / increment;
  plan->rotor_frequency_in_Hz = rotor_frequency_in_Hz;
  plan->is_static = (rotor_frequency_in_Hz < 1.0e-3) ? true : false;

  plan->vr_freq = get_FFT_order_freq(plan->number_of_sidebands, rotor_frequency_in_Hz);
  // cblas_dscal(plan->number_of_sidebands, increment_inverse, plan->vr_freq,
  // 1);

  /**
   * calculating the sideband phase multiplier.
   *    pre_phase(m, t) =  I 2π [(exp(I m wr t) - 1)/(I m wr)].
   * for m = [-4, -3, -2, -1]
   * @see get_sideband_phase_components()
   */
  size_4 = 4 * plan->number_of_sidebands;
  plan->pre_phase = malloc_complex128(size_4);
  get_sideband_phase_components(plan->number_of_sidebands, rotor_frequency_in_Hz,
                                (double *)plan->pre_phase);
}

/**
 * Update the MRS plan for the given rotor angle in radians.
 *
 * The method updates the wigner d^l_{m,0}(rotor_angle_in_rad) vectors used in
 * transforming the l-rank tensors from the rotor frame to lab frame. Here l is either 2
 * or 4.
 */
void MRS_plan_update_from_rotor_angle_in_rad(MRS_plan *plan, double rotor_angle_in_rad,
                                             bool allow_4th_rank) {
  unsigned int size_2, size_4, i, j;
  plan->rotor_angle_in_rad = rotor_angle_in_rad;
  /**
   * Calculate wigner-2j d^2_{m,0} vector where m ∈ [-2, 2]. This vector is used to
   * rotate the second-rank tensors from the rotor frame to the lab frame.
   * @see wigner_dm0_vector()
   */
  plan->wigner_d2m0_vector = malloc_double(5);
  wigner_dm0_vector(2, rotor_angle_in_rad, plan->wigner_d2m0_vector);

  plan->wigner_d4m0_vector = NULL;
  if (allow_4th_rank) {
    /**
     * Calculate wigner-4j d^4_{m,0} vector where m ∈ [-4, 4]. This vector is used to
     * rotate the fourth-rank tensors from the rotor frame to the lab frame.
     * @see wigner_dm0_vector()
     */
    plan->wigner_d4m0_vector = malloc_double(9);
    wigner_dm0_vector(4, rotor_angle_in_rad, plan->wigner_d4m0_vector);
  }

  // pre_phase_2 is only calculated for m=-2 and -1 for l=2 rank tensor calculation.
  size_2 = 2 * plan->number_of_sidebands;
  plan->pre_phase_2 = malloc_complex128(size_2);

  /* Copy the pre_phase[m=-2 to 2] to pre_phase2 */
  cblas_zcopy(size_2, (double *)(plan->pre_phase[2 * plan->number_of_sidebands]), 1,
              (double *)(plan->pre_phase_2), 1);
  /**
   * Multiply the wigner-2j d^2_{m,0}(rotor_angle_in_rad) vector to the sideband phase
   * multiplier, pre_phase2. This multiplication accounts for the rotation of the
   * second-rank tensors from the rotor-frame to the lab-frame, thereby, reducing the
   * number of calculations involved per site. This step assumes that the Euler angles
   * involved in the rotation of the 2nd-rank tensors to the lab frame is (0,
   * rotor_angle_in_rad, 0).
   */

  j = 0;
  for (i = 0; i < 2; i++) {
    cblas_zdscal(plan->number_of_sidebands, plan->wigner_d2m0_vector[i],
                 (double *)(plan->pre_phase_2[j]), 1);
    j += plan->number_of_sidebands;
  }

  plan->pre_phase_4 = NULL;

  /* Setup for processing the fourth rank tensors. */
  if (allow_4th_rank) {
    /* pre_phase_4 is only calculated for m=-4, -3, -2, and -1 for l=4 rank tensor
     * calculation. */
    size_4 = 4 * plan->number_of_sidebands;
    plan->pre_phase_4 = malloc_complex128(size_4);
    /* Copy the pre_phase[m=-4 to 4] to pre_phase4 */
    cblas_zcopy(size_4, (double *)(plan->pre_phase), 1, (double *)(plan->pre_phase_4),
                1);

    /**
     * Multiply the wigner-4j d^4_{m,0} vector to the sideband phase multiplier,
     * pre_phase4. This multiplication accounts for the rotation of the fourth rank
     * tensors from the-rotor frame to the lab-frame, therefore, reducing the number of
     * calculations involved per site. This step assumes that the Euler angles involved
     * in the rotation of the 4th rank tensors to the lab frame is (0,
     * rotor_angle_in_rad, 0).
     */

    j = 0;
    for (i = 0; i < 4; i++) {
      cblas_zdscal(plan->number_of_sidebands, plan->wigner_d4m0_vector[i],
                   (double *)(plan->pre_phase_4[j]), 1);
      j += plan->number_of_sidebands;
    }
  }
}

/**
 * Returns a copy of the mrsimulator plan.
 */
MRS_plan *MRS_copy_plan(MRS_plan *plan) {
  MRS_plan *new_plan = malloc(sizeof(MRS_plan));
  new_plan->number_of_sidebands = plan->number_of_sidebands;
  new_plan->rotor_frequency_in_Hz = plan->rotor_frequency_in_Hz;
  new_plan->rotor_angle_in_rad = plan->rotor_angle_in_rad;
  new_plan->vr_freq = plan->vr_freq;
  new_plan->allow_4th_rank = plan->allow_4th_rank;
  new_plan->size = plan->size;
  new_plan->n_octants = plan->n_octants;
  new_plan->norm_amplitudes = plan->norm_amplitudes;
  new_plan->wigner_d2m0_vector = plan->wigner_d2m0_vector;
  new_plan->wigner_d4m0_vector = plan->wigner_d4m0_vector;
  new_plan->pre_phase = plan->pre_phase;
  new_plan->pre_phase_2 = plan->pre_phase_2;
  new_plan->pre_phase_4 = plan->pre_phase_4;
  new_plan->buffer = plan->buffer;
  new_plan->copy = plan->copy;
  new_plan->copy_for_rotor_angle = plan->copy_for_rotor_angle;
  new_plan->copy_for_rotor_freq = plan->copy_for_rotor_freq;
  new_plan->is_static = plan->is_static;
  return new_plan;
}

/**
 * The function evaluates the amplitudes at every orientation and at every sideband per
 * orientation. This is done in two steps.
 * 1) Rotate F2 and F4, given in the crystal or common frame to w2 and w4 in the lab
 *    frame using wigner 2j and 4j rotation matrices, respectively, at all orientations.
 * 2) Evaluate the sideband amplitudes using equation [39] of the reference
 *    https://doi.org/10.1006/jmre.1998.1427.
 */
void MRS_get_amplitudes_from_plan(MRS_averaging_scheme *scheme, MRS_plan *plan,
                                  MRS_fftw_scheme *fftw_scheme,
                                  complex128 *event_amplitudes, bool reset) {
  /* If the number of sidebands is 1, the sideband amplitude at every sideband order is
   * one. In this case, return null,
   */
  if (plan->number_of_sidebands == 1) return;

  /* ================ Calculate the spinning sideband amplitude. ==================== */

  // if (reset) {
  //   cblas_dscal(2 * plan->size, 0.0, (double *)(fftw_scheme->vector), 1);
  // }

  /**
   * Evaluate the exponent of the sideband phase w.r.t the second-rank tensor
   * components. The exponent is given as,
   *
   * w2(Θ) * d^2_{m,0}(rotor_angle_in_rad) * 2πI [(exp(I m ωr t) - 1)/(I m ωr)]
   * |-----lab frame 2nd-rank tensors----|
   *         |------------------------- pre_phase_2 --------------------------|
   *
   * A given element of this product is given as the summation,
   *
   *           res[i, j] = \sum_{m=-2}^2 w2[i, m] * pre_phase_2[m, j],              (1)
   *
   * where the following symmetry holds,
   *
   *    w2[i, m] * pre_phase_2[m, j] = conj(w2[i, -m] * pre_phase_2[-m, j]).
   *
   * The above symmetry simplifies Eq (1) to
   *
   *         res[i, j] = \sum_{m=1}^2 2*imag(w2[i, m] * pre_phase_2[m, j]).         (2)
   *
   * From Eq(2), we find that evaluating half the calculations is sufficient. Since
   * pre_phase_2[0, j] is zero, the m=0 term is dropped from Eq. (2). Notice the scaling
   * factor 2 in Eq. (2). For computation efficiency, this factor is added to the
   * `pre_phase_2` term in the one-time computation step.
   *
   * Here, `pre_phase_2` is pre-calculated and stored in the plan. The calculated
   * product is stored in the fftw_scheme as a complex double array under the variable
   * name `vector`, which is interpreted as a row major matrix of shape
   * `number_of_sidebands` x `total_orientations` with `total_orientations` as the
   * leading dimension.
   */
  cblas_zgemm(CblasRowMajor, CblasTrans, CblasTrans, plan->number_of_sidebands,
              scheme->total_orientations, 2, ONE, (double *)(plan->pre_phase_2),
              plan->number_of_sidebands, (double *)(scheme->w2), 3, ZERO,
              (double *)(fftw_scheme->vector), scheme->total_orientations);

  if (scheme->w4 != NULL) {
    /**
     * Similarly, evaluate the exponent of the sideband phase w.r.t the fourth-rank
     * tensor components. The exponent is given as,
     *
     * w4(Θ) * d^4_{m, 0}(rotor_angle_in_rad) * 2πI[(exp(I m ωr t) - 1)/(I m ωr)]
     * |-----lab frame 4th rank tensors-----|
     *         |-------------------------- pre_phase_4--------------------------|
     *
     * * A given element of this product is given as the summation,
     *
     *           res[i, j] = \sum_{m=-4}^4 w4[i, m] * pre_phase_4[m, j],            (3)
     *
     * where the following symmetry holds,
     *
     *    w2[i, m] * pre_phase_2[m, j] = conj(w4[i, -m] * pre_phase_4[-m, j]).
     *
     * The above symmetry simplifies Eq (3) to
     *
     *         res[i, j] = \sum_{m=1}^4 2*imag(w4[i, m] * pre_phase_4[m, j]).       (4)
     *
     * From Eq(2), we find that evaluating half the calculations is sufficient. Since
     * pre_phase_4[0, j] is zero, the m=0 term is dropped from Eq. (4). Notice the
     * scaling factor 2 in Eq. (4). For computation efficiency, this factor is added to
     * the `pre_phase_4` term in the one-time computation step.
     *
     * where `pre_phase_4` is pre-calculated and stored in the plan. This operation will
     * add and update the values stored in the variable `vector`.
     */
    cblas_zgemm(CblasRowMajor, CblasTrans, CblasTrans, plan->number_of_sidebands,
                scheme->total_orientations, 4, ONE, (double *)(plan->pre_phase_4),
                plan->number_of_sidebands, (double *)(scheme->w4), 5, ONE,
                (double *)(fftw_scheme->vector), scheme->total_orientations);
  }

  /**
   * Evaluate the sideband phase -> exp(vector). Since the real part of the complex data
   * is zero, evaluate the exponential for only the imaginary part. The evaluated value
   * is overwritten on the variable `vector`. */
  vm_double_complex_exp_imag_only(plan->size, fftw_scheme->vector, fftw_scheme->vector);

  /**
   * Evaluate the Fourier transform of the variable, `vector`, -> fft(vector). The fft
   * operation again updates the values of the array, `vector`. */
  fftw_execute(fftw_scheme->the_fftw_plan);

  cblas_dcopy(2 * plan->size, (double *)fftw_scheme->vector, 1,
              (double *)event_amplitudes, 1);

  /** ONLY VALID FOR SINGLE EVENT **/
  /**
   * Evaluate the absolute value square of the `vector` array. The absolute value square
   * is stores as the real part of the `vector` array. The imaginary part is now
   * garbage. This method avoids creating new arrays. */
  vm_double_square_inplace(2 * plan->size, (double *)fftw_scheme->vector);
  cblas_daxpy(plan->size, 1.0, (double *)fftw_scheme->vector + 1, 2,
              (double *)fftw_scheme->vector, 2);
}

void MRS_get_sideband_intensity(MRS_plan *plan, MRS_fftw_scheme *fftw_scheme) {
  /** ONLY VALID FOR SINGLE EVENT **/
  /**
   * Evaluate the absolute value square of the `vector` array. The absolute value square
   * is stores as the real part of the `vector` array. The imaginary part is now
   * garbage. This method avoids creating new arrays. */
  vm_double_square_inplace(2 * plan->size, (double *)fftw_scheme->vector);
  cblas_daxpy(plan->size, 1.0, (double *)fftw_scheme->vector + 1, 2,
              (double *)fftw_scheme->vector, 2);
}

/**
 * Get the lab-frame normalized frequency contributions from the zeroth, second,
 * fourth-rank tensors. Here, normalization refers to dividing the calculated
 * frequencies by the increment of the respective spectral dimension. Normalization
 * makes binning of frequencies on the spectrum faster as bins can then be of 1 unit
 * increments.
 */
void MRS_get_normalized_frequencies_from_plan(MRS_averaging_scheme *scheme,
                                              MRS_plan *plan, double F0, complex128 *F2,
                                              complex128 *F4, MRS_dimension *dim,
                                              double fraction,
                                              unsigned char is_spectral,
                                              double duration) {
  unsigned int i, g_idx, ptr;
  unsigned int freq_size = scheme->n_gamma * scheme->total_orientations;
  double temp, fraction_duration;
  double *f_complex, *local_frequency;

  if (is_spectral) {
    local_frequency = dim->local_frequency;
    fraction_duration = dim->inverse_increment * fraction;
  } else {
    local_frequency = dim->local_phase;
    vm_double_zeros(freq_size, local_frequency);
    fraction_duration = CONST_2PI * duration;
  }

  /**
   * Rotate the F2 and F4 components from the common frame to the rotor frame over all
   * the orientations (alpha, beta). The components are stored in w2 and w4 of the
   * averaging scheme, respectively.
   */
  __batch_wigner_rotation(scheme->octant_orientations, plan->n_octants,
                          scheme->wigner_2j_matrices, F2, scheme->wigner_4j_matrices,
                          F4, scheme->exp_Im_alpha, scheme->w2, scheme->w4);

  /* Normalized the isotropic frequency contribution from the zeroth-rank tensor. */
  if (is_spectral) dim->R0_offset += F0 * fraction_duration;

  /**
   * Rotate the w2 and w4 components from the rotor-frame to the lab-frame. Since only
   * the zeroth-order is relevant in the lab-frame, only evaluate the R20 and R40
   * components. This is equivalent to scaling the w2(0) term by
   * `wigner_d2m0_vector[2]`, that is, d^2(0,0)(rotor_angle).
   */

  /* Normalized local anisotropic frequency contributions from the 2nd-rank tensor. */
  for (g_idx = 0; g_idx < scheme->n_gamma; g_idx++) {
    ptr = scheme->total_orientations * g_idx;
    temp = 2 * fraction_duration;

    if (plan->is_static) {
      for (i = 0; i < 2; i++) {
        // [w2] (a + ib) * [gamma] (c + id) = (ac - bd) + i(ad + bc) -- (1)
        // [w2*] (a - ib) * [gamma*] (c - id) = (ac - bd) - i(ad + bc) -- (2)
        // (1) + (2) = 2 * (ac - bd)
        f_complex =
            (double *)&(scheme->exp_Im_gamma[(2 + i) * scheme->n_gamma + g_idx]);
        plan->buffer = temp * plan->wigner_d2m0_vector[i] * f_complex[0];
        cblas_daxpy(scheme->total_orientations, plan->buffer,
                    (double *)&(scheme->w2[i]), 6, &local_frequency[ptr], 1);
        plan->buffer = -temp * plan->wigner_d2m0_vector[i] * f_complex[1];
        cblas_daxpy(scheme->total_orientations, plan->buffer,
                    (double *)&(scheme->w2[i]) + 1, 6, &local_frequency[ptr], 1);
      }
    }
    plan->buffer = plan->wigner_d2m0_vector[2] * fraction_duration;
    cblas_daxpy(scheme->total_orientations, plan->buffer, (double *)&(scheme->w2[2]), 6,
                &local_frequency[ptr], 1);
  }
  if (plan->allow_4th_rank) {
    /**
     * Similarly, calculate the normalized local anisotropic frequency contributions
     * from the fourth-rank tensor. `wigner_d2m0_vector[4] = d^4(0,0)(rotor_angle)`.
     */
    for (g_idx = 0; g_idx < scheme->n_gamma; g_idx++) {
      ptr = scheme->total_orientations * g_idx;
      temp = 2 * fraction_duration;
      if (plan->is_static) {
        for (i = 0; i < 4; i++) {
          f_complex = (double *)&(scheme->exp_Im_gamma[i * scheme->n_gamma + g_idx]);
          plan->buffer = temp * plan->wigner_d4m0_vector[i] * f_complex[0];
          cblas_daxpy(scheme->total_orientations, plan->buffer,
                      (double *)&scheme->w4[i], 10, &local_frequency[ptr], 1);
          plan->buffer = -temp * plan->wigner_d4m0_vector[i] * f_complex[1];
          cblas_daxpy(scheme->total_orientations, plan->buffer,
                      (double *)&scheme->w4[i] + 1, 10, &local_frequency[ptr], 1);
        }
      }
      plan->buffer = plan->wigner_d4m0_vector[4] * fraction_duration;
      cblas_daxpy(scheme->total_orientations, plan->buffer, (double *)&scheme->w4[4],
                  10, &local_frequency[ptr], 1);
    }
  }

  if (!is_spectral) {  // compute phase.
    plan->buffer = F0 * fraction_duration;
    vm_double_add_vector_offset_inplace(freq_size, local_frequency, plan->buffer,
                                        scheme->phase);
  }
}

static inline void MRS_rotate_single_site_interaction_components(
    site_struct *sites,          // Pointer to a list of sites within a spin system.
    float *transition,           // The spin transition.
    bool allow_4th_rank,         // if true, prep for 4th rank computation.
    double *F0,                  // The F0 components.
    complex128 *F2,              // The F2 components.
    complex128 *F4,              // The F4 components.
    double *F0_temp,             // The temporary F0 components.
    complex128 *F2_temp,         // The temporary F2 components.
    complex128 *F4_temp,         // The temporary R3 components.
    double B0_in_T,              // Magnetic flux density in T.
    unsigned char *freq_contrib  // The pointer to freq contribs boolean.
) {
  unsigned int i, n_sites = sites->number_of_sites;
  double larmor_freq_in_MHz, larmor_freq_in_Hz;
  float *mf = &transition[n_sites], *mi = transition;
  double F2_shield[10];
  double F2_quad[10];

  /* Frequency computation for sites */
  for (i = 0; i < n_sites; i++) {
    if (*mi == *mf) {
      mi++;
      mf++;
      continue;
    }
    larmor_freq_in_MHz = -B0_in_T * sites->gyromagnetic_ratio[i];
    larmor_freq_in_Hz = larmor_freq_in_MHz * 1.0e6;

    /* Nuclear shielding components ================================================= */
    /*  Upto the first order */
    FT_1st_order_nuclear_shielding_tensor_components(
        F0_temp, F2_shield,
        sites->isotropic_chemical_shift_in_ppm[i] * larmor_freq_in_MHz *
            sites->one_minus_sigma_iso_ref[i],
        sites->shielding_symmetric_zeta_in_ppm[i] * larmor_freq_in_MHz,
        sites->shielding_symmetric_eta[i], &sites->shielding_orientation[3 * i], *mf,
        *mi);

    // Note:
    // In PAS, F2_shield = [1/√6 v_0ζη, 0, -v_0ζ, 0 , 1/√6 v_0ζη] * p(mf, mi)

    // in-place update the F0 and F2 components.
    if (freq_contrib[0]) *F0 += *F0_temp;
    if (freq_contrib[1]) vm_double_add_inplace(10, (double *)F2_shield, (double *)F2);
    /* ============================================================================== */
    if (sites->spin[i] == 0.5) {
      mi++;
      mf++;
      continue;
    }
    /* Electric quadrupolar components ============================================== */
    /*  Upto the first order */
    FT_1st_order_electric_quadrupole_tensor_components(
        F2_quad, sites->spin[i], sites->quadrupolar_Cq_in_Hz[i],
        sites->quadrupolar_eta[i], &sites->quadrupolar_orientation[3 * i], *mf, *mi);

    // Note:
    // In PAS, F2_quad = [-1/6 v_q η, 0, 1/√6 v_q, 0 , -1/6 v_q η] * d(mf, mi)

    // in-place update the F2 components.
    if (freq_contrib[2]) vm_double_add_inplace(10, (double *)F2_quad, (double *)F2);

    /*  Upto the second order */
    if (allow_4th_rank) {
      FT_2nd_order_electric_quadrupole_tensor_components(
          F0_temp, F2_temp, F4_temp, sites->spin[i], larmor_freq_in_Hz,
          sites->quadrupolar_Cq_in_Hz[i], sites->quadrupolar_eta[i],
          &sites->quadrupolar_orientation[3 * i], *mf, *mi);

      // in-place update the F0, F2, and F4 components.
      if (freq_contrib[3]) *F0 += *F0_temp;
      if (freq_contrib[4]) vm_double_add_inplace(10, (double *)F2_temp, (double *)F2);
      if (freq_contrib[5]) vm_double_add_inplace(18, (double *)F4_temp, (double *)F4);
    }
    /* ============================================================================== */

    /* Shielding Quad cross term components ========================================= */
    FT_NS_EQ_cross_tensor_components(F0_temp, F2_temp, F4_temp, F2_shield, F2_quad,
                                     larmor_freq_in_Hz, *mf, *mi);
    if (freq_contrib[9]) *F0 += *F0_temp;
    if (freq_contrib[10]) vm_double_add_inplace(10, (double *)F2_temp, (double *)F2);
    if (freq_contrib[11]) vm_double_add_inplace(18, (double *)F4_temp, (double *)F4);

    mi++;
    mf++;
  }
}

static inline void quad_coupling_cross_terms(
    coupling_struct *couplings,  // Pointer to a list of couplings within a spin system.
    site_struct *sites,          // Pointer to a list of sites within a spin system.
    int site_index,              // site index
    unsigned int coupling_index,  // coupled index
    double *F0,                   // The F0 components.
    complex128 *F2,               // The F2 components.
    complex128 *F4,               // The F4 components.
    double *F0_temp,              // The temporary F0 components.
    complex128 *F2_temp,          // The temporary F2 components.
    complex128 *F4_temp,          // The temporary F4 components.
    double B0_in_T,               // Magnetic flux density in T.
    unsigned char *freq_contrib,  // The pointer to freq contribs boolean.
    float mIi,   // inital transition state of site coupled to the quad site (p)
    float mSqi,  // inital transition state of the quad site (d)
    float mIf,   // final transition state of site coupled to the quad site (p)
    float mSqf   // final transition state of the quad site (d)
) {
  double larmor_freq_in_Hz, spin, Delta_2q[10];

  spin = sites->spin[site_index];
  if (spin == 0.5) return;

  fsSST_1st_order_electric_quadrupole_tensor_components(
      Delta_2q, spin, sites->quadrupolar_Cq_in_Hz[site_index],
      sites->quadrupolar_eta[site_index],
      &sites->quadrupolar_orientation[3 * site_index]);

  larmor_freq_in_Hz = -B0_in_T * sites->gyromagnetic_ratio[site_index] * 1.0e6;

  // Site S cross J-coupling
  FT_Quad_coupling_cross_tensor_components(
      F0_temp, F2_temp, F4_temp, spin, Delta_2q,
      couplings->j_symmetric_zeta_in_Hz[coupling_index],
      couplings->j_symmetric_eta[coupling_index],
      &couplings->j_orientation[3 * coupling_index], larmor_freq_in_Hz, mIf, mIi, mSqf,
      mSqi);  // transition function is pd (p for spin I, and d for spin S)

  if (freq_contrib[12]) *F0 += *F0_temp;
  if (freq_contrib[13]) vm_double_add_inplace(10, (double *)F2_temp, (double *)F2);
  if (freq_contrib[14]) vm_double_add_inplace(18, (double *)F4_temp, (double *)F4);

  // Site S cross dipolar-coupling
  FT_Quad_coupling_cross_tensor_components(
      F0_temp, F2_temp, F4_temp, spin, Delta_2q,
      2.0 * couplings->dipolar_coupling_in_Hz[coupling_index], 0,
      &couplings->dipolar_orientation[3 * coupling_index], larmor_freq_in_Hz, mIf, mIi,
      mSqf, mSqi);  // transition function is pd (p for spin I, and d for spin S)

  if (freq_contrib[15]) *F0 += *F0_temp;
  if (freq_contrib[16]) vm_double_add_inplace(10, (double *)F2_temp, (double *)F2);
  if (freq_contrib[17]) vm_double_add_inplace(18, (double *)F4_temp, (double *)F4);
}

static inline void MRS_rotate_coupled_site_interaction_components(
    coupling_struct *couplings,  // Pointer to a list of couplings within a spin system.
    site_struct *sites,          // Pointer to a list of sites within a spin system.
    float *transition,           // The spin transition.
    unsigned int n_sites,        // The number of sites.
    double *F0,                  // The F0 components.
    complex128 *F2,              // The F2 components.
    complex128 *F4,              // The F4 components.
    double *F0_temp,             // The temporary F0 components.
    complex128 *F2_temp,         // The temporary F2 components.
    complex128 *F4_temp,         // The temporary F4 components.
    double B0_in_T,              // Magnetic flux density in T.
    unsigned char *freq_contrib  // The pointer to freq contribs boolean.
) {
  unsigned int i, j = 0, n_couplings = couplings->number_of_couplings;
  int site_index_I, site_index_S;
  float mIf, mSf, mIi, mSi;

  /* Frequency computation for couplings */
  for (i = 0; i < n_couplings; i++) {
    site_index_I = couplings->site_index[j++];
    site_index_S = couplings->site_index[j++];

    mIi = transition[site_index_I];
    mSi = transition[site_index_S];
    mIf = transition[site_index_I + n_sites];
    mSf = transition[site_index_S + n_sites];

    // Weakly coupled J-couplings
    FT_1st_order_weak_J_coupling_tensor_components(
        F0_temp, F2_temp, couplings->isotropic_j_in_Hz[i],
        couplings->j_symmetric_zeta_in_Hz[i], couplings->j_symmetric_eta[i],
        &couplings->j_orientation[3 * i], mIf, mIi, mSf, mSi);

    // in-place update the F0 and F2 components.
    if (freq_contrib[6]) *F0 += *F0_temp;
    if (freq_contrib[7]) vm_double_add_inplace(10, (double *)F2_temp, (double *)F2);

    // Weakly coupled dipolar-couplings
    FT_1st_order_weak_dipolar_coupling_tensor_components(
        F2_temp, couplings->dipolar_coupling_in_Hz[i],
        &couplings->dipolar_orientation[3 * i], mIf, mIi, mSf, mSi);

    // in-place update the F2 components.
    if (freq_contrib[8]) vm_double_add_inplace(10, (double *)F2_temp, (double *)F2);

    // Quad J-couplings cross components
    // Site S
    quad_coupling_cross_terms(couplings, sites, site_index_S, i, F0, F2, F4, F0_temp,
                              F2_temp, F4_temp, B0_in_T, freq_contrib, mIi, mSi, mIf,
                              mSf);
    // Site I
    quad_coupling_cross_terms(couplings, sites, site_index_I, i, F0, F2, F4, F0_temp,
                              F2_temp, F4_temp, B0_in_T, freq_contrib, mSi, mIi, mSf,
                              mIf);
  }
}

/**
 * The function evaluates the tensor components from the principal axis system (PAS) to
 * the common frame of the spin system.
 */
void MRS_rotate_components_from_PAS_to_common_frame(
    site_struct *sites,          // Pointer to a list of sites within a spin system.
    coupling_struct *couplings,  // Pointer to a list of couplings within a spin system.
    float *transition,           // The spin transition.
    bool allow_4th_rank,         // If true, prep for 4th rank computation.
    double *F0,                  // The F0 components.
    complex128 *F2,              // The F2 components.
    complex128 *F4,              // The F4 components.
    double *F0_temp,             // The temporary F0 components.
    complex128 *F2_temp,         // The temporary F2 components.
    complex128 *F4_temp,         // The temporary R3 components.
    double B0_in_T,              // Magnetic flux density in T.
    unsigned char *freq_contrib  // The pointer to freq contribs boolean.
) {
  /* The following codeblock populates the product of spatial part, Rlm, of the tensor
   * and the spin transition function, T(mf, mi) for
   *      zeroth rank, F0 = [ R00 ] * T(mf, mi)
   *      second rank, F2 = [ R2m ] * T(mf, mi) where m ∈ [-2, 2].
   *      fourth rank, F4 = [ R4m ] * T(mf, mi) where m ∈ [-4, 4].
   * Here, mf, mi are the spin quantum numbers of the final and initial energy state of
   * the spin transition. The term `Rlm` is the coefficient of the irreducible spherical
   * tensor of rank `l` and order `m`. For more information, see reference
   *
   *   Symmetry pathways in solid-state NMR. PNMRS 2011 59(2):12 1-96.
   *   https://doi.org/10.1016/j.pnmrs.2010.11.003
   *
   */
  MRS_rotate_single_site_interaction_components(sites, transition, allow_4th_rank, F0,
                                                F2, F4, F0_temp, F2_temp, F4_temp,
                                                B0_in_T, freq_contrib);

  if (couplings->number_of_couplings == 0) return;
  MRS_rotate_coupled_site_interaction_components(
      couplings, sites, transition, sites->number_of_sites, F0, F2, F4, F0_temp,
      F2_temp, F4_temp, B0_in_T, freq_contrib);
}

/**
 * The function calculates the following.
 *
 *   pre_phase(m, t) = I 2π [(exp(I m ωr t) - 1)/(I m ωr)]
 *                   = (2π / m ωr) (exp(I m ωr t) - 1)
 *                     |--scale--|
 *                   = scale (exp(I m ωr t) - 1)
 *                   = scale [[cos(m ωr t) -1] + I sin(m ωr t)],
 *
 * where ωr is the sample spinning frequency in Hz, m goes from -4 to 4, t is a vector
 * of length `number_of_sidebands` given as
 *
 *    t = [0, 1, ... number_of_sidebands-1]/(ωr*number_of_sidebands),
 *
 * and `pre_phase` is a matrix of shape, `9 x number_of_sidebands`.
 *
 * Also,
 *   pre_phase(-m, t) = (-2π / m ωr) (exp(-I m ωr t) - 1)
 *                    = -scale [[cos(m ωr t) -1] -Isin(m ωr t)]
 *                    = scale [-[cos(m ωr t) -1] +Isin(m ωr t)]
 * That is, pre_phase[-m] = -Re(pre_phase[m]) + Im(pre_phase[m])
 */
void get_sideband_phase_components_2(unsigned int number_of_sidebands,
                                     double rotor_frequency_in_Hz,
                                     complex128 *pre_phase) {
  int m, i;
  double spin_angular_freq, tau, scale;

  double *input = malloc_double(number_of_sidebands);
  double *ones = malloc_double(number_of_sidebands);
  double *phase = malloc_double(number_of_sidebands);

  vm_double_ones(number_of_sidebands, ones);
  vm_double_arange(number_of_sidebands, input);

  // Calculate the spin angular frequency
  spin_angular_freq = rotor_frequency_in_Hz * CONST_2PI;

  // Calculate tau, where tau = (rotor period / number of phase steps)
  tau = 1.0 / ((double)number_of_sidebands * rotor_frequency_in_Hz);

  // pre-calculate the m omega spinning frequencies
  double m_wr[9] = {-4., -3., -2., -1., 0., 1., 2., 3., 4.};
  cblas_dscal(9, spin_angular_freq, m_wr, 1);

  for (m = 0; m <= 3; m++) {
    /**
     * Evaluate pre_phase = scale * (cexp(I * phase) - 1.0), where
     *    phase = m_wr[m] * tau * [0 .. number_of_sidebands-1] and
     *    scale = 2π/m_wr[m].
     */
    i = m * number_of_sidebands;
    scale = CONST_2PI / m_wr[m];

    // step 1. calculate phase
    vm_double_ramp(number_of_sidebands, input, m_wr[m] * tau, 0.0, phase);

    // step 2. evaluate cexp(I * phase) = cos(phase) + I sin(phase)
    vm_cosine_I_sine(number_of_sidebands, phase, &pre_phase[i]);

    // step 3. subtract 1.0 from pre_phase
    cblas_daxpy(number_of_sidebands, -1.0, ones, 1, (double *)(pre_phase[i]), 2);

    // step 4. scale pre_phase with factor `scale`
    cblas_zdscal(number_of_sidebands, scale, (double *)(pre_phase[i]), 1);

    /**
     * The expression pre_phase[m] = scale * (cexp(I * phase) - 1.0) given above for
     * positive m is related to -m as
     *
     * pre_phase[-m] = -Re(pre_phase[m]) + Im(pre_phase[m])
     */
    cblas_zcopy(number_of_sidebands, (double *)(pre_phase[i]), 1,
                (double *)(pre_phase[(8 - m) * number_of_sidebands]), 1);
    cblas_dscal(number_of_sidebands, -1.0, (double *)(pre_phase[i]), 2);
  }
  vm_double_zeros(2 * number_of_sidebands,
                  (double *)(pre_phase[4 * number_of_sidebands]));

  free(input);
  free(phase);
  free(ones);
}

/**
 * The function calculates the following.
 *   pre_phase(m, t) = I 2π [(exp(I m ωr t) - 1)/(I m ωr)]
 *                   = (2π / m ωr) (exp(I m ωr t) - 1)
 *                     |--scale--|
 *                   = scale * (exp(I m ωr t) - 1)
 * where ωr is the sample spinning frequency in Hz, m goes from -4 to -1, and t is a
 * vector of length `number_of_sidebands` given as
 *    t = [0, 1, ... number_of_sidebands-1]/(ωr*number_of_sidebands).
 *
 * `pre_phase` is a matrix of shape, `9 x number_of_sidebands`, with number_of_sidebands
 * as the leading dimension. The first number_of_sidebands entries corresponds to
 * m_wr=-4.
 */
void get_sideband_phase_components(unsigned int number_of_sidebands,
                                   double sample_rotation_frequency,
                                   double *restrict pre_phase) {
  double spin_angular_freq, tau, wrt, pht, scale;
  unsigned int step, m;

  // Calculate the spin angular frequency
  spin_angular_freq = sample_rotation_frequency * CONST_2PI;

  // Calculate tau increments, where tau = (rotor period / number of phase steps)
  tau = 1.0 / ((double)number_of_sidebands * sample_rotation_frequency);

  double m_wr[4] = {-4., -3., -2., -1.};
  cblas_dscal(4, spin_angular_freq, m_wr, 1);

  for (m = 0; m < 4; m++) {
    wrt = m_wr[m] * tau;
    pht = 0.0;
    // scale = 2 * CONST_2PI / m_wr[m]. See Eq.(2) and (4) for reason for the factor 2.
    scale = CONST_4PI / m_wr[m];
    for (step = 0; step < number_of_sidebands; step++) {
      *pre_phase++ = scale * (cos(pht) - 1.0);
      *pre_phase++ = scale * sin(pht);
      pht += wrt;
    }
  }
}
