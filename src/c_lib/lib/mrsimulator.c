// -*- coding: utf-8 -*-
//
//  mrsimulator.c
//
//  @copyright Deepansh J. Srivastava, 2019-2020.
//  Created by Deepansh J. Srivastava, Jun 9, 2019
//  Contact email = srivastava.89@osu.edu
//

#include "mrsimulator.h"

/**
 * @func MRS_free_plan
 *
 * Free the buffers and pre-calculated tables from the mrsimulator plan.
 */
void MRS_free_plan(MRS_plan *the_plan) {
  if (!the_plan->vr_freq) {
    free(the_plan->vr_freq);
  }

  if (!the_plan->wigner_d2m0_vector) {
    free(the_plan->wigner_d2m0_vector);
  }

  if (!the_plan->wigner_d4m0_vector) {
    free(the_plan->wigner_d4m0_vector);
  }

  if (!the_plan->norm_amplitudes) {
    free(the_plan->norm_amplitudes);
  }

  if (!the_plan->pre_phase) {
    free(the_plan->pre_phase);
  }

  if (!the_plan->pre_phase_2) {
    free(the_plan->pre_phase_2);
  }

  if (!the_plan->pre_phase_4) {
    free(the_plan->pre_phase_4);
  }
}

/**
 * @func MRS_plan_free_rotor_angle_in_rad
 *
 * Free the memory from the mrsimulator plan associated with the wigner
 * d^l_{m,0}(rotor_angle_in_rad) vectors. Here, l=2 or 4.
 */
void MRS_plan_free_rotor_angle_in_rad(MRS_plan *plan) {
  free(plan->wigner_d2m0_vector);
  free(plan->wigner_d4m0_vector);
  plan->wigner_d2m0_vector = NULL;
  plan->wigner_d4m0_vector = NULL;
}

/**
 * @func MRS_create_plan
 *
 * Create a new mrsimulator plan.
 *
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
MRS_plan *MRS_create_plan(MRS_averaging_scheme *scheme,
                          unsigned int number_of_sidebands,
                          double sample_rotation_frequency_in_Hz,
                          double rotor_angle_in_rad, double increment,
                          bool allow_fourth_rank) {
  MRS_plan *plan = malloc(sizeof(MRS_plan));
  plan->number_of_sidebands = number_of_sidebands;
  plan->sample_rotation_frequency_in_Hz = sample_rotation_frequency_in_Hz;
  plan->rotor_angle_in_rad = rotor_angle_in_rad;

  plan->allow_fourth_rank = allow_fourth_rank;
  plan->one[0] = 1.0;
  plan->one[1] = 0.0;
  plan->zero[0] = 0.0;
  plan->zero[1] = 0.0;

  /**
   * Update the mrsimulator plan with the given spherical averaging scheme. We
   * create the coordinates on the surface of the unit sphere by projecting the
   * points on the face of the octahedron to a unit sphere. Usually, before
   * updating the averaging scheme, the memory allocated by the previous scheme
   * must be freed. Since, we are creating the scheme for this plan for the very
   * first time, there is no need to call MRS_free_averaging_plan() method.
   */

  plan->n_octants = 1;
  if (scheme->integration_volume == 1) plan->n_octants = 4;
  if (scheme->integration_volume == 2) plan->n_octants = 8;

  /**
   * Normalizing amplitudes from the spherical averaging scheme by the number
   * of sidebands square times the number of octants.
   */
  plan->norm_amplitudes = malloc_double(scheme->octant_orientations);
  cblas_dcopy(scheme->octant_orientations, scheme->amplitudes, 1,
              plan->norm_amplitudes, 1);
  double scale = (1.0 / (double)(plan->number_of_sidebands *
                                 plan->number_of_sidebands * plan->n_octants));
  cblas_dscal(scheme->octant_orientations, scale, plan->norm_amplitudes, 1);

  plan->size = scheme->total_orientations * plan->number_of_sidebands;

  MRS_plan_update_from_sample_rotation_frequency_in_Hz(
      plan, increment, sample_rotation_frequency_in_Hz);

  return plan;
}

/**
 * @func MRS_plan_update_from_sample_rotation_frequency_in_Hz
 *
 * Update the MRS plan for the given sample rotation frequency in Hz.
 */
void MRS_plan_update_from_sample_rotation_frequency_in_Hz(
    MRS_plan *plan, double increment, double sample_rotation_frequency_in_Hz) {
  unsigned int size_4;
  // double increment_inverse = 1.0 / increment;
  plan->sample_rotation_frequency_in_Hz = sample_rotation_frequency_in_Hz;

  plan->vr_freq = __get_frequency_in_FFT_order(plan->number_of_sidebands,
                                               sample_rotation_frequency_in_Hz);
  // cblas_dscal(plan->number_of_sidebands, increment_inverse, plan->vr_freq,
  // 1);

  /**
   * calculating the sideband phase multiplier.
   *    pre_phase(m, t) =  I 2π [(exp(I m wr t) - 1)/(I m wr)].
   * for m = [-4, -3, .. 3, 4]
   * @see __get_components()
   */
  size_4 = 9 * plan->number_of_sidebands;
  plan->pre_phase = malloc_complex128(size_4);
  __get_components(plan->number_of_sidebands, sample_rotation_frequency_in_Hz,
                   (double *)plan->pre_phase);

  /**
   * Update the mrsimulator plan with the given rotor angle in radian.
   * This method updates the wigner d^l_{m,0}(rotor_angle_in_rad) vectors used
   * in tranforming the l-rank tensors from the rotor frame to lab frame. Here l
   * is either 2 or 4.
   */
  MRS_plan_update_from_rotor_angle_in_rad(plan, plan->rotor_angle_in_rad,
                                          plan->allow_fourth_rank);
}

/**
 * @func MRS_plan_update_from_rotor_angle_in_rad
 *
 * Update the MRS plan for the given rotor angle in radians.
 */
void MRS_plan_update_from_rotor_angle_in_rad(MRS_plan *plan,
                                             double rotor_angle_in_rad,
                                             bool allow_fourth_rank) {
  unsigned int size_2, size_4, i, j;
  plan->rotor_angle_in_rad = rotor_angle_in_rad;
  /**
   * Calculate wigner-2j d^2_{m,0} vector where m ∈ [-2, 2]. This vector is
   * used to rotate the second-rank tensors from the rotor frame to the lab
   * frame.
   * @see wigner_dm0_vector()
   */
  plan->wigner_d2m0_vector = malloc_double(5);
  wigner_dm0_vector(2, rotor_angle_in_rad, plan->wigner_d2m0_vector);

  plan->wigner_d4m0_vector = NULL;
  if (allow_fourth_rank) {
    /**
     * Calculate wigner-4j d^4_{m,0} vector where m ∈ [-4, 4]. This vector is
     * used to rotate the fourth-rank tensors from the rotor frame to the lab
     * frame.
     * @see wigner_dm0_vector()
     */
    plan->wigner_d4m0_vector = malloc_double(9);
    wigner_dm0_vector(4, rotor_angle_in_rad, plan->wigner_d4m0_vector);
  }

  size_2 = 5 * plan->number_of_sidebands;
  plan->pre_phase_2 = malloc_complex128(size_2);

  /* Copy the pre_phase[m=-2 to 2] to pre_phase2 */
  cblas_zcopy(size_2,
              (double *)(plan->pre_phase[2 * plan->number_of_sidebands]), 1,
              (double *)(plan->pre_phase_2), 1);
  /**
   * Multiply the wigner-2j d^2_{m,0}(rotor_angle_in_rad) vector to the sideband
   * phase multiplier, pre_phase2. This multiplication accounts for the rotation
   * of the second-rank tensors from the rotor-frame to the lab-frame, thereby,
   * reducing the number of calculations involved per site. This step assumes
   * that the Euler angles invloved in the rotation of the 2nd-rank tensors to
   * the lab frame is (0, rotor_angle_in_rad, 0).
   */

  j = 0;
  for (i = 0; i < 5; i++) {
    cblas_zdscal(plan->number_of_sidebands, plan->wigner_d2m0_vector[i],
                 (double *)(plan->pre_phase_2[j]), 1);
    j += plan->number_of_sidebands;
  }

  plan->pre_phase_4 = NULL;

  /* Setup for processing the fourth rank tensors. */
  if (allow_fourth_rank) {
    /* Copy the pre_phase[m=-4 to 4] to pre_phase4 */
    size_4 = 9 * plan->number_of_sidebands;
    plan->pre_phase_4 = malloc_complex128(size_4);
    cblas_zcopy(size_4, (double *)(plan->pre_phase), 1,
                (double *)(plan->pre_phase_4), 1);

    /**
     * Multiply the wigner-4j d^4_{m,0} vector to the sideband phase multiplier,
     * pre_phase4. This multiplication accounts for the rotation of the fourth
     * rank tensors from the-rotor frame to the lab-frame, thereby, reducing the
     * number of calculations involved per site. This step assumes that the
     * Euler angles involved in the rotation of the 4th rank tensors to the lab
     * frame is (0, rotor_angle_in_rad, 0).
     */

    j = 0;
    for (i = 0; i < 9; i++) {
      cblas_zdscal(plan->number_of_sidebands, plan->wigner_d4m0_vector[i],
                   (double *)(plan->pre_phase_4[j]), 1);
      j += plan->number_of_sidebands;
    }
  }
}

/**
 * @func MRS_copy_plan
 *
 * Returns a copy of the mrsimulator plan.
 */
MRS_plan *MRS_copy_plan(MRS_plan *plan) {
  MRS_plan *new_plan = malloc(sizeof(MRS_plan));
  new_plan->averaging_scheme = plan->averaging_scheme;
  new_plan->number_of_sidebands = plan->number_of_sidebands;
  new_plan->sample_rotation_frequency_in_Hz =
      plan->sample_rotation_frequency_in_Hz;
  new_plan->rotor_angle_in_rad = plan->rotor_angle_in_rad;
  new_plan->vr_freq = plan->vr_freq;
  new_plan->allow_fourth_rank = plan->allow_fourth_rank;
  new_plan->size = plan->size;
  new_plan->n_octants = plan->n_octants;
  new_plan->norm_amplitudes = plan->norm_amplitudes;
  new_plan->wigner_d2m0_vector = plan->wigner_d2m0_vector;
  new_plan->wigner_d4m0_vector = plan->wigner_d4m0_vector;
  new_plan->pre_phase = plan->pre_phase;
  new_plan->pre_phase_2 = plan->pre_phase_2;
  new_plan->pre_phase_4 = plan->pre_phase_4;
  new_plan->one[0] = plan->one[0];
  new_plan->one[1] = plan->one[1];
  new_plan->zero[0] = plan->zero[0];
  new_plan->zero[1] = plan->zero[1];
  new_plan->buffer = plan->buffer;
  return new_plan;
}

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
void MRS_get_amplitudes_from_plan(MRS_averaging_scheme *scheme, MRS_plan *plan,
                                  MRS_fftw_scheme *fftw_scheme, bool refresh) {
  /* If the number of sidebands is 1, the sideband amplitude at every sideband
   * order is one. In this case, return null,
   */
  if (plan->number_of_sidebands == 1) {
    return;
  }

  /* =========== Calculate the spinning sideband amplitude. ================= */

  // if (refresh) {
  //   cblas_dscal(2 * plan->size, 0.0, (double *)(fftw_scheme->vector), 1);
  // }

  /**
   * Evaluate the exponent of the sideband phase w.r.t the second-rank tensors.
   * The exponent is given as,
   *
   * w2(Θ)*d^2_{m,0}(rotor_angle_in_rad) * 2πI [(exp(I m ωr t) - 1)/(I m ωr)]
   * |----lab frame 2nd-rank tensors---|
   *       |------------------------- pre_phase_2 ---------------------------|
   *
   * where `pre_phase_2` is pre-calculated and stored in the plan. The product
   * is stored in the fftw_scheme as a complex double array under the variable
   * `vector`, which is interpreted as a row major matrix of shape
   * `number_of_sidebands` x `total_orientations` with `total_orientations`
   * as the leading dimension.
   */
  cblas_zgemm(CblasRowMajor, CblasTrans, CblasTrans, plan->number_of_sidebands,
              scheme->total_orientations, 5, (double *)(plan->one),
              (double *)(plan->pre_phase_2), plan->number_of_sidebands,
              (double *)(scheme->w2), 5, (double *)(plan->zero),
              (double *)(fftw_scheme->vector), scheme->total_orientations);

  if (scheme->w4 != NULL) {
    /**
     * Similarly, evaluate the exponent of the sideband phase w.r.t the fourth
     * rank tensors. The exponent is given as,
     *
     * w4(Θ)*d^4_{m, 0}(rotor_angle_in_rad) * 2πI[(exp(I m ωr t) - 1)/(I m ωr)]
     * |----lab frame 4th rank tensors----|
     *       |-------------------------- pre_phase_4--------------------------|
     *
     * where `pre_phase_4` is pre-calculated and stored in the plan. This
     * operation will add and update the values stored in the variable `vector`.
     */
    cblas_zgemm(CblasRowMajor, CblasTrans, CblasTrans,
                plan->number_of_sidebands, scheme->total_orientations, 9,
                (double *)(plan->one), (double *)(plan->pre_phase_4),
                plan->number_of_sidebands, (double *)(scheme->w4), 9,
                (double *)(plan->one), (double *)(fftw_scheme->vector),
                scheme->total_orientations);
  }

  /**
   * Evaluate the sideband phase -> exp(vector). Since the real part of the
   * complex data is zero, evaluate the exponential for only imaginary part.
   * The evaluated value is overwritten on the variable `vector`.
   */
  vm_double_complex_exp_imag_only(plan->size, fftw_scheme->vector,
                                  fftw_scheme->vector);

  /**
   * Evaluate the Fourier transform of the variable, `vector`, -> fft(vector).
   * The fft operation again updates the values of the array, `vector`.
   */
  fftw_execute(fftw_scheme->the_fftw_plan);

  /**
   * Evaluate the absolute value square of the vector array. The absolute value
   * square is stores as the real part of the `vector` array. The imaginary
   * part is now garbage. This method avoids creating new arrays.
   */
  vm_double_square_inplace(2 * plan->size, (double *)fftw_scheme->vector);
  cblas_daxpy(plan->size, 1.0, (double *)fftw_scheme->vector + 1, 2,
              (double *)fftw_scheme->vector, 2);

  /* Scaling the absolute value square with the powder scheme weights. Only
   * the real part is scaled and the imaginary part is left as is.
   */
  // for (i = 0; i < scheme->octant_orientations; i++) {
  //   cblas_dscal(plan->n_octants * plan->number_of_sidebands,
  //               plan->norm_amplitudes[i], (double *)&fftw_scheme->vector[i],
  //               2 * scheme->octant_orientations);
  // }
}

/**
 * @func MRS_get_frequencies_from_plan
 *
 * Get the lab-frame frequency contributions from the zeroth, second,
 * fourth-rank tensors.
 */
void MRS_get_frequencies_from_plan(MRS_averaging_scheme *scheme, MRS_plan *plan,
                                   double R0, complex128 *R2, complex128 *R4,
                                   bool refresh, MRS_sequence *seq) {
  /**
   * Rotate the R2 and R4 components from the common frame to the rotor frame
   * over all the orientations. The componets are stored in w2 and w4 of the
   * averaging scheme, respectively.
   */
  __batch_wigner_rotation(scheme->octant_orientations, plan->n_octants,
                          scheme->wigner_2j_matrices, R2,
                          scheme->wigner_4j_matrices, R4, scheme->exp_Im_alpha,
                          scheme->w2, scheme->w4);

  /* If refresh is true, zero the local_frequencies before update. */
  if (refresh) {
    cblas_dscal(scheme->total_orientations, 0.0, seq->local_frequency, 1);
    seq->R0_offset = 0.0;
  }

  /* Add the isotropic frequency contribution from the zeroth-rank tensor. */
  seq->R0_offset += R0;
  // vm_double_add_offset_inplace(scheme->total_orientations, plan->R0_offset,
  //                              seq->local_frequency);

  /**
   * Calculate the local anisotropic frequency contributions from the 2nd-rank
   * tensor. The w2 and w4 frequencies from the plan are in the rotor-frame. Use
   * the wigner-2j and 4j rotations to transform the frequencies in the
   * lab-frame.
   */
  /* Wigner 2j rotation for the second-rank tensor frequency contributions. */
  plan->buffer = plan->wigner_d2m0_vector[2];
  cblas_daxpy(scheme->total_orientations, plan->buffer,
              (double *)&scheme->w2[2], 10, seq->local_frequency, 1);
  if (plan->allow_fourth_rank) {
    /* Wigner 4j rotation for the fourth-rank tensor frequency contributions. */
    plan->buffer = plan->wigner_d4m0_vector[4];
    cblas_daxpy(scheme->total_orientations, plan->buffer,
                (double *)&scheme->w4[4], 18, seq->local_frequency, 1);
  }
}

/**
 * @func MRS_get_normalized_frequencies_from_plan
 *
 * Get the lab-frame normalized frequency contributions from the zeroth, second,
 * fourth-rank tensors. Here, normalization refers to dividing the frequencies
 * by the increment of the spectral dimension.
 */
void MRS_get_normalized_frequencies_from_plan(MRS_averaging_scheme *scheme,
                                              MRS_plan *plan, double R0,
                                              complex128 *R2, complex128 *R4,
                                              bool refresh, MRS_sequence *seq,
                                              double fraction) {
  /**
   * Rotate the R2 and R4 components from the common frame to the rotor frame
   * over all the orientations. The componets are stored in w2 and w4 of the
   * averaging scheme, respectively.
   */
  __batch_wigner_rotation(scheme->octant_orientations, plan->n_octants,
                          scheme->wigner_2j_matrices, R2,
                          scheme->wigner_4j_matrices, R4, scheme->exp_Im_alpha,
                          scheme->w2, scheme->w4);

  /* If refresh is true, zero the local_frequencies before update */
  if (refresh) {
    cblas_dscal(scheme->total_orientations, 0.0, seq->local_frequency, 1);
    seq->R0_offset = 0.0;
  }

  /* Normalized isotropic frequency contribution from the zeroth-rank tensor. */
  seq->R0_offset += R0 * seq->inverse_increment * fraction;
  // vm_double_add_offset_inplace(scheme->total_orientations, seq->R0_offset,
  //                              seq->local_frequency);

  /**
   * Rotate the w2 and w4 components from the rotor-frame to the lab-frame.
   * Since only the zeroth-order is relevent in the lab-frame, only evalute the
   * R20 and R40 components. This is equivalent to scaling the w2(0) term by
   * `wigner_d2m0_vector[2]`, that is, d^2(0,0)(rotor_angle).
   */

  /**
   * Normalized local anisotropic frequency contributions from the 2nd-rank
   * tensor.
   */
  plan->buffer =
      seq->inverse_increment * plan->wigner_d2m0_vector[2] * fraction;
  cblas_daxpy(scheme->total_orientations, plan->buffer,
              (double *)&scheme->w2[2], 10, seq->local_frequency, 1);
  if (plan->allow_fourth_rank) {
    /**
     * Similarly, calculate the normalized local anisotropic frequency
     * contributions from the fourth-rank tensor. The `wigner_d2m0_vector[4]` is
     * d^4(0,0)(rotor_angle).
     */
    plan->buffer =
        seq->inverse_increment * plan->wigner_d4m0_vector[4] * fraction;
    cblas_daxpy(scheme->total_orientations, plan->buffer,
                (double *)&scheme->w4[4], 18, seq->local_frequency, 1);
  }
}

/**
 * @func MRS_rotate_components_from_PAS_to_common_frame
 *
 * The function evaluates the tensor components from the principal axis system
 * (PAS) to the common frame of the isotopomer.
 */
void MRS_rotate_components_from_PAS_to_common_frame(
    isotopomer_ravel *ravel_isotopomer,  // isotopomer structure
    float *transition,                   // The spin transition
    bool allow_fourth_rank,  // if true, prep for 4th rank computation
    double *R0,              // the R0 components
    complex128 *R2,          // the R2 components
    complex128 *R4,          // the R4 components
    double *R0_temp,         // the temporary R0 components
    complex128 *R2_temp,     // the temporary R2 components
    complex128 *R4_temp,     // the temporary R3 components
    double B0_in_T,          // magnetic flux density in T
    bool *freq_contrib       // the pointer to freq contribs boolean
) {
  /* The following codeblock populates the product of spatial part, Rlm, of
   * the tensor and the spin transition function, T(mf, mi) for
   *      zeroth rank, R0 = [ R00 ] * T(mf, mi)
   *      second rank, R2 = [ R2m ] * T(mf, mi) where m ∈ [-2, 2].
   *      fourth rank, R4 = [ R4m ] * T(mf, mi) where m ∈ [-4, 4].
   * Here, mf, mi are the spin quantum numbers of the final and initial
   * energy state of the spin transition. The term `Rlm` is the coefficient
   * of the irreducible spherical tensor of rank `l` and order `m`. For more
   * information, see reference
   *
   *   Symmetry pathways in solid-state NMR. PNMRS 2011 59(2):12 1-96.
   *   https://doi.org/10.1016/j.pnmrs.2010.11.003
   *
   */

  unsigned int site, n_sites = ravel_isotopomer->number_of_sites;
  double larmor_freq_in_MHz;
  float *mf = &transition[n_sites], *mi = transition;

  for (site = 0; site < n_sites; site++) {
    if (*mi == *mf) {
      mi++;
      mf++;
      continue;
    }
    larmor_freq_in_MHz = -B0_in_T * ravel_isotopomer->gyromagnetic_ratio[site];
    /* Nuclear shielding components ======================================== */
    /*  Upto the first order */
    FCF_1st_order_nuclear_shielding_tensor_components(
        R0_temp, R2_temp,
        ravel_isotopomer->isotropic_chemical_shift_in_ppm[site] *
            larmor_freq_in_MHz,
        ravel_isotopomer->shielding_symmetric_zeta_in_ppm[site] *
            larmor_freq_in_MHz,
        ravel_isotopomer->shielding_symmetric_eta[site],
        &ravel_isotopomer->shielding_orientation[3 * site], *mf, *mi);

    // in-place update the R0 and R2 components.
    if (*freq_contrib++) *R0 += *R0_temp;
    if (*freq_contrib++) {
      vm_double_add_inplace(10, (double *)R2_temp, (double *)R2);
    }
    /* ===================================================================== */

    /* Electric quadrupolar components ===================================== */
    if (ravel_isotopomer->spin[site] > 0.5) {
      /*  Upto the first order */
      if (*freq_contrib++) {
        FCF_1st_order_electric_quadrupole_tensor_components(
            R2_temp, ravel_isotopomer->spin[site],
            ravel_isotopomer->quadrupolar_Cq_in_Hz[site],
            ravel_isotopomer->quadrupolar_eta[site],
            &ravel_isotopomer->quadrupolar_orientation[3 * site], *mf, *mi);

        // in-place update the R2 components.
        vm_double_add_inplace(10, (double *)R2_temp, (double *)R2);
      }

      /*  Upto the second order */
      if (allow_fourth_rank) {
        FCF_2nd_order_electric_quadrupole_tensor_components(
            R0_temp, R2_temp, R4_temp, ravel_isotopomer->spin[site],
            larmor_freq_in_MHz * 1e6,
            ravel_isotopomer->quadrupolar_Cq_in_Hz[site],
            ravel_isotopomer->quadrupolar_eta[site],
            &ravel_isotopomer->quadrupolar_orientation[3 * site], *mf, *mi);

        // in-place update the R0 component.
        if (*freq_contrib++) *R0 += *R0_temp;

        // in-place update the R2 components.
        if (*freq_contrib++) {
          vm_double_add_inplace(10, (double *)R2_temp, (double *)R2);
        }
        // in-place update the R4 components.
        if (*freq_contrib++) {
          vm_double_add_inplace(18, (double *)R4_temp, (double *)R4);
        }
      }
    }

    if (n_sites > 1) {
      /* Weakly coupled direct-dipole components =========================== */
      /*      Upto the first order (to do..-> add orientation dependence)    */
      weakly_coupled_direct_dipole_frequencies_to_first_order(
          R0, R2_temp, ravel_isotopomer->dipolar_couplings[site], *mf, *mi, 0.5,
          0.5);

      // in-place update the R2 components.
      vm_double_add_inplace(10, (double *)R2_temp, (double *)R2);
      /* =================================================================== */
    }
    mi++;
    mf++;
  }
}

/**
 * @func __get_components_2
 *
 * The function calculates the following.
 *
 *   pre_phase(m, t) = I 2π [(exp(I m ωr t) - 1)/(I m ωr)]
 *                   = (2π / m ωr) (exp(I m ωr t) - 1)
 *                     |--scale--|
 *                   = scale (exp(I m ωr t) - 1)
 *                   = scale [[cos(m ωr t) -1] +Isin(m ωr t)],
 *
 * where ωr is the sample spinning frequency in Hz, m goes from -4 to 4,
 * t is a vector of length `number_of_sidebands` given as
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
void __get_components_2(unsigned int number_of_sidebands,
                        double sample_rotation_frequency_in_Hz,
                        complex128 *pre_phase) {
  int m, i;
  double spin_angular_freq, tau, scale;

  double *input = malloc_double(number_of_sidebands);
  double *ones = malloc_double(number_of_sidebands);
  double *phase = malloc_double(number_of_sidebands);

  vm_double_ones(number_of_sidebands, ones);
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
     * Evaluate pre_phase = scale * (cexp(I * phase) - 1.0)
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

    /**
     * The expression pre_phase[m] = scale * (cexp(I * phase) - 1.0) given
     * above for positive m is related to -m as
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
 * @func __get_components
 *
 * The function calculates the following.
 *   pre_phase(m, t) = I 2π [(exp(I m ωr t) - 1)/(I m ωr)]
 *                   = (2π / m ωr) (exp(I m ωr t) - 1)
 *                     |--scale--|
 *                   = scale * (exp(I m ωr t) - 1)
 * where ωr is the sample spinning frequency in Hz, m goes from -4 to 4, and
 * t is a vector of length `number_of_sidebands` given as
 *    t = [0, 1, ... number_of_sidebands-1]/(ωr*number_of_sidebands)
 * `pre_phase` is a matrix of shape, `9 x number_of_sidebands`, with
 * number_of_sidebands as the leading dimension. The first number_of_sidebands
 * entries corresponds to m_wr=-4.
 */
void __get_components(unsigned int number_of_sidebands,
                      double sample_rotation_frequency,
                      double *restrict pre_phase) {
  double spin_angular_freq, tau, wrt, pht, scale;
  unsigned int step, m;
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
