// -*- coding: utf-8 -*-
//
//  schemes.c
//
//  @copyright Deepansh J. Srivastava, 2019-2021.
//  Created by Deepansh J. Srivastava.
//  Contact email = srivastava.89@osu.edu
//

#include "schemes.h"

static inline void averaging_scheme_setup(MRS_averaging_scheme *scheme,
                                          complex128 *exp_I_beta, bool allow_4th_rank) {
  unsigned int allocate_size_2, allocate_size_4;
  double delta_alpha = 0.0;
  scheme->total_orientations = scheme->octant_orientations;

  switch (scheme->integration_volume) {
  case 0:  // positive octant
    break;
  case 1:  // positive hemisphere
    scheme->total_orientations *= 4;
    break;
  case 2:  // full sphere
    scheme->total_orientations *= 8;
    break;
  }

  /**
   * Calculating exp(-Imα) at every orientation angle α form m=-4 to -1, where α is the
   * azimuthal angle over the positive octant. The input α is in the form of phase
   * exp(Iα). */
  if (scheme->integration_volume == 2) {
    delta_alpha = CONST_PI / (4.0 * scheme->integration_density);
  }
  get_exp_Im_angle(scheme->octant_orientations, allow_4th_rank, scheme->exp_Im_alpha,
                   delta_alpha);
  /* -------------------------------------------------------------------------------- */

  /**
   * Wigner matrices corresponding to the upper hemisphere.
   *
   * The wigner matrices are evaluated at every β orientation over the positive upper
   * octant. Note, the β angles from this octant repeat for the other three octant in
   * the upper hemisphere, therfore, only one set of second and fourth rank wigner
   * matrices should suffice.
   */

  // calculating the required space for storing wigner matrices.
  allocate_size_2 = 15 * scheme->octant_orientations;  // (5 x 3) half-matrix
  allocate_size_4 = 45 * scheme->octant_orientations;  // (9 x 5) half-matrix
  if (scheme->integration_volume == 2) {
    allocate_size_2 *= 2;
    allocate_size_4 *= 2;
  }

  /* Second-rank reduced wigner matrices at every β orientation from the positive upper
   * octant. */
  scheme->wigner_2j_matrices = malloc_double(allocate_size_2);
  wigner_d_matrices_from_exp_I_beta(2, scheme->octant_orientations, true, exp_I_beta,
                                    scheme->wigner_2j_matrices);

  scheme->wigner_4j_matrices = NULL;
  if (allow_4th_rank) {
    /* Fourt-rank reduced wigner matrices at every β orientation from the positive upper
     * octant. */
    scheme->wigner_4j_matrices = malloc_double(allocate_size_4);
    wigner_d_matrices_from_exp_I_beta(4, scheme->octant_orientations, true, exp_I_beta,
                                      scheme->wigner_4j_matrices);
  }

  /**
   * If averaging over a sphere is selected, then calculate the wigner matrices
   * corresponding to the lower hemisphere.
   *
   * The wigner matrices only dependents on the β angles. Going from upper to the lower
   * hemisphere, β -> β+π/2. This implies,
   *      cos(β+π/2) -> -cos(β)
   *      sin(β+π/2) -> sin(β)
   *
   * For evaluating the reduced wigner matrices from the lower hemisphere, the sign of
   * cosine beta is changed. As before, the β angles from any octant from the lower
   * hemisphere repeat for the other three octant in the lower hemisphere, therfore,
   * only one set of second rank and fourth rank reduced  wigner matrices should
   * suffice. */
  if (scheme->integration_volume == 2) {
    allocate_size_2 /= 2;
    allocate_size_4 /= 2;
    /* cos(beta) is negative in the lower hemisphere */
    cblas_dscal(scheme->octant_orientations, -1.0, (double *)exp_I_beta, 2);

    /* Second-rank reduced wigner matrices at every β orientation over an octant from
     * the lower hemisphere */
    wigner_d_matrices_from_exp_I_beta(2, scheme->octant_orientations, true, exp_I_beta,
                                      &scheme->wigner_2j_matrices[allocate_size_2]);
    if (allow_4th_rank) {
      /* Fourth-rank reduced wigner matrices at every β orientation. */
      wigner_d_matrices_from_exp_I_beta(4, scheme->octant_orientations, true,
                                        exp_I_beta,
                                        &scheme->wigner_4j_matrices[allocate_size_4]);
    }
  }
  /* -------------------------------------------------------------------------------- */

  /* Setting up buffers and tables for processing the second rank tensors. . */
  /* ................................................................................ */
  /* w2 is the buffer for storing the frequencies calculated from the second-rank
   * tensors. Only calcuate the -2, -1, and 0 tensor components.*/
  scheme->w2 = malloc_complex128(3 * scheme->total_orientations);

  scheme->w4 = NULL;
  if (allow_4th_rank) {
    /* w4 is the buffer for storing the frequencies calculated from the fourth-rank
     * tensors. Only calcuate the -4, -3, -2, -1, and 0 tensor components.*/
    scheme->w4 = malloc_complex128(5 * scheme->total_orientations);
  }
}

/* Free the memory from the mrsimulator plan associated with the spherical averaging
 * scheme */
void MRS_free_averaging_scheme(MRS_averaging_scheme *scheme) {
  free(scheme->amplitudes);
  free(scheme->exp_Im_alpha);
  free(scheme->exp_Im_gamma);
  free(scheme->w2);
  free(scheme->w4);
  free(scheme->wigner_2j_matrices);
  free(scheme->wigner_4j_matrices);
  free(scheme->scratch);
  free(scheme->amps_real);
  free(scheme->amps_imag);
  free(scheme->phase);
  free(scheme->exp_I_phase);
  free(scheme);
}

/* Create a new orientation averaging scheme. */
MRS_averaging_scheme *MRS_create_averaging_scheme(unsigned int integration_density,
                                                  bool allow_4th_rank,
                                                  unsigned int n_gamma,
                                                  unsigned int integration_volume,
                                                  bool interpolation, bool is_complex) {
  int alpha_size;
  MRS_averaging_scheme *scheme = malloc(sizeof(MRS_averaging_scheme));

  scheme->interpolation = interpolation;
  scheme->user_defined = false;
  scheme->positions = NULL;
  scheme->position_size = 0;
  scheme->n_gamma = n_gamma;
  scheme->integration_density = integration_density;
  scheme->integration_volume = integration_volume;
  scheme->allow_4th_rank = allow_4th_rank;
  scheme->is_complex = is_complex;

  scheme->octant_orientations =
      ((integration_density + 1) * (integration_density + 2)) / 2;

  /* Calculate α, β, and weights over the positive octant. .......................... */
  /* ................................................................................ */
  // The 4 * octant_orientations memory allocation is for m=4, 3, 2, and 1
  alpha_size = 4 * scheme->octant_orientations;
  alpha_size *= (integration_volume == 2) ? 2 : 1;
  scheme->exp_Im_alpha = malloc_complex128(alpha_size);
  scheme->exp_Im_gamma = malloc_complex128(4 * n_gamma);
  complex128 *exp_I_beta = malloc_complex128(scheme->octant_orientations);
  scheme->amplitudes = malloc_double(scheme->octant_orientations);

  averaging_setup(integration_density,
                  &scheme->exp_Im_alpha[3 * scheme->octant_orientations], exp_I_beta,
                  scheme->amplitudes, interpolation);

  averaging_scheme_setup(scheme, exp_I_beta, allow_4th_rank);

  // exp(-im gamma) for m=[-4,-1], and gamma=[0..8]*2pi/9
  double *gamma = malloc_double(n_gamma);
  double *temp = malloc_double(n_gamma);
  double factor = CONST_2PI / (double)n_gamma;
  vm_double_arange(n_gamma, gamma);
  cblas_dscal(n_gamma, factor, gamma, 1);
  vm_cosine_I_sine(n_gamma, gamma, &scheme->exp_Im_gamma[3 * n_gamma]);
  get_exp_Im_angle(n_gamma, allow_4th_rank, scheme->exp_Im_gamma, 0.0);  // for gamma

  // reallocate exp_I_beta memory as scratch.
  scheme->scratch = (double *)exp_I_beta;
  scheme->amps_real = malloc_double(scheme->total_orientations);
  scheme->amps_imag = malloc_double(scheme->total_orientations);
  scheme->phase = malloc_double(scheme->n_gamma * scheme->total_orientations);
  scheme->exp_I_phase = malloc_complex128(scheme->n_gamma * scheme->total_orientations);

  free(gamma);
  free(temp);
  return scheme;
}

/* Create a new orientation averaging scheme. */
MRS_averaging_scheme *MRS_create_averaging_scheme_from_alpha_beta(
    double *alpha, double *beta, double *weight, unsigned int n_angles,
    bool allow_4th_rank, unsigned int n_gamma, const unsigned int position_size,
    int32_t *positions, bool interpolation, bool is_complex) {
  double scale;
  MRS_averaging_scheme *scheme = malloc(sizeof(MRS_averaging_scheme));

  scheme->interpolation = interpolation;
  scheme->user_defined = true;
  scheme->positions = positions;
  scheme->position_size = position_size;
  scheme->n_gamma = n_gamma;
  scheme->integration_density = 0;
  scheme->integration_volume = 0;
  scheme->allow_4th_rank = allow_4th_rank;
  scheme->is_complex = is_complex;

  scheme->octant_orientations = n_angles;

  /* Calculate α, β, and weights over the positive octant. .......................... */
  /* ................................................................................ */
  // The 4 * octant_orientations memory allocation is for m=4, 3, 2, and 1
  scheme->exp_Im_alpha = malloc_complex128(4 * scheme->octant_orientations);
  scheme->exp_Im_gamma = malloc_complex128(4 * n_gamma);
  complex128 *exp_I_beta = malloc_complex128(scheme->octant_orientations);
  scheme->amplitudes = malloc_double(scheme->octant_orientations);

  // scale = (1 / 6) factor for 1 point to 6 triangle ratio.
  scale = (interpolation) ? 0.1666666667 : 1.0;
  vm_double_ramp(scheme->octant_orientations, weight, scale, 0.0, scheme->amplitudes);

  /* Calculate cos(α) + isin(α) from α. ............................................. */
  vm_cosine_I_sine(n_angles, alpha,
                   scheme->exp_Im_alpha[3 * scheme->octant_orientations]);

  /* Calculate cos(β) + isin(β) from β. ............................................. */
  vm_cosine_I_sine(n_angles, beta, exp_I_beta);

  averaging_scheme_setup(scheme, exp_I_beta, allow_4th_rank);

  // exp(-im gamma) for m=[-4,-1], and gamma=[0..8]*2pi/9
  double *gamma = malloc_double(n_gamma);
  double *temp = malloc_double(n_gamma);
  double factor = CONST_2PI / (double)n_gamma;
  vm_double_arange(n_gamma, gamma);
  cblas_dscal(n_gamma, factor, gamma, 1);
  vm_cosine_I_sine(n_gamma, gamma, &scheme->exp_Im_gamma[3 * n_gamma]);
  get_exp_Im_angle(n_gamma, allow_4th_rank, scheme->exp_Im_gamma, 0.0);  // for gamma

  // reallocate exp_I_beta memory as scratch.
  scheme->scratch = (double *)exp_I_beta;
  scheme->amps_real = malloc_double(scheme->total_orientations);
  scheme->amps_imag = malloc_double(scheme->total_orientations);
  scheme->phase = malloc_double(scheme->n_gamma * scheme->total_orientations);
  scheme->exp_I_phase = malloc_complex128(scheme->n_gamma * scheme->total_orientations);

  free(gamma);
  free(temp);
  return scheme;
}

/* ---------------------------------------------------------------------------------- */
/* fftw routine setup ............................................................... */
/* .................................................................................. */
MRS_fftw_scheme *create_fftw_scheme(unsigned int total_orientations,
                                    unsigned int number_of_sidebands) {
  unsigned int size = total_orientations * number_of_sidebands;
  int nssb = (int)number_of_sidebands;
  MRS_fftw_scheme *fftw_scheme = malloc(sizeof(MRS_fftw_scheme));

  // fftw_scheme->vector = fftw_alloc_complex(size);
  fftw_scheme->vector = (fftw_complex *)fftw_malloc(sizeof(fftw_complex) * size);
  // malloc_complex128(plan->size);
  // gettimeofday(&fft_setup_time, NULL);

  // int fftw_thread = fftw_init_threads();
  // if (fftw_thread == 0) {
  //   printf("failed to initialize fftw threading");
  // }
  // fftw_plan_with_nthreads(2);

  fftw_scheme->the_fftw_plan = fftw_plan_many_dft(
      1, &nssb, total_orientations, fftw_scheme->vector, NULL, total_orientations, 1,
      fftw_scheme->vector, NULL, total_orientations, 1, FFTW_FORWARD, FFTW_ESTIMATE);

  // char *filename = "128_sidebands.wisdom";
  // int status = fftw_export_wisdom_to_filename(filename);
  // printf("file save status %i \n", status);
  /* ----------------------------------------------------------------------- */
  return fftw_scheme;
}

void MRS_free_fftw_scheme(MRS_fftw_scheme *fftw_scheme) {
  fftw_destroy_plan(fftw_scheme->the_fftw_plan);
  fftw_free(fftw_scheme->vector);
  // fftw_cleanup();
  free(fftw_scheme);
}
