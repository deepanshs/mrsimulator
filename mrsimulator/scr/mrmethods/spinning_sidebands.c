
//
//  sideband_simulator.c
//
//  Created by Deepansh J. Srivastava, Apr 11, 2019
//  Copyright © 2019 Deepansh J. Srivastava. All rights reserved.
//  Contact email = srivastava.89@osu.edu, deepansh2012@gmail.com
//

#include "spinning_sidebands.h"

static inline void __zero_components(double *R0, double complex *R2,
                                     double complex *R4) {
  int i;
  R0[0] = 0.0;
  for (i = 0; i < 5; i++) {
    R2[i] = 0.0;
  }
  for (i = 0; i < 9; i++) {
    R4[i] = 0.0;
  }
}

static inline void __spinning_sideband_core(
    // spectrum information and related amplitude
    double *spec,             // amplitude vector representing the spectrum.
    double coordinate_offset, // starting coordinate of the dimension.
    double increment,         // increment of coordinates along the dimension.
    int count,                // number of coordinates along the dimension.

    isotopomer_ravel *ravel_isotopomer, // isotopomer structure

    bool ALLOW_FOURTH_RANK,           // Quad theory for second order,
    int remove_second_order_quad_iso, // remove the isotropic contribution from
                                      // the second order quad Hamiltonian.

    // spin rate, spin angle and number spinning sidebands
    int number_of_sidebands, // The number of spinning sidebands to evaluate
    double sample_rotation_frequency, // The rotor spin frequency
    double rotor_angle, // The rotor angle relative to lab-frame z-axis

    // Pointer to the transitions. transition[0] = mi and transition[1] = mf
    double *transition, mrsimulator_plan *plan) {

  /*
  The computation of the spinning sidebands is based on the method described by
  Eden and Levitt et. al.
    `Computation of Orientational Averages in Solid-State NMR by Gaussian
     Spherical Quadrature`
      JMR, 132, 1998. https://doi.org/10.1006/jmre.1998.1427
  */

  // Sampled over an octant
  int i, site;
  unsigned int orientation;

  double inverse_increment = 1.0 / increment;
  double iso_n_, zeta_n_, eta_n_, Cq_e_, eta_e_, d_, offset, scale;

  double R0 = 0.0;
  double complex *R2 = malloc_double_complex(5);
  double complex *R4 = malloc_double_complex(9);

  double complex one = 1.0, zero = 0.0;
  double shift_half_bin = 0.5;

  int spec_site;
  double *spec_site_ptr;
  double local_frequency_offset;

  // Per site base calculation
  for (site = 0; site < ravel_isotopomer->number_of_sites; site++) {
    // gettimeofday(&start_site_time, NULL);

    spec_site = site * count;
    spec_site_ptr = &spec[spec_site];

    /* Nuclear shielding terms                                               */
    iso_n_ = ravel_isotopomer->isotropic_chemical_shift_in_Hz[site];
    zeta_n_ = ravel_isotopomer->shielding_anisotropy_in_Hz[site];
    eta_n_ = ravel_isotopomer->shielding_asymmetry[site];

    /* Electric quadrupolar terms                                            */
    Cq_e_ = ravel_isotopomer->quadrupolar_constant_in_Hz[site];
    eta_e_ = ravel_isotopomer->quadrupolar_asymmetry[site];

    /* Magnetic dipole terms                                                 */
    d_ = ravel_isotopomer->dipolar_couplings[site];

    /* The following codeblock populates the spatial part of the tensor for   /
    /   zeroth rank, R0 = [ R00 ] * transition element                        /
    /   second rank, R2 = [ R2m ] * transition element where m ∈ [-2, 2].     /
    /   fourth rank, R4 = [ R4m ] * transition element where m ∈ [-4, 4].     /
    / Here, transition element is derived from symmetry pathways. See ref.    /
    /   Symmetry pathways in solid-state NMR. PNMRS 2011 59(2):12 1-96.       /
    /   https://doi.org/10.1016/j.pnmrs.2010.11.003                          */

    /* Initialize with zeroing all spatial components                        */
    __zero_components(&R0, R2, R4);

    /* get nuclear shielding components upto first order                     */
    get_nuclear_shielding_hamiltonian_to_first_order(&R0, R2, iso_n_, zeta_n_,
                                                     eta_n_, transition);

    /* get weakly coupled direct dipole components upto first order          */
    get_weakly_coupled_direct_dipole_hamiltonian_to_first_order(&R0, R2, d_,
                                                                transition);

    if (ravel_isotopomer->spin > 0.5) {
      /* get electric quadrupolar components upto first order                */
      get_quadrupole_hamiltonian_to_first_order(&R0, R2, ravel_isotopomer->spin,
                                                Cq_e_, eta_e_, transition);

      /* get electric quadrupolar components upto second order               */
      if (ALLOW_FOURTH_RANK) {
        get_quadrupole_hamiltonian_to_second_order(
            &R0, R2, R4, ravel_isotopomer->spin, Cq_e_, eta_e_, transition,
            ravel_isotopomer->larmor_frequency, remove_second_order_quad_iso);
      }
    }

    /* Calculating the local isotropic offset with respect to frequency bins */
    local_frequency_offset =
        shift_half_bin + (R0 - coordinate_offset) * inverse_increment;

    // ------------------------------------------------------------------- //
    // Equation [39] in the refernce https://doi.org/10.1006/jmre.1998.1427.
    //
    // w_cs^{m}(O_MR) = iso delta(m,0) + sum_{m', m" =-2}^{2} A[m"]
    // D^2_{m"m'}(O_PM) D^2_{m'm}(O_MR) d^2_{m'm}(b_RL)
    //

    /* Wigner second rank rotation */
    __wigner_rotation_2(2, plan->n_orientations, plan->wigner_2j_matrices,
                        plan->exp_Im_alpha, R2, plan->w2);

    /* Evaluating the exponent of the sideband phase w.r.t fourth rank       *
     * tensor given as,                                                      *
     * w_cs(Θ) * I 2π [(exp(I m wr t) - 1)/(I m wr)] * d^2m0(rotor_angle)    *
     *           |------------------- pre_phase_2 ----------------------|    *
     * The result, `vector` is interpreted as a matrix of dimensions         *
     * `number_of_sidebands` x `n_orientations` with `n_orientations` as the *
     * leading dimension.                                                    */
    cblas_zgemm(CblasRowMajor, CblasTrans, CblasTrans,
                plan->number_of_sidebands, plan->n_orientations, 5, &one,
                plan->pre_phase_2, plan->number_of_sidebands, plan->w2, 5,
                &zero, plan->vector, plan->n_orientations);

    if (ALLOW_FOURTH_RANK) {
      /* Wigner fourth rank rotation */
      __wigner_rotation_2(4, plan->n_orientations, plan->wigner_4j_matrices,
                          plan->exp_Im_alpha, R4, plan->w4);

      /* Evaluating the exponent of the sideband phase w.r.t fourth rank     *
       *  tensor given as,                                                   *
       * w_cs(Θ) * I 2π [(exp(I m wr t) - 1)/(I m wr)] * d^4m0(rotor_angle)  *
       *            |------------------- pre_phase_4 ----------------------| *
       * The result, `vector` is interpreted as a matrix of dimensions       *
       * `number_of_sidebands` x `n_orientations` with `n_orientations` as   *
       * the leading dimension.                                              */
      cblas_zgemm(CblasRowMajor, CblasTrans, CblasTrans,
                  plan->number_of_sidebands, plan->n_orientations, 9, &one,
                  plan->pre_phase_4, plan->number_of_sidebands, plan->w4, 9,
                  &one, plan->vector, plan->n_orientations);
    }

    /* Evaluating the sideband phase exp(vector)                             */
    vmzExp(plan->size, plan->vector, plan->vector, VML_EP);

    /* Evaluate the Fourier transform of vector, fft(vector).                *
     * The fft operation updates the value of the array, `vector`            */
    fftw_execute(plan->the_fftw_plan);
    // DftiComputeForward(plan, vector);

    /* Taking the absolute value square of the vector array. The absolute    *
     *  value square is stores as the real part of the vector array. The     *
     *  imaginary part is garbage. This method avoids creating new arrays.   */
    vmdSqr(2 * plan->size, (double *)plan->vector, (double *)plan->vector,
           VML_EP);
    cblas_daxpby(plan->size, 1.0, (double *)plan->vector + 1, 2, 1.0,
                 (double *)plan->vector, 2);

    /* Scaling the absolute square value with the power scheme weights.      *
     *  Only the real part is scaled and the imaginary part is discarded.    */
    for (orientation = 0; orientation < plan->n_orientations; orientation++) {
      cblas_dscal(number_of_sidebands, plan->amplitudes[orientation],
                  (double *)&plan->vector[orientation],
                  2 * plan->n_orientations);
    }

    /* Calculating local anisotropic frequency contributions                 *
     * contribution from the second rank tensor.                             */
    scale = inverse_increment * plan->rotor_lab_2[2];
    cblas_daxpby(plan->n_orientations, scale, (double *)&plan->w2[2], 10, 0.0,
                 plan->local_frequency, 1);

    if (ALLOW_FOURTH_RANK) {
      /* contribution from the fourth rank tensor.                           */
      scale = inverse_increment * plan->rotor_lab_4[4];
      cblas_daxpby(plan->n_orientations, scale, (double *)&plan->w4[4], 18, 1.0,
                   plan->local_frequency, 1);
    }

    // ---------------------------------------------------------------------
    //              Calculating the tent for every sideband
    // Allowing only sidebands that are within the spectral bandwidth
    //
    for (i = 0; i < number_of_sidebands; i++) {
      offset = plan->vr_freq[i] + local_frequency_offset;
      if ((int)offset >= 0 && (int)offset <= count) {

        vdLinearFrac(plan->n_orientations, plan->local_frequency,
                     plan->local_frequency, 1.0, offset, 0.0, 1.0,
                     plan->freq_offset);
        octahedronInterpolation(
            spec_site_ptr, plan->freq_offset,
            plan->geodesic_polyhedron_frequency,
            (double *)&plan->vector[i * plan->n_orientations], 2, count);
      }
    }
    // gettimeofday(&end_site_time, NULL);
    // clock_time =
    //     (double)(end_site_time.tv_usec - start_site_time.tv_usec) / 1000000.
    //     + (double)(end_site_time.tv_sec - start_site_time.tv_sec);
    // printf("Total time per site %f \n", clock_time);
  }
}

void spinning_sideband_core(
    // spectrum information and related amplitude
    double *spec,             // The amplitude of the spectrum.
    double *cpu_time_,        // Execution time
    double coordinate_offset, // The start of the frequency spectrum.
    double increment,         // The increment of the frequency spectrum.
    int count,                // Number of points on the frequency spectrum.

    isotopomer_ravel *ravel_isotopomer, // isotopomer structure

    int quadSecondOrder,              // Quad theory for second order,
    int remove_second_order_quad_iso, // remove the isotropic contribution
                                      // from the second order quad
                                      // Hamiltonian.

    // spin rate, spin angle and number spinning sidebands
    int number_of_sidebands,          // The number of sidebands
    double sample_rotation_frequency, // The rotor spin frequency
    double rotor_angle, // The rotor angle relative to lab-frame z-axis

    // Pointer to the transitions. transition[0] = mi and transition[1] = mf
    double *transition,

    // powder orientation average
    int geodesic_polyhedron_frequency // The number of triangle along the edge
                                      // of octahedron
) {
  bool ALLOW_FOURTH_RANK = false;
  if (ravel_isotopomer[0].spin > 0.5 && quadSecondOrder == 1) {
    ALLOW_FOURTH_RANK = true;
  }
  // mkl_set_threading_layer(MKL_THREADING_INTEL);
  int max_threads = mkl_get_max_threads();
  mkl_set_num_threads(max_threads);
  printf("Using %d threads\n", max_threads);

  // check for spinning speed
  if (sample_rotation_frequency < 1.0e-3) {
    sample_rotation_frequency = 1.0e9;
    rotor_angle = 0.0;
    number_of_sidebands = 1;
  }

  struct timeval begin, end; // all_site_time, all_c_time;
  double clock_time;

  gettimeofday(&begin, NULL);
  mrsimulator_plan *plan = mrsimulator_create_plan(
      geodesic_polyhedron_frequency, number_of_sidebands,
      sample_rotation_frequency, rotor_angle, increment, ALLOW_FOURTH_RANK);
  gettimeofday(&end, NULL);
  clock_time = (double)(end.tv_usec - begin.tv_usec) / 1000000. +
               (double)(end.tv_sec - begin.tv_sec);
  printf("mrsimulator plan time %f \n", clock_time);

  // gettimeofday(&all_site_time, NULL);
  __spinning_sideband_core(
      // spectrum information and related amplitude
      spec,              // The amplitude of the spectrum.
      coordinate_offset, // The start of the frequency spectrum.
      increment,         // The increment of the frequency spectrum.
      count,             // Number of points on the frequency spectrum.

      ravel_isotopomer, // isotopomer structure

      ALLOW_FOURTH_RANK,            // Quad theory for second order,
      remove_second_order_quad_iso, // remove the isotropic contribution
                                    // from the second order quad
                                    // Hamiltonian.

      // spin rate, spin angle and number spinning sidebands
      number_of_sidebands,       // The number of spinning sidebands to evaluate
      sample_rotation_frequency, // The rotor spin frequency
      rotor_angle,               // The rotor angle relative to lab-frame z-axis

      // Pointer to the transitions. transition[0] = mi and transition[1] = mf
      transition,

      plan);

  // gettimeofday(&all_site_time_end, NULL);
  // clock_time =
  //     (double)(all_site_time_end.tv_usec - all_site_time.tv_usec) / 1000000.
  //     + (double)(all_site_time_end.tv_sec - all_site_time.tv_sec);
  // printf("all site time %f \n", clock_time);

  // gettimeofday(&all_c_time_end, NULL);
  // clock_time =
  //     (double)(all_c_time_end.tv_usec - all_c_time.tv_usec) / 1000000. +
  //     (double)(all_c_time_end.tv_sec - all_c_time.tv_sec);
  // printf("all c time %f \n", clock_time);

  // gettimeofday(&end, NULL);
  // clock_time = (double)(end.tv_usec - begin.tv_usec) / 1000000. +
  //              (double)(end.tv_sec - begin.tv_sec);
  // printf("time %f s\n", clock_time);
  // cpu_time_[0] += clock_time;

  mrsimulator_free_plan(plan); /* clean up */
}
