
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
    double *spec, // amplitude vector representing the spectrum.

    isotopomer_ravel *ravel_isotopomer, // isotopomer structure

    int remove_second_order_quad_iso, // remove the isotropic contribution from
                                      // the second order quad Hamiltonian.

    // Pointer to the transitions. transition[0] = mi and transition[1] = mf
    double *transition, MRS_plan *plan, MRS_dimension *dimension) {

  /*
  The computation of the spinning sidebands is based on the method described by
  Eden and Levitt et. al.
    `Computation of Orientational Averages in Solid-State NMR by Gaussian
     Spherical Quadrature`
      JMR, 132, 1998. https://doi.org/10.1006/jmre.1998.1427
  */

  int i, site;
  double iso_n_, zeta_n_, eta_n_, Cq_e_, eta_e_, d_, offset;
  double R0 = 0.0;
  double complex *R2 = malloc_double_complex(5);
  double complex *R4 = malloc_double_complex(9);

  int spec_site;
  double *spec_site_ptr;

  // Per site base calculation
  for (site = 0; site < ravel_isotopomer->number_of_sites; site++) {
    // gettimeofday(&start_site_time, NULL);

    spec_site = site * dimension->count;
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

    /* The following codeblock populates the product of spatial part, Rlm, of
     * the tensor and the spin transition function, T(mf, mi) for
     *      zeroth rank, R0 = [ R00 ] * T(mf, mi)
     *      second rank, R2 = [ R2m ] * T(mf, mi) where m ∈ [-2, 2].
     *      fourth rank, R4 = [ R4m ] * T(mf, mi) where m ∈ [-4, 4].
     * Here, mf, mi are the spin quantum numbers of the final and initial
     * energy state of the spin transition. The term `Rlm` is the coefficient
     * of the irreducible spherical tensor of rank `l` and order `m`. For more
     * information, see reference
     *   Symmetry pathways in solid-state NMR. PNMRS 2011 59(2):12 1-96.
     *   https://doi.org/10.1016/j.pnmrs.2010.11.003                         */

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
      if (plan->ALLOW_FOURTH_RANK) {
        get_quadrupole_hamiltonian_to_second_order(
            &R0, R2, R4, ravel_isotopomer->spin, Cq_e_, eta_e_, transition,
            ravel_isotopomer->larmor_frequency, remove_second_order_quad_iso);
      }
    }

    /*  */
    MRS_get_amplitudes_from_plan(plan, R2, R4);
    MRS_get_normalized_frequencies_from_plan(plan, dimension, R0, R2, R4);

    // ---------------------------------------------------------------------
    //              Calculating the tent for every sideband
    // Allowing only sidebands that are within the spectral bandwidth
    //
    for (i = 0; i < plan->number_of_sidebands; i++) {
      offset = plan->vr_freq[i] + plan->isotropic_offset;
      if ((int)offset >= 0 && (int)offset <= dimension->count) {

        vdLinearFrac(plan->n_orientations, plan->local_frequency,
                     plan->local_frequency, 1.0, offset, 0.0, 1.0,
                     plan->freq_offset);
        octahedronInterpolation(
            spec_site_ptr, plan->freq_offset,
            plan->geodesic_polyhedron_frequency,
            (double *)&plan->vector[i * plan->n_orientations], 2,
            dimension->count);
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
    double *spec,              // The amplitude of the spectrum.
    double *cpu_time_,         // Execution time
    double coordinates_offset, // The start of the frequency spectrum.
    double increment,          // The increment of the frequency spectrum.
    int count,                 // Number of points on the frequency spectrum.

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
  printf("Using upto %d threads for simulation.\n", max_threads);

  // check for spinning speed
  if (sample_rotation_frequency < 1.0e-3) {
    sample_rotation_frequency = 1.0e9;
    rotor_angle = 0.0;
    number_of_sidebands = 1;
  }

  MRS_dimension *dimension =
      MRS_create_dimension(count, coordinates_offset, increment);

  struct timeval begin, end; // all_site_time, all_c_time;
  double clock_time;

  gettimeofday(&begin, NULL);
  MRS_plan *plan = MRS_create_plan(
      geodesic_polyhedron_frequency, number_of_sidebands,
      sample_rotation_frequency, rotor_angle, increment, ALLOW_FOURTH_RANK);
  gettimeofday(&end, NULL);
  clock_time = (double)(end.tv_usec - begin.tv_usec) / 1000000. +
               (double)(end.tv_sec - begin.tv_sec);
  printf("Created mrsimulator plan in %f s.\n", clock_time);

  // gettimeofday(&all_site_time, NULL);
  __spinning_sideband_core(
      // spectrum information and related amplitude
      spec, // The amplitude of the spectrum.

      ravel_isotopomer, // isotopomer structure

      remove_second_order_quad_iso, // remove the isotropic contribution
                                    // from the second order quad
                                    // Hamiltonian.

      // Pointer to the transitions. transition[0] = mi and transition[1] = mf
      transition,

      plan, dimension);

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

  MRS_free_plan(plan); /* clean up */
}
