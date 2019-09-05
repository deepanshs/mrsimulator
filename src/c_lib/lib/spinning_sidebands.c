
//
//  sideband_simulator.c
//
//  Created by Deepansh J. Srivastava, Apr 11, 2019
//  Copyright © 2019 Deepansh J. Srivastava. All rights reserved.
//  Contact email = srivastava.89@osu.edu, deepansh2012@gmail.com
//

#include "spinning_sidebands.h"

static inline void __zero_components(double *R0, complex128 *R2,
                                     complex128 *R4) {
  R0[0] = 0.0;
  vm_double_zeros(10, (double *)R2);
  vm_double_zeros(18, (double *)R4);
}

static inline void __spinning_sideband_core(
    // spectrum information and related amplitude
    double *spec, // amplitude vector representing the spectrum.

    isotopomer_ravel *ravel_isotopomer, // isotopomer structure

    int remove_second_order_quad_isotropic, // remove the isotropic contribution
                                            // from the second order quad
                                            // Hamiltonian.

    // Pointer to the transitions. transition[0] = mi and transition[1] = mf
    double *transition, MRS_plan *plan, MRS_dimension *dimension) {

  /*
  The computation of the spinning sidebands is based on the method described by
  Eden and Levitt et. al.
    `Computation of Orientational Averages in Solid-State NMR by Gaussian
     Spherical Quadrature`
      JMR, 132, 1998. https://doi.org/10.1006/jmre.1998.1427
  */

  int i;
  unsigned int site, j;
  double iso_n_, zeta_n_, eta_n_, Cq_e_, eta_e_, d_, offset;
  double *shielding_orientation, *quadrupole_orientation;

  double R0 = 0.0;
  double R0_temp = 0.0;
  complex128 *R2_temp = malloc_complex128(5);
  complex128 *R4_temp = malloc_complex128(9);

  complex128 *R2 = malloc_complex128(5);
  complex128 *R4 = malloc_complex128(9);

  int spec_site;
  double *spec_site_ptr;
  MRS_averaging_scheme *scheme = plan->averaging_scheme;

  // Per site base calculation
  for (site = 0; site < ravel_isotopomer->number_of_sites; site++) {
    // gettimeofday(&start_site_time, NULL);

    spec_site = site * dimension->count;
    spec_site_ptr = &spec[spec_site];

    /* Nuclear shielding terms                                               */
    iso_n_ = ravel_isotopomer->isotropic_chemical_shift_in_Hz[site];
    zeta_n_ = ravel_isotopomer->shielding_anisotropy_in_Hz[site];
    eta_n_ = ravel_isotopomer->shielding_asymmetry[site];
    shielding_orientation = &ravel_isotopomer->shielding_orientation[3 * site];

    /* Electric quadrupole terms                                            */
    Cq_e_ = ravel_isotopomer->quadrupole_coupling_constant_in_Hz[site];
    eta_e_ = ravel_isotopomer->quadrupole_asymmetry[site];
    quadrupole_orientation =
        &ravel_isotopomer->quadrupole_orientation[3 * site];

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

    /* get nuclear shielding components upto first order ................... */
    FCF_1st_order_nuclear_shielding_Hamiltonian(
        &R0_temp, R2_temp, iso_n_, zeta_n_, eta_n_, shielding_orientation,
        transition);
    R0 += R0_temp;
    vm_double_add_inplace(10, (double *)R2_temp, (double *)R2);

    /* get weakly coupled direct dipole components upto first order ........ */
    weakly_coupled_direct_dipole_frequencies_to_first_order(&R0, R2_temp, d_,
                                                            transition);
    vm_double_add_inplace(10, (double *)R2_temp, (double *)R2);
    // add orientation dependence

    if (ravel_isotopomer->spin > 0.5) {
      /* get electric quadrupole frequency tensors upto first order .............. */
      FCF_1st_order_electric_quadrupole_Hamiltonian(
          R2_temp, ravel_isotopomer->spin, Cq_e_, eta_e_,
          quadrupole_orientation, transition);
      vm_double_add_inplace(10, (double *)R2_temp, (double *)R2);

      /* get electric quadrupole frequency tensors upto second order ............. */
      if (plan->allow_fourth_rank) {
        FCF_2nd_order_electric_quadrupole_Hamiltonian(
            &R0_temp, R2_temp, R4_temp, ravel_isotopomer->spin,
            ravel_isotopomer->larmor_frequency, Cq_e_, eta_e_,
            quadrupole_orientation, transition);
        if (remove_second_order_quad_isotropic == 0) {
          R0 += R0_temp;
        }

        vm_double_add_inplace(10, (double *)R2_temp, (double *)R2);
        vm_double_add_inplace(18, (double *)R4_temp, (double *)R4);
      }
    }

    /* Get frequencies and amplitudes per octant ........................... */
    MRS_get_amplitudes_from_plan(plan, R2, R4);
    MRS_get_normalized_frequencies_from_plan(plan, dimension, R0);

    // ---------------------------------------------------------------------
    //              Calculating the tent for every sideband
    // Allowing only sidebands that are within the spectral bandwidth
    //
    // for (i = 0; i < plan->number_of_sidebands; i++) {
    //   offset = plan->vr_freq[i] + plan->isotropic_offset;
    //   if ((int)offset >= 0 && (int)offset <= dimension->count) {

    //     vm_double_ramp(plan->octant_orientations, plan->local_frequency, 1.0,
    //                    offset, plan->freq_offset);
    //     octahedronInterpolation(
    //         spec_site_ptr, plan->freq_offset,
    //         plan->geodesic_polyhedron_frequency,
    //         (double *)&plan->vector[i * plan->octant_orientations], 2,
    //         dimension->count);
    //   }
    // }

    unsigned int step_vector = 0, address;
    for (i = 0; i < plan->number_of_sidebands; i++) {
      offset = plan->vr_freq[i] + plan->isotropic_offset;
      if ((int)offset >= 0 && (int)offset <= dimension->count) {
        step_vector = i * scheme->total_orientations;
        for (j = 0; j < plan->n_octants; j++) {
          address = j * scheme->octant_orientations;

          vm_double_ramp(scheme->octant_orientations,
                         &scheme->local_frequency[address], 1.0, offset,
                         scheme->freq_offset);
          octahedronInterpolation(spec_site_ptr, scheme->freq_offset,
                                  scheme->geodesic_polyhedron_frequency,
                                  (double *)&plan->vector[step_vector], 2,
                                  dimension->count);
          step_vector += scheme->octant_orientations;
        }
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
    double coordinates_offset, // The start of the frequency spectrum.
    double increment,          // The increment of the frequency spectrum.
    int count,                 // Number of points on the frequency spectrum.
    isotopomer_ravel *ravel_isotopomer,     // Isotopomer structure
    int quad_second_order,                  // Quad theory for second order,
    int remove_second_order_quad_isotropic, // remove the isotropic contribution
                                            // from the second order quad
                                            // interaction.

    // spin rate, spin angle and number spinning sidebands
    int number_of_sidebands,                // The number of sidebands
    double sample_rotation_frequency_in_Hz, // The rotor spin frequency
    double rotor_angle_in_rad, // The rotor angle relative to lab-frame z-axis

    // Pointer to the transitions. transition[0] = mi and transition[1] = mf
    double *transition,

    // powder orientation average
    int geodesic_polyhedron_frequency // The number of triangle along the edge
                                      // of octahedron
) {

  // int num_process = openblas_get_num_procs();
  // int num_threads = openblas_get_num_threads();
  // // openblas_set_num_threads(1);
  // printf("%d processors", num_process);
  // printf("%d threads", num_threads);
  // int parallel = openblas_get_parallel();
  // printf("%d parallel", parallel);

  bool allow_fourth_rank = false;
  if (ravel_isotopomer[0].spin > 0.5 && quad_second_order == 1) {
    allow_fourth_rank = true;
  }

  // check for spinning speed
  if (sample_rotation_frequency_in_Hz < 1.0e-3) {
    sample_rotation_frequency_in_Hz = 1.0e9;
    rotor_angle_in_rad = 0.0;
    number_of_sidebands = 1;
  }

  MRS_averaging_scheme *scheme = MRS_create_averaging_scheme(
      geodesic_polyhedron_frequency, allow_fourth_rank);

  MRS_dimension *dimension =
      MRS_create_dimension(count, coordinates_offset, increment);

  // gettimeofday(&begin, NULL);
  MRS_plan *plan =
      MRS_create_plan(scheme, number_of_sidebands,
                      sample_rotation_frequency_in_Hz, rotor_angle_in_rad,
                      increment, allow_fourth_rank);

  // gettimeofday(&all_site_time, NULL);
  __spinning_sideband_core(
      // spectrum information and related amplitude
      spec, // The amplitude of the spectrum.

      ravel_isotopomer, // isotopomer structure

      remove_second_order_quad_isotropic, // remove the isotropic contribution
                                          // from the second order quad
                                          // Hamiltonian.

      // Pointer to the transitions. transition[0] = mi and transition[1] = mf
      transition,

      plan, dimension);

  // gettimeofday(&end, NULL);
  // clock_time = (double)(end.tv_usec - begin.tv_usec) / 1000000. +
  //              (double)(end.tv_sec - begin.tv_sec);
  // printf("time %f s\n", clock_time);
  // cpu_time_[0] += clock_time;

  MRS_free_plan(plan); /* clean up */
}
