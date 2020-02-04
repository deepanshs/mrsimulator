// -*- coding: utf-8 -*-
//
//  sideband_simulator.c
//
//  @copyright Deepansh J. Srivastava, 2019-2020.
//  Created by Deepansh J. Srivastava, Apr 11, 2019
//  Contact email = deepansh2012@gmail.com
//

#include "simulation.h"

static inline void __zero_components(double *R0, complex128 *R2,
                                     complex128 *R4) {
  R0[0] = 0.0;
  vm_double_zeros(10, (double *)R2);
  vm_double_zeros(18, (double *)R4);
}

/**
 * @func __MRS_rotate_components_from_PAS_to_common_frame
 *
 * The function evaluates the tensor components, at every orientation, from the
 * principal axis system (PAS) to the common frame of the isotopomer.
 */
static inline void __MRS_rotate_components_from_PAS_to_common_frame(
    isotopomer_ravel *ravel_isotopomer,  // isotopomer structure
    int site,                            // the site index
    double *transition,                  // the  transition
    bool allow_fourth_rank,  // if true, pre for 4th rank computation
    double *R0,              // the R0 components
    complex128 *R2,          // the R2 components
    complex128 *R4,          // the R4 components
    double *R0_temp,         // the temporary R0 components
    complex128 *R2_temp,     // the temporary R2 components
    complex128 *R4_temp,     // the temporary R3 components
    int remove_second_order_quad_isotropic  // if true, remove second order quad
                                            // isotropic shift
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
   *   Symmetry pathways in solid-state NMR. PNMRS 2011 59(2):12 1-96.
   *   https://doi.org/10.1016/j.pnmrs.2010.11.003                         */

  /* Initialize with zeroing all spatial components                        */
  __zero_components(R0, R2, R4);

  /* Nuclear shielding components ======================================== */
  /*      Upto the first order                                             */
  FCF_1st_order_nuclear_shielding_Hamiltonian(
      R0_temp, R2_temp, ravel_isotopomer->isotropic_chemical_shift_in_Hz[site],
      ravel_isotopomer->shielding_anisotropy_in_Hz[site],
      ravel_isotopomer->shielding_asymmetry[site],
      &ravel_isotopomer->shielding_orientation[3 * site], transition);

  // in-place update the R0 and R2 components.
  *R0 += *R0_temp;
  vm_double_add_inplace(10, (double *)R2_temp, (double *)R2);
  /* ===================================================================== */

  /* Weakly coupled direct-dipole components ============================= */
  /*      Upto the first order (to do..-> add orientation dependence)      */
  weakly_coupled_direct_dipole_frequencies_to_first_order(
      R0, R2_temp, ravel_isotopomer->dipolar_couplings[site], transition);

  // in-place update the R2 components.
  vm_double_add_inplace(10, (double *)R2_temp, (double *)R2);
  /* ===================================================================== */

  /* Electric quadrupolar components ===================================== */
  if (ravel_isotopomer->spin > 0.5) {
    /*   Upto the first order                                              */
    FCF_1st_order_electric_quadrupole_Hamiltonian(
        R2_temp, ravel_isotopomer->spin,
        ravel_isotopomer->quadrupole_coupling_constant_in_Hz[site],
        ravel_isotopomer->quadrupole_asymmetry[site],
        &ravel_isotopomer->quadrupole_orientation[3 * site], transition);

    // in-place update the R2 components.
    vm_double_add_inplace(10, (double *)R2_temp, (double *)R2);

    /*  Upto the second order                                              */
    if (allow_fourth_rank) {
      FCF_2nd_order_electric_quadrupole_Hamiltonian(
          R0_temp, R2_temp, R4_temp, ravel_isotopomer->spin,
          ravel_isotopomer->larmor_frequency,
          ravel_isotopomer->quadrupole_coupling_constant_in_Hz[site],
          ravel_isotopomer->quadrupole_asymmetry[site],
          &ravel_isotopomer->quadrupole_orientation[3 * site], transition);

      // in-place update the R0 component.
      if (remove_second_order_quad_isotropic == 0) *R0 += *R0_temp;

      // in-place update the R2 and R4 components.
      vm_double_add_inplace(10, (double *)R2_temp, (double *)R2);
      vm_double_add_inplace(18, (double *)R4_temp, (double *)R4);
    }
  }
}

void __mrsimulator_core(
    // spectrum information and related amplitude
    double *spec,  // amplitude vector representing the spectrum.

    isotopomer_ravel *ravel_isotopomer,  // isotopomer structure

    int remove_second_order_quad_isotropic,  // remove the isotropic
                                             // contribution from the second
                                             // order quad Hamiltonian.

    // Pointer to the transitions. transition[0] = mi and transition[1] = mf
    double *transition, MRS_plan *plan, MRS_dimension *dimension,
    bool interpolation  // if true, perform a 1D interpolation
) {
  /*
  The sideband computation is based on the method described by Eden and Levitt
  et. al. `Computation of Orientational Averages in Solid-State NMR by Gaussian
  Spherical Quadrature` JMR, 132, 1998. https://doi.org/10.1006/jmre.1998.1427
  */

  unsigned int site, i, j;
  double offset;                          // for 1D interpolation
  unsigned int step_vector = 0, address;  // for 1D interpolation

  double R0 = 0.0;
  complex128 *R2 = malloc_complex128(5);
  complex128 *R4 = malloc_complex128(9);

  double R0_temp = 0.0;
  complex128 *R2_temp = malloc_complex128(5);
  complex128 *R4_temp = malloc_complex128(9);

  int spec_site;
  double *spec_site_ptr;
  MRS_averaging_scheme *scheme = plan->averaging_scheme;

  // Per site base calculation
  for (site = 0; site < ravel_isotopomer->number_of_sites; site++) {
    // gettimeofday(&start_site_time, NULL);

    spec_site = site * dimension->count;
    spec_site_ptr = &spec[spec_site];

    /* Rotate all frequency components from PAS to a common frame */
    __MRS_rotate_components_from_PAS_to_common_frame(
        ravel_isotopomer,         // isotopomer structure
        site,                     // the site index
        transition,               // the  transition
        plan->allow_fourth_rank,  // if true, prepare for 4th rank computation
        &R0,                      // the R0 components
        R2,                       // the R2 components
        R4,                       // the R4 components
        &R0_temp,                 // the temporary R0 components
        R2_temp,                  // the temporary R2 components
        R4_temp,                  // the temporary R4 components
        remove_second_order_quad_isotropic  // if true, remove second order quad
                                            // isotropic shift
    );

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
    //         plan->integration_density,
    //         (double *)&plan->vector[i * plan->octant_orientations], 2,
    //         dimension->count);
    //   }
    // }
    if (interpolation) {
      for (i = 0; i < plan->number_of_sidebands; i++) {
        offset = plan->vr_freq[i] + plan->isotropic_offset;
        if ((int)offset >= 0 && (int)offset <= dimension->count) {
          step_vector = i * scheme->total_orientations;
          for (j = 0; j < plan->n_octants; j++) {
            address = j * scheme->octant_orientations;

            vm_double_ramp(scheme->octant_orientations,
                           &scheme->local_frequency[address], 1.0, offset,
                           scheme->freq_offset);
            octahedronInterpolation(
                spec_site_ptr, scheme->freq_offset, scheme->integration_density,
                (double *)&plan->vector[step_vector], 2, dimension->count);
            step_vector += scheme->octant_orientations;
          }
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

void mrsimulator_core(
    // spectrum information and related amplitude
    double *spec,               // The amplitude of the spectrum.
    double coordinates_offset,  // The start of the frequency spectrum.
    double increment,           // The increment of the frequency spectrum.
    int count,                  // Number of points on the frequency spectrum.
    isotopomer_ravel *ravel_isotopomer,      // Isotopomer structure
    int quad_second_order,                   // Quad theory for second order,
    int remove_second_order_quad_isotropic,  // remove the isotropic
                                             // contribution from the second
                                             // order quad interaction.

    // spin rate, spin angle and number spinning sidebands
    int number_of_sidebands,                 // The number of sidebands
    double sample_rotation_frequency_in_Hz,  // The rotor spin frequency
    double rotor_angle_in_rad,  // The rotor angle relative to lab-frame z-axis

    // Pointer to the transitions. transition[0] = mi and transition[1] = mf
    double *transition,

    // powder orientation average
    int integration_density,          // The number of triangle along the edge
                                      // of octahedron
    unsigned int integration_volume,  // 0-octant, 1-hemisphere, 2-sphere.
    bool interpolation) {
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
      integration_density, allow_fourth_rank, integration_volume);

  MRS_dimension *dimension =
      MRS_create_dimension(count, coordinates_offset, increment);

  // gettimeofday(&begin, NULL);
  MRS_plan *plan = MRS_create_plan(
      scheme, number_of_sidebands, sample_rotation_frequency_in_Hz,
      rotor_angle_in_rad, increment, allow_fourth_rank);

  // gettimeofday(&all_site_time, NULL);
  __mrsimulator_core(
      // spectrum information and related amplitude
      spec,  // The amplitude of the spectrum.

      ravel_isotopomer,  // isotopomer structure

      remove_second_order_quad_isotropic,  // remove the isotropic contribution
                                           // from the second order quad
                                           // Hamiltonian.

      // Pointer to the transitions. transition[0] = mi and transition[1] = mf
      transition,

      plan, dimension, interpolation);

  // gettimeofday(&end, NULL);
  // clock_time = (double)(end.tv_usec - begin.tv_usec) / 1000000. +
  //              (double)(end.tv_sec - begin.tv_sec);
  // printf("time %f s\n", clock_time);
  // cpu_time_[0] += clock_time;

  MRS_free_plan(plan); /* clean up */
}
