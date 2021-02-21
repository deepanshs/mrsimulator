// -*- coding: utf-8 -*-
//
//  object_struct.h
//
//  @copyright Deepansh J. Srivastava, 2019-2021.
//  Created by Deepansh J. Srivastava, Apr 11, 2019.
//  Contact email = srivastava.89@osu.edu
//

#ifndef object_struct_h
#define object_struct_h

/**
 * =====================================================================================
 *                                    Site structure
 * =====================================================================================
 * Site structure is a collection of site interaction parameters.
 **/
struct __site_struct {
  unsigned int number_of_sites; /**< Number of sites */

  /* Pointer to an array of spin quantum numbers for each site within a spin system. */
  float *spin;

  /* Pointer to an array of gyromagnetic ratio in MHz/T for each site within a spin
   * system. */
  double *gyromagnetic_ratio;

  /* Pointer to an array of isotropic chemical shifts in ppm for each site within a spin
   * system. */
  double *isotropic_chemical_shift_in_ppm;

  /* Pointer to an array of symmetric nuclear shielding tensor anisotropies (zeta) in
  ppm of each site within a spin system. */
  double *shielding_symmetric_zeta_in_ppm;

  /* Pointer to an array of symmetric nuclear shielding tensor asymmetry parameters
   * (eta) of each site within a spin system. */
  double *shielding_symmetric_eta;

  /* Pointer to an array of symmetric nuclear shielding tensor Euler angles in radians,
     packed as alpha, beta, and gamma. The array size is 3*number_of_sites. The Euler
     angle rotate the tensor from the principal axis system (PAS) to a common frame. */
  double *shielding_orientation;

  /* Pointer to an array of quadrupole coupling constants in Hz of each site within a
   * spin system. If the site is not a quadrupole, a value of zero is set instead. */
  double *quadrupolar_Cq_in_Hz;

  /* Pointer to an array of symmetric electric field gradient tensor asymmetry
   * parameters (eta) of each site within a spin system. If the site is not a
   * quadrupole, a value of zero is set instead. */
  double *quadrupolar_eta;

  /* Pointer to an array of symmetric electric field gradient tensor Euler angles in
     radians, packed as alpha, beta, and gamma. The array size is 3*number_of_sites. The
     Euler angle rotate the tensor from the principal axis system (PAS) to a common
     frame. If the site is not a quadrupole, a value of zero is set instead. */
  double *quadrupolar_orientation;
};
typedef struct __site_struct site_struct;

/**
 * =====================================================================================
 *                                   Coupling structure
 * =====================================================================================
 * Coupling structure is a collection of coupled site interaction parameters.
 **/
struct __coupling_struct {
  unsigned int number_of_couplings; /**< Number of couplings */

  /* Pointer to an array of the site indexes for each coupled pair within a spin system,
   * incremented in stride of 2/pair. The array size is 2*number_of_couplings. */
  int *site_index;

  /* Pointer to an array of isotropic J-couplings in Hz for each pair within the spin
   * system. */
  double *isotropic_j_in_Hz;

  /* Pointer to an array of J anisotropies (zeta) in Hz for each pair within the spin
   * system. */
  double *j_symmetric_zeta_in_Hz;

  /* Pointer to an array of J asymmetry parameters (eta) for each pair within the spin
   * system. */
  double *j_symmetric_eta;

  /* Pointer to an array of symmetric J tensor Euler angles in radians, packed as alpha,
     beta, and gamma. The array size is 3*number_of_couplings. The Euler angle rotate
     the J tensor from the principal axis system (PAS) to a common frame. */
  double *j_orientation;

  /* Pointer to an array of direct-dipolar coupling constants in Hz for each pair within
   * the spin system. */
  double *dipolar_coupling_in_Hz;

  /* Pointer to an array of direct-dipolar asymmetry parameters (eta) for each pair
   * within the spin system. */
  double *dipolar_eta;

  /* Pointer to an array of direct-dipolar tensor Euler angles in radians, packed as
     alpha, beta, and gamma. The array size is 3*number_of_couplings. The Euler angle
     rotate the tensor from the principal axis system (PAS) to a common frame. */
  double *dipolar_orientation;
};
typedef struct __coupling_struct coupling_struct;

#endif /* object_struct_h */
