// -*- coding: utf-8 -*-
//
//  isotopomer_ravel.h
//
//  @copyright Deepansh J. Srivastava, 2019-2020.
//  Created by Deepansh J. Srivastava, Apr 11, 2019.
//  Contact email = deepansh2012@gmail.com
//

#ifndef isotopomer_ravel_h
#define isotopomer_ravel_h

// isotopomer like structure
struct __isotopomer_ravel {
  unsigned int number_of_sites; /**< Number of sites */
  /* Pointer to an array of the spin quantum numbers for each site. */
  float *spin;

  /* Pointer to an array of the gyromagnetic ratio in MHz/T for each site. */
  double *gyromagnetic_ratio;

  /* Pointer to an array of the Isotropic chemical shifts of each site given in
   * units of ppm. */
  double *isotropic_chemical_shift_in_ppm;

  /* Pointer to an array of symmetric Nuclear shielding tensor anisotropy (zeta)
   * in ppm for each site. */
  double *shielding_symmetric_zeta_in_ppm;

  /* Pointer to an array of symmetric Nuclear shielding tensor asymmetry
   * parameter (eta) for each site. */
  double *shielding_symmetric_eta;

  /* Pointer to an array of symmetric Nuclear shielding tensor euler angles
     given in radians and packed as alpha, beta, gamma. The array size is
     3*number_of_sites. The euler angle rotate the tensor from the principal
     axis system (PAS) to a common frame. */
  double *shielding_orientation;

  /* Pointer to an array of Quadrupole coupling constant for each site given in
   * units of Hz. If the site is not a quadrupole, a value of zero is applied
   * instead.*/
  double *quadrupolar_Cq_in_Hz;

  /* Pointer to an array of symmetric Electric field gradient tensor asymmetry
   * parameter (eta) for each site. */
  double *quadrupolar_eta;

  /* Pointer to an array of symmetric Electric field gradient tensor euler
     angles given in radians and packed as alpha, beta, gamma. The array size is
     3*number_of_sites. The euler angle rotate the tensor from the principal
     axis system (PAS) to a common frame. */
  double *quadrupolar_orientation;

  double *dipolar_couplings; /**< Pointer to an array of dipolar couplings */
};

typedef struct __isotopomer_ravel isotopomer_ravel;

#endif /* isotopomer_h */
