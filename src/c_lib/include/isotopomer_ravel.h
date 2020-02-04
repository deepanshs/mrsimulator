// -*- coding: utf-8 -*-
//
//  isotopomer_ravel.h
//
//  @copyright Deepansh J. Srivastava, 2019-2020.
//  Created by Deepansh J. Srivastava, Apr 11, 2019
//  Contact email = deepansh2012@gmail.com
//

#ifndef isotopomer_ravel_h
#define isotopomer_ravel_h

// isotopomer like structure
struct __isotopomer_ravel {
  unsigned int number_of_sites;           /**< Number of sites */
  float spin;                             /**< The spin quantum number */
  double larmor_frequency;                /**< Larmor frequency (MHz) */
  double *isotropic_chemical_shift_in_Hz; /**< Isotropic chemical shift (Hz) */
  double *shielding_anisotropy_in_Hz; /**< Nuclear shielding anisotropy (Hz) */
  double *shielding_asymmetry;   /**< Nuclear shielding asymmetry parameter */
  double *shielding_orientation; /**< Nuclear shielding PAS to CRS euler angles
                                    (rad.) */

  double *quadrupole_coupling_constant_in_Hz; /**< Quadrupole coupling constant
                                                 (Hz) */
  double *quadrupole_asymmetry; /**< Quadrupole asymmetry parameter */
  double
      *quadrupole_orientation; /**< Quadrupole PAS to CRS euler angles (rad.) */
  double *dipolar_couplings;   /**< dipolar coupling stored as list of lists */
};

typedef struct __isotopomer_ravel isotopomer_ravel;

#endif /* isotopomer_h */
