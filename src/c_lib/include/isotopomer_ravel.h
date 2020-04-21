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
  unsigned int number_of_sites; /**< Number of sites */
  float *spin; /**< Pointer to an array of spin quantum numbers. */
  double *gyromagnetic_ratio; /**< Pointer to an array of gyromagnetic ratio
                                 (MHz/T). */
  double *isotropic_chemical_shift_in_Hz; /**< Pointer to an array of Isotropic
                                             chemical shift (Hz). */
  double *shielding_anisotropy_in_Hz;     /**< Pointer to an array of Nuclear
                                             shielding anisotropy (Hz). */
  double *shielding_asymmetry;   /**< Pointer to an array of Nuclear shielding
                                    asymmetry parameter. */
  double *shielding_orientation; /**< Pointer to an array of Nuclear shielding
                                    PAS to CRS euler angles (rad.). The array
                                    size is 3*number_of_sites */

  double *
      quadrupole_coupling_constant_in_Hz; /**< Pointer to an array of Quadrupole
                                             coupling constant (Hz). */
  double *quadrupole_asymmetry; /**< Pointer to an array of Quadrupole asymmetry
                                   parameter. */
  double *quadrupole_orientation; /**< Pointer to an array of Quadrupole PAS to
                                     CRS euler angles (rad.).  The array size is
                                     3*number_of_sites. */
  double *dipolar_couplings; /**< Pointer to an array of dipolar couplings */
};

typedef struct __isotopomer_ravel isotopomer_ravel;

#endif /* isotopomer_h */
