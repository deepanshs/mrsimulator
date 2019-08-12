//
//  isotopomer_ravel.h
//
//  Created by Deepansh J. Srivastava, Apr 11, 2019
//  Copyright © 2019 Deepansh J. Srivastava. All rights reserved.
//  Contact email = srivastava.89@osu.edu, deepansh2012@gmail.com
//

#ifndef isotopomer_ravel_h
#define isotopomer_ravel_h

// isotopomer like structure
struct __isotopomer_ravel {
  int number_of_sites;                    /* Number of sites */
  float spin;                             /* The spin quantum number */
  double larmor_frequency;                /* Larmor frequency (MHz) */
  double *isotropic_chemical_shift_in_Hz; /* Isotropic chemical shift (Hz) */
  double *shielding_anisotropy_in_Hz; /* Nuclear shielding anisotropy (Hz) */
  double *shielding_asymmetry;   /* Nuclear shielding asymmetry parameter */
  double *shielding_orientation; /* Nuclear shielding PAS to CRS euler angles
                                    (rad.) */
  double *quadrupolar_constant_in_Hz; /* Quadrupolar coupling constant (Hz) */
  double *quadrupolar_asymmetry;      /* Quadrupolar asymmetry parameter */
  double
      *quadrupolar_orientation; /* Quadrupolar PAS to CRS euler angles (rad.) */
  double *dipolar_couplings;    /* dipolar coupling stored as list of lists */
};

typedef struct __isotopomer_ravel isotopomer_ravel;

struct __isotopomers_list {
  isotopomer_ravel *isotopomers;
};

typedef struct __isotopomers_list isotopomers_list;

#endif /* isotopomer_h */
