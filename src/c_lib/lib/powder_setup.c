
//
//  powder_setup.c
//
//  Created by Deepansh J. Srivastava, Apr 11, 2019
//  Copyright © 2019 Deepansh J. Srivastava. All rights reserved.
//  Contact email = srivastava.89@osu.edu, deepansh2012@gmail.com
//

#include "powder_setup.h"
#include "octahedron.h"

void __powder_averaging_setup(int nt, double *cos_alpha, double *cos_beta,
                              double *amp) {
  octahedronGetPolarAngleTrigOverAnOctant(nt, cos_alpha, cos_beta, amp);
}
