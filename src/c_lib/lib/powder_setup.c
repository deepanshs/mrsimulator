// -*- coding: utf-8 -*-
//
//  powder_setup.c
//
//  @copyright Deepansh J. Srivastava, 2019-2020.
//  Created by Deepansh J. Srivastava, Apr 11, 2019
//  Contact email = srivastava.89@osu.edu
//

#include "powder_setup.h"

#include "octahedron.h"

void octahedron_averaging_setup(int nt, void *exp_I_alpha, void *exp_I_beta,
                                double *amp) {
  // octahedronGetPolarAngleTrigOverAnOctant(nt, cos_alpha, cos_beta, amp);
  octahedronGetPolarAngleCosineAzimuthalAnglePhaseOverOctant(nt, exp_I_alpha,
                                                             exp_I_beta, amp);
}
