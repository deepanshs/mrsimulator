// -*- coding: utf-8 -*-
//
//  powder_setup.c
//
//  @copyright Deepansh J. Srivastava, 2019-2021.
//  Created by Deepansh J. Srivastava, Apr 11, 2019
//  Contact email = srivastava.89@osu.edu
//

#include "powder_setup.h"

#include "octahedron.h"

void averaging_setup(unsigned int nt, void *exp_I_alpha, void *exp_I_beta,
                     double *amp) {
  // octahedronGetPolarAngleTrigOverOctant(nt, cos_alpha, cos_beta, amp);
  octahedronGetComplexExpOfPolarAngleOverOctant(nt, exp_I_alpha, exp_I_beta, amp);
}
