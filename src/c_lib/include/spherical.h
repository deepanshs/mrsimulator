// -*- coding: utf-8 -*-
//
//  spherical.h
//
//  @copyright Deepansh J. Srivastava, 2019-2021.
//  Created by Deepansh J. Srivastava, Feb 12, 2021.
//  Contact email = srivastava.89@osu.edu
//

#include "config.h"
#include "interpolation.h"

extern void sphericalGetCoordinates(int nt, double *restrict x,
                                    double *restrict y, double *restrict z);

extern void sphericalGetComplexExpOfPolarAngleOverOctant(int nt,
                                                         void *exp_I_alpha,
                                                         void *exp_I_beta,
                                                         double *amp);
