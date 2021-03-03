// -*- coding: utf-8 -*-
//
//  octahedron.h
//
//  @copyright Deepansh J. Srivastava, 2019-2021.
//  Created by Deepansh J. Srivastava, Apr 11, 2019.
//  Contact email = srivastava.89@osu.edu
//

#include "config.h"
#include "interpolation.h"

extern void octahedronGetDirectionCosineSquareAndWeightsOverOctant(
    const unsigned int nt, double *restrict xr, double *restrict yr,
    double *restrict zr, double *restrict amp);

extern void octahedronGetPolarAngleTrigOverOctant(const unsigned int nt,
                                                  double *restrict cos_alpha,
                                                  double *restrict cos_beta,
                                                  double *restrict amp);

extern void octahedronGetComplexExpOfPolarAngleOverOctant(const unsigned int nt,
                                                          void *restrict exp_I_alpha,
                                                          void *restrict exp_I_beta,
                                                          double *restrict amp);

void get_total_amplitude(const unsigned int nt, double *amp, double *amp_sum);

extern void averaging_setup(unsigned int nt, void *exp_I_alpha, void *exp_I_beta,
                            double *amp);
