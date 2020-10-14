// -*- coding: utf-8 -*-
//
//  powder_setup.h
//
//  @copyright Deepansh J. Srivastava, 2019-2020.
//  Created by Deepansh J. Srivastava, Apr 11, 2019.
//  Contact email = srivastava.89@osu.edu
//

#include "config.h"
#include "interpolation.h"

extern void octahedronGetDirectionCosineSquareOverOctantAndWeights(
    int nt, double *xr, double *yr, double *zr, double *amp);

extern void octahedronGetPolarAngleTrigOverAnOctant(int nt, double *cos_alpha,
                                                    double *cos_beta,
                                                    double *amp);

extern void octahedronGetPolarAngleCosineAzimuthalAnglePhaseOverOctant(
    int nt, void *exp_I_alpha, void *exp_I_beta, double *amp);

extern void octahedronInterpolation(double *spec, double *freq, int nt,
                                    double *amp, int stride, int m);

extern void octahedronInterpolation2D(double *spec, double *freq1,
                                      double *freq2, int nt, double *amp,
                                      int stride, int m0, int m1);
