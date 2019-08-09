
//
//  powder_setup.h
//
//  Created by Deepansh J. Srivastava, Apr 11, 2019
//  Copyright © 2019 Deepansh J. Srivastava. All rights reserved.
//  Contact email = srivastava.89@osu.edu, deepansh2012@gmail.com
//

#include "mrsimulator.h"

extern void octahedronGetDirectionCosineSquareOverOctantAndWeights(
    int nt, double *xr, double *yr, double *zr, double *amp);

extern void octahedronGetPolarAngleTrigOverAnOctant(int nt, double *cos_alpha,
                                                    double *cos_beta,
                                                    double *amp);

extern void octahedronInterpolation(double *spec, double *freq, int nt,
                                    double *amp, int stride, int m);
