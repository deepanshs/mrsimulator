// -*- coding: utf-8 -*-
//
//  frequency_averaging.h
//
//  @copyright Deepansh J. Srivastava, 2019-2021.
//  Created by Deepansh J. Srivastava, Mar 2, 2021.
//  Contact email = srivastava.89@osu.edu
//

#include "simulation.h"

void one_dimensional_averaging(MRS_dimension *dimensions, MRS_averaging_scheme *scheme,
                               double *spec, unsigned int iso_intrp,
                               complex128 *exp_I_phase);

void two_dimensional_averaging(MRS_dimension *dimensions, MRS_averaging_scheme *scheme,
                               double *spec, double *affine_matrix,
                               unsigned int iso_intrp, complex128 *exp_I_phase);
