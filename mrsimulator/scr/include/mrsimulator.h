//
//  mrsimulator.h
//
//  Created by Deepansh J. Srivastava, Jun 30, 2019
//  Copyright © 2019 Deepansh J. Srivastava. All rights reserved.
//  Contact email = srivastava.89@osu.edu, deepansh2012@gmail.com
//

#ifndef mrsimulator_h
#define mrsimulator_h

// library definition
#include <complex.h>
#define MKL_Complex16 double complex

#include "mkl.h"
#include "mkl_dfti.h"
// #include "omp.h"
#include <math.h>
// to use calloc, malloc, and free methods
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>

// user definition
#define PI2 6.2831853072
#define PI2I PI2 *I

#include "Hamiltonian.h"
#include "angular_momentum.h"
#include "array.h"
#include "fftw/fftw3.h"
#include "interpolation.h"
#include "isotopomer_ravel.h"
#include "octahedron.h"
#include "powder_setup.h"

#endif /* mrsimulator_h */
