
//
//  array.h
//
//  Created by Deepansh J. Srivastava, Apr 11, 2019
//  Copyright © 2019 Deepansh J. Srivastava. All rights reserved.
//  Contact email = srivastava.89@osu.edu, deepansh2012@gmail.com
//

#include "mrsimulator.h"

// allocate memory for a float array of size m .
extern float *malloc_float(int m);
// free allocated memory of a float array of size m.
extern void free_float(float *arr);

// allocate memory for a float complex array of size m.
// extern float complex *malloc_float_complex(int m);
// // free allocated memory of a float complex array of size m.
// extern void free_float_complex(float complex *arr);

// allocate memory for a double array of size m.
extern double *malloc_double(int m);
// free allocated memory of a double array of size m.
extern void free_double(double *arr);

// allocate memory for a complex128 array of size m.
extern complex128 *malloc_double_complex(int m);
// free allocated memory of a complex128 array of size m.
extern void free_double_complex(complex128 *arr);
