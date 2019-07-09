//
//  array.c
//
//  Created by Deepansh J. Srivastava, Apr 11, 2019
//  Copyright © 2019 Deepansh J. Srivastava. All rights reserved.
//  Contact email = srivastava.89@osu.edu, deepansh2012@gmail.com
//

#include "array.h"
#include <complex.h>

/* float array */
float *malloc_float(int m) {
  float *values = malloc(m * sizeof(float));
  return values;
}

/* free float array */
void free_float(float *arr) { free(arr); }

/* float complex array */
float complex *malloc_float_complex(int m) {
  float complex *values = malloc(m * sizeof(float complex));
  return values;
}

/* free float complex array */
void free_float_complex(float complex *arr) { free(arr); }

/* double array */
double *malloc_double(int m) {
  double *values = malloc(m * sizeof(double));
  return values;
}

/* free double array */
void free_double(double *arr) { free(arr); }

/* complex double array */
double complex *malloc_double_complex(int m) {
  double complex *values = malloc(m * sizeof(double complex));
  return values;
}

/* free complex double array */
void free_double_complex(double complex *arr) { free(arr); }
