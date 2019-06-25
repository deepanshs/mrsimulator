
//
//  c_array.h
//
//  Created by Deepansh J. Srivastava, Apr 11, 2019
//  Copyright © 2019 Deepansh J. Srivastava. All rights reserved.
//  Contact email = srivastava.89@osu.edu, deepansh2012@gmail.com
//

#include "complex.h"
// to use calloc, malloc, and free methods
#include <stdlib.h>

extern float *createFloat1DArray(int m);
extern void destroyFloat1DArray(float *arr);

extern float complex *createFloatComplex1DArray(int m);
extern void destroyFloatComplex1DArray(float complex *arr);

extern double *createDouble1DArray(int m);
extern void destroyDouble1DArray(double *arr);

extern double complex *createDoubleComplex1DArray(int m);
extern void destroyDoubleComplex1DArray(double complex *arr);

extern float **createFloat2DMatrix(int m, int n);
extern void destroyFloat2DMatrix(float **arr);

extern double **createDouble2DMatrix(int n, int m);
extern void destroyDouble2DMatrix(double **arr);

extern double ***createDouble3DArray(int n, int m, int o);
extern void destroyDouble3DArray(double ***arr);

// extern OCPolarAngleTrig** create2DOCPolarAngleTrigArray(int m, int n);
// extern void destroy2DOCPolarAngleTrigArray(OCPolarAngleTrig** arr);

// extern OCDirectionCosineSquare** create2DOCDirectionCosineSquareArray(int m,
// int n);
// extern void destroy2DOCDirectionCosineSquareArray(OCDirectionCosineSquare**
// arr);
