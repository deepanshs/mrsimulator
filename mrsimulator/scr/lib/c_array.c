//
//  c_array.c
//
//  Created by Deepansh J. Srivastava, Apr 11, 2019
//  Copyright © 2019 Deepansh J. Srivastava. All rights reserved.
//  Contact email = srivastava.89@osu.edu, deepansh2012@gmail.com
//

#include "c_array.h"
#include <complex.h>

// float 1d array
float *createFloat1DArray(int m)
{
    float *values = calloc(m, sizeof(float));
    return values;
}

void destroyFloat1DArray(float *arr)
{
    free(arr);
}

// complex float 1d array
float complex *createFloatComplex1DArray(int m)
{
    float complex *values = calloc(m, sizeof(float complex));
    return values;
}

void destroyFloatComplex1DArray(float complex *arr)
{
    free(arr);
}

// double 1d array
double *createDouble1DArray(int m)
{
    double *values = calloc(m, sizeof(double));
    return values;
}

void destroyDouble1DArray(double *arr)
{
    free(arr);
}

// complex double 1d array
double complex *createDoubleComplex1DArray(int m)
{
    double complex *values = calloc(m, sizeof(double complex));
    return values;
}

void destroyDoubleComplex1DArray(double complex *arr)
{
    free(arr);
}

// float 2d matrix
float **createFloat2DMatrix(int n, int m)
{
    float *values = calloc(m * n, sizeof(float));
    float **rows = malloc(n * sizeof(float *));
    for (int i = 0; i < n; ++i)
    {
        rows[i] = values + i * m;
    }
    return rows;
}

void destroyFloat2DMatrix(float **arr)
{
    free(*arr);
    free(arr);
}

// double array
double **createDouble2DMatrix(int n, int m)
{
    double *values = calloc(m * n, sizeof(double));
    double **rows = malloc(n * sizeof(double *));
    for (int i = 0; i < n; ++i)
    {
        rows[i] = values + i * m;
    }
    return rows;
}

void destroyDouble2DMatrix(double **arr)
{
    free(*arr);
    free(arr);
}

// double 3d array
double ***createDouble3DArray(int n, int m, int o)
{
    double ***p = malloc(n * sizeof(double ***));
    for (int i = 0; i < n; i++)
    {
        p[i] = (double **)malloc(m * sizeof(double *));

        for (int j = 0; j < m; j++)
        {
            p[i][j] = (double *)malloc(o * sizeof(double));
        }
    }
    return p;
}

void destroyDouble3DArray(double ***arr)
{
    free(**arr);
    free(*arr);
    free(arr);
}
