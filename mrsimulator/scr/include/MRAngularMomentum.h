
#include <complex.h>
#define MKL_Complex16 double complex

#include "mkl.h"
#include <math.h>
#include <stdio.h>

extern void full_DLM(double complex *wigner, int l, double *omega);

extern double wigner_d(int l, int m1, int m2, double beta);

// extern double complex DLM(int l, int  m1, int m2, OCEulerAngle omega);

void full_DLM_trig(double complex *wigner, int l, double cosAlpha,
                   double sinAlpha, double cosBeta, double sinBeta);

void get_even_DLM_4_from_2(double complex *wigner, double cosBeta);

double wigner_d_trig(int l, int m1, int m2, double cx, double sx);

// extern double wigner4(double beta, int m1, int m2);

extern void wigner_d_matrix(double *wigner, int l, double *value, int trig);
