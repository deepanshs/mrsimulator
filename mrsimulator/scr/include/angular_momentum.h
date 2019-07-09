//
//  angular_momentum.c
//
//  Created by Philip Grandinetti on 4/12/17.
//  Copyright © 2017 Philip Grandinetti. All rights reserved.
//  Contribution: Deepansh J. Srivatava. contact: srivastava.89@osu.edu
//      Contact email = srivastava.89@osu.edu, deepansh2012@gmail.com

#include "mrsimulator.h"

extern void full_DLM(double complex *wigner, int l, double *omega);

extern double wigner_d(int l, int m1, int m2, double beta);

// extern double complex DLM(int l, int  m1, int m2, OCEulerAngle omega);

void full_DLM_trig(double complex *wigner, int l, double cosAlpha,
                   double sinAlpha, double cosBeta, double sinBeta);

void get_even_DLM_4_from_2(double complex *wigner, double cosBeta);

double wigner_d_trig(int l, int m1, int m2, double cx, double sx);

/*!
@function __wigner_d_matrix
@abstract Evaluates nx(2l+1)x(2l+1) wigner-d matrices of rank `l` for `n`
angles given in radians. Here `angle` is a 1D array of size `n`.
@param l: The rank of the matrix of type int.
@param n: The number of angles of type int.
@param angle: Pointer to a 1D-array of angles `beta` expressed in radian of
              type double.
@param wigner: Pointer to the wigner-d matrix of type double.
 */
extern void __wigner_d_matrix(int l, int n, double *angle, double *wigner);

/*!
@function __wigner_d_matrix_cosine
@abstract Evaluates nx(2l+1)x(2l+1) wigner-d matrices of rank `l` for `n`
angles given as cosine of angles. Here `cos_angle` is a 1D array of size `n`.
@param l: The rank of the matrix of type int.
@param n: The number of cosine angles of type int.
@param cos_angle: Pointer to a 1D-array of angles `beta` in radian of type
                 double.
@param wigner: Pointer to the wigner-d matrix of type double.
 */
extern void __wigner_d_matrix_cosine(int l, int n, double *cos_angle,
                                     double *wigner);

/*!
@function __wigner_rotation
@abstract Evaluates the wigner rotation of a vector of length `l` over
`n` wigner matrices. The result is a stack of output vector of size `nxl`
evaluated at all `n` wigner matrices.
@param l: The rank `l` of the wigner-d matrix of type int.
@param n: The number of cosine alpha angles and wigner-lj matrices of type int.
@param cos_alpha: Pointer to the 1d array of cosine alpha of type double.
@param wigner: Pointer to the nx(2l+1)x(2l+1) wigner-d matrices of type double.
               The wigner matrices are stacked in a row major order.
@param R_in: Pointer to a 1D-array of initial vector of type double. The length
             of this vector is `l`.
@param R_out: Pointer to a 1D-array of final vectors of type double. The length
              of this vector is `n*l`. The final vectors are stacked in a row
              major order similar to the ordering of the wigner-matrices.
 */
extern void __wigner_rotation(int l, int n, double *wigner, double *cos_alpha,
                              double complex *R_in, double complex *R_out);

extern void __wigner_rotation_2(int l, int n, double *wigner,
                                double complex *exp_Im_alpha,
                                double complex *R_in, double complex *R_out);

extern void __wigner_dm0_vector(int l, double beta, double *R_out);
