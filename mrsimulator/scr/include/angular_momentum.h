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

/**
 * @brief Evaluates `n` wigner-d matrices of rank `l` at every given angles
 *          @f$\beta@f$ (in radians).
 *
 * Each wigner-d matrix is `(2l+1) x (2l+1)` in dimension, where @p l is the
 * rank. When evaluating @p n wigner-d matrices, the matrices are stored such
 * that the wigner-d matrix corresponding to the angle `alpha[i]` starts at the
 * index `i*(2*l+1)*(2*l+1)`.
 *
 * @param l The rank of the wigner-d matrix.
 * @param n The number of wigner-d matrix to evaluate.
 * @param angle A pointer to a 1D-array of angles @f$\beta@f$ of length @p n,
 *          given in radians.
 * @param wigner A pointer to the start of wigner-d matrices.
 */
extern void AMT_wigner_d_matrix(int l, int n, double *angle, double *wigner);

/**
 * @brief Evaluates @p n wigner-d matrices of rank @p l at every given
 *          **cosine** of angles, @f$\beta@f$.
 *
 * Each wigner-d matrix is `(2l+1) x (2l+1)` in dimension, where @p l is the
 * rank. When evaluating @p n wigner-d matrices, the matrices are stored such
 * that the wigner-d matrix corresponding to the angle `alpha[i]` starts at the
 * index `i*(2*l+1)*(2*l+1)`.
 *
 * @param l The rank of the wigner-d matrix.
 * @param n The number of wigner-d matrix to evaluate.
 * @param cos_angle A pointer to a 1D-array of cosine of angles @f$\beta@f$ of
 *          length @p n.
 * @param wigner A pointer to the start of wigner-d matrices.
 */
extern void __wigner_d_matrix_cosine(int l, int n, double *cos_angle,
                                     double *wigner);

/**
 * @brief Rotate @p R_in of length @p l using wigner-d matrices of rank @p l.
 * wigner matrices. The result is a stack of output vector of size `nxl`
evaluated
 * at all `n` wigner matrices.
 * @param l The rank of the wigner-d matrices.
 * @param n The number of cosine alpha angles and wigner-lj matrices.
 * @param cos_alpha A pointer to the 1d array of @f$\cos\alphas@f$.
 * @param wigner A pointer to nx(2l+1)x(2l+1) wigner-d matrices of rank @p l.
 *          The wigner matrices are stacked in a row major order.
 * @param R_in A pointer to a 1D-array of initial vector of length `2l+1`.
 * @param R_out A pointer to a 1D-array of final vectors after rotation. The
 *          length of this vector is `n*(2*l+1)`, where the vector of index
 *          `i*(2*l+1)` is rotated with the wigner-lj matrix at index
 *          `i*(2*l+1)*(2*l+1)`.
 */
extern void __wigner_rotation(int l, int n, double *wigner, double *cos_alpha,
                              double complex *R_in, double complex *R_out);

extern void __wigner_rotation_2(int l, int n, double *wigner,
                                double complex *exp_Im_alpha,
                                double complex *R_in, double complex *R_out);

extern void __wigner_dm0_vector(int l, double beta, double *R_out);
