//
//  angular_momentum.c
//
//  Created by Philip Grandinetti on 4/12/17.
//  Copyright © 2017 Philip Grandinetti. All rights reserved.
//  Contribution: Deepansh J. Srivatava. contact: srivastava.89@osu.edu
//      Contact email = srivastava.89@osu.edu, deepansh2012@gmail.com

#include "mrsimulator.h"

extern void full_DLM(complex128 *wigner, int l, double *omega);

extern double wigner_d(int l, int m1, int m2, double beta);

// extern complex128 DLM(int l, int  m1, int m2, OCEulerAngle omega);

void full_DLM_trig(complex128 *wigner, int l, double cosAlpha, double sinAlpha,
                   double cosBeta, double sinBeta);

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
extern void __wigner_d_matrix_cosine(const int l, const int n,
                                     const double *cos_angle, double *wigner);

/**
 * @brief Rotate @p R_in of length @p l using wigner-d matrices of rank @p l.
 * wigner matrices. The result is a stack of output vector of size `nxl`
 * evaluated at all `n` wigner matrices.
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
extern void __wigner_rotation(const int l, const int n, const double *wigner,
                              const double *cos_alpha, const complex128 *R_in,
                              complex128 *R_out);

extern void __wigner_rotation_2(const int l, const int n, const double *wigner,
                                const complex128 *exp_Im_alpha,
                                const complex128 *R_in, complex128 *R_out);

extern void __wigner_dm0_vector(const int l, const double beta, double *R_out);

/**
 * ✅
 * @brief Performs a rank l wigner rotation of the coefficients from the l rank
 * spherical tensors.
 *
 * @param l The rank of the tensor.
 * @param euler_angles A pointer to the array of three euler angles.
 * @param R_in A pointer to the array of coefficients from the l rank tensors of
 *          length 2xl+1 before rotation.
 * @param R_out A pointer to the array of coefficients from the l rank tensors
 *          of length 2xl+1 after rotation.
 */
extern void single_wigner_rotation(const int l, const double *euler_angles,
                                   const complex128 *R_in, complex128 *R_out);

extern void wigner_d_matrix(const int l, const int n, const double *angle,
                            double *wigner);

extern void __batch_wigner_rotation(
    const unsigned int octant_orientations, const unsigned int n_octants,
    const double *wigner_2j_matrices, const complex128 *R2,
    const double *wigner_4j_matrices, const complex128 *R4,
    complex128 *exp_Im_alpha, complex128 *w2, complex128 *w4);

extern void get_exp_Im_alpha(const unsigned int n, const double *cos_alpha,
                             const bool allow_fourth_rank,
                             complex128 *exp_Im_alpha);
