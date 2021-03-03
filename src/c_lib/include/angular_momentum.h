// -*- coding: utf-8 -*-
//
//  angular_momentum.c
//
//  @copyright Deepansh J. Srivastava, 2019-2021.
//  Created by Deepansh J. Srivastava.
//  Contact email = srivastava.89@osu.edu

#include "config.h"

// Generic function ................................................................. //

/**
 * @brief Evaluates @f$d^{l}_{m_1, m_2}(\beta)@f$ wigner-d element of the given angle
 * @f$\beta@f$.
 *
 * @param l The rank of the wigner-d matrix element.
 * @param m1 The quantum number @f$m_1@f$.
 * @param m2 The quantum number @f$m_2@f$.
 * @param beta The angle @f$\beta@f$ given in radian.
 * @return The wigner-d element, @f$d^{l}_{m_1, m_2}(\beta)@f$.
 */
extern double wigner_d_element(const int l, const int m1, const int m2,
                               const double beta);

/**
 * @brief Evaluates @f$n@f$ wigner-d matrices of rank @f$l@f$ at @f$n@f$ angles given in
 * radians.
 *
 * At a given angle, @f$\beta@f$, the wigner-d matrix, @f$d^{l}\left(m_1, m_2 |
 * \beta\right)@f$, is a @f$(2l+1) \times (2l+1)@f$ square matrix where @f$m_1@f$ and
 * @f$m_2@f$ range from @f$-l@f$ to @f$l@f$. The wigner-d elements, @f$d^{l}_{m_1,
 * m_2}(\beta)@f$, are ordered with @f$m_1@f$ as the leading dimension. For example,
 * when @f$l=2@f$, the wigner-d elements are ordered according to
 *
 * @f[
 *    \left[d^{2}_{-2, -2}(\beta),~d^{2}_{-1, -2}(\beta),~d^{2}_{0, -2}(\beta),
 *    ~~...~~d^{2}_{1, 2}(\beta),~d^{2}_{2, 2}(\beta) \right].
 * @f]
 *
 * The @f$n@f$ matrices are stored such that the wigner-d matrix from angle at index
 * @f$i@f$ is stacked after the wigner-d matrix from angle at index @f$i-1@f$, that is,
 * the wigner-d matrix corresponding to `angle[i]` starts at the index
 * `i*(2*l+1)*(2*l+1)`.
 *
 * @param l The rank of the wigner-d matrix.
 * @param n The number of wigner-d matrix and the size of `beta` array.
 * @param beta Pointer to an array of length @f$n@f$ of @f$\beta@f$ values (in rads).
 * @param wigner A pointer to an array of length `n*(2*l+1)*(2*l+1)`.
 */
extern void wigner_d_matrices(const int l, const int n, const double *beta,
                              double *wigner);

// Specialized methods .............................................................. //

/**
 * @brief Evaluates @f$d^{l}_{m_1, m_2}(\beta)@f$ wigner-d element from
 * @f$\exp(i\beta)@f$ of the angle @f$\beta@f$.
 *
 * @sa wigner_d_element(const int, const int, const int, const double)
 *
 * @param l The rank of the wigner-d matrix element.
 * @param m1 The quantum number @f$m_1@f$.
 * @param m2 The quantum number @f$m_2@f$.
 * @param exp_I_beta A pointer to an array of size two, where cosine and sine of the
 *      angle @f$\beta@f$ is stored.
 * @return The wigner-d element, @f$d^{l}_{m_1, m_2}(\beta)@f$.
 */
double wigner_d_element_from_exp_I_beta(const int l, const int m1, const int m2,
                                        const void *exp_I_beta);

/**
 * @brief Evaluates @f$n@f$ wigner-d matrices of rank @f$l@f$ from an array of angles
 * expressed as @f$\exp(i\beta)@f$.
 *
 * @see `wigner_d_matrices()` for details.
 *
 * @param l The rank of the wigner-d matrix.
 * @param n The number of wigner-d matrix to evaluate.
 * @param half Boolean. If true, calculate half of the wigner matrix.
 * @param exp_I_beta A pointer to an array of length `2*n`, where the cosine and sine of
 *      the angles are stored.
 * @param wigner A pointer to an array of length `n*(2*l+1)*(2*l+1)`.
 */
extern void wigner_d_matrices_from_exp_I_beta(const int l, const int n, const bool half,
                                              const void *restrict exp_I_beta,
                                              double *wigner);

/**
 * @brief Rotate @p R_in of length @p l using wigner-d matrices of rank @p l.
 *
 * The result is a stack of output vector of size `nxl` evaluated at all `n` wigner
 * matrices.
 *
 * @param l The rank of the wigner-d matrices.
 * @param n The number of cosine alpha angles and wigner-lj matrices.
 * @param wigner A pointer to nx(2l+1)x(2l+1) wigner-d matrices of rank @p l. The wigner
 *      matrices are stacked in a row major order.
 * @param exp_Im_alpha A pointer to the 1d array of @f$\exp(-i\alpha)@f$.
 * @param R_in A pointer to a 1D-array of initial vector of length `2l+1`.
 * @param R_out A pointer to a 1D-array of final vectors after rotation. The length of
 *      this vector is `n*(2*l+1)`, where the vector of index `i*(2*l+1)` is rotated
 *      with the wigner-lj matrix at index `i*(2*l+1)*(2*l+1)`.
 */
extern void __wigner_rotation_2(const int l, const int n, const double *wigner,
                                const void *exp_Im_alpha, const void *R_in,
                                void *R_out);

extern void wigner_dm0_vector(const int l, const double beta, double *R_out);

/**
 * ✅
 * @brief Performs a rank l wigner rotation of the coefficients from the l rank
 * spherical tensors.
 *
 * @param l The rank of the tensor.
 * @param euler_angles A pointer to the array of three euler angles.
 * @param R_in A pointer to the array of coefficients from the l rank tensors of length
 *      2xl+1 before rotation.
 * @param R_out A pointer to the array of coefficients from the l rank tensors of length
 *      2xl+1 after rotation.
 */
extern void single_wigner_rotation(const int l, const double *euler_angles,
                                   const void *R_in, void *R_out);

/**
 * ❌ Performs wigner rotations on a batch of wigner matrices and initial tensor
 * orientation. The wigner matrices corresponds to the beta orientations. The
 * orientations cover either an octant, hemisphere, or the sphere surface. This is
 * specified by the value of the `n_octant` variables.
 *
 *    n_octant = 1 for octant
 *    n_octant = 4 for hemisphere
 *    n_octant = 8 for sphere
 *
 * The function performs both second rank and fourth rank wigner rotations.
 *
 * @param octant_orientations Number of orientations on an octant.
 * @param n_octants Number of octants.
 * @param wigner_2j_matrices A pointer to a stack of 5x5 second rank wigner matrices.
 * @param R2 A pointer to the second rank tensor coefficients of length 5 to be rotated.
 * @param wigner_4j_matrices A pointer to a stack of 9x9 fourth rank wigner matrices.
 * @param R4 A pointer to the fourth rank tensor coefficients of length 9 to be rotated.
 * @param exp_Im_alpha A pointer to a `4 x octant_orientations` array with the exp(-Imα)
 *      with `octant_orientations` as the leading dimension, ordered as m=[-4,-3,-2,-1].
 * @param w2 A pointer to a stack of second rank tensor coefficients after rotation with
 *      second-rank wigner matrices. The length of w2 is `octant_orientations x
 *      n_octants x 5` with 5 as the leading dimension.
 * @param w4 A pointer to a stack of fourth rank tensor coefficients after rotation with
 *      fourth-rank wigner matrices. The length of w4 is `octant_orientations x
 *      n_octants x 9` with 9 as the leading dimension.
 */
extern void __batch_wigner_rotation(const unsigned int octant_orientations,
                                    const unsigned int n_octants,
                                    double *wigner_2j_matrices, complex128 *R2,
                                    double *wigner_4j_matrices, complex128 *R4,
                                    complex128 *exp_Im_alpha, complex128 *w2,
                                    complex128 *w4);

/**
 * ✅ Calculates exp(-Im alpha) where alpha is an array of size n.
 * The function accepts cos_alpha = cos(alpha).
 * The result is stored in exp_Im_alpha as m x n matrix where m = [-4,-3,-2,-1]
 */
extern void get_exp_Im_alpha(const unsigned int n, const bool allow_fourth_rank,
                             void *exp_Im_alpha);
