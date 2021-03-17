// -*- coding: utf-8 -*-
//
//  spin_transition_function.h
//
//  @copyright Deepansh J. Srivastava, 2019-2021.
//  Created by Deepansh J. Srivastava, Apr 11, 2019.
//  Contact email = srivastava.89@osu.edu
//

// =====================================================================================
//                      Single nucleus spin transition functions
// =====================================================================================

/**
 * @brief Single nucleus spin transition function from the irreducible spherical tensor
 * of rank @f$L=1@f$ is given as
 * @f[
 *    \mathbb{p}(m_f, m_i) &= \left< m_f | \hat{T}_{10} | m_f \right> -
 *                            \left< m_i | \hat{T}_{10} | m_i \right> \\
 *                         &= m_f - m_i,
 * @f]
 * where @f$\hat{T}_{10}@f$ is the irreducible 1st-rank spherical tensor operator in the
 * rotating tilted frame.
 *
 * @param mi The quantum number associated with the quantized initial energy level.
 * @param mf The quantum number associated with the quantized final energy level.
 * @returns The spin transition function @f$\mathbb{p}@f$.
 */
static inline double STF_p(const double mf, const double mi) { return (mf - mi); }

/**
 * Single nucleus spin transition function from the irreducible spherical tensor of rank
 * @f$L=2@f$ is given as
 * @f[
 *    \mathbb{d}(m_f, m_i) &= \left< m_f | \hat{T}_{20} | m_f \right> -
 *                            \left< m_i | \hat{T}_{20} | m_i \right> \\
 *    &= \sqrt{\frac{3}{2}} \left(m_f^2 - m_i^2 \right),
 * @f]
 * where @f$\hat{T}_{20}@f$ is the irreducible 2nd-rank spherical tensor operator in the
 * rotating tilted frame.
 *
 * @param mi The quantum number associated with the quantized initial energy level.
 * @param mf The quantum number associated with the quantized final energy level.
 * @returns The spin transition function @f$\mathbb{d}@f$.
 */
static inline double STF_d(const double mf, const double mi) {
  return 1.2247448714 * (mf * mf - mi * mi);
}

/**
 * Single nucleus spin transition function from the irreducible spherical tensor of rank
 * @f$L=3@f$ is given as
 * @f[
 *    \mathbb{f}(m_f, m_i) &= \left< m_f | \hat{T}_{30} | m_f \right> -
 *                            \left< m_i | \hat{T}_{30} | m_i \right> \\
 *    &= \frac{1}{\sqrt{10}} [5(m_f^3 - m_i^3) + (1 - 3I(I+1))(m_f-m_i)],
 * @f]
 * where @f$\hat{T}_{30}@f$ is the irreducible 3rd-rank spherical tensor operator in the
 * rotating tilted frame.
 *
 * @param mi The quantum number associated with the quantized initial energy level.
 * @param mf The quantum number associated with the quantized final energy level.
 * @param spin The spin quantum angular momentum number.
 * @return The spin transition function @f$\mathbb{f}@f$.
 */
static inline double STF_f(const double mf, const double mi, const double spin) {
  double f_value = 1.0 - 3.0 * spin * (spin + 1.0);
  f_value *= (mf - mi);
  f_value += 5.0 * (mf * mf * mf - mi * mi * mi);
  f_value *= 0.316227766;
  return f_value;
}

/**
 * The following single nucleus composite spin transition functions corresponding to
 * rank @f$L=[0,2,4]@f$ irreducible tensors results from the second-order corrections to
 * the quadrupole frequency. The functions are defined as
 * @f[
 *   \mathbb{c}_{0}(m_f, m_i) &= \frac{4}{\sqrt{125}} \left[I(I+1) -
 *          \frac{3}{4}\right] \mathbb{p}(m_f, m_i) +
 *          \sqrt{\frac{18}{25}} \mathbb{f}(m_f, m_i), \\
 *   \mathbb{c}_{2}(m_f, m_i) &= \sqrt{\frac{2}{175}} \left[I(I+1) -
 *          \frac{3}{4}\right] \mathbb{p}(m_f, m_i) -
 *          \frac{6}{\sqrt{35}} \mathbb{f}(m_f, m_i), \\
 *   \mathbb{c}_{4}(m_f, m_i) &= -\sqrt{\frac{18}{875}} \left[I(I+1) -
 *          \frac{3}{4}\right] \mathbb{p}(m_f, m_i) -
 *          \frac{17}{\sqrt{175}} \mathbb{f}(m_f, m_i),
 * @f]
 * where @f$\mathbb{p}(m_f, m_i)@f$ and @f$\mathbb{f}(m_f, m_i)@f$ are single nucleus
 * spin transition functions described before, and @f$I@f$ is the spin quantum number.
 *
 * @param mi The quantum number associated with the quantized initial energy level.
 * @param mf The quantum number associated with the quantized final energy level.
 * @param spin The spin quantum number, @f$I@f$.
 * @param cl_value A pointer to an array of size 3 where the spin transition functions,
 *      @f$\mathbb{c}_{L}@f$, will be stored ordered according to @f$L=[0,2,4]@f$.
 */
static inline void STF_cL(double *restrict cl_value, const double mf, const double mi,
                          const double spin) {
  double f_value = STF_f(mf, mi, spin);
  double p_value = STF_p(mf, mi);
  double temp = spin * (spin + 1.0) - 0.75;
  temp *= p_value;

  *cl_value++ = 0.3577708764 * temp + 0.8485281374 * f_value;
  *cl_value++ = 0.1069044968 * temp + -1.0141851057 * f_value;
  *cl_value++ = -0.1434274331 * temp + -1.2850792082 * f_value;
}

// =====================================================================================
//              Two weakly coupled nuclei spin transition functions
// =====================================================================================
/**
 * @brief The @f$\mathbb{d}_{IS}@f$ spin transition symmetry function.
 *
 * Two weakly coupled nuclei spin transition function from the irreducible spherical
 * tensors is given as
 * @f[
 *   \mathbb{d}_{IS}(m_{f_I}, m_{f_S}, m_{i_I}, m_{i_S}) &=
 *   \left<m_{f_I}m_{f_S}|\hat{T}_{10}(I)~\hat{T}_{10}(S)|m_{f_I} m_{f_S}\right>
 * -\left<m_{i_I}m_{i_S}|\hat{T}_{10}(I)~\hat{T}_{10}(S)|m_{i_I} m_{i_S}\right>
 * \\
 *   &= m_{f_I} m_{f_S} - m_{i_I} m_{i_S},
 * @f]
 * where @f$\hat{T}_{10}(I)@f$ and @f$\hat{T}_{10}(S)@f$ are the irreducible first-rank
 * spherical tensor operators in the rotating tilted frame for spin I and S,
 * respectively.
 *
 * @param mIf The quantum number associated with the quantized final energy state
 *      corresponding to spin I.
 * @param mSf The quantum number associated with the quantized final energy state
 *      corresponding to spin S.
 * @param mIi The quantum number associated with the quantized initial energy state
 *      corresponding to spin I.
 * @param mSi The quantum number associated with the quantized initial energy state
 *      corresponding to spin S.
 * @return The spin transition symmetry function @f$\mathbb{d}_{IS}@f$.
 */
static inline double STF_dIS(const double mIf, const double mIi, const double mSf,
                             const double mSi) {
  return mIf * mSf - mIi * mSi;
}
