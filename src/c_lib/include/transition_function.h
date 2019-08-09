//
//  transition_function.h
//
//  Created by Deepansh J. Srivastava, Apr 11, 2019
//  Copyright © 2019 Deepansh J. Srivastava. All rights reserved.
//  Contact email = srivastava.89@osu.edu, deepansh2012@gmail.com
//

#include "mrsimulator.h"

// Definning pi^{2,2}_{L,J} as piLJ //
// #define pi01 = 0.3577708764
// #define pi21 = 0.1069044968
// #define pi41 = -0.1434274331

// #define pi03 = 0.8485281374
// #define pi23 = -1.0141851057
// #define pi43 = -1.2850792082
// --------------------------------- //

/**
 * @brief The @f$\mathbb{p}@f$ spin symmetry transition function.
 *
 * The transition symmetry function from irreducible 1st rank tensor, given as
 * @f[
 *    \mathbb{p}(m_f, m_i) &= \left< m_f | T_{10} | m_f \right> -
 *                            \left< m_i | T_{10} | m_i \right> \\
 *    &= m_f - m_i
 * @f]
 * where @f$T_{10}@f$ is the irreducible 1st rank tensor operator in the
 * rotating tilted frame.
 *
 * @param mf The quantum number associated with the final energy state.
 * @param mi The quantum number associated with the initial energy state.
 * @return The spin transition symmetry function @f$\mathbb{p}@f$.
 */
static inline double p(double mf, double mi) { return (mf - mi); }


/**
 * @brief The @f$\mathbb{d}@f$ spin transition symmetry function.
 *
 * The transition symmetry function from irreducible 2nd rank tensor, given as
 * @f[
 *    \mathbb{d}(m_f, m_i) &= \left< m_f | T_{20} | m_f \right> -
 *                            \left< m_i | T_{20} | m_i \right> \\
 *    &= \sqrt{\frac{3}{2}} \left(m_f^2 - m_i^2 \right)
 * @f]
 * where @f$T_{20}@f$ is the irreducible 2nd rank tensor operator in the
 * rotating tilted frame.
 *
 * @param mf The quantum number associated with the final energy state.
 * @param mi The quantum number associated with the initial energy state.
 * @return The spin transition symmetry function @f$\mathbb{d}@f$.
 */
static inline double d(double mf, double mi)
{
    return 1.2247448714 * (mf * mf - mi * mi);
}


/**
 * @brief The @f$\mathbb{f}@f$ spin transition symmetry function.
 *
 * The transition symmetry function from irreducible 3rd rank tensor, given as
 * @f[
 *    \mathbb{f}(m_f, m_i) &= \left< m_f | T_{30} | m_f \right> -
 *                            \left< m_i | T_{30} | m_i \right> \\
 *    &= \frac{1}{\sqrt{10}} [5(m_f^3 - m_i^3) + (1 - 3I(I+1))(m_f-m_i)]
 * @f]
 * where @f$T_{30}@f$ is the irreducible 3rd rank tensor operator in the
 * rotating tilted frame.
 *
 * @param mf The quantum number associated with the final energy state.
 * @param mi The quantum number associated with the initial energy state.
 * @return The spin transition symmetry function @f$\mathbb{f}@f$.
 */
static inline double f(double mf, double mi, double spin)
{
    double f_ = 1.0 - 3.0 * spin * (spin + 1.0);
    f_ *= (mf - mi);
    f_ += 5.0 * (mf * mf * mf - mi * mi * mi);
    f_ *= 0.316227766;
    return f_;
}

/**
 * @brief The @f$\mathbb{d_{IS}}@f$ spin transition symmetry function.
 *
 * The transition symmetry function from irreducible tensors for a coupled
 * spin system is defined as
 * @f[
 *   \mathbb{d}_{IS}(m_{If}, m_{Sf}, m_{Ii}, m_{Si}) &=
 *        \left< m_{If} m_{Sf} | T_{10}(I)~T_{10}(S) | m_{If} m_{Sf} \right> \\
 *   &~~- \left< m_{Ii} m_{Si} | T_{10}(I)~T_{10}(S) | m_{Ii} m_{Si} \right> \\
 *   &= m_{If} m_{Sf} - m_{Ii} m_{Si}
 * @f]
 * where @f$T_{10}(I)@f$ and @f$T_{10}(S)@f$ are the irreducible 1st rank tensor
 * operators in the rotating tilted frame for spin I and S, respectively.
 *
 * @param mIi The quantum number associated with the initial state of spin I.
 * @param mIf The quantum number associated with the final state of spin I.
 * @param mSi The quantum number associated with the initial state of spin S.
 * @param mSf The quantum number associated with the final state of spin S.
 * @return The spin transition symmetry function @f$\mathbb{d_{IS}}@f$.
 */
static inline double dIS(double mIf, double mIi, double mSf, double mSi)
{
    return mIf * mSf - mIi * mSi;
}

static inline void quad_ci(double *c0, double *c2, double *c4, double mf, double mi, double spin)
{
    double f_ = f(mf, mi, spin);
    double p_ = p(mf, mi);

    double temp = spin * (spin + 1.0) - 0.75;
    c0[0] = 0.3577708764 * temp * p_ + 0.8485281374 * f_;
    c2[0] = 0.1069044968 * temp * p_ + -1.0141851057 * f_;
    c4[0] = -0.1434274331 * temp * p_ + -1.2850792082 * f_;
}
