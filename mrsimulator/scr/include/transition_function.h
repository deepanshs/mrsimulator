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

/*
Transition symmetry function from irreducible 1st rank tensor.
The function evaluates,

.. math::

    \mathbb{p}(m_f, m_i) &= \left< m_f | T_{10} | m_f \right> -
                            \left< m_i | T_{10} | m_i \right> \\
    &= m_f - m_i

where :math:`T_{10}` is the irreducible 1st rank tensor operator in the
rotating tilted frame.

@param double mf = The quantum number associated with the final energy state.
@param double mi = The quantum number associated with the initial energy state.
@rtype = float.
@return double p = The transition symmetry function p.
*/
static inline double p(double mf, double mi) { return (mf - mi); }

/*
Transition symmetry function from irreducible 2nd rank tensor.
The function evaluates,

.. math::

    \mathbb{p}(m_f, m_i) &= \left< m_f | T_{20} | m_f \right> -
                            \left< m_i | T_{20} | m_i \right> \\
    &= \sqrt{\frac{3}{2}} \left(m_f^2 - m_i^2 \right)

where :math:`T_{20}` is the irreducible 2nd rank tensor operator in the
rotating tilted frame.

@param float mf = The quantum number associated with the final energy state.
@param float mi = The quantum number associated with the initial energy state.
@rtype = float.
@return = The transition symmetry function d.
*/
static inline double d(double mf, double mi)
{
    return 1.2247448714 * (mf * mf - mi * mi);
}

/*
Transition symmetry function from irreducible 3rd rank tensor.
The function evaluates,

.. math::

    \mathbb{p}(m_f, m_i) &= \left< m_f | T_{30} | m_f \right> -
                            \left< m_i | T_{30} | m_i \right> \\
    &= \frac{1}{\sqrt{10}} [5(m_f^3 - m_i^3) + (1 - 3I(I+1))(m_f-m_i)]

where :math:`T_{30}` is the irreducible 3rd rank tensor operator in the
rotating tilted frame.

@param float mf = The quantum number associated with the final energy state.
@param float mi = The quantum number associated with the initial energy state.
@rtype = float
@return = The transition symmetry function f.
*/
static inline double f(double mf, double mi, double spin)
{
    double f_ = 1.0 - 3.0 * spin * (spin + 1.0);
    f_ *= (mf - mi);
    f_ += 5.0 * (mf * mf * mf - mi * mi * mi);
    f_ *= 0.316227766;
    return f_;
}

/* Returns the dIS(mIi, mIf, mSi, mSf) transition element.                   *
 * The expression follows,                                                   *
 *     d(mf, mi) = < mIf mSf | T10(I) T10(S) | mIf mSf > -                   *
 *                 < mIi mSi | T10(I) T10(S) | mTi mSi >                     *
 *               = (mIf * mSf - mIi * mSi)                                   *
 *                                                                           *
 * @param double mIi = quantum number for the initial state of spin I.      *
 * @param double mIf = quantum number for the final state of spin I.        *
 * @param double mSi = quantum number for the initial state of spin S.      *
 * @param double mSf = quantum number for the final state of spin S.        *
 * @return double dIS = The value.                                          */
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
