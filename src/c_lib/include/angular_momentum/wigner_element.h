// -*- coding: utf-8 -*-
//
//  wigner_element.h
//
//  @copyright Deepansh J. Srivastava, 2019-2021.
//  Created by Deepansh J. Srivastava, Mar 8, 2021.
//  Contact email = srivastava.89@osu.edu
//

#include "config.h"

/** ✅  Tested with pytest
 * @brief Evaluates @f$d^{l}_{m_1, m_2}(\beta)@f$ wigner-d element of the given angle
 * @f$\beta@f$.
 *
 * @param l The rank of the wigner-d matrix element.
 * @param m1 The quantum number @f$m_1@f$.
 * @param m2 The quantum number @f$m_2@f$.
 * @param beta The angle @f$\beta@f$ given in radian.
 * @return The wigner-d element, @f$d^{l}_{m_1, m_2}(\beta)@f$.
 */
extern double wigner_d_element(const float l, const float m1, const float m2,
                               const double beta);

// /**
//  * @brief Evaluate the probability of connecting two transitions driven by an
//  external
//  *      rf pulse. The connected transitions are |m1_f >< m1_i | --> | m2_f > < m2_i
//  |.
//  *
//  * @param l The angular momentum quantum number of the spin involved in the
//  transition.
//  * @param m1_f Final quantum number of the starting transition.
//  * @param m1_i Initial quantum number of the starting transition.
//  * @param m2_f Final quantum number of the connecting transition.
//  * @param m2_i Initial quantum number of the connecting transition.
//  * @param theta The tip-angle of the rf pulse.
//  * @param phi The phase of the rf pulse.
//  * @param factor The complex probability of connection.
//  */
// extern void transition_connect_factor(const float l, const float m1_f, const float
// m1_i,
//                                       const float m2_f, const float m2_i,
//                                       const double theta, const double phi,
//                                       double *restrict factor);

/**
 * @brief Evaluate the probability of connecting two transitions by a rotation
 *      defined by the euler angles alpha, beta, and gamma in the ZYZ convention.
 *      The connected transitions are | m1_f >< m1_i | --> | m2_f > < m2_i |.
 *
 * @param l The angular momentum quantum number of the spin involved in the transition.
 * @param m1_f Final quantum number of the starting transition.
 * @param m1_i Initial quantum number of the starting transition.
 * @param m2_f Final quantum number of the connecting transition.
 * @param m2_i Initial quantum number of the connecting transition.
 * @param alpha First euler angle.
 * @param beta Second euler angle.
 * @param gamma Third euler angle.
 * @param factor The complex probability of connection.
 */
extern void transition_connect_factor(const float l, const float m1_f, const float m1_i,
                                      const float m2_f, const float m2_i,
                                      const double alpha, const double beta,
                                      const double gamma, double *restrict factor);
