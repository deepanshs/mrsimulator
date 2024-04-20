// -*- coding: utf-8 -*-
//
//  wigner_element.c
//
//  @copyright Deepansh J. Srivastava, 2019-2021.
//  Created by Deepansh J. Srivastava, Mar 8, 2021.
//  Contact email = srivastava.89@osu.edu
//

#include "angular_momentum/wigner_element.h"

double SQRT_2_INVERSE = 0.7071067811865475;
double SQRT_3 = 1.732050807568877;

/* This routine calculates the factorial of x */
static inline double fac(double x) {
  if (x < 0) exit(1);

  int ix = (int)x;
  double sum = 1.0;
  while (ix > 1) sum *= ix--;
  return sum;
}

/* power function */
static inline double my_power(double x, int n) {
  if (n == 0) return 1.0;
  double temp = 1.0;
  while (n-- >= 1) temp *= x;
  return temp;
}

// Generic wigner d element for a given angular momentum l, m1, m2.
static inline double __generic_wigner_d_element(const float l, const float m1,
                                                const float m2, const double beta) {
  double sx = sin(beta / 2.);
  double cx = cos(beta / 2.);
  double sum = 0.0, k1, k2, k3, x, y;
  int sign = 1;
  int k;

  for (k = 0; k <= (int)(l - m1); k++) {
    k1 = (int)(l - m1 - k);
    k2 = (int)(l + m2 - k);
    k3 = (int)(k + m1 - m2);

    if (k1 >= 0 && k2 >= 0 && k3 >= 0) {
      int n1 = (int)(2 * l + m2 - m1 - 2 * k);
      int n2 = (int)(m1 - m2 + 2 * k);
      x = my_power(cx, n1);
      y = my_power(sx, n2);
      sum += sign * x * y /
             (fac((double)k1) * fac((double)k2) * fac((double)k3) * fac((double)k));
    }
    sign = -sign;
  }
  double f = fac(l + m2) * fac(l - m2) * fac(l + m1) * fac(l - m1);
  sign = ((int)(m1 - m2) % 2 == 0) ? 1 : -1;
  f = sqrt(f);
  return (sign * sum * f);
}

// Wigner d^{1/2} (m1, m2) elements.
static inline double wigner_one_half_d_elements(const float m1, const float m2,
                                                const double beta) {
  short m1_ = (short)(2 * m1), m2_ = (short)(2 * m2);
  double b2 = 0.5 * beta;
  switch (m1_) {
  case -1:             // -1/2
    switch (m2_) {     //
    case -1:           // d(-1/2 -1/2)
      return cos(b2);  //
    case 1:            // d(-1/2 +1/2)
      return sin(b2);
    }
  case 1:               // +1/2
    switch (m2_) {      //
    case -1:            // d(+1/2 -1/2)
      return -sin(b2);  //
    case 1:             // d(+1/2 +1/2)
      return cos(b2);
    }
  }
  return 0;
}
// ---------------------------------------------------------------------------------- //

// Wigner d^{1} (m1, m2) elements.
static inline double __wigner_one_m1_1(short m1, short m2, const double beta) {
  double m = 1.0;
  if (m1 > 0) {
    m2 *= -1;
    m = ((short)(m1 - m2) % 2 == 0) ? 1.0 : -1.0;
  }
  switch (m2) {                             //
  case -1:                                  // d(-1 -1)
    return m * 0.5 * (1.0 + cos(beta));     //
  case 0:                                   // d(-1  0)
    return m * SQRT_2_INVERSE * sin(beta);  //
  case 1:                                   // d(-1  1)
    return m * 0.5 * (1.0 - cos(beta));     //
  }
  return 0;
}

static inline double wigner_one_d_elements(const float m1, const float m2,
                                           const double beta) {
  short m1_ = (short)m1, m2_ = (short)m2;
  switch (m1_) {
  case -1:  // -1
  case 1:
    return __wigner_one_m1_1(m1_, m2_, beta);
  case 0:                                  // 0
    switch (m2_) {                         //
    case -1:                               // d(0 -1)
      return -SQRT_2_INVERSE * sin(beta);  //
    case 0:                                // d(0  0)
      return cos(beta);                    //
    case 1:                                // d(0  1)
      return SQRT_2_INVERSE * sin(beta);   //
    }
  }
  return 0;
}
// ---------------------------------------------------------------------------------- //

// Wigner d^{3/2} (m1, m2) elements.
static inline double __wigner_three_half_m1_3(short m1, short m2, const double b2,
                                              const double b3) {
  double m = 1.0;
  if (m1 > 0) {
    m2 *= -1;
    m = ((short)((m1 - m2) / 2) % 2 == 0) ? -1.0 : 1.0;
  }
  switch (m2) {                                      //
  case -3:                                           // d(-3/2 -3/2)
    return m * 0.25 * (3.0 * cos(b2) + cos(b3));     //
  case -1:                                           // d(-3/2 -1/2)
    return m * 0.25 * SQRT_3 * (sin(b2) + sin(b3));  //
  case 1:                                            // d(-3/2 +1/2)
    return m * 0.25 * SQRT_3 * (cos(b2) - cos(b3));  //
  case 3:                                            // d(-3/2 +3/2)
    return m * 0.25 * (3.0 * sin(b2) - sin(b3));     //
  }
  return 0;
}

static inline double __wigner_three_half_m1_1(short m1, short m2, const double b2,
                                              const double b3) {
  double m = 1.0;
  if (m1 > 0) {
    m2 *= -1;
    m = ((short)((m1 - m2) / 2) % 2 == 0) ? -1.0 : 1.0;
  }
  switch (m2) {                                       //
  case -3:                                            // d(-1/2 -3/2)
    return m * -0.25 * SQRT_3 * (sin(b2) + sin(b3));  //
  case -1:                                            // d(-1/2  -1/2)
    return m * 0.25 * (cos(b2) + 3.0 * cos(b3));      //
  case 1:                                             // d(-1/2  +1/2)
    return m * -0.25 * (sin(b2) - 3.0 * sin(b3));     //
  case 3:                                             // d(-1/2  +3/2)
    return m * 0.25 * SQRT_3 * (cos(b2) - cos(b3));   //
  }
  return 0;
}

static inline double wigner_three_half_d_elements(const float m1, const float m2,
                                                  const double beta) {
  short m1_ = (short)(2 * m1), m2_ = (short)(2 * m2);
  double b2 = 0.5 * beta, b3 = 1.5 * beta;
  switch (m1_) {
  case -3:  // -3/2
  case 3:   // 3/2
    return __wigner_three_half_m1_3(m1_, m2_, b2, b3);
  case -1:  // -1/2
  case 1:   // 1/2
    return __wigner_three_half_m1_1(m1_, m2_, b2, b3);
  }
  return 0;
}
// ---------------------------------------------------------------------------------- //

// Wigner d^{2} (m1, m2) elements.
static inline double __wigner_two_m1_2(short m1, short m2, const double cx,
                                       const double sx) {
  double m = 1.0;
  if (m1 > 0) {
    m2 *= -1;
    m = ((short)(m1 - m2) % 2 == 0) ? 1.0 : -1.0;
  }
  switch (m2) {                            // -2
  case -2:                                 //
    return m * (1 + cx) * (1. + cx) / 4.;  // d(-2 -2)
  case -1:                                 //
    return m * sx * (1. + cx) / 2.;        // d(-2 -1)
  case 0:                                  //
    return m * 0.6123724355 * sx * sx;     // d(-2 +0)
  case 1:                                  //
    return m * sx * (1. - cx) / 2.;        // d(-2 +1)
  case 2:                                  //
    return m * (1 - cx) * (1. - cx) / 4.;  // d(-2 +2)
  }
  return 0;
}

static inline double __wigner_two_m1_1(short m1, short m2, const double cx,
                                       const double sx) {
  double m = 1.0;
  if (m1 > 0) {
    m2 *= -1;
    m = ((short)(m1 - m2) % 2 == 0) ? 1.0 : -1.0;
  }
  switch (m2) {                                // -1
  case -2:                                     //
    return m * -sx * (1 + cx) / 2.;            // d(-1 -2)
  case -1:                                     //
    return m * (2 * cx * cx + cx - 1.) / 2.;   // d(-1 -1)
  case 0:                                      //
    return m * 1.224744871 * sx * cx;          // d(-1 +0)
  case 1:                                      //
    return m * -(2 * cx * cx - cx - 1.) / 2.;  // d(-1 +1)
  case 2:                                      //
    return m * sx * (1 - cx) / 2.;             // d(-1 +2)
  }
  return 0;
}

static inline double wigner_two_d_elements(const float m1, const float m2,
                                           const double beta) {
  short m1_ = (short)m1, m2_ = (short)m2;
  double cx = cos(beta), sx = sin(beta);
  switch (m1_) {
  case -2:  // -2
  case 2:   // 2
    return __wigner_two_m1_2(m1_, m2_, cx, sx);
  case -1:  // -1
  case 1:   // 1
    return __wigner_two_m1_1(m1_, m2_, cx, sx);
  case 0:
    switch (m2_) {                    // +0
    case -2:                          //
      return 0.6123724355 * sx * sx;  // d(+0 -2)
    case -1:                          //
      return -1.224744871 * sx * cx;  // d(+0 -1)
    case 0:                           //
      return 1.5 * cx * cx - .5;      // d(+0 +0)
    case 1:                           //
      return 1.224744871 * sx * cx;   // d(+0 +1)
    case 2:                           //
      return 0.6123724355 * sx * sx;  // d(+0 +2)
    }
  }
  return 0;
}
// ---------------------------------------------------------------------------------- //

static inline double __wigner_d_element(const float l, const float m1, const float m2,
                                        const double beta) {
  short l1 = (short)(2 * l);
  switch (l1) {
  case 1:  // wigner j=1/2 elements
    return wigner_one_half_d_elements(m1, m2, beta);
  case 2:  // wigner j=1 elements
    return wigner_one_d_elements(m1, m2, beta);
  case 3:  // wigner j=3/2 elements
    return wigner_three_half_d_elements(m1, m2, beta);
  case 4:  // wigner j=2 elements
    return wigner_two_d_elements(m1, m2, beta);
  default:  // remaining
    return __generic_wigner_d_element(l, m1, m2, beta);
  }
}

double wigner_d_element(const float l, const float m1, const float m2,
                        const double beta) {
  return __wigner_d_element(l, m1, m2, beta);
}

// Function has been rewritten using general euler angles instead of angle and phase
// /** Computes the connect factor of transitions |m1_f >< m1_i | --> | m2_f > < m2_i |
//  * when the rotation vector lies in the XY-plane. Said equivalently, the Euler
//  * angles defining the SO(3) rotation has the relation `alpha = -gamma`
//  */
// void transition_connect_factor(const float l, const float m1_f, const float m1_i,
//                                const float m2_f, const float m2_i, const double
//                                theta, const double phi, double *restrict factor) {
//   short delta_p = (short)((m2_f - m2_i) - (m1_f - m1_i));
//   double scale, phase = (double)delta_p * phi;
//   double cx = cos(phase), sx = sin(phase), re = 0.0, im = 0.0, temp;
//   delta_p %= 4;
//   delta_p += 4;
//   delta_p %= 4;

//   scale = __wigner_d_element(l, m2_f, m1_f, theta);
//   scale *= __wigner_d_element(
//       l, m2_i, m1_i,
//       theta);  // Why is this not negative theta or m2_i and m1_i not reversed?

//   switch (delta_p) {
//   case 0:  // * 1
//     re = scale * cx;
//     im = scale * sx;
//     break;
//   case 1:  // * -I
//     re = scale * sx;
//     im = -scale * cx;
//     break;
//   case 2:  // * -1
//     re = -scale * cx;
//     im = -scale * sx;
//     break;
//   case 3:  // * I
//     re = -scale * sx;
//     im = scale * cx;
//     break;
//   }
//   temp = factor[0] * re - factor[1] * im;
//   factor[1] = factor[0] * im + factor[1] * re;
//   factor[0] = temp;
// }

/** Computes the connect factor of transitions |m1_f >< m1_i | --> | m2_f > < m2_i |
 * when the rotation vector does not lie in the XY-plane (i.e. `alpha != -gamma).
 * Alpha, beta, and gamma defined in the ZYZ convention.
 *
 * @param l Spin quantum number for total angular momentum
 * @param m1_f Final angular momentum in the first transition
 * @param m1_i Initial angular momentum in the first transition
 * @param m2_f Final angular momentum in the second transition
 * @param m2_i Initial angular momentum in the second transition
 * @param alpha
 * @param beta
 * @param gamma
 * @param factor
 */
void transition_connect_factor(const float l, const float m1_f, const float m1_i,
                               const float m2_f, const float m2_i, const double alpha,
                               const double beta, const double gamma,
                               double *restrict factor) {
  // Calculate the exponent factors
  double alpha_fac = (double)(m2_i - m2_f) * alpha;
  double gamma_fac = (double)(m1_i - m1_f) * gamma;

  // Calculate the real and imaginary parts
  double re = cos(alpha_fac + gamma_fac), im = sin(alpha_fac + gamma_fac);

  // Calculate the magnitude of the transition factor
  double scale = __wigner_d_element(l, m2_f, m1_f, beta);
  scale *= __wigner_d_element(l, m2_i, m1_i, beta);

  factor[0] = scale * re;
  factor[1] = scale * im;
}

// void transition_connect_factor(const float *l, const float *transition_inital,
//                                const float *transition_final, const double theta,
//                                const double phi, double *restrict factor, int
//                                n_sites) {
//   double m1_f, m1_i, m2_f, m2_i;
//   complex128 *weight = malloc_complex128(n_sites);
//   for (int i = 1; i < n_sites; i++) {
//     m1_i = transition_inital[i];
//     m1_f = (transition_inital + n_sites)[i];
//     m2_i = transition_final[i];
//     m2_f = (transition_final + n_sites)[i];
//     __transition_connect_factor(l[i], m1_f, m1_i, m2_f, m2_i, theta, phi,
//     &weight[i]);
//   }
// }
