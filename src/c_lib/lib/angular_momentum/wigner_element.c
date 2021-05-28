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
  while (ix-- > 1) sum *= ix;
  return sum;
}

/* power function */
static inline double my_power(double x, int n) {
  if (n == 0) return 1.0;
  double temp = 1.0;
  while (n-- >= 1) temp *= x;
  return temp;
}

// Generic wigner d element for a given angulae momentum l, m1, m2.
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
  f = sqrt(f);
  return (sum * f);
}

// Wigner d^{1/2} (m1, m2) elements.
static inline double __wigner_one_half_d_elements(const float m1, const float m2,
                                                  const double beta) {
  unsigned int m1_ = (unsigned int)2 * m1, m2_ = (unsigned int)2 * m2;
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

// Wigner d^{1} (m1, m2) elements.
static inline double __wigner_one_d_elements(const float m1, const float m2,
                                             const double beta) {
  unsigned int m1_ = (unsigned int)m1, m2_ = (unsigned int)m2;
  switch (m1_) {
  case -1:                                // -1
    switch (m2_) {                        //
    case -1:                              // d(-1 -1)
      return 0.5 * (1.0 + cos(beta));     //
    case 0:                               // d(-1  0)
      return SQRT_2_INVERSE * sin(beta);  //
    case 1:                               // d(-1  1)
      return 0.5 * (1.0 - cos(beta));     //
    }
  case 0:                                  // 0
    switch (m2_) {                         //
    case -1:                               // d(0 -1)
      return -SQRT_2_INVERSE * sin(beta);  //
    case 0:                                // d(0  0)
      return cos(beta);                    //
    case 1:                                // d(0  1)
      return SQRT_2_INVERSE * sin(beta);   //
    }
  case 1:                                  // 1
    switch (m2_) {                         //
    case -1:                               // d(1 -1)
      return 0.5 * (1.0 - cos(beta));      //
    case 0:                                // d(1  0)
      return -SQRT_2_INVERSE * sin(beta);  //
    case 1:                                // d(1  1)
      return 0.5 * (1.0 + cos(beta));      //
    }
  }
  return 0;
}

// Wigner d^{3/2} (m1, m2) elements.
static inline double __wigner_three_half_d_elements(const float m1, const float m2,
                                                    const double beta) {
  unsigned int m1_ = (unsigned int)2 * m1, m2_ = (unsigned int)2 * m2;
  double b2 = 0.5 * beta, b3 = 1.5 * beta;
  switch (m1_) {
  case -3:                                          // -3/2
    switch (m2_) {                                  //
    case -3:                                        // d(-3/2 -3/2)
      return 0.25 * (3.0 * cos(b2) + cos(b3));      //
    case -1:                                        // d(-3/2 -1/2)
      return 0.25 * SQRT_3 * (sin(b2) + sin(b3));   //
    case 1:                                         // d(-3/2 +1/2)
      return 0.25 * SQRT_3 * (cos(b2) - cos(b3));   //
    case 3:                                         // d(-3/2 +3/2)
      return 0.25 * (3.0 * sin(b2) - sin(b3));      //
    }                                               //
  case -1:                                          // -1/2
    switch (m2_) {                                  //
    case -3:                                        // d(-1/2 -3/2)
      return -0.25 * SQRT_3 * (sin(b2) + sin(b3));  //
    case -1:                                        // d(-1/2  -1/2)
      return 0.25 * (cos(b2) + 3.0 * cos(b3));      //
    case 1:                                         // d(-1/2  +1/2)
      return 0.25 * (-sin(b2) + 3.0 * sin(b3));     //
    case 3:                                         // d(-1/2  +3/2)
      return 0.25 * SQRT_3 * (cos(b2) - cos(b3));   //
    }                                               //
  case 1:                                           // +1/2
    switch (m2_) {                                  //
    case -3:                                        // d(+1/2 -3/2)
      return 0.25 * SQRT_3 * (cos(b2) - cos(b3));   //
    case -1:                                        // d(+1/2  -1/2)
      return 0.25 * (sin(b2) - 3.0 * sin(b3));      //
    case 1:                                         // d(+1/2  +1/2)
      return 0.25 * (cos(b2) + 3.0 * cos(b3));      //
    case 3:                                         // d(+1/2  +3/2)
      return 0.25 * SQRT_3 * (sin(b2) + sin(b3));   //
    }                                               //
  case 3:                                           // +3/2
    switch (m2_) {                                  //
    case -3:                                        // d(+3/2 -3/2)
      return 0.25 * (-3.0 * sin(b2) + sin(b3));     //
    case -1:                                        // d(+3/2 -1/2)
      return 0.25 * SQRT_3 * (cos(b2) - cos(b3));   //
    case 1:                                         // d(+3/2 +1/2)
      return -0.25 * SQRT_3 * (sin(b2) + sin(b3));  //
    case 3:                                         // d(+3/2 +3/2)
      return 0.25 * (3.0 * cos(b2) + cos(b3));      //
    }
  }
  return 0;
}

// Wigner d^{2} (m1, m2) elements.
static inline double __wigner_two_d_elements(const float m1, const float m2,
                                             const double beta) {
  unsigned int m1_ = (unsigned int)m1, m2_ = (unsigned int)m2;
  double cx = cos(beta), sx = sin(beta);
  switch (m1_) {
  case -2:
    switch (m2_) {                           // -2
    case -2:                                 //
      return (1 + cx) * (1. + cx) / 4.;      // d(-2 -2)
    case -1:                                 //
      return sx * (1. + cx) / 2.;            // d(-2 -1)
    case 0:                                  //
      return 0.6123724355 * sx * sx;         // d(-2 +0)
    case 1:                                  //
      return sx * (1. - cx) / 2.;            // d(-2 +1)
    case 2:                                  //
      return (1 - cx) * (1. - cx) / 4.;      // d(-2 +2)
    }                                        //
  case -1:                                   //
    switch (m2_) {                           // -1
    case -2:                                 //
      return -sx * (1 + cx) / 2.;            // d(-1 -2)
    case -1:                                 //
      return (2 * cx * cx + cx - 1.) / 2.;   // d(-1 -1)
    case 0:                                  //
      return 1.224744871 * sx * cx;          // d(-1 +0)
    case 1:                                  //
      return -(2 * cx * cx - cx - 1.) / 2.;  // d(-1 +1)
    case 2:                                  //
      return sx * (1 - cx) / 2.;             // d(-1 +2)
    }                                        //
  case 0:                                    //
    switch (m2_) {                           // +0
    case -2:                                 //
      return 0.6123724355 * sx * sx;         // d(+0 -2)
    case -1:                                 //
      return -1.224744871 * sx * cx;         // d(+0 -1)
    case 0:                                  //
      return 1.5 * cx * cx - .5;             // d(+0 +0)
    case 1:                                  //
      return 1.224744871 * sx * cx;          // d(+0 +1)
    case 2:                                  //
      return 0.6123724355 * sx * sx;         // d(+0 +2)
    }                                        //
  case 1:                                    //
    switch (m2_) {                           // +1
    case -2:                                 //
      return -sx * (1 - cx) / 2.;            // d(+1 -2)
    case -1:                                 //
      return -(2 * cx * cx - cx - 1.) / 2.;  // d(+1 -1)
    case 0:                                  //
      return -1.224744871 * sx * cx;         // d(+1 +0)
    case 1:                                  //
      return (2 * cx * cx + cx - 1.) / 2.;   // d(+1 +1)
    case 2:                                  //
      return sx * (1 + cx) / 2.;             // d(+1 +2)
    }                                        //
  case 2:                                    //
    switch (m2_) {                           // +2
    case -2:                                 //
      return (1 - cx) * (1. - cx) / 4.;      // d(+2 -2)
    case -1:                                 //
      return -sx * (1. - cx) / 2.;           // d(+2 -1)
    case 0:                                  //
      return 0.6123724355 * sx * sx;         // d(+2 +0)
    case 1:                                  //
      return -sx * (1. + cx) / 2.;           // d(+2 +1)
    case 2:                                  //
      return (1 + cx) * (1. + cx) / 4.;      // d(+2 +2)
    }                                        //
  }
  return 0;
}

static inline double __wigner_d_element(const float l, const float m1, const float m2,
                                        const double beta) {
  unsigned int l1 = (unsigned int)(2 * l);
  switch (l1) {
  case 1:  // wigner j=1/2 elements
    return __wigner_one_half_d_elements(m1, m2, beta);
  case 2:  // wigner j=1 elements
    return __wigner_one_d_elements(m1, m2, beta);
  case 3:  // wigner j=3/2 elements
    return __wigner_three_half_d_elements(m1, m2, beta);
  case 4:  // wigner j=2 elements
    return __wigner_two_d_elements(m1, m2, beta);
  default:  // remaining
    return __generic_wigner_d_element(l, m1, m2, beta);
  }
}

double wigner_d_element(const float l, const float m1, const float m2,
                        const double beta) {
  return __wigner_d_element(l, m1, m2, beta);
}

// For transition |m1_f >< m1_i | --> | m2_f > < m2_i |
void transition_connect_factor(const float l, const float m1_f, const float m1_i,
                               const float m2_f, const float m2_i, const double theta,
                               const double phi, double *restrict factor) {
  unsigned int delta_p = (unsigned int)((m2_f - m2_i) - (m1_f - m1_i));
  double scale, phase = (double)delta_p * phi;
  double cx = cos(phase), sx = sin(phase);
  delta_p %= 4;

  scale = __wigner_d_element(l, m2_f, m1_f, theta);
  scale *= __wigner_d_element(l, m2_i, m1_i, theta);

  switch (delta_p) {
  case 0:  // * 1
    factor[0] = scale * cx;
    factor[1] = scale * sx;
    break;
  case 1:  // * -I
    factor[0] = scale * sx;
    factor[1] = -scale * cx;
    break;
  case 2:  // * -1
    factor[0] = -scale * cx;
    factor[1] = -scale * sx;
    break;
  case 3:  // * I
    factor[0] = -scale * sx;
    factor[1] = scale * cx;
    break;
  }
}
