//
//  my_complex.h
//
//  Created by Deepansh J. Srivastava, Jun 30, 2019
//  Copyright © 2019 Deepansh J. Srivastava. All rights reserved.
//  Contact email = srivastava.89@osu.edu, deepansh2012@gmail.com
//

#include <math.h>

// inline complex128 cconj(complex128 a) {
//   complex128 z = {a.real, -a.imag};
//   return z;
// }
static inline complex128 cadd(complex128 a, complex128 b) {
  complex128 z = {a.real + b.real, a.imag + b.imag};
  return z;
}
static inline complex128 csub(complex128 a, complex128 b) {
  complex128 z = {a.real - b.real, a.imag - b.imag};
  return z;
}
static inline complex128 cmult(complex128 a, complex128 b) {
  complex128 z;
  z.real = a.real * b.real - a.imag * b.imag;
  z.imag = a.real * b.imag + a.imag * b.real;
  return z;
}
static inline double ccabs(complex128 a) {
  double res = sqrt(a.real * a.real + a.imag * a.imag);
  return res;
}
static inline complex128 cdiv(complex128 a, complex128 b) {
  complex128 t;
  t.real = b.real;
  t.imag = -b.imag;
  complex128 z = cmult(a, t);
  double res = b.real * b.real + b.imag * b.imag;
  z.real /= res;
  z.imag /= res;
  return z;
}
