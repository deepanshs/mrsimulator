

#include "core.h"

inline double fastPow(double a, double b) {
  union {
    double d;
    int x[2];
  } u = { a };
  u.x[0] = (int)(b * (u.x[1] - 1072632447) + 1072632447);
  u.x[1] = 0;
  return u.d;
}

double complex complexMult(double complex *A1, double complex *B1){
  double a, b, c, d, k1, k2, k3;
  a = creal(A1[0]);
  b = cimag(A1[0]);
  c = creal(B1[0]);
  d = cimag(B1[0]);

  k1 = a*(c + d);
  k2 = d*(a + b);   
  k3 = c*(b - a);
  return (k1 - k2) + I * (k1 + k3);
}