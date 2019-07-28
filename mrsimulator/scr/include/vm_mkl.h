//
//  vm_mkl.h
//
//  Created by Deepansh J. Srivastava, Jul 26, 2019
//  Copyright © 2019 Deepansh J. Srivastava. All rights reserved.
//  Contact email = srivastava.89@osu.edu, deepansh2012@gmail.com
//

// static inline complex128 cmult(complex128 x, complex128 y) {
//   complex128 res;
//   res.real = x.real * y.real - x.imag * y.imag;
//   res.imag = x.real * y.imag + x.imag * y.real;
//   return res;
// }

/**
 * Add the elements of vector x and y and store in res of type double.
 * res = x + y
 */
static inline void vmd_add(int count, double *x, double *y, double *res) {
  vdAdd(count, x, y, res);
}

/**
 * Subtract the elements of vector x from y and store in res of type double.
 * res = x - y
 */
static inline void vmd_sub(int count, double *x, double *y, double *res) {
  vdSub(count, x, y, res);
}

/**
 * Multiply the elements of vector x and y and store in res of type double.
 * res = x * y
 */
static inline void vmd_mul(int count, double *x, double *y, double *res) {
  vdMul(count, x, y, res);
}

/**
 * Divide the elements of vector x by y and store in res of type double.
 * res = x / y
 */
static inline void vmd_div(int count, double *x, double *y, double *res) {
  vdDiv(count, x, y, res);
}

/**
 * Square the elements of vector x and store in res of type double.
 * res = x * x
 */
static inline void vmd_sqr(int count, double *x, double *res) {
  vmdSqr(count, x, res, VML_EP);
}

/**
 * Square root of the elements of vector x stored in res of type double.
 * res = sqrt(x)
 */
static inline void vmd_sqrt(int count, double *x, double *res) {
  vmdSqrt(count, x, res, VML_EP);
}

/**
 * Multiply the elements of vector x and y and store in res of type double
 * complex.
 * res = x * y
 */
static inline void vmz_mul(int count, complex128 *x, complex128 *y,
                           complex128 *res) {
  vmzMul(count, x, y, res, VML_EP);
}

// Trignometry

/**
 * Cosine of the elements of vector x stored in res of type double.
 * res = cos(x)
 */
static inline void vmd_cos(int count, double *x, double *res) {
  vmdCos(count, x, res, VML_EP);
}

/**
 * Sine of the elements of vector x stored in res of type double.
 * res = sin(x)
 */
static inline void vmd_sin(int count, double *x, double *res) {
  vmdSin(count, x, res, VML_EP);
}

/**
 * Cosine + I Sine of the elements of vector x in rad and stored in
 * res of type complex128.
 * res = cos(x) + I sin(x)
 */
static inline void vmz_CIS(int count, double *x, complex128 *res) {
  vmzCIS(count, x, res, VML_EP);
}

static inline void vm_dlinear(int count, double *x, double scale, double offset,
                              double *res) {
  vmdLinearFrac(count, x, x, scale, offset, 0.0, 1.0, res, VML_EP);
}
// Exponent

/**
 * Exponent of the elements of vector x stored in res of type double.
 * res = exp(x)
 */
static inline void vm_dexp(int count, double *x, double *res) {
  vmdExp(count, x, res, VML_EP);
}

/**
 * Exponent of the elements of vector x stored in res of type complex128.
 * res = exp(x)
 */
static inline void vm_zexp(int count, complex128 *x, complex128 *res) {
  vmzExp(count, x, res, VML_EP);
}

/**
 * Copy elements of vector x to vector y of type complex128.
 * y = a*x + b*y
 * Equivalent to cblas_daxpby.
 */
static inline void catlas_daxpby(int count, double a, double *x, int stride_x,
                                 double b, double *y, int stride_y) {
  cblas_daxpby(count, a, x, stride_x, b, y, stride_y);
}
