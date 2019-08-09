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
static inline void vm_double_add(int count, const double *x, const double *y,
                                 double *res) {
  vdAdd(count, x, y, res);
}

/**
 * Subtract the elements of vector x from y and store in res of type double.
 * res = x - y
 */
static inline void vm_double_subtract(int count, const double *x,
                                      const double *y, double *res) {
  vdSub(count, x, y, res);
}

/**
 * Multiply the elements of vector x and y and store in res of type double.
 * res = x * y
 */
static inline void vm_double_multiply(int count, const double *x,
                                      const double *y, double *res) {
  vdMul(count, x, y, res);
}

/**
 * Divide the elements of vector x by y and store in res of type double.
 * res = x / y
 */
static inline void vm_double_divide(int count, const double *x, const double *y,
                                    double *res) {
  vdDiv(count, x, y, res);
}

/**
 * Square the elements of vector x and store in res of type double.
 * res = x * x
 */
static inline void vm_double_square(int count, const double *x, double *res) {
  vdSqr(count, x, res);
}

/**
 * Square root of the elements of vector x stored in res of type double.
 * res = sqrt(x)
 */
static inline void vm_double_square_root(int count, const double *x,
                                         double *res) {
  vdSqrt(count, x, res);
}

/**
 * Multiply the elements of vector x and y and store in res of type double
 * complex.
 * res = x * y
 */
static inline void vm_double_complex_multiply(int count, const complex128 *x,
                                              const complex128 *y,
                                              complex128 *res) {
  vmzMul(count, x, y, res, VML_EP);
}

// Trignometry

/**
 * Cosine of the elements of vector x stored in res of type double.
 * res = cos(x)
 */
static inline void vm_double_cosine(int count, const double *x, double *res) {
  vmdCos(count, x, res, VML_EP);
}

/**
 * Sine of the elements of vector x stored in res of type double.
 * res = sin(x)
 */
static inline void vm_double_sine(int count, const double *x, double *res) {
  vmdSin(count, x, res, VML_EP);
}

/**
 * Cosine + I Sine of the elements of vector x in rad and stored in
 * res of type complex128.
 * res = cos(x) + I sin(x)
 */
static inline void vm_cosine_I_sine(int count, const double *x,
                                    complex128 *res) {
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
static inline void vm_dexp(int count, const double *x, double *res) {
  vmdExp(count, x, res, VML_EP);
}

/**
 * Exponent of the elements of vector x stored in res of type complex128.
 * res = exp(x)
 */
static inline void vm_double_complex_exp(int count, const complex128 *x,
                                         complex128 *res) {
  vmzExp(count, x, res, VML_EP);
}

/**
 * Copy elements of vector x to vector y of type complex128.
 * y = a*x + b*y
 * Equivalent to cblas_daxpby.
 */
static inline void catlas_daxpby(int count, const double a, const double *x,
                                 const int stride_x, const double b, double *y,
                                 const int stride_y) {
  cblas_daxpby(count, a, x, stride_x, b, y, stride_y);
}
