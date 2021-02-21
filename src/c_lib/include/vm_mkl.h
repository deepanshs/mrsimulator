// -*- coding: utf-8 -*-
//
//  vm_mkl.h
//
//  @copyright Deepansh J. Srivastava, 2019-2021.
//  Created by Deepansh J. Srivastava, Jul 26, 2019.
//  Contact email = srivastava.89@osu.edu
//

#ifndef DOXYGEN_SHOULD_SKIP_THIS
/**
 * Add the elements of vector x and y and store in res of type double.
 * res = x + y
 */
static inline void vm_double_add(int count, const double *x, const double *y,
                                 double *res) {
  vdAdd(count, x, y, res);
}

/**
 * Add the elements of vector y inplace with the elements from vector x.
 * y += x
 */
static inline void vm_double_add_inplace(int count, const double *x, double *y) {
  vdAdd(count, x, y, y);
}

/**
 * Subtract the elements of vector x from y and store in res of type double.
 * res = x - y
 */
static inline void vm_double_subtract(int count, const double *x, const double *y,
                                      double *res) {
  vdSub(count, x, y, res);
}

/**
 * Multiply the elements of vector x and y and store in res of type double.
 * res = x * y
 */
static inline void vm_double_multiply(int count, const double *x, const double *y,
                                      double *res) {
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
 * Square the elements of vector y inplace.
 * x *= x
 */
static inline void vm_double_square_inplace(int count, double *x) {
  vmdSqr(count, x, x, VML_EP);
}

/**
 * Square root of the elements of vector x stored in res of type double.
 * res = sqrt(x)
 */
static inline void vm_double_square_root(int count, const double *x, double *res) {
  vdSqrt(count, x, res);
}

/**
 * Multiply the elements of vector x and y and store in res of type double
 * complex.
 * res = x * y
 */
static inline void vm_double_complex_multiply(int count, const void *x, const void *y,
                                              void *res) {
  vzMul(count, x, y, res);
}

// Trignometry

/**
 * Cosine of the elements of vector x stored in res of type double.
 * res = cos(x)
 */
static inline void vm_double_cosine(int count, const double *x, double *res) {
  vdCos(count, x, res);
}

/**
 * Sine of the elements of vector x stored in res of type double.
 * res = sin(x)
 */
static inline void vm_double_sine(int count, const double *x, double *res) {
  vdSin(count, x, res);
}

/**
 * Cosine + I Sine of the elements of vector x in rad and stored in
 * res of type complex128.
 * res = cos(x) + I sin(x)
 */
static inline void vm_cosine_I_sine(int count, const double *x, void *res) {
  vzCIS(count, x, res);
}

static inline void vm_dlinear(int count, double *x, double scale, double offset,
                              double *res) {
  vdLinearFrac(count, x, x, scale, offset, 0.0, 1.0, res);
}
// Exponent

/**
 * Exponent of the elements of vector x stored in res of type double.
 * res = exp(x)
 */
static inline void vm_double_exp(int count, const double *x, double *res) {
  vdExp(count, x, res);
}

/**
 * Exponent of the elements of vector x stored in res of type complex128.
 * res = exp(x)
 */
static inline void vm_double_complex_exp(int count, const void *x, void *res) {
  vmzExp(count, x, res, VML_EP);
}
#endif /* DOXYGEN_SHOULD_SKIP_THIS */
