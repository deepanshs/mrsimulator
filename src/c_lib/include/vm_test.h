// -*- coding: utf-8 -*-
//
//  vm.h
//
//  @copyright Deepansh J. Srivastava, 2019-2021.
//  Created by Deepansh J. Srivastava, Jul 26, 2019.
//  Contact email = srivastava.89@osu.edu
//

#include "vm.h"

#ifndef DOXYGEN_SHOULD_SKIP_THIS

/** Arithmetic suit ======================================================== */

/**
 * Absolute double test
 * res = |x|
 */
double test_vm_absd(double a) { return absd(a); }

/**
 * Add the elements of vector x and y and store in res of type double.
 *      res = x + y
 */
void test_vm_double_add(int count, const double *restrict x, const double *restrict y,
                        double *restrict res) {
  vm_double_add(count, x, y, res);
}

/**
 * Add the elements of vector y inplace with the elements from vector x.
 *      y += x
 */
void test_vm_double_add_inplace(int count, const double *restrict x,
                                double *restrict y) {
  vm_double_add_inplace(count, x, y);
}

/**
 * Subtract the elements of vector x from y and store in res of type double.
 *      res = x - y
 */
void test_vm_double_subtract(int count, const double *restrict x,
                             const double *restrict y, double *restrict res) {
  vm_double_subtract(count, x, y, res);
}

/**
 * Subtract the elements of vector y inplace with the elements from vector x.
 *      y -= x
 */
void test_vm_double_subtract_inplace(int count, const double *restrict x,
                                     double *restrict y) {
  vm_double_subtract_inplace(count, x, y);
}

/**
 * Multiply the elements of vector x and y and store in res of type double.
 *      res = x * y
 */
void test_vm_double_multiply(int count, const double *restrict x, const int stride_x,
                             const double *restrict y, double *restrict res) {
  vm_double_multiply(count, x, stride_x, y, res);
}

/**
 * Multiply the elements of vector y inplace with the elements from vector x.
 *      y *= x
 */
void test_vm_double_multiply_inplace(int count, const double *restrict x,
                                     const int stride_x, double *restrict y,
                                     const int stride_y) {
  vm_double_multiply_inplace(count, x, stride_x, y, stride_y);
}

/**
 * Divide the elements of vector x by y and store in res of type double.
 *      res = x / y
 */
void test_vm_double_divide(int count, const double *restrict x,
                           const double *restrict y, double *restrict res) {
  vm_double_divide(count, x, y, res);
}

/**
 * Divide the elements of vector y inplace with the elements from vector x.
 *      y /= x
 */
void test_vm_double_divide_inplace(int count, const double *restrict x,
                                   double *restrict y) {
  vm_double_divide_inplace(count, x, y);
}

/**
 * Square the elements of vector x and store in res of type double.
 *      res = x * x
 */
void test_vm_double_square(int count, const double *restrict x, double *restrict res) {
  vm_double_square(count, x, res);
}

/**
 * Square the elements of vector y inplace.
 *      x *= x
 */
void test_vm_double_square_inplace(int count, double *restrict x) {
  vm_double_square_inplace(count, x);
}

/** Power and roots suit =================================================== */
/**
 * Square root of the elements of vector x stored in res of type double.
 *    res = sqrt(x)
 */
void test_vm_double_square_root(int count, const double *restrict x,
                                double *restrict res) {
  vm_double_square_root(count, x, res);
}

/**
 * Square root of the elements of vector x inplace.
 *      x = sqrt(x)
 */
void test_vm_double_square_root_inplace(int count, double *restrict x) {
  vm_double_square_root_inplace(count, x);
}

/**
 * Multiply the elements of vector x and y and store in res of type double
 * complex.
 *    res = x * y
 */
void test_vm_double_complex_multiply(int count, const void *restrict x,
                                     const void *restrict y, void *restrict res) {
  vm_double_complex_multiply(count, x, y, res);
}

/**
 * Multiply the elements of vector x and conj(y) and store in res of type double
 * complex.
 *    res = x * conj(y)
 */
void test_vm_double_complex_conj_multiply(int count, const void *restrict x,
                                          const void *restrict y, void *restrict res) {
  vm_double_complex_conj_multiply(count, x, y, res);
}

// Trignometry

/**
 * Cosine of the elements of vector x stored in res of type double.
 *    res = cos(x)
 */
void test_vm_double_cosine(int count, const double *restrict x, double *restrict res) {
  vm_double_cosine(count, x, res);
}

/**
 * Sine of the elements of vector x stored in res of type double.
 *      res = sin(x)
 */
void test_vm_double_sine(int count, const double *restrict x, double *restrict res) {
  vm_double_sine(count, x, res);
}

/**
 * Elementwise (Cosine + I Sine) of vector x in rad. Stores result in 'res' with
 * complex128 type. res = cos(x) + I sin(x)
 */
void test_vm_cosine_I_sine(int count, const double *restrict x, void *restrict res) {
  vm_cosine_I_sine(count, x, res);
}

// Exponent
/**
 * Exponent of the elements of vector x stored in res of type double.
 *      res = exp(x)
 */
void test_vm_double_exp(int count, double *restrict x, double *restrict res,
                        const int ix) {
  vm_double_exp(count, x, res, ix);
}

/**
 * Exponent of the elements of vector x stored in res of type complex128.
 *      res = exp(x)
 */
void test_vm_double_complex_exp(int count, const void *restrict x, void *restrict res) {
  vm_double_complex_exp(count, x, res);
}

/**
 * Exponent of the elements of vector x stored in res of type complex128.
 *      res = exp(x(imag))
 */
void test_vm_double_complex_exp_imag_only(int count, const void *restrict x,
                                          void *restrict res) {
  vm_double_complex_exp_imag_only(count, x, res);
}

#endif /* DOXYGEN_SHOULD_SKIP_THIS */
