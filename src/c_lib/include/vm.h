// -*- coding: utf-8 -*-
//
//  vm.h
//
//  @copyright Deepansh J. Srivastava, 2019-2021.
//  Created by Deepansh J. Srivastava, Jul 26, 2019.
//  Contact email = srivastava.89@osu.edu
//

#include <string.h>

#include "config.h"

#ifndef DOXYGEN_SHOULD_SKIP_THIS

static inline double absd(double a) {
  *((unsigned __int64_ *)&a) &= ~(1ULL << 63);
  return a;
}

/** Arithmetic suit ======================================================== */

/**
 * Add the elements of vector x and y and store in res of type double.
 *      res = x + y
 */
static inline void vm_double_add(int count, const double *restrict x,
                                 const double *restrict y, double *restrict res) {
  // x = __builtin_assume_aligned(x, 32);
  // y = __builtin_assume_aligned(y, 32);
  // res = __builtin_assume_aligned(res, 32);
  while (count-- > 0) {
    *res++ = *x++ + *y++;
    // x += stride_x;
    // y += stride_y;
    // res += stride_res;
  }
}

/**
 * Add the elements of vector y inplace with the elements from vector x.
 *      y += x
 */
static inline void vm_double_add_inplace(int count, const double *restrict x,
                                         double *restrict y) {
  // x = __builtin_assume_aligned(x, 32);
  // y = __builtin_assume_aligned(y, 32);
  while (count-- > 0) {
    *y++ += *x++;
    // x += stride_x;
    // y += stride_y;
  }
}

/**
 * Subtract the elements of vector x from y and store in res of type double.
 *      res = x - y
 */
static inline void vm_double_subtract(int count, const double *restrict x,
                                      const double *restrict y, double *restrict res) {
  // x = __builtin_assume_aligned(x, 32);
  // y = __builtin_assume_aligned(y, 32);
  // res = __builtin_assume_aligned(res, 32);
  while (count-- > 0) {
    *res++ = *x++ - *y++;
    // x += stride_x;
    // y += stride_y;
    // res += stride_res;
  }
}

/**
 * Subtract the elements of vector y inplace with the elements from vector x.
 *      y -= x
 */
static inline void vm_double_subtract_inplace(int count, const double *restrict x,
                                              double *restrict y) {
  // x = __builtin_assume_aligned(x, 32);
  // y = __builtin_assume_aligned(y, 32);
  while (count-- > 0) {
    *y++ -= *x++;
    // x += stride_x;
    // y += stride_y;
  }
}

/**
 * Multiply the elements of vector x and y and store in res of type double.
 *      res = x * y
 */
static inline void vm_double_multiply(int count, const double *restrict x,
                                      const double *restrict y, double *restrict res) {
  // x = __builtin_assume_aligned(x, 32);
  // y = __builtin_assume_aligned(y, 32);
  // res = __builtin_assume_aligned(res, 32);
  while (count-- > 0) {
    *res++ = *x++ * *y++;
    // x += stride_x;
    // y += stride_y;
    // res += stride_res;
  }
}

/**
 * Multiply the elements of vector y inplace with the elements from vector x.
 *      y *= x
 */
static inline void vm_double_multiply_inplace(int count, const double *restrict x,
                                              const int stride_x, double *restrict y,
                                              const int stride_y) {
  // x = __builtin_assume_aligned(x, 32);
  // y = __builtin_assume_aligned(y, 32);
  while (count-- > 0) {
    *y *= *x;
    x += stride_x;
    y += stride_y;
  }
}

/**
 * Divide the elements of vector x by y and store in res of type double.
 *      res = x / y
 */
static inline void vm_double_divide(int count, const double *restrict x,
                                    const double *restrict y, double *restrict res) {
  // x = __builtin_assume_aligned(x, 32);
  // y = __builtin_assume_aligned(y, 32);
  // res = __builtin_assume_aligned(res, 32);
  while (count-- > 0) {
    *res++ = *x++ / *y++;
    // x += stride_x;
    // y += stride_y;
    // res += stride_res;
  }
}

/**
 * Divide the elements of vector y inplace with the elements from vector x.
 *      y /= x
 */
static inline void vm_double_divide_inplace(int count, const double *restrict x,
                                            double *restrict y) {
  // x = __builtin_assume_aligned(x, 32);
  // y = __builtin_assume_aligned(y, 32);
  while (count-- > 0) {
    *y++ /= *x++;
    // x += stride_x;
    // y += stride_y;
  }
}

/**
 * Square the elements of vector x and store in res of type double.
 *      res = x * x
 */
static inline void vm_double_square(int count, const double *restrict x,
                                    double *restrict res) {
  // x = __builtin_assume_aligned(x, 32);
  // res = __builtin_assume_aligned(res, 32);
  while (count-- > 0) {
    *res++ = *x * *x;
    x++;
  }
}

/**
 * Square the elements of vector y inplace.
 *      x *= x
 */
static inline void vm_double_square_inplace(int count, double *restrict x) {
  // x = __builtin_assume_aligned(x, 32);
  // y = __builtin_assume_aligned(y, 32);
  while (count-- > 0) {
    *x *= *x;
    x++;
    // x += stride_x;
    // y += stride_y;
  }
}

/** Power and roots suit =================================================== */
/**
 * Square root of the elements of vector x stored in res of type double.
 *    res = sqrt(x)
 */
static inline void vm_double_square_root(int count, const double *restrict x,
                                         double *restrict res) {
  // x = __builtin_assume_aligned(x, 32);
  // res = __builtin_assume_aligned(res, 32);
  while (count-- > 0) {
    *res++ = sqrt(*x++);
  }
}

/**
 * Square root of the elements of vector x inplace.
 *      x = sqrt(x)
 */
static inline void vm_double_square_root_inplace(int count, double *restrict x) {
  // x = __builtin_assume_aligned(x, 32);
  // res = __builtin_assume_aligned(res, 32);
  while (count-- > 0) {
    *x = sqrt(*x);
    x++;
  }
}

/**
 * Multiply the elements of vector x and y and store in res of type double
 * complex.
 *    res = x * y
 */
static inline void vm_double_complex_multiply(int count, const void *restrict x,
                                              const void *restrict y,
                                              void *restrict res) {
  // x = __builtin_assume_aligned(x, 32);
  // y = __builtin_assume_aligned(y, 32);
  // res = __builtin_assume_aligned(res, 32);
  double *res_ = (double *)res;
  double *x_ = (double *)x;
  double *y_ = (double *)y;
  double real, imag, a, b, c, d;

  while (count-- > 0) {
    real = *x_++;
    imag = *x_++;
    a = real * *y_;    // real real
    c = imag * *y_++;  // imag real
    b = imag * *y_;    // imag imag
    d = real * *y_++;  // real imag
    *res_++ = a - b;
    *res_++ = c + d;

    // *res++ = *x++ * *y++;
  }
}

// Trignometry

/**
 * Cosine of the elements of vector x stored in res of type double.
 *    res = cos(x)
 */
static inline double get_cos_from_table(double x) {
  x = absd(x);
  x = modd(x, CONST_2PI);
  double res = x * table_precision_inverse;
  int index = (int)res;
  return lerp(res - index, cos_table[index], cos_table[index + 1]);
}

static inline void vm_double_cosine(int count, const double *restrict x,
                                    double *restrict res) {
  // x = __builtin_assume_aligned(x, 32);
  // res = __builtin_assume_aligned(res, 32);
  while (count-- > 0) {
    *res++ = get_cos_from_table(*x++);
    // *res++ = cos(*x++);
  }
}

/**
 * Sine of the elements of vector x stored in res of type double.
 *      res = sin(x)
 */
static inline double get_sin_from_table(double x) {
  x = absd(x);
  x = modd(x, CONST_2PI);
  double res = x * table_precision_inverse;
  int index = (int)res;
  res = lerp(res - index, sin_table[index], sin_table[index + 1]);
  res *= sign(x);
  return res;
}

static inline void vm_double_sine(int count, const double *restrict x,
                                  double *restrict res) {
  // x = __builtin_assume_aligned(x, 32);
  // res = __builtin_assume_aligned(res, 32);
  while (count-- > 0) {
    *res++ = get_sin_from_table(*x++);
    // *res++ = sin(*x++);
  }
}

/**
 * Cosine + I Sine of the elements of vector x in rad and stored in res of type
 * complex128. res = cos(x) + I sin(x)
 */
static inline void get_cos_sin_from_table(double x, double *restrict res_) {
  x = absd(x);
  x = modd(x, CONST_2PI);
  x *= table_precision_inverse;
  int i = (int)x;
  double wt = x - i;
  *res_++ = lerp(wt, cos_table[i], cos_table[i + 1]);
  *res_ = lerp(wt, sin_table[i], sin_table[i + 1]);
  *res_ *= sign(x);
}

static inline void vm_cosine_I_sine(int count, const double *restrict x,
                                    void *restrict res) {
  // x = __builtin_assume_aligned(x, 32);
  // res = __builtin_assume_aligned(res, 32);
  double *res_ = (double *)res;
  while (count-- > 0) {
    get_cos_sin_from_table(*x++, res_);
    res_ += 2;
    // *res_++ = cos(*x);
    // *res_++ = sin(*x++);
  }
}

// Exponent

static inline double my_exp(double x) {
  x /= 1024.0;
  x += 1.0;
  x *= x;  // 1
  x *= x;  // 2
  x *= x;  // 3
  x *= x;  // 4
  x *= x;  // 5
  x *= x;  // 6
  x *= x;  // 7
  x *= x;  // 8
  x *= x;  // 9
  x *= x;  // 10
  return x;
}

/**
 * Exponent of the elements of vector x stored in res of type double.
 *      res = exp(x)
 */
static inline void vm_double_exp(int count, double *restrict x, double *restrict res,
                                 const int ix) {
  // x = __builtin_assume_aligned(x, 32);
  // res = __builtin_assume_aligned(res, 32);
  while (count-- > 0) {
    *res = my_exp(*x);
    x += ix;
    res += ix;
  }
}

/**
 * Exponent of the elements of vector x stored in res of type complex128.
 *      res = exp(x)
 */
static inline void vm_double_complex_exp(int count, const void *restrict x,
                                         void *restrict res) {
  // x = __builtin_assume_aligned(x, 32);
  // res = __builtin_assume_aligned(res, 32);
  double *x_ = (double *)x;
  double *res_ = (double *)res;
  double temp;

  while (count-- > 0) {
    // *res++ = cexp(*x++);
    temp = my_exp(*x_++);
    *res_++ = cos(*x_) * temp;
    *res_++ = sin(*x_++) * temp;
  }
}

/**
 * Exponent of the elements of vector x stored in res of type complex128.
 *      res = exp(x(imag))
 */
static inline void vm_double_complex_exp_imag_only(int count, const void *restrict x,
                                                   void *restrict res) {
  // x = __builtin_assume_aligned(x, 32);
  // res = __builtin_assume_aligned(res, 32);
  double *x_ = (double *)x;
  double *res_ = (double *)res;
  double y, wt;
  int i;

  while (count-- > 0) {
    x_++;
    y = absd(*x_);
    y = modd(y, CONST_2PI);
    y *= table_precision_inverse;
    i = (int)y;
    wt = y - i;
    *res_++ = lerp(wt, cos_table[i], cos_table[i + 1]);
    *res_ = lerp(wt, sin_table[i], sin_table[i + 1]);
    *res_ *= sign(*x_);

    res_++;
    x_++;

    // x_++;
    // *res_++ = cos(*x_);
    // *res_++ = sin(*x_++);
  }
}

#ifndef __blas_activate
//========================================================================== //
//                  Wrapper for blas and blas like functions                 //
//========================================================================== //

/**
 * Scale the elements of vector x by a of type double.
 * x *= a
 * Equivalent to cblas_dscal.
 */
static inline void cblas_dscal(int count, const double a, double *restrict x,
                               const int stride_x) {
  // x = __builtin_assume_aligned(x, 32);
  while (count-- > 0) {
    *x *= a;
    x += stride_x;
  }
}

/**
 * Scale the elements of complex128 vector x by a double a.
 * x *= a
 * Equivalent to cblas_zdscal.
 */
static inline void cblas_zdscal(int count, const double a, void *restrict x,
                                const int stride_x) {
  // x = __builtin_assume_aligned(x, 32);
  double *x_ = (double *)x;
  int stride_xp = stride_x * 2 - 1;
  while (count-- > 0) {
    *x_++ *= a;
    *x_++ *= a;
    x_ += stride_xp;
  }
}

/**
 * Copy elements of vector x to vector y of type double.
 * y = x
 * Equivalent to cblas_dcopy.
 */
static inline void cblas_dcopy(int count, const double *restrict x, const int stride_x,
                               double *restrict y, const int stride_y) {
  // x = __builtin_assume_aligned(x, 32);
  // y = __builtin_assume_aligned(y, 32);
  while (count-- > 0) {
    *y = *x;
    x += stride_x;
    y += stride_y;
  }
}

/**
 * Copy elements of vector x to vector y of type complex128.
 * y = x
 * Equivalent to cblas_zcopy.
 */
static inline void cblas_zcopy(int count, const void *restrict x, const int stride_x,
                               void *restrict y, const int stride_y) {
  // x = __builtin_assume_aligned(x, 32);
  // y = __builtin_assume_aligned(y, 32);
  double *x_ = (double *)x;
  double *y_ = (double *)y;
  int stride_xp = stride_x * 2 - 1;
  int stride_yp = stride_y * 2 - 1;
  while (count-- > 0) {
    *y_++ = *x_++;
    *y_++ = *x_++;
    x_ += stride_xp;
    y_ += stride_yp;
  }
}

/**
 * Compute the following
 * y = a*x + y
 * Equivalent to cblas_daxpy.
 */
static inline void cblas_daxpy(int count, const double a, const double *restrict x,
                               const int stride_x, double *restrict y,
                               const int stride_y) {
  // x = __builtin_assume_aligned(x, 32);
  // y = __builtin_assume_aligned(y, 32);
  while (count-- > 0) {
    *y += a * *x;
    x += stride_x;
    y += stride_y;
  }
}

/**
 * Compute the following
 * y = a*x + b*y
 * Equivalent to catlas_daxpby.
 */
static inline void catlas_daxpby(int count, const double a, const double *restrict x,
                                 const int stride_x, const double b, double *restrict y,
                                 const int stride_y) {
  // x = __builtin_assume_aligned(x, 32);
  // y = __builtin_assume_aligned(y, 32);
  while (count-- > 0) {
    *y *= b;
    *y += a * *x;
    x += stride_x;
    y += stride_y;
  }
}
#endif /* __blas_activate */

#endif /* DOXYGEN_SHOULD_SKIP_THIS */
