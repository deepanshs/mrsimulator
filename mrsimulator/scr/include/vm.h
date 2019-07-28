//
//  vm.h
//
//  Created by Deepansh J. Srivastava, Jul 26, 2019
//  Copyright © 2019 Deepansh J. Srivastava. All rights reserved.
//  Contact email = srivastava.89@osu.edu, deepansh2012@gmail.com
//

// #include <math.h>

/**
 * Multiply a vector of type double by `scale` and add an `offset` to its
 * elements. res = scale*x + offset
 */
static inline void vm_dlinear(int count, double *restrict x, double scale,
                              double offset, double *restrict res) {
  // x = __builtin_assume_aligned(x, 32);
  // res = __builtin_assume_aligned(res, 32);
  while (count-- > 0) {
    *res++ = scale * *x++ + offset;
  }
}

/**
 * Add the elements of vector x and y and store in res of type double.
 * res = x + y
 */
static inline void vmd_add(int count, double *restrict x, double *restrict y,
                           double *restrict res) {
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
 * Subtract the elements of vector x from y and store in res of type double.
 * res = x - y
 */
static inline void vmd_sub(int count, double *restrict x, double *restrict y,
                           double *restrict res) {
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
 * Multiply the elements of vector x and y and store in res of type double.
 * res = x * y
 */
static inline void vmd_mul(int count, double *restrict x, double *restrict y,
                           double *restrict res) {
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
 * Divide the elements of vector x by y and store in res of type double.
 * res = x / y
 */
static inline void vmd_div(int count, double *restrict x, double *restrict y,
                           double *restrict res) {
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
 * Square the elements of vector x and store in res of type double.
 * res = x * x
 */
static inline void vmd_sqr(int count, double *restrict x,
                           double *restrict res) {
  // x = __builtin_assume_aligned(x, 32);
  // res = __builtin_assume_aligned(res, 32);
  while (count-- > 0) {
    *res++ = *x * *x;
    x++;
  }
}

/**
 * Square root of the elements of vector x stored in res of type double.
 * res = sqrt(x)
 */
static inline void vmd_sqrt(int count, double *restrict x,
                            double *restrict res) {
  // x = __builtin_assume_aligned(x, 32);
  // res = __builtin_assume_aligned(res, 32);
  while (count-- > 0) {
    *res++ = sqrt(*x++);
  }
}

/**
 * Multiply the elements of vector x and y and store in res of type double
 * complex.
 * res = x * y
 */
static inline void vmz_mul(int count, complex128 *restrict x,
                           complex128 *restrict y, complex128 *restrict res) {
  // x = __builtin_assume_aligned(x, 32);
  // y = __builtin_assume_aligned(y, 32);
  // res = __builtin_assume_aligned(res, 32);
  // double *res_ = (double *)res;
  // double *x_ = (double *)x;
  // double *y_ = (double *)y;
  // int count_ = 2 * count;

  while (count-- > 0) {
    *res++ = *x++ * *y++;
    // *(res+1) -=
  }
}

// Trignometry

/**
 * Cosine of the elements of vector x stored in res of type double.
 * res = cos(x)
 */
static inline void vmd_cos(int count, double *restrict x,
                           double *restrict res) {
  // x = __builtin_assume_aligned(x, 32);
  // res = __builtin_assume_aligned(res, 32);
  while (count-- > 0) {
    *res++ = cos(*x++);
    // x += stride_x;
    // res += stride_res;
  }
}

/**
 * Sine of the elements of vector x stored in res of type double.
 * res = sin(x)
 */
static inline void vmd_sin(int count, double *restrict x,
                           double *restrict res) {
  // x = __builtin_assume_aligned(x, 32);
  // res = __builtin_assume_aligned(res, 32);
  while (count-- > 0) {
    *res++ = sin(*x++);
    // x += stride_x;
    // res += stride_res;
  }
}

/**
 * Cosine + I Sine of the elements of vector x in rad and stored in
 * res of type complex128.
 * res = cos(x) + I sin(x)
 */
static inline void vmd_cospisin(int count, double *restrict x,
                                complex128 *restrict res) {
  // x = __builtin_assume_aligned(x, 32);
  // res = __builtin_assume_aligned(res, 32);
  double *res_ = (double *)res;
  *res_++ = cos(*x);
  *res_++ = sin(*x++);
}

// Exponent

static inline double my_exp(double x) {
  x /= 1024.0;
  x += 1.0;
  x *= x; // 1
  x *= x; // 2
  x *= x; // 3
  x *= x; // 4
  x *= x; // 5
  x *= x; // 6
  x *= x; // 7
  x *= x; // 8
  x *= x; // 9
  x *= x; // 10
  return x;
}
/**
 * Exponent of the elements of vector x stored in res of type double.
 * res = exp(x)
 */
static inline void vm_dexp(int count, double *restrict x,
                           double *restrict res) {
  // x = __builtin_assume_aligned(x, 32);
  // res = __builtin_assume_aligned(res, 32);
  while (count-- > 0) {
    *res++ = my_exp(*x++);
  }
}

/**
 * Exponent of the elements of vector x stored in res of type complex128.
 * res = exp(x)
 */
static inline void vm_zexp(int count, complex128 *restrict x,
                           complex128 *restrict res) {
  // x = __builtin_assume_aligned(x, 32);
  // res = __builtin_assume_aligned(res, 32);
  // double *x_ = (double *)x;
  // double *res_ = (double *)res;
  // double *res_1 = (double *)res + 1;
  // int count_ = 2 * count;
  // int i = 1;
  // double factor, num_;

  // // double temp;
  while (count-- > 0) {
    //   i = 1;
    //   factor = *++x_;
    //   num_ = factor;
    //   *res_ = 1;
    //   *res_1 = factor;
    //   while (i < 10)
    //   {
    //     factor *= num_ / ++i;
    //     *res_ -= factor;
    //     factor *= num_ / ++i;
    //     *res_1 -= factor;
    //     factor *= num_ / ++i;
    //     *res_ += factor;
    //     factor *= num_ / ++i;
    //     *res_1 += factor;
    //   }
    //   res_ += 2;
    //   res_1 += 2;
    //   x_++;

    *res++ = cexp(*x++);
    // temp = my_exp(*x_++);
    // *res_++ = cos(*++x_) * temp;
    // *res_++ = sin(*x_++) * temp;
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
static inline void cblas_dscal(int count, double a, double *restrict x,
                               int stride_x) {
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
static inline void cblas_zdscal(int count, double a, complex128 *restrict x,
                                int stride_x) {
  // x = __builtin_assume_aligned(x, 32);
  while (count-- > 0) {
    *x *= a;
    x += stride_x;
  }
}

/**
 * Copy elements of vector x to vector y of type double.
 * y = x
 * Equivalent to cblas_dcopy.
 */
static inline void cblas_dcopy(int count, double *restrict x, int stride_x,
                               double *restrict y, int stride_y) {
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
static inline void cblas_zcopy(int count, complex128 *restrict x, int stride_x,
                               complex128 *restrict y, int stride_y) {
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
 * y = a*x + b*y
 * Equivalent to catlas_daxpby.
 */
static inline void catlas_daxpby(int count, double a, double *restrict x,
                                 int stride_x, double b, double *restrict y,
                                 int stride_y) {
  // x = __builtin_assume_aligned(x, 32);
  // y = __builtin_assume_aligned(y, 32);
  while (count-- > 0) {
    *y *= b;
    *y += a * *x;
    x += stride_x;
    y += stride_y;
  }
}

// static inline void vm_zgemm(const CBLAS_TRANSPOSE transa, const
// CBLAS_TRANSPOSE transb, int m,
//                             int n, int k, complex128 *alpha,
//                             complex128 *restrict a, int lda, complex128
//                             *restrict b, int ldb, complex128 *beta,
//                             complex128 *restrict c, int ldc)
// {
//   int i, j;
// }
#endif
