//
//  vm.h
//
//  Created by Deepansh J. Srivastava, Jul 26, 2019
//  Copyright © 2019 Deepansh J. Srivastava. All rights reserved.
//  Contact email = srivastava.89@osu.edu, deepansh2012@gmail.com
//

#ifndef DOXYGEN_SHOULD_SKIP_THIS
/** Arithmetic suit ======================================================== */

/**
 * Add the elements of vector x and y and store in res of type double.
 * res = x + y
 */
static inline void vm_double_add(int count, const double *restrict x,
                                 const double *restrict y,
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
 * Add the elements of vector y inplace with the elements from vector x.
 * y += x
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
 * res = x - y
 */
static inline void vm_double_subtract(int count, const double *restrict x,
                                      const double *restrict y,
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
 * Subtract the elements of vector y inplace with the elements from vector x.
 * y -= x
 */
static inline void vm_double_subtract_inplace(int count,
                                              const double *restrict x,
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
 * res = x * y
 */
static inline void vm_double_multiply(int count, const double *restrict x,
                                      const double *restrict y,
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
 * Multiply the elements of vector y inplace with the elements from vector x.
 * y *= x
 */
static inline void vm_double_multiply_inplace(int count,
                                              const double *restrict x,
                                              double *restrict y) {
  // x = __builtin_assume_aligned(x, 32);
  // y = __builtin_assume_aligned(y, 32);
  while (count-- > 0) {
    *y++ *= *x++;
    // x += stride_x;
    // y += stride_y;
  }
}

/**
 * Divide the elements of vector x by y and store in res of type double.
 * res = x / y
 */
static inline void vm_double_divide(int count, const double *restrict x,
                                    const double *restrict y,
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
 * Divide the elements of vector y inplace with the elements from vector x.
 * y /= x
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
 * res = x * x
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
 * x *= x
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
 * res = sqrt(x)
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
 * x = sqrt(x)
 */
static inline void vm_double_square_root_inplace(int count,
                                                 double *restrict x) {
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
 * res = x * y
 */
static inline void vm_double_complex_multiply(int count,
                                              const complex128 *restrict x,
                                              const complex128 *restrict y,
                                              complex128 *restrict res) {
  // x = __builtin_assume_aligned(x, 32);
  // y = __builtin_assume_aligned(y, 32);
  // res = __builtin_assume_aligned(res, 32);
  // double *res_ = (double *)res;
  // double *x_ = (double *)x;
  // double *y_ = (double *)y;
  // int count_ = 2 * count;

  while (count-- > 0) {
    *res++ = *x++ * *y++;
  }
}

// Trignometry

/**
 * Cosine of the elements of vector x stored in res of type double.
 * res = cos(x)
 */
static inline void vm_double_cosine(int count, const double *restrict x,
                                    double *restrict res) {
  // x = __builtin_assume_aligned(x, 32);
  // res = __builtin_assume_aligned(res, 32);
  while (count-- > 0) {
    *res = cos(*x);
    res++;
    x++;
    // x += stride_x;
    // res += stride_res;
  }
}

/**
 * Sine of the elements of vector x stored in res of type double.
 * res = sin(x)
 */
static inline void vm_double_sine(int count, const double *restrict x,
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
static inline void vm_cosine_I_sine(int count, const double *restrict x,
                                    complex128 *restrict res) {
  // x = __builtin_assume_aligned(x, 32);
  // res = __builtin_assume_aligned(res, 32);
  double *res_ = (double *)res;
  while (count-- > 0) {
    *res_++ = cos(*x);
    *res_++ = sin(*x++);
  }
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
static inline void vm_double_exp(int count, double *restrict x,
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
static inline void vm_double_complex_exp(int count,
                                         const complex128 *restrict x,
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
static inline void cblas_zdscal(int count, const double a,
                                complex128 *restrict x, const int stride_x) {
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
static inline void cblas_dcopy(int count, const double *restrict x,
                               const int stride_x, double *restrict y,
                               const int stride_y) {
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
static inline void cblas_zcopy(int count, const complex128 *restrict x,
                               const int stride_x, complex128 *restrict y,
                               const int stride_y) {
  // x = __builtin_assume_aligned(x, 32);
  // y = __builtin_assume_aligned(y, 32);
  while (count-- > 0) {
    *y = *x;
    x += stride_x;
    y += stride_y;
  }
}

/**
 * Compute the following
 * y = a*x + y
 * Equivalent to cblas_daxpy.
 */
static inline void cblas_daxpy(int count, const double a,
                               const double *restrict x, const int stride_x,
                               double *restrict y, const int stride_y) {
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
static inline void catlas_daxpby(int count, const double a,
                                 const double *restrict x, const int stride_x,
                                 const double b, double *restrict y,
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
