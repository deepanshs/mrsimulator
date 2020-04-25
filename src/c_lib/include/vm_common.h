// -*- coding: utf-8 -*-
//
//  vm_common.h
//
//  @copyright Deepansh J. Srivastava, 2019-2020.
//  Created by Deepansh J. Srivastava, Jul 26, 2019.
//  Contact email = deepansh2012@gmail.com
//

#include <math.h>
#include <string.h>

/**
 * List of function -
 *
 * vm_double_ramp -> creates a linear ramp, y = scale*x + offset.
 * vm_double_arrange -> Create a vector x = [0 .. count-1].
 * vm_double_zeros -> Create a vector of zeros.
 * vm_double_ones -> Create a vector of ones.
 *
 **/

/**
 * Multiply a vector of type double by `scale` and add an `offset` to its
 * elements. res = scale*x + offset
 */
static inline void vm_double_ramp(int count, const double *restrict x,
                                  const double scale, const double offset,
                                  double *restrict res) {
  // x = __builtin_assume_aligned(x, 32);
  // res = __builtin_assume_aligned(res, 32);
  while (count-- > 0) {
    *res++ = scale * *x++ + offset;
  }
}

/**
 * Create a vector x = [0 .. count-1]
 * res = 0 .. count-1
 */
static inline void vm_double_arrange(int count, double *restrict res) {
  //   x = __builtin_assume_aligned(x, 32);
  //   res = __builtin_assume_aligned(res, 32);
  double i = 0.0;
  while (count-- > 0) {
    *res++ = i++;
  }
}

/**
 * Create a vector of length count with all zero entries
 * res = [0.0, 0.0, 0.0, ... ]
 */
static inline void vm_double_zeros(int count, double *restrict res) {
  //   x = __builtin_assume_aligned(x, 32);
  //   res = __builtin_assume_aligned(res, 32);
  while (count-- > 0) {
    *res++ = 0.0;
  }
  // memset(res, 0, count * sizeof(double));
}

/**
 * Create a vector of length count with all one entries
 * res = [1.0, 1.0, 1.0, ... ]
 */
static inline void vm_double_ones(int count, double *restrict res) {
  //   x = __builtin_assume_aligned(x, 32);
  //   res = __builtin_assume_aligned(res, 32);
  while (count-- > 0) {
    *res++ = 1.0;
  }
  // memset(res, 0, count * sizeof(double));
}

/**
 * @brief Return a vector ordered according to the fft output order.
 *
 * @param n The number of points.
 * @param increment The increment along the dimension axis (sampling interval).
 * @returns values A pointer to the fft output order vector of size @p n.
 */
static inline double *__get_frequency_in_FFT_order(int n, double increment) {
  double *vr_freq = malloc_double(n);
  int i = 0, m, positive_limit, negative_limit;

  if (n % 2 == 0) {
    negative_limit = (int)(-n / 2);
    positive_limit = -negative_limit - 1;
  } else {
    negative_limit = (int)(-(n - 1) / 2);
    positive_limit = -negative_limit;
  }

  for (m = 0; m <= positive_limit; m++) {
    vr_freq[i] = (double)m * increment;
    i++;
  }
  for (m = negative_limit; m < 0; m++) {
    vr_freq[i] = (double)m * increment;
    i++;
  }
  return vr_freq;
};
