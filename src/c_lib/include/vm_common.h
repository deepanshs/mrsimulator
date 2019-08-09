//
//  vm_common.h
//
//  Created by Deepansh J. Srivastava, Jul 26, 2019
//  Copyright © 2019 Deepansh J. Srivastava. All rights reserved.
//  Contact email = srivastava.89@osu.edu, deepansh2012@gmail.com
//

#include <math.h>
#include <string.h>

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
static inline void vm_double_arange(int count, double *restrict res) {
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
  // while (count-- > 0) {
  //   *res++ = 0.0;
  // }
  memset(res, 0, count * sizeof(double));
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
