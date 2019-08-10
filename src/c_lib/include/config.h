

//
//  config.h
//
//  Created by Deepansh J. Srivastava, Aug 10, 2019
//  Copyright © 2019 Deepansh J. Srivastava. All rights reserved.
//  Contact email = srivastava.89@osu.edu, deepansh2012@gmail.com
//

#ifndef __config__

#define __config__

#include "vm_common.h"
#include <complex.h>

#if __STDC_VERSION__ >= 199901L
#define complex128 double _Complex
#define MKL_Complex16 double _Complex
#define complex64 _Complex
#define MKL_Complex6 _Complex
#define complex128_add_inplace(a, b) (a += b)
// similarly for other operations

#else // not C99
typedef struct complex128_ {
  double real;
  double imag;
} complex128;
typedef struct complex64_ {
  float real;
  float imag;
} complex64;
#define restrict __restrict
inline complex128 complex128_add_inplace(complex128 a, double b) {
  a.real += b;
  a.imag += b;
  return a;
}
#endif

// library definition
#if __has_include("mkl.h")
#include "mkl.h"
#define __blas_activate
#include "vm_mkl.h"

#elif __has_include("cblas.h")
#include "cblas.h"
#define __blas_activate
#include "vm.h"
#endif

// user definition
#define PI2 6.2831853072
#define PI2I PI2 *I

// #ifdef __APPLE__
// #include <Accelerate/Accelerate.h>
// #define __blas_activate
// #include "vm.h"
// #include "mkl.h"
// #include "vm_mkl.h"
// #endif

// #ifdef linux
// #include "mkl.h"
// #include "vm_mkl.h"
// // mkl_set_threading_layer(MKL_THREADING_INTEL);
// // int max_threads = mkl_get_max_threads();
// // mkl_set_num_threads(max_threads);
// // printf("Using upto %d threads for simulation.\n", max_threads);
// #endif

// #ifdef _WIN32
// #include "mkl.h"
// #include "vm_mkl.h"
// // mkl_set_threading_layer(MKL_THREADING_INTEL);
// // int max_threads = mkl_get_max_threads();
// // mkl_set_num_threads(max_threads);
// // printf("Using upto %d threads for simulation.\n", max_threads);
// #endif



#endif
