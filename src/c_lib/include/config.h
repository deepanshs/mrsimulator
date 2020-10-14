// -*- coding: utf-8 -*-
//
//  config.h
//
//  @copyright Deepansh J. Srivastava, 2019-2020.
//  Created by Deepansh J. Srivastava, Aug 10, 2019
//  Contact email = srivastava.89@osu.edu
//

#ifndef __config__

#define __config__

typedef double complex128[2];
typedef float complex64[2];

#if __STDC_VERSION__ >= 199901L
#define MKL_Complex16 complex128
#else  // not C99
#define restrict __restrict
#endif

// for windows msvc compiler
#if _MSC_VER && !__INTEL_COMPILER
#include "cblas.h"
#else

// #ifdef __APPLE__

#if __has_include(<Accelerate/Accelerate.h>)
#include <Accelerate/Accelerate.h>

#elif __has_include("cblas.h")
#include "cblas.h"

// installing openblas on cent-os with
//      yum install -y openblas-devel
// include cblas-openblas.h header file
#elif __has_include("cblas-openblas.h")
#include "cblas-openblas.h"

// installing openblas on macos with
//      port install openblas
// include cblas_openblas.h header file
#elif __has_include("cblas_openblas.h")
#include "cblas_openblas.h"

#endif
#endif

#define __blas_activate
#include "vm.h"

// user definition
#define PI2 6.2831853072
#define PI2I PI2 *I

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

#include <stdbool.h>
#include <stdio.h>
#include <time.h>

#include "array.h"
#include "vm_common.h"

#endif
