// -*- coding: utf-8 -*-
//
//  config.h
//
//  @copyright Deepansh J. Srivastava, 2019-2021.
//  Created by Deepansh J. Srivastava, Aug 10, 2019.
//  Contact email = srivastava.89@osu.edu
//
#include <math.h>
#include <stdbool.h>
#include <stdio.h>
#include <time.h>

#ifndef __config__

#define __config__

// user definitions ----------------------------------------------------------- //
typedef double complex128[2];
typedef float complex64[2];

#include "tables/trig.h"

#define CONST_PI 3.14159265358979323846264338327950288419716939937510
#define CONST_2PI 6.28318530717958623199592693708837032318115234375000
#define CONST_4PI 2 * CONST_2PI
#define CONST_iPI CONST_2PI *I
#define TOL 1.0e-6

#define modd(x, y) ((x) - (int)((x) / (y)) * (y))  // fold x within range y
#define lerp(w, v1, v2) ((1.0 - (w)) * (v1) + (w) * (v2))
#define sign(x) (int)(((x) > 0) - ((x) < 0))  // return sign of x

// Compiler version check
#if __STDC_VERSION__ >= 199901L
#define MKL_Complex16 complex128
#else  // not C99
#define restrict __restrict
#endif
// ---------------------------------------------------------------------------- //

// Blas definitions ----------------------------------------------------------- //
#ifdef USE_OPENBLAS  // openblas header
// when openblas is installed on cent-os with yum
//      yum install -y openblas-devel
#if __has_include("cblas-openblas.h")  // openblas header
#include "cblas-openblas.h"
#else
#include "cblas.h"
#endif
#endif

#ifdef USE_ACCELERATE  // accelerate header
#include <Accelerate/Accelerate.h>
#endif  // end blas headers

#ifdef USE_MKL  // mkl header
#include "mkl.h"
#endif  // mkl header
// ---------------------------------------------------------------------------- //

// OS base definitions -------------------------------------------------------- //
#define __int64_ long  // for both mac-os and linux (unix system)

#ifdef _WIN32  // windows 32-bit or 64-bit
#define __int64_ long long
// #if _MSC_VER && !__INTEL_COMPILER  // windows msvc compiler
// #endif  // end windows 32-bit or 64-bit
#endif  // end windows msvc compiler
// ---------------------------------------------------------------------------- //

#define __blas_activate

#include "array.h"
#include "vm.h"
#include "vm_common.h"

#endif /* __config__ */
