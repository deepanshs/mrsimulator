//
//  mrsimulator.h
//
//  Created by Deepansh J. Srivastava, Jun 30, 2019
//  Copyright © 2019 Deepansh J. Srivastava. All rights reserved.
//  Contact email = srivastava.89@osu.edu, deepansh2012@gmail.com
//

#ifndef mrsimulator_h
#define mrsimulator_h

// library definition
#include <complex.h>

#define MKL_Complex16 double complex

#include "fftw/fftw3.h"
#include "mkl.h"
#include "omp.h"
#include <math.h>
#include <stdbool.h>
#include <stdio.h> // to use calloc, malloc, and free methods
#include <stdlib.h>
#include <string.h>
#include <time.h>

// user definition
#define PI2 6.2831853072
#define PI2I PI2 *I

#include "Hamiltonian.h"
#include "angular_momentum.h"
#include "array.h"
#include "interpolation.h"
#include "isotopomer_ravel.h"
#include "octahedron.h"
#include "powder_setup.h"

/* A plan for mrsimulator containing buffers and tabulated values for faster *
 * simulation. The plan includes pre-calculating 1) wigner-2j and wigner-4j  *
 * matrices at every orientation (β), 2) pre-calculating exp(-Imα), and the  *
 * exponent of the sideband phase.                                           */
typedef struct mrsimulator_plan_t {
  fftw_plan the_fftw_plan;
  fftw_complex *vector;
  unsigned short geodesic_polyhedron_frequency; // number of triangles along
                                                // the edge of octahedron.
  int number_of_sidebands;                      // number of sidebands.
  unsigned int n_orientations; // number of unique orientations.
  unsigned int size;           // number of orientations * number of sizebands.
  double *vr_freq;    // sideband order frequencies in fft output order.
  double *amplitudes; // array of amplitude per orientation.
  double complex *exp_Im_alpha; // array of cos_alpha per orientation.
  double complex *w2;           // buffer for 2nd rank frequency calculation.
  double *wigner_2j_matrices;   // wigner-d 2j matrix per orientation.
  double *rotor_lab_2;          // wigner-2j dm0 vector, n ∈ [-2, 2].
  double complex *pre_phase_2; // buffer for 2nk rank sideband phase calculation
  double complex *w4;          // buffer for 4nd rank frequency calculation.
  double *wigner_4j_matrices;  // wigner-d 4j matrix per orientation.
  double *rotor_lab_4;         // wigner-4j dm0 vector, n ∈ [-4, 4].
  double complex *pre_phase_4; // buffer for 4th rank sideband phase calculation
  double *local_frequency;     // buffer for local frequencies.
  double *freq_offset;         // buffer for local + sideband frequencies.
} mrsimulator_plan;

typedef struct dimension_t {
  int count;                // number of coordinates along the dimension.
  double coordinate_offset; // starting coordinate of the dimension.
  double increment;         // increment of coordinates along the dimension.

  /* private attributes */
  double normalize_offset; // fixed value = 0.5 - coordinate_offset/increment
} dimension;

mrsimulator_plan *
mrsimulator_create_plan(unsigned int geodesic_polyhedron_frequency,
                        int number_of_sidebands,
                        double sample_rotation_frequency, double rotor_angle,
                        double spectral_increment, bool ALLOW_FOURTH_RANK);

void mrsimulator_free_plan(mrsimulator_plan *the_plan);

mrsimulator_plan *
mrsimulator_update_plan(unsigned int *geodesic_polyhedron_frequency,
                        int *number_of_sidebands,
                        double *sample_rotation_frequency, double *rotor_angle,
                        double *spectral_increment, bool ALLOW_FOURTH_RANK);

/* Return a vector ordered according to the fft output order.                *
 * @params int n - The number of points                                      *
 * @params double increment - The increment (sampling interval)              *
 * @returns *double values = The pointer to the fft output order vector      */
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

extern void __get_pre_phase_components(int number_of_sidebands,
                                       double spin_frequency,
                                       double complex *pre_phase);

#endif /* mrsimulator_h */
