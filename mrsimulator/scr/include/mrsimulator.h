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
// #include <complex.h>

// #if __STDC_VERSION__ >= 199901L
// // using a C99 compiler
// #include <complex.h>
// typedef double _Complex complex128;
// #else
// typedef struct {
//   double real, imag;
// } complex128;
// #endif

// #if __STDC_VERSION__ >= 199901L
// // creal, cimag already defined in complex.h
// inline complex128 make_complex_double(double real, double imag) {
//   return real + imag * I;
// }
// #else
// #define creal(z) ((z).real)
// #define cimag(z) ((z).imag)

// // extern const complex128 complex_i; // put in a translation unit somewhere
// // #define I complex_i

// inline complex128 make_complex_double(double real, double imag) {
//   complex128 z = {real, imag};
//   return z;
// }
// #endif

// #if __STDC_VERSION__ >= 199901L
// #define cadd(a, b) ((a) + (b))
// #define csub(a, b) ((a) - (b))
// #define cmult(a, b) ((a) * (b))
// #define cdiv(a, b) ((a) / (b))
// // similarly for other operations
// #else // not C99
// inline complex128 conj(complex128 a) {
//   complex128 z = {a.real, -a.imag};
//   return z;
// }
// inline complex128 cadd(complex128 a, complex128 b) {
//   complex128 z = {a.real + b.real, a.imag + b.imag};
//   return z;
// }
// inline complex128 csub(complex128 a, complex128 b) {
//   complex128 z = {a.real - b.real, a.imag - b.imag};
//   return z;
// }
// inline complex128 cmult(complex128 a, complex128 b) {
//   complex128 z;
//   z.real = a.real * b.real - a.imag * b.imag;
//   z.imag = a.real * b.imag + a.imag * b.real;
//   return z;
// }
// inline double cabs(complex a) {
//   double res = sqrt(a.real * a.real + b.imag * b.imag)
// }
// inline complex128 cdiv(complex128 a, complex128 b) {
//   complex128 z = cmult(a, conj(b));
//   double res = b.real * b.real + b.imag * b.imag;
//   z.real /= res;
//   z.imag /= res;
//   return z;
// }
// #endif

// #include "mkl.h"
// #define complex128 MKL_Complex16
// #include "my_complex.h"

#ifdef _WIN32
#include "mkl.h"
#define complex128 MKL_Complex16
#include "my_complex.h"
// #if !(__STDC_VERSION__ >= 199901L)
// #define restrict __restrict
// #endif
#include "vm_mkl.h"
#endif

#ifdef unix
// #include "mkl.h"
// #define complex128 MKL_Complex16
#include "vm_mkl.h"
#endif

#ifdef __APPLE__
#include "mkl.h"
#define complex128 MKL_Complex16
#include "my_complex.h"
#include "vm_mkl.h"
// int max_threads = mkl_get_max_threads();
// mkl_set_num_threads(max_threads);
// #include "vm_accelerate.h"
#endif

#include "fftw/fftw3.h"

// #include "mkl.h"
// #include "omp.h"
// #include "vm.h"

#include <stdbool.h>
#include <stdio.h>
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

/**
 * @struct MRS_plan_t
 * @brief Create a mrsimulator plan for faster lineshape simulation.
 */
struct MRS_plan_t {
  /**
   * The number of triangles along the edge of polyhedra. This value is a
   * positive integer which represents the frequency of class I geodesic
   * polyhedra. These polyhedra are used in calculating the spherical
   * average. Currently, we are using octahedral as the frequency 1 polyhedra.
   *
   * Generally, as the frequency of the geodesic polyhedron increases, the
   * polyhedra geometry approach to that of a sphere. For line-shape simulation,
   * a higher geodesic polyhedron frequency results in a better spherical
   * averaging. The default value is 72. Read more on <a
   * href="https://en.wikipedia.org/wiki/Geodesic_polyhedron">Geodesic
   * polyhedron</a>.
   */
  unsigned short geodesic_polyhedron_frequency;

  int number_of_sidebands; /**< The number of sidebands to compute. */
  /** The sample rotation frequency in Hz. */
  double sample_rotation_frequency_in_Hz;
  /**
   * The polar angle, in radians, describing the axis of rotation of the sample
   * with respect to the lab-frame z-axis.
   */
  double rotor_angle_in_rad;
  bool allow_fourth_rank; /**< Allow buffer for fourth rank tensors. */
  /**
   * The pointer to sideband frequency ratio array stored in the fft output
   * order. The sideband frequency ratio is defined as
   *    @f[\frac{\mathbf{J} \omega_r}{n_i}@f]
   * where
   *    @f$\mathbf{J}@f = [0, 1, 2, ... number_of_sideband/2,
   * -number_of_sideband/2+1, -2, -1]$ is an integer, @f$\omega_r@f$ is the
   * spinning frequency frequency in Hz, and @f$n_i@f$ is the `increment` along
   * the spectroscopic grid dimension.
   */
  double *vr_freq;
  /** Isotropic frequency offset ratio. The ratio is similarly defined as
   *  before.
   */
  double isotropic_offset;
  /** The buffer to hold the sideband amplitudes as stride 2 array after
   *  mrsimulator processing.
   */
  complex128 *vector;

  /* private attributes */
  fftw_plan the_fftw_plan;     /**< The plan for fftw routine. */
  unsigned int n_orientations; /**<  Number of unique orientations. */
  unsigned int size; /**<  Number of orientations * number of sidebands. */

  /**
   *  Array of amplitudes of size `n_orientations` representing the spherical
   * average.
   */
  double *amplitudes;
  /**
   * Array of @f$\exp(Im\alpha)@r$, where @f$m \in [-4,0] @f$ and @f$\alpha@f$
   * is the azimuthal angle per orientation. The shape of this array is
   * 5 x n_orientations.
   */
  complex128 *exp_Im_alpha;
  /**
   * Buffer for processing the frequencies from the second rank tensor. The
   * shape of this array is 5 * n_orientations.
   */
  complex128 *w2;
  /**
   * Wigner @f$d^2{m,n}(\beta)@f$ matrices at every orientation angle
   * @f$\beta@f$. Here, `m` and `n` ranges from -2 to 2. The shape of this
   * array is 25 x n_orientations.
   */
  double *wigner_2j_matrices;
  /**
   * Wigner @f$d^2_{m,0}@f$ vector where @f$m \in [-2, 2]@f$. This vector is
   * used to rotate the second rank tensor from the rotor-frame to the
   * lab-frame.
   */
  double *rotor_lab_2;

  complex128
      *pre_phase_2; /**<  buffer for 2nk rank sideband phase calculation */

  /**
   * Buffer for processing the frequencies from the fourth rank tensors. The
   * shape of this array is 5 * n_orientations.
   */
  complex128 *w4;

  /**
   * Wigner @f$d^4{m,n}(\beta)@f$ matrices at every orientation angle
   * @f$\beta@f$. Here, `m` and `n` ranges from -4 to 4. The shape of this
   * array is 81 x n_orientations.
   */
  double *wigner_4j_matrices;

  /**
   * Wigner @f$d^4_{m,0}@f$ vector where @f$m \in [-2, 2]@f$. This vector is
   * used to rotate the second rank tensor from the rotor-frame to the
   * lab-frame.
   */
  double *rotor_lab_4; /**<  wigner-4j dm0 vector, n ∈ [-4, 4]. */
  complex128
      *pre_phase_4; /**<  buffer for 4th rank sideband phase calculation */
  double *local_frequency; /**<  buffer for local frequencies. */
  double *freq_offset;     /**<  buffer for local + sideband frequencies. */
  complex128 one;          /**<  holds complex value 1. */
  complex128 zero;         /**<  hold complex value 0. */
  double buffer;           /**<  buffer for temporary storage. */
};

typedef struct MRS_plan_t MRS_plan;

typedef struct MRS_dimension_t {
  int count;                 // the number of coordinates along the dimension.
  double coordinates_offset; // starting coordinate of the dimension.
  double increment;          // increment of coordinates along the dimension.

  /* private attributes */
  double normalize_offset; // fixed value = 0.5 - coordinate_offset/increment
  double inverse_increment;
} MRS_dimension;

MRS_dimension *MRS_create_dimension(int count, double coordinates_offset,
                                    double increment);

/**
 * @brief Release the allocated memory from a mrsimulator plan.
 *
 * @param	MRS_plan *plan The pointer to the mrsimulator plan to be freed.
 */
void MRS_free_plan(MRS_plan *plan);

/**
 * @brief Create a new mrsimulator plan.
 *
 * @param	geodesic_polyhedron_frequency The number of triangles along the
 *            edge of the octahedron.
 * @param number_of_sidebands The number of sideband to compute.
 * @param sample_rotation_frequency_in_Hz The sample rotation frequency in Hz.
 * @param rotor_angle_in_rad The polar angle in radians with respect to z-axis
 *            describing the axis of rotation.
 * @param increment The increment along the spectroscopic dimension.
 * @param allow_fourth_rank When true, the plan calculates matrices for
 *            processing the fourth rank tensor.
 */
MRS_plan *MRS_create_plan(unsigned int geodesic_polyhedron_frequency,
                          int number_of_sidebands,
                          double sample_rotation_frequency_in_Hz,
                          double rotor_angle_in_rad, double increment,
                          bool allow_fourth_rank);

/**
 * @brief Update the mrsimulator plan.
 *
 * @param	geodesic_polyhedron_frequency The number of triangles along the
 *            edge of the octahedron.
 * @param number_of_sidebands The number of sideband to compute.
 * @param sample_rotation_frequency_in_Hz The sample rotation frequency in Hz.
 * @param rotor_angle_in_rad The polar angle in radians with respect to z-axis
 *            describing the axis of rotation.
 * @param increment The increment along the spectroscopic dimension.
 * @param allow_fourth_rank When true, the plan calculates matrices for
 *            processing the fourth rank tensor.
 */
// MRS_plan *MRS_update_plan(unsigned int *geodesic_polyhedron_frequency,
//                           int *number_of_sidebands,
//                           double *sample_rotation_frequency_in_Hz,
//                           double *rotor_angle_in_rad, double *increment,
//                           bool allow_fourth_rank);

/**
 * @brief Return a copy of the mrsimulator plan.
 *
 * @param *plan The pointer to the plan to be copied.
 * @return MRS_plan = A pointer to the copied plan.
 * This function is incomplete.
 */
MRS_plan *MRS_copy_plan(MRS_plan *plan);

/**
 * @brief Process the plan for the amplitudes at every orientation.
 *
 * The method takes the arguments @p R2 and @p R4 vectors defined in a crystal /
 * commmon frame and evaluates the amplitudes corresponding to the @p R2 and @p
 * R4 vectors in the lab frame. The transformation from the crystal / commmon
 * frame to the lab frame is done using the wigner 2j and 4j rotation matrices
 * over all orientations. The sideband amplitudes are evaluated using equation
 * [39] of the reference https://doi.org/10.1006/jmre.1998.1427.
 *
 * @param	plan A pointer to the mrsimulator plan of type MRS_plan.
 * @param R2 A pointer to the product of the spatial part coefficients of the
 *            second rank tensor and the spin transition functions. The vector
 *            @p R2 is a complex128 array of length 5 with the first element
 *            corresponding to the product of the spin transition function and
 *            the coefficient of the @f$T_{2,-2}@f$ spatial irreducible tensor.
 * @param R4 A pointer to the product of the spatial part coefficients of the
 *            fourth rank tensor and the spin transition functions. The vector
 *            @p R4 is a complex128 array of length 9 with the first element
 *            corresponding to the product of the spin transition function and
 *            the coefficient of the @f$T_{4,-4}@f$ spatial irreducible tensor.
 */
void MRS_get_amplitudes_from_plan(MRS_plan *plan, complex128 *R2,
                                  complex128 *R4);

/**
 * @brief Process the plan for normalized frequencies at every orientation.
 *
 * @param	plan
 * @param	dim
 * @param	R0
 * @param	R2
 * @param	R4
 */
void MRS_get_normalized_frequencies_from_plan(MRS_plan *plan,
                                              MRS_dimension *dim, double R0);

void MRS_get_frequencies_from_plan(MRS_plan *plan, double R0);

/**
 * @brief Return a vector ordered according to the fft output order.
 *
 * @params n The number of points.
 * @params increment The increment along the dimension axis (sampling interval).
 * @returns values A pointer to the fft output order vector of size @p p.
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

extern void __get_components(int number_of_sidebands, double spin_frequency,
                             complex128 *pre_phase);

#endif /* mrsimulator_h */
