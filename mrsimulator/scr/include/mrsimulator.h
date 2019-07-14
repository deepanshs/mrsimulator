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

/**
 * @struct MRS_plan_t
 * @brief Create a mrsimulator plan for faster lineshape simulation.
 */
struct MRS_plan_t {
  /**
   * The number of triangles along the edge of octahedron. This value is a
   * positive integer which represents the frequency of class I geodesic
   * polyhedra. These polyhedra may be used in calculating the spherical
   * average. Currently, we only use octahedral as the frequency 1 polyhedra. As
   * the frequency of the geodesic polyhedron increases, the polyhedra approach
   * to a sphere geometry. For line-shape simulation, a higher geodesic
   * polyhedron frequency will result in a better spherical averaging. The
   * default value is 72. Read more on the <a
   * href="https://en.wikipedia.org/wiki/Geodesic_polyhedron">Geodesic
   * polyhedron</a>.
   */
  unsigned short geodesic_polyhedron_frequency;

  int number_of_sidebands; /**< The number of sidebands to compute. */

  double sample_rotation_frequency; /**< The sample rotation frequency in Hz. */

  /**
   * The polar angle, in radians, describing the axis of rotation of the sample
   * with respect to the lab-frame z-axis.
   */
  double rotor_angle;

  bool allow_fourth_rank; /**< Allow buffer for fourth rank tensors. */

  /**
   * The sideband frequency ratio stored in the fft output order. The sideband
   * frequency ratio is defined as the ratio -
   *    @f[\frac{n \omega_r}{n_i}@f]
   * where `n` is an integer, @f$\omega_r@f$ is the spinning frequency frequency
   * in Hz, and @f$n_i@f$ is the `increment` along the spectroscopic grid
   * dimension.
   */
  double *vr_freq;

  /** The isotropic frequency offset ratio. The ratio is similarly defined as
   * before. */
  double isotropic_offset;

  // The buffer to hold the sideband amplitudes as stride 2 array after
  // mrsimulator processing.
  fftw_complex *vector;

  /* private attributes */

  fftw_plan the_fftw_plan;      // The plan for fftw routine.
  unsigned int n_orientations;  // number of unique orientations.
  unsigned int size;            // number of orientations * number of sizebands.
  double *amplitudes;           // array of amplitude scaling per orientation.
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
  double complex one;          // holds complex value 1.
  double complex zero;         // hold complex value 0.
  double buffer;               // buffer for temporary storage.
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
 * @param sample_rotation_frequency The sample rotation frequency in Hz.
 * @param rotor_angle The polar angle in radians with respect to z-axis
 *            describing the axis of rotation.
 * @param increment The increment along the spectroscopic dimension.
 * @param allow_fourth_rank When true, the plan calculates matrices for
 *            processing the fourth rank tensor.
 */
MRS_plan *MRS_create_plan(unsigned int geodesic_polyhedron_frequency,
                          int number_of_sidebands,
                          double sample_rotation_frequency, double rotor_angle,
                          double increment, bool allow_fourth_rank);

/**
 * @brief Update the mrsimulator plan.
 *
 * @param	geodesic_polyhedron_frequency The number of triangles along the
 *            edge of the octahedron.
 * @param number_of_sidebands The number of sideband to compute.
 * @param sample_rotation_frequency The sample rotation frequency in Hz.
 * @param rotor_angle The polar angle in radians with respect to z-axis
 *            describing the axis of rotation.
 * @param increment The increment along the spectroscopic dimension.
 * @param allow_fourth_rank When true, the plan calculates matrices for
 *            processing the fourth rank tensor.
 */
// MRS_plan *MRS_update_plan(unsigned int *geodesic_polyhedron_frequency,
//                           int *number_of_sidebands,
//                           double *sample_rotation_frequency,
//                           double *rotor_angle, double *increment,
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
 *            @p R2 is a double complex array of length 5 with the first element
 *            corresponding to the product of the spin transition function and
 *            the coefficient of the @f$T_{2,-2}@f$ spatial irreducible tensor.
 * @param R4 A pointer to the product of the spatial part coefficients of the
 *            fourth rank tensor and the spin transition functions. The vector
 *            @p R4 is a double complex array of length 9 with the first element
 *            corresponding to the product of the spin transition function and
 *            the coefficient of the @f$T_{4,-4}@f$ spatial irreducible tensor.
 */
void MRS_get_amplitudes_from_plan(MRS_plan *plan, double complex *R2,
                                  double complex *R4);

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
                                              MRS_dimension *dim, double R0,
                                              double complex *R2,
                                              double complex *R4);

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
                             double complex *pre_phase);

#endif /* mrsimulator_h */
