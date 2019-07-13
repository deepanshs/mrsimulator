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


typedef struct MRS_plan_t {
  // number of triangles along  the edge of octahedron.
  unsigned short geodesic_polyhedron_frequency;
  int number_of_sidebands;          // The number of sidebands to compute.
  double sample_rotation_frequency; // sample rotation frequency in Hz.
  double rotor_angle;     // polar angle describing the axis of rotation in rad.
  bool ALLOW_FOURTH_RANK; // Allow buffer for fourth rank tensors.
  double *vr_freq;        // sideband order frequencies in fft output order.
  double isotropic_offset; // isotropic frequency offset.
  fftw_complex *vector; // array of sideband amplitudes stored as strides of 2.

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
} MRS_plan;

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

void MRS_free_plan(MRS_plan *the_plan);

/* Create a new mrsimulator plan.
 * @func MRS_create_plan.
 *
 * @author	Deepansh J. Srivastava
 * @since	v0.1.1
 * @version	v0.1.1	Friday, July 12th, 2019.
 *
 * @param	unsigned int geodesic_polyhedron_frequency = The number of
 *      triangles along the edge of the octahedron.
 * @param	int number_of_sidebands = The number of sideband to compute.
 * @param	double sample_rotation_frequency = The sample rotation
 *      frequency in Hz.
 * @param	double rotor_angle = The polar angle in radians with respect to
 *      z-axis describing the axis of rotation.
 * @param	double increment = The increment along the spectroscopic
 *      dimension.
 * @param	bool ALLOW_FOURTH_RANK = When true, the plan calculates
 *      matrices for processing the fourth rank tensor.
 * @return	void
 */
MRS_plan *MRS_create_plan(unsigned int geodesic_polyhedron_frequency,
                          int number_of_sidebands,
                          double sample_rotation_frequency, double rotor_angle,
                          double increment, bool ALLOW_FOURTH_RANK);

/* Update the mrsimulator plan.
//  * @func MRS_update_plan.
//  *
//  * @author	Deepansh J. Srivastava
//  * @since	v0.1.1
//  * @version	v0.1.1	Friday, July 12th, 2019.
//  *
//  * @param	int   	geodesic_polyhedron_frequency = The number of triangles
//  *      along the edge of the octahedron.
//  * @param	int   	number_of_sidebands = The number of sideband to compute.
//  * @param	double	sample_rotation_frequency = The sample rotation
//  *      frequency in Hz.
//  * @param	double	rotor_angle = The polar angle in radians with respect to
//  *      z-axis describing the axis of rotation.
//  * @param	double	increment = The increment along the spectroscopic
//  *      dimension.
//  * @param	bool  	ALLOW_FOURTH_RANK = When true, the plan calculates
//  *      matrices for processing the fourth rank tensor.
//  * @return	void
//  */
// MRS_plan *MRS_update_plan(unsigned int *geodesic_polyhedron_frequency,
//                           int *number_of_sidebands,
//                           double *sample_rotation_frequency,
//                           double *rotor_angle, double *increment,
//                           bool ALLOW_FOURTH_RANK);

/* Return a copy of thr mrsimulator plan.
 * @func MRS_copy_plan.
 *
 * @author	Deepansh J. Srivastava
 * @since	v0.1.1
 * @version	v0.1.1	Friday, July 12th, 2019.
 * @param MRS_plan *plan = The pointer to the plan to be copied.
 * @return MRS_plan = A pointer to the copied plan.
 * This function is incomplete.
 */
MRS_plan *MRS_copy_plan(MRS_plan *plan);

/* Compute the R2 and R4 in the lab frame using wigner 2j and 4j rotation
 * matrices at all orientations. Evalute the sideband amplitudes
 * Equation [39] in the refernce https://doi.org/10.1006/jmre.1998.1427.
 *
 * @func MRS_get_amplitudes_from_plan.
 *
 * @author	Deepansh J. Srivastava
 * @since	v0.1.1
 * @version	v0.1.1	Friday, July 12th, 2019.
 *
 * @param	MRS_plan       *plan = The pointer to the mrsimulator plan.
 * @param	double complex *R2 = The pointer to the product of the spatial
 *      part coefficient of the second rank tensor and the spin transition
 *      function. R2 is an array of length 5 with the first element
 *      corresponding to the product of the spin transition function and the
 *      coefficient of the T_{2,-2} irreducible tensor.
 * @param double complex *R4 = The pointer to the product of the spatial part
 *      coefficient of the fourth rank tensor and the spin transition function.
 *      R4 is an array of length 9 with the first element corresponding to the
 *      product of the spin transition function and the coefficient of the
 *      T_{4,-4} irreducible tensor.
 * @return void
 */
void MRS_get_amplitudes_from_plan(MRS_plan *plan, double complex *R2,
                                  double complex *R4);

/**
 * MRS_get_normalized_frequencies_from_plan.
 *
 * @author	Deepansh J. Srivastava
 * @since	v0.0.1
 * @version	v1.0.0	Friday, July 12th, 2019.
 * @global
 * @param	mrs_plan     	*plan
 * @param	mrs_dimension	*dim
 * @param	double       	R0
 * @param	complex      	*R2
 * @param	complex      	*R4
 * @return	mixed
 */
void MRS_get_normalized_frequencies_from_plan(MRS_plan *plan,
                                              MRS_dimension *dim, double R0,
                                              double complex *R2,
                                              double complex *R4);

/* Return a vector ordered according to the fft output order. *
 * @params int n - The number of points *
 * @params double increment - The increment (sampling interval) *
 * @returns *double values = The pointer to the fft output order vector */
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
