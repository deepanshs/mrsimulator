// -*- coding: utf-8 -*-
//
//  mrsimulator.h
//
//  @copyright Deepansh J. Srivastava, 2019-2020.
//  Created by Deepansh J. Srivastava, Jun 30, 2019
//  Contact email = deepansh2012@gmail.com
//

#ifndef mrsimulator_h
#define mrsimulator_h

#include "angular_momentum.h"
#include "averaging_scheme.h"
#include "config.h"
#include "fftw3.h"
#include "frequency_tensor.h"
#include "isotopomer_ravel.h"
#include "octahedron.h"

typedef struct MRS_dimension {
  int count; /**<  The number of coordinates along the dimension. */
  double coordinates_offset; /**<  Start coordinate of the dimension. */
  double increment; /**<  Increment of coordinates along the dimension. */

  /* private attributes */
  double normalize_offset;  // fixed value = 0.5 - coordinate_offset/increment
  double inverse_increment;
} MRS_dimension;

MRS_dimension *MRS_create_dimension(int count, double coordinates_offset,
                                    double increment);

/**
 * @struct MRS_plan
 * An mrsimulator plan for computing lineshapes lineshape. An mrsimulator plan,
 * MRS_plan includes,
 *    - a pre-calculated MRS_averaging_scheme.
 *    - pre-calculating stacked arrays of irreducible second rank,
 * wigner-2j(β), and fourth rank, wigner-4j(β), matrices at every
 * orientation angle β,
 *    - pre-calculating the exponent of the sideband order phase,
 *      @f$\exp(-im\alpha)@f$, at every orientation angle α,
 *    - creating the fftw plan, and
 *    - allocating buffer for storing the evaluated frequencies and their
 *      respective amplitudes.
 *
 * Creating a plan adds an overhead to the lineshape simulation. We suggest
 * creating a plan at the start and re-using it as necessary. This is
 * especially efficient when performing a batch simulation, such as,
 * simulating lineshapes from thousands of sites.
 */

struct MRS_plan {
  /**
   * A pointer to the MRS_averaging_scheme orientation averaging scheme.
   */
  MRS_averaging_scheme *averaging_scheme;

  int number_of_sidebands; /**< The number of sidebands to compute. */

  double sample_rotation_frequency_in_Hz; /**< The sample rotation frequency in
                                             Hz. */

  /**
   * The angle, in radians, describing the sample axis-rotation with respect to
   * the lab-frame z-axis.
   */
  double rotor_angle_in_rad;

  /** \privatesection */
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

  /** The buffer to hold the sideband amplitudes as stride 2 array after
   * mrsimulator processing.
   */
  fftw_complex *vector;

  bool allow_fourth_rank;   // If true, creates buffer and tables for processing
                            // fourth rank tensors.
  fftw_plan the_fftw_plan;  //  The plan for fftw routine.
  unsigned int size;        //  number of orientations * number of sizebands.
  unsigned int n_octants;   //  number of octants used in the simulation.
  double *norm_amplitudes;  //  array of normalized amplitudes per orientation.
  double *wigner_d2m0_vector;  //  wigner-2j dm0 vector, n ∈ [-2, 2].
  double *wigner_d4m0_vector;  //  wigner-4j dm0 vector, n ∈ [-4, 4].
  complex128 *pre_phase_2;  //  buffer for 2nk rank sideband phase calculation.
  complex128 *pre_phase_4;  //  buffer for 4th rank sideband phase calculation.
  complex128 one;           //  holds complex value 1.
  complex128 zero;          //  hold complex value 0.
  double buffer;            //  buffer for temporary storage.
};

typedef struct MRS_plan MRS_plan;

/**
 * Create a new mrsimulator plan.
 *
 * @param	scheme The MRS_averaging_scheme.
 * @param number_of_sidebands The number of sideband to compute.
 * @param sample_rotation_frequency_in_Hz The sample rotation frequency in Hz.
 * @param rotor_angle_in_rad The polar angle in radians with respect to z-axis
 *            describing the axis of rotation.
 * @param increment The increment along the spectroscopic dimension.
 * @param allow_fourth_rank When true, the plan calculates matrices for
 *            processing the fourth rank tensor.
 * @return A pointer to the MRS_plan.
 */
MRS_plan *MRS_create_plan(MRS_averaging_scheme *scheme, int number_of_sidebands,
                          double sample_rotation_frequency_in_Hz,
                          double rotor_angle_in_rad, double increment,
                          bool allow_fourth_rank);

/**
 * @brief Release the memory allocated for the given mrsimulator plan.
 *
 * @param	plan The pointer to the MRS_plan.
 */
void MRS_free_plan(MRS_plan *plan);

/* Update the MRS plan for the given rotor angle in radians. */
void MRS_plan_update_rotor_angle_in_rad(MRS_plan *plan,
                                        double rotor_angle_in_rad,
                                        bool allow_fourth_rank);

/**
 * Free the memory from the mrsimulator plan associated with the wigner
 * d^l_{m,0}(rotor_angle_in_rad) vectors. Here, l=2 or 4.
 */
void MRS_plan_free_rotor_angle_in_rad(MRS_plan *plan);

/**
 * @brief Return a copy of the mrsimulator plan.
 *
 * @param plan The pointer to the plan to be copied.
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
 * @param	plan The pointer to the MRS_plan.
 * @param	dim The pointer to the MRS_dimension.
 * @param	R0 The irreducible zeroth rank frequency component.
 */
void MRS_get_normalized_frequencies_from_plan(MRS_plan *plan,
                                              MRS_dimension *dim, double R0);

void MRS_get_frequencies_from_plan(MRS_plan *plan, double R0);

extern void __get_components(int number_of_sidebands, double spin_frequency,
                             double *restrict pre_phase);

#endif /* mrsimulator_h */
