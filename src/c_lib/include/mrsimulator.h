// -*- coding: utf-8 -*-
//
//  mrsimulator.h
//
//  @copyright Deepansh J. Srivastava, 2019-2021.
//  Created by Deepansh J. Srivastava, Jun 30, 2019.
//  Contact email = srivastava.89@osu.edu
//

#ifndef mrsimulator_h
#define mrsimulator_h

#include "angular_momentum.h"
#include "config.h"
#include "fftw3.h"
#include "frequency/frequency_tensor_components.h"
#include "object_struct.h"
#include "schemes.h"
/**
 * @struct MRS_plan
 * A mrsimulator plan for computing spectra. The plan, MRS_plan includes,
 *    - a pre-calculated MRS_averaging_scheme,
 *    - pre-calculates stacked arrays of irreducible second-rank, wigner-2j(β), and
 *      fourth-rank, wigner-4j(β), matrices at every orientation angle β,
 *    - pre-calculates the exponent of the sideband order phase, @f$\exp(-im\alpha)@f$,
 *      at every orientation angle α,
 *    - create the fftw plan, and
 *    - allocating buffer for storing the evaluated frequencies and their respective
 *      amplitudes.
 *
 * Creating a plan adds an overhead to the simulation. Recommendation is to create the
 * plan at the start and re-using it as necessary. This is especially efficient when
 * performing a batch simulation, such as, simulating spectra from thousands of spin
 * systems.
 */

struct MRS_plan {
  /**
   * A pointer to the MRS_averaging_scheme orientation averaging scheme.
   */
  MRS_averaging_scheme *averaging_scheme;

  unsigned int number_of_sidebands; /**< The number of sidebands to compute. */
  double rotor_frequency_in_Hz;     /**< The sample rotation frequency in  Hz. */

  /**
   * The angle, in radians, describing the sample axis-rotation with respect to
   * the lab-frame z-axis.
   */
  double rotor_angle_in_rad;

  /** \privatesection */
  /**
   * A pointer to an array of sideband frequency ratio stored in the fft output order.
   * The sideband frequency ratio is defined as the ratio -
   *    @f[\frac{n \omega_r}{n_i}@f]
   * where `n` is an integer, @f$\omega_r@f$ is the spinning frequency in Hz, and
   * @f$n_i@f$ is the `increment` along the spectroscopic grid dimension.
   */
  double *vr_freq;

  /** The buffer to hold the sideband amplitudes as stride 2 array after mrsimulator
   * processing.
   */

  bool allow_fourth_rank;      // If true, creates buffer/tables for 4th-rank tensors.
  unsigned int size;           // # of angular orientations * number of sizebands.
  unsigned int n_octants;      // # of octants used in the orientational averaging.
  double *norm_amplitudes;     // array of normalized amplitudes per orientation.
  double *wigner_d2m0_vector;  // wigner-2j dm0 vector, n ∈ [-2, 2].
  double *wigner_d4m0_vector;  // wigner-4j dm0 vector, n ∈ [-4, 4].
  complex128 *pre_phase;       // temp buffer to hold sideband phase calculation.
  complex128 *pre_phase_2;     // buffer for 2nk rank sideband phase calculation.
  complex128 *pre_phase_4;     // buffer for 4th rank sideband phase calculation.
  double buffer;               // buffer for temporary storage.
};

typedef struct MRS_plan MRS_plan;

/**
 * @brief Create a new mrsimulator plan.
 *
 * @param scheme The MRS_averaging_scheme.
 * @param number_of_sidebands The number of sidebands.
 * @param rotor_frequency_in_Hz The sample rotation frequency in Hz.
 * @param rotor_angle_in_rad The polar angle in radians with respect to the
 *          z-axis describing the axis of rotation.
 * @param increment The increment along the spectroscopic dimension in Hz.
 * @param allow_fourth_rank When true, the plan calculates matrices for
 *          processing the fourth-rank tensors.
 * @return A pointer to the MRS_plan.
 */
MRS_plan *MRS_create_plan(MRS_averaging_scheme *scheme,
                          unsigned int number_of_sidebands,
                          double rotor_frequency_in_Hz, double rotor_angle_in_rad,
                          double increment, bool allow_fourth_rank);

/**
 * @brief Release the memory allocated for the given mrsimulator plan.
 *
 * @param plan The pointer to the MRS_plan.
 */
void MRS_free_plan(MRS_plan *plan);

/* Update the MRS plan when sample rotation frequency is changed. */
void MRS_plan_update_from_rotor_frequency_in_Hz(MRS_plan *plan, double increment,
                                                double rotor_frequency_in_Hz);

/* Update the MRS plan when the rotor angle is changed. */
void MRS_plan_update_from_rotor_angle_in_rad(MRS_plan *plan, double rotor_angle_in_rad,
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
 * @brief Process the plan and evaluates the amplitudes at every orientation.
 *
 * The method takes the arguments @p R2 and @p R4 vectors defined in a crystal / commmon
 * frame and evaluates the amplitudes corresponding to the @p R2 and @p R4 vectors in
 * the lab frame. The transformation from the crystal / commmon frame to the lab frame
 * is done using the wigner 2j and 4j rotation matrices over all orientations. The
 * sideband amplitudes are evaluated using equation [39] of the reference
 * https://doi.org/10.1006/jmre.1998.1427.
 *
 * @param scheme The pointer to the powder averaging scheme of type
 *            MRS_averaging_scheme.
 * @param plan A pointer to the mrsimulator plan of type MRS_plan.
 * @param fftw_scheme A pointer to the fftw scheme of type MRS_fftw_scheme.
 * @param refresh If true, zero the output array before proceeding, else add to
 *            the existing array.
 */
void MRS_get_amplitudes_from_plan(MRS_averaging_scheme *scheme, MRS_plan *plan,
                                  MRS_fftw_scheme *fftw_scheme, bool refresh);

// Important: `method.h` header file must be included after defining MRS_plan.
#include "method.h"

/**
 * @brief Process the plan and evaluates normalized frequencies at every orientation.
 * Normalization refers to scaling the frequencies by the corresponding spectroscopic
 * dimension increment.
 *
 * @param scheme The pointer to the powder averaging scheme of type
 *      MRS_averaging_scheme.
 * @param plan A pointer to the mrsimulator plan of type MRS_plan.
 * @param R0 The irreducible zeroth-rank frequency component.
 * @param R2 A pointer to an array of second-rank frequency components. The frequency
 *      components are the product of the size of interaction, spatial symmetry
 *      functions, and the spin transition functions. The vector @p R2 is a complex128
 *      array of length 5, ordered as m = [-2, -1, 0, 1, 2].
 * @param R4 A pointer to an array of fourth-rank frequency components. The frequency
 *      components are the product of the size of interaction, spatial symmetry
 *      functions, and the spin transition functions. The vector @p R4 is a complex128
 *      array of length 9, ordered as m = [-4, -3, -2, -1, 0, 1, 2, 3, 4].
 * @param refresh If true, zero the frequencies before update, else self update.
 * @param dim The pointer to the dimension of type MRS_dimension.
 * @param fraction A float representing the fraction of dimension during an event.
 */
void MRS_get_normalized_frequencies_from_plan(MRS_averaging_scheme *scheme,
                                              MRS_plan *plan, double R0, complex128 *R2,
                                              complex128 *R4, bool refresh,
                                              MRS_dimension *dim, double fraction);

void MRS_get_frequencies_from_plan(MRS_averaging_scheme *scheme, MRS_plan *plan,
                                   double R0, complex128 *R2, complex128 *R4,
                                   bool refresh, MRS_dimension *dim);

/**
 * @brief The function rotates the tensor components from the principal axis system
 * (PAS) to the common frame of the spin system.
 *
 * @param sites A pointer to the site_struct structure.
 * @param couplings A pointer to the coupling_struct structure.
 * @param transition A pointer to the spin quantum numbers from the inital and final
 *      states of the spin transition packed as initial quantum numbers followed by the
 *      final quantum numbers.
 * @param allow_fourth_rank A boolean, if true, evalutes the frequency contributions
 *      from the fourth-rank tensor.
 * @param R0 A pointer to an array where the frequency contribution from the zeroth-rank
 *      tensor is stored.
 * @param R2 A pointer to a complex array where the frequency contributions from the
 *      second-rank tensor are stored.
 * @param R4 A pointer to a complex array where the frequency contributions from the
 *      fourth-rank tensor are stored.
 * @param R0_temp A pointer to an array where the frequency contribution from the
 *      zeroth-rank tensor is temporarily stored.
 * @param R2_temp A pointer to a complex array where the frequency contributions from
 *      the second-rank tensor are temporarily stored.
 * @param R4_temp A pointer to a complex array where the frequency contributions from
 *      the fourth-rank tensor are temporarily stored.
 * @param B0_in_T The magnetic flux density of the macroscopic external magnetic field
 *      in T.
 * @param freq_contrib A pointer to a stack of boolean frequency contribs repeated
 *      number_of_events times. The order of the freq contribs follow
 *          	1. Shielding 1st order 0th rank
 *          	2. Shielding 1st order 2th rank
 *          	3. Quad 1st order 2th rank
 *          	4. Quad 2st order 0th rank
 *          	5. Quad 2st order 2th rank
 *          	6. Quad 2st order 4th rank
 */
void MRS_rotate_components_from_PAS_to_common_frame(
    site_struct *sites,          // Pointer to a list of sites in the spin system.
    coupling_struct *couplings,  // Pointer to a list of couplings within a spin system.
    float *transition,           // The pointer to the spin transition.
    bool allow_fourth_rank,      // If true, pre for 4th rank computation.
    double *R0,                  // The R0 components.
    complex128 *R2,              // The R2 components.
    complex128 *R4,              // The R4 components.
    double *R0_temp,             // The temporary R0 components.
    complex128 *R2_temp,         // The temporary R2 components.
    complex128 *R4_temp,         // The temporary R3 components.
    double B0_in_T,              // Magnetic flux density in T.
    bool *freq_contrib           // The pointer to freq contribs boolean.
);

extern void get_sideband_phase_components(unsigned int number_of_sidebands,
                                          double spin_frequency,
                                          double *restrict pre_phase);

#endif /* mrsimulator_h */
