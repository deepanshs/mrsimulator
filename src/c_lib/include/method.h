
// -*- coding: utf-8 -*-
//
//  method.h
//
//  @copyright Deepansh J. Srivastava, 2019-2020.
//  Created by Deepansh J. Srivastava, Apr 15, 2020.
//  Contact email = srivastava.89@osu.edu
//
#include "mrsimulator.h"

#ifndef method_h
#define method_h

typedef struct MRS_event {
  double fraction; /**< The weighted frequency contribution from the event. */
  double magnetic_flux_density_in_T; /**<  The magnetic flux density in T. */
  double rotor_angle_in_rad;         /**<  The rotor angle in radians. */
  double sample_rotation_frequency_in_Hz; /**<  The sample rotation frequency in
                                             Hz. */

  MRS_plan *plan;          /**< The plan for every event. */
  double *freq_amplitude;  // buffer for event amplitude
} MRS_event;

typedef struct MRS_sequence {
  int count;        /**<  The number of coordinates along the dimension. */
  double increment; /**<  Increment of coordinates along the dimension. */
  double coordinates_offset; /**<  Start coordinate of the dimension. */

  MRS_event *events;     /**< Holds a list of events. */
  unsigned int n_events; /**< The number of events. */

  /* private attributes */
  double R0_offset;  // holds the isotropic offset. This is used in determining
                     // if or not to bin the frequencies, especially for
                     // sideband order.
  double *local_frequency;  //  buffer for local frequencies.
  double *freq_offset;      //  buffer for local + sideband frequencies.
  double normalize_offset;  // fixed value = 0.5 - coordinate_offset/increment
  double inverse_increment;
} MRS_sequence;

/**
 * @brief Free the memory allocation for the MRS event.
 *
 * @param the_event The pointer to the MRS_event structs.
 */
void MRS_free_event(MRS_event *the_event);

/**
 * @brief Allocate memory for the MRS sequences.
 *
 * @param	n An interger defining the number of MRS_sequence structs for
 *     which the memory is allocated.
 */
MRS_sequence *MRS_sequence_malloc(int n);

/**
 * @brief Create plans for every event with the array of sequences.
 *
 * @param scheme A pointer to the MRS_averaging_scheme.
 * @param count A pointer to an array of number of points along each sequence.
 * @param coordinates_offset A pointer to an array of coordinates_offsets along
 *      each sequence.
 * @param increment A pointer to an array of increments along each sequence.
 * @param magnetic_flux_density_in_T A pointer to an array of magnetic flux
 *      density given in T along each sequence.
 * @param sample_rotation_frequency_in_Hz A pointer to an array of the sample
 *      rotation frequency given in Hz along each sequence.
 * @param rotor_angle_in_rad A pointer to an array of the the rotor angle given
 *      in radians along each sequence.
 * @param	n_events A pointer to a list of number of events within each
 *      sequence.
 * @param n_seq An unsigned int with the number of sequences.
 * @param number_of_sidebands An int with the number of sidebands.
 */
MRS_sequence *MRS_create_sequences(
    MRS_averaging_scheme *scheme, int *count, double *coordinates_offset,
    double *increment, double *fractions, double *magnetic_flux_density_in_T,
    double *sample_rotation_frequency_in_Hz, double *rotor_angle_in_rad,
    int *n_events, unsigned int n_seq, unsigned int number_of_sidebands);

/**
 * @brief Free the memory allocation for the MRS sequences.
 *
 * @param the_sequence The pointer to an array of MRS_sequence structs.
 * @param	n An interger defining the number of MRS_sequence structs in
 *     `the_sequence`.
 */
void MRS_free_sequence(MRS_sequence *the_sequence, unsigned int n);

#endif /* method_h */
