
// -*- coding: utf-8 -*-
//
//  method.h
//
//  @copyright Deepansh J. Srivastava, 2019-2021.
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
  double rotor_frequency_in_Hz;      /**<  The sample rotation frequency in Hz. */
  MRS_plan *plan;                    /**< The plan for every event. */
  double *freq_amplitude;            // buffer for event amplitude
} MRS_event;

typedef struct MRS_dimension {
  int count;                 /**<  The number of coordinates along the dimension. */
  double increment;          /**<  Increment of coordinates along the dimension. */
  double coordinates_offset; /**<  Start coordinate of the dimension. */
  MRS_event *events;         /**< Holds a list of events. */
  unsigned int n_events;     /**< The number of events. */

  /* private attributes */
  double R0_offset;  // holds the isotropic offset. This is used in determining if or
                     // not to bin the frequencies, especially for sideband order.
  double *local_frequency;  // buffer for local frequencies.
  double *freq_offset;      // buffer for local + sideband frequencies.
  double normalize_offset;  // fixed value = 0.5 - coordinate_offset/increment
  double inverse_increment;
} MRS_dimension;

/**
 * @brief Free the memory allocation for the MRS event.
 *
 * @param the_event The pointer to the MRS_event structs.
 */
void MRS_free_event(MRS_event *the_event);

/**
 * @brief Create plans for every event with the array of dimensions.
 *
 * @param scheme Pointer to the MRS_averaging_scheme.
 * @param count Pointer to number of points array along each dimension.
 * @param coordinates_offset Pointer to coordinates_offsets array along each dimension.
 * @param increment Pointer to increment array along each dimension.
 * @param magnetic_flux_density_in_T Pointer to magnetic flux density array per event.
 * @param rotor_frequency_in_Hz Pointer to rotation frequency array per event.
 * @param rotor_angle_in_rad Pointer to rotor angle array per event.
 * @param	n_events Pointer to number of events list within each dimension.
 * @param n_dim Unsigned int with the number of dimensions.
 * @param number_of_sidebands An int with the number of sidebands.
 * @return MRS_dimension*
 */
MRS_dimension *MRS_create_dimensions(
    MRS_averaging_scheme *scheme, int *count, double *coordinates_offset,
    double *increment, double *fractions, double *magnetic_flux_density_in_T,
    double *rotor_frequency_in_Hz, double *rotor_angle_in_rad, int *n_events,
    unsigned int n_dim, unsigned int number_of_sidebands);

/**
 * @brief Free the memory allocation for the MRS dimensions.
 *
 * @param dimensions The pointer to an array of MRS_dimension structs.
 * @param	n An interger defining the number of MRS_dimension structs in `dimensions`.
 */
void MRS_free_dimension(MRS_dimension *dimensions, unsigned int n);

#endif /* method_h */
