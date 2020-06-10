
// -*- coding: utf-8 -*-
//
//  method.h
//
//  @copyright Deepansh J. Srivastava, 2019-2020.
//  Created by Deepansh J. Srivastava, Apr 15, 2020.
//  Contact email = deepansh2012@gmail.com
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

  MRS_plan *plan; /**< The plan for every event. */
} MRS_event;

typedef struct MRS_sequence {
  int count;        /**<  The number of coordinates along the dimension. */
  double increment; /**<  Increment of coordinates along the dimension. */
  double coordinates_offset; /**<  Start coordinate of the dimension. */

  MRS_event *events;     /**< Holds a list of events. */
  unsigned int n_events; /**< The number of events. */

  /* private attributes */
  double normalize_offset;  // fixed value = 0.5 - coordinate_offset/increment
  double inverse_increment;
} MRS_sequence;

void MRS_set_event(MRS_event *event, double magnetic_flux_density_in_T,
                   double sample_rotation_frequency_in_Hz,
                   double rotor_angle_in_rad, int increment, MRS_plan *plan);

void MRS_free_event(MRS_event *the_event);

MRS_sequence *MRS_sequence_malloc(int n);

MRS_sequence *MRS_create_plans_for_sequence(
    MRS_averaging_scheme *scheme, int *count, double *coordinates_offset,
    double *increment, double *magnetic_flux_density_in_T,
    double *sample_rotation_frequency_in_Hz, double *rotor_angle_in_rad,
    int *n_events, unsigned int n_seq, int number_of_sidebands);

#endif /* method_h */
