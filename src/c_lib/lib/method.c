// -*- coding: utf-8 -*-
//
//  method.c
//
//  @copyright Deepansh J. Srivastava, 2019-2021.
//  Created by Deepansh J. Srivastava, Apr 15, 2020.
//  Contact email = srivastava.89@osu.edu
//

#include "method.h"

/* free the buffer and pre-calculated tables from the mrsimulator plan. */
void MRS_free_event(MRS_event *the_event) {
  if (!the_event->plan) {
    MRS_free_plan(the_event->plan);
  }
  free(the_event->freq_amplitude);
}

static inline void MRS_set_event(MRS_event *event, double fraction,
                                 double magnetic_flux_density_in_T,
                                 double sample_rotation_frequency_in_Hz,
                                 double rotor_angle_in_rad, double increment,
                                 MRS_plan *plan) {
  event->fraction = fraction;
  event->sample_rotation_frequency_in_Hz = sample_rotation_frequency_in_Hz;
  event->rotor_angle_in_rad = rotor_angle_in_rad;
  event->magnetic_flux_density_in_T = magnetic_flux_density_in_T;

  if (sample_rotation_frequency_in_Hz == plan->sample_rotation_frequency_in_Hz &&
      rotor_angle_in_rad == plan->rotor_angle_in_rad) {
    event->plan = plan;
    return;
  }

  if (sample_rotation_frequency_in_Hz != plan->sample_rotation_frequency_in_Hz &&
      rotor_angle_in_rad == plan->rotor_angle_in_rad) {
    MRS_plan *new_plan = MRS_copy_plan(plan);
    MRS_plan_update_from_sample_rotation_frequency_in_Hz(
        new_plan, increment, sample_rotation_frequency_in_Hz);
    event->plan = new_plan;
    return;
  }

  if (sample_rotation_frequency_in_Hz == plan->sample_rotation_frequency_in_Hz &&
      rotor_angle_in_rad != plan->rotor_angle_in_rad) {
    MRS_plan *new_plan = MRS_copy_plan(plan);
    MRS_plan_update_from_rotor_angle_in_rad(new_plan, rotor_angle_in_rad,
                                            plan->allow_fourth_rank);
    event->plan = new_plan;
    return;
  }
  return;
}

MRS_dimension *MRS_dimension_malloc(int n) {
  MRS_dimension *dimensions = (MRS_dimension *)malloc(n * sizeof(MRS_dimension));
  return dimensions;
}

static inline void create_plans_for_events_in_dimension(
    MRS_dimension *dimension, MRS_averaging_scheme *scheme, int count, double increment,
    double coordinates_offset, int n_events, double *fraction,
    double *sample_rotation_frequency_in_Hz, double *rotor_angle_in_rad,
    double *magnetic_flux_density_in_T, unsigned int number_of_sidebands) {
  int i;
  dimension->count = count;
  dimension->coordinates_offset = coordinates_offset;
  dimension->increment = increment;
  dimension->n_events = n_events;
  dimension->events = (MRS_event *)malloc(n_events * sizeof(MRS_event));

  MRS_plan *the_plan =
      MRS_create_plan(scheme, number_of_sidebands, sample_rotation_frequency_in_Hz[0],
                      rotor_angle_in_rad[0], increment, scheme->allow_fourth_rank);

  for (i = 0; i < n_events; i++) {
    dimension->events[i].freq_amplitude = malloc_double(the_plan->size);
    vm_double_ones(the_plan->size, dimension->events[i].freq_amplitude);
    MRS_set_event(&(dimension->events[i]), *fraction++, *magnetic_flux_density_in_T++,
                  *sample_rotation_frequency_in_Hz++, *rotor_angle_in_rad++, increment,
                  the_plan);
  }
  dimension->inverse_increment = 1.0 / increment;
  dimension->normalize_offset =
      0.5 - (coordinates_offset * dimension->inverse_increment);
  dimension->R0_offset = 0.0;
  /* buffer to hold the local frequencies and frequency offset. The buffer   *
   * is useful when the rotor angle is off magic angle (54.735 deg). */
  dimension->local_frequency = malloc_double(scheme->total_orientations);
  dimension->freq_offset = malloc_double(scheme->octant_orientations);
}

MRS_dimension *MRS_create_dimensions(
    MRS_averaging_scheme *scheme, int *count, double *coordinates_offset,
    double *increment, double *fractions, double *magnetic_flux_density_in_T,
    double *sample_rotation_frequency_in_Hz, double *rotor_angle_in_rad, int *n_events,
    unsigned int n_dim, unsigned int number_of_sidebands) {
  unsigned int i;
  MRS_dimension *dimension = (MRS_dimension *)malloc(n_dim * sizeof(MRS_dimension));

  for (i = 0; i < n_dim; i++) {
    create_plans_for_events_in_dimension(
        &(dimension[i]), scheme, count[i], increment[i], coordinates_offset[i],
        n_events[i], fractions, sample_rotation_frequency_in_Hz, rotor_angle_in_rad,
        magnetic_flux_density_in_T, number_of_sidebands);

    fractions += n_events[i];
    sample_rotation_frequency_in_Hz += n_events[i];
    rotor_angle_in_rad += n_events[i];
    magnetic_flux_density_in_T += n_events[i];

    // printf("dimension %d\n", i);
    // printf("\tcount %d\n", dimension[i].count);
    // printf("\tincrement %f Hz\n", dimension[i].increment);
    // printf("\tcoordinates offset %f Hz\n", dimension[i].coordinates_offset);
    // for (int j = 0; j < n_events[i]; j++) {
    //   printf("\tEvent %d\n", j);
    //   printf("\t\tfraction %f\n", dimension[i].events[j].fraction);
    //   printf("\t\trotor frequency %f Hz\n",
    //          dimension[i].events[j].sample_rotation_frequency_in_Hz);
    //   printf("\t\trotor angle %f rad\n",
    //          dimension[i].events[j].rotor_angle_in_rad);
    //   printf("\t\tmagnetic flux density %f T\n",
    //          dimension[i].events[j].magnetic_flux_density_in_T);
    //   printf("\n");
    // }
  }
  return dimension;
}

void MRS_free_dimension(MRS_dimension *dimensions, unsigned int n) {
  unsigned int dim, evt;
  MRS_dimension *dimension;
  for (dim = 0; dim < n; dim++) {
    dimension = &dimensions[dim];
    for (evt = 0; evt < dimension->n_events; evt++) {
      MRS_free_event(&dimension->events[evt]);
    }
    free(dimension->local_frequency);
    free(dimension->freq_offset);
  }
}
