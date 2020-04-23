// -*- coding: utf-8 -*-
//
//  method.c
//
//  @copyright Deepansh J. Srivastava, 2019-2020.
//  Created by Deepansh J. Srivastava, Apr 15, 2020
//  Contact email = deepansh2012@gmail.com
//

#include "method.h"

/* free the buffer and pre-calculated tables from the mrsimulator plan. */
void MRS_free_event(MRS_event *the_event) { MRS_free_plan(the_event->plan); }

void MRS_set_event(MRS_event *event, double magnetic_flux_density_in_T,
                   double sample_rotation_frequency_in_Hz,
                   double rotor_angle_in_rad, int increment, MRS_plan *plan) {
  event->sample_rotation_frequency_in_Hz = sample_rotation_frequency_in_Hz;
  event->rotor_angle_in_rad = rotor_angle_in_rad;
  event->magnetic_flux_density_in_T = magnetic_flux_density_in_T;

  if (sample_rotation_frequency_in_Hz ==
          plan->sample_rotation_frequency_in_Hz &&
      rotor_angle_in_rad == plan->rotor_angle_in_rad) {
    event->plan = plan;
    return;
  }

  if (sample_rotation_frequency_in_Hz !=
          plan->sample_rotation_frequency_in_Hz &&
      rotor_angle_in_rad == plan->rotor_angle_in_rad) {
    MRS_plan *new_plan = MRS_copy_plan(plan);
    MRS_plan_update_from_sample_rotation_frequency_in_Hz(
        new_plan, increment, sample_rotation_frequency_in_Hz);
    event->plan = new_plan;
    return;
  }

  if (sample_rotation_frequency_in_Hz ==
          plan->sample_rotation_frequency_in_Hz &&
      rotor_angle_in_rad != plan->rotor_angle_in_rad) {
    MRS_plan *new_plan = MRS_copy_plan(plan);
    MRS_plan_update_from_rotor_angle_in_rad(new_plan, rotor_angle_in_rad,
                                            plan->allow_fourth_rank);
    event->plan = new_plan;
    return;
  }
  return;
}

MRS_sequence *MRS_sequence_malloc(int n) {
  MRS_sequence *the_sequence = (MRS_sequence *)malloc(n * sizeof(MRS_sequence));
  return the_sequence;
}

static inline void create_plans_for_events_in_sequence(
    MRS_sequence *sequence, MRS_averaging_scheme *scheme, int count,
    double increment, double coordinates_offset, int n_events,
    double *sample_rotation_frequency_in_Hz, double *rotor_angle_in_rad,
    double *magnetic_flux_density_in_T, int number_of_sidebands) {
  int i;
  sequence->count = count;
  sequence->coordinates_offset = coordinates_offset;
  sequence->increment = increment;
  sequence->n_events = n_events;
  sequence->events = (MRS_event *)malloc(n_events * sizeof(MRS_event));

  MRS_plan *the_plan = MRS_create_plan(
      scheme, number_of_sidebands, sample_rotation_frequency_in_Hz[0],
      rotor_angle_in_rad[0], increment, scheme->allow_fourth_rank);

  for (i = 0; i < n_events; i++) {
    MRS_set_event(&(sequence->events[i]), *magnetic_flux_density_in_T++,
                  *sample_rotation_frequency_in_Hz++, *rotor_angle_in_rad++,
                  increment, the_plan);
  }
  sequence->inverse_increment = 1.0 / increment;
  sequence->normalize_offset =
      0.5 - (coordinates_offset * sequence->inverse_increment);
}

MRS_sequence *MRS_create_plans_for_sequence(
    MRS_averaging_scheme *scheme, int *count, double *coordinates_offset,
    double *increment, double *magnetic_flux_density_in_T,
    double *sample_rotation_frequency_in_Hz, double *rotor_angle_in_rad,
    int *n_events, unsigned int n_seq, int number_of_sidebands) {
  unsigned int i;
  MRS_sequence *sequence = (MRS_sequence *)malloc(n_seq * sizeof(MRS_sequence));

  for (i = 0; i < n_seq; i++) {
    create_plans_for_events_in_sequence(
        &(sequence[i]), scheme, count[i], increment[i], coordinates_offset[i],
        n_events[i], sample_rotation_frequency_in_Hz, rotor_angle_in_rad,
        magnetic_flux_density_in_T, number_of_sidebands);

    sample_rotation_frequency_in_Hz += n_events[i];
    rotor_angle_in_rad += n_events[i];
    magnetic_flux_density_in_T += n_events[i];

    // printf("Sequence %d\n", i);
    // printf("\tcount %d\n", sequence[i].count);
    // printf("\tincrement %f Hz\n", sequence[i].increment);
    // printf("\tcoordinates offset %f Hz\n", sequence[i].coordinates_offset);
    // for (j = 0; j < n_events[i]; j++) {
    //   printf("\tEvent %d\n", j);
    //   printf("\t\trotor frequency %f Hz\n",
    //          sequence[i].events[j].sample_rotation_frequency_in_Hz);
    //   printf("\t\trotor angle %f rad\n",
    //          sequence[i].events[j].rotor_angle_in_rad);
    //   printf("\t\tmagnetic flux density %f T\n",
    //          sequence[i].events[j].magnetic_flux_density_in_T);
    //   printf("\n");
    // }
  }
  return sequence;
}
