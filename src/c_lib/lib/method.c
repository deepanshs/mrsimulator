// -*- coding: utf-8 -*-
//
//  method.c
//
//  @copyright Deepansh J. Srivastava, 2019-2020.
//  Created by Deepansh J. Srivastava, Apr 15, 2020
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
    double increment, double coordinates_offset, int n_events, double *fraction,
    double *sample_rotation_frequency_in_Hz, double *rotor_angle_in_rad,
    double *magnetic_flux_density_in_T, unsigned int number_of_sidebands) {
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
    sequence->events[i].freq_amplitude = malloc_double(the_plan->size);
    vm_double_ones(the_plan->size, sequence->events[i].freq_amplitude);
    MRS_set_event(&(sequence->events[i]), *fraction++,
                  *magnetic_flux_density_in_T++,
                  *sample_rotation_frequency_in_Hz++, *rotor_angle_in_rad++,
                  increment, the_plan);
  }
  sequence->inverse_increment = 1.0 / increment;
  sequence->normalize_offset =
      0.5 - (coordinates_offset * sequence->inverse_increment);
  sequence->R0_offset = 0.0;
  /* buffer to hold the local frequencies and frequency offset. The buffer   *
   * is useful when the rotor angle is off magic angle (54.735 deg). */
  sequence->local_frequency = malloc_double(scheme->total_orientations);
  sequence->freq_offset = malloc_double(scheme->octant_orientations);
}

MRS_sequence *MRS_create_sequences(
    MRS_averaging_scheme *scheme, int *count, double *coordinates_offset,
    double *increment, double *fractions, double *magnetic_flux_density_in_T,
    double *sample_rotation_frequency_in_Hz, double *rotor_angle_in_rad,
    int *n_events, unsigned int n_seq, unsigned int number_of_sidebands) {
  unsigned int i;
  MRS_sequence *sequence = (MRS_sequence *)malloc(n_seq * sizeof(MRS_sequence));

  for (i = 0; i < n_seq; i++) {
    create_plans_for_events_in_sequence(
        &(sequence[i]), scheme, count[i], increment[i], coordinates_offset[i],
        n_events[i], fractions, sample_rotation_frequency_in_Hz,
        rotor_angle_in_rad, magnetic_flux_density_in_T, number_of_sidebands);

    fractions += n_events[i];
    sample_rotation_frequency_in_Hz += n_events[i];
    rotor_angle_in_rad += n_events[i];
    magnetic_flux_density_in_T += n_events[i];

    // printf("Sequence %d\n", i);
    // printf("\tcount %d\n", sequence[i].count);
    // printf("\tincrement %f Hz\n", sequence[i].increment);
    // printf("\tcoordinates offset %f Hz\n", sequence[i].coordinates_offset);
    // for (int j = 0; j < n_events[i]; j++) {
    //   printf("\tEvent %d\n", j);
    //   printf("\t\tfraction %f\n", sequence[i].events[j].fraction);
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

void MRS_free_sequence(MRS_sequence *the_sequence, unsigned int n) {
  unsigned int seq, evt;
  MRS_sequence *sequence;
  for (seq = 0; seq < n; seq++) {
    sequence = &the_sequence[seq];
    for (evt = 0; evt < sequence->n_events; evt++) {
      MRS_free_event(&sequence->events[evt]);
    }
    free(sequence->local_frequency);
    free(sequence->freq_offset);
  }
}
