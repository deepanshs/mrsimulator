// -*- coding: utf-8 -*-
//
//  method.c
//
//  @copyright Deepansh J. Srivastava, 2019-2021.
//  Created by Deepansh J. Srivastava, Apr 15, 2020.
//  Contact email = srivastava.89@osu.edu
//

#include "method.h"

/* Free the buffer and pre-calculated tables from the mrsimulator plan. */
void MRS_free_event(MRS_event *the_event) {
  if (!the_event->plan) MRS_free_plan(the_event->plan);
  free(the_event->freq_amplitude);
}

/** Free the memory/buffer allocation for the MRS dimensions and events within. **/
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

/**
 * @brief Polulates event structs and creates/updates new plans based on the plan from
 * the first event along the dimension.
 *
 * @param event The pointer to the event.
 * @param fraction The fraction/weight of the event.
 * @param magnetic_flux_density_in_T The external field flux density at the event.
 * @param rotor_frequency_in_Hz The rotor frequency at the event.
 * @param rotor_angle_in_rad The rotor angle at the event.
 * @param increment The increment of the dimension to which the event belongs.
 * @param plan The MRS plan of the first event along the dimension.
 */
static inline void MRS_set_event(MRS_event *event, double fraction,
                                 double magnetic_flux_density_in_T,
                                 double rotor_frequency_in_Hz,
                                 double rotor_angle_in_rad, double increment,
                                 MRS_plan *plan) {
  event->fraction = fraction;
  event->rotor_frequency_in_Hz = rotor_frequency_in_Hz;
  event->rotor_angle_in_rad = rotor_angle_in_rad;
  event->magnetic_flux_density_in_T = magnetic_flux_density_in_T;

  bool rotor_frequency_equal = rotor_frequency_in_Hz == plan->rotor_frequency_in_Hz;
  bool rotor_angle_equal = rotor_angle_in_rad == plan->rotor_angle_in_rad;

  /* When both rotor angle and rotor freq is the same as the plan, return plan. */
  if (rotor_frequency_equal && rotor_angle_equal) {
    event->plan = plan;
    return;
  }

  /* When only rotor freq is different, update the plan accordingly and return. */
  if (!rotor_frequency_equal && rotor_angle_equal) {
    MRS_plan *new_plan = MRS_copy_plan(plan);
    MRS_plan_update_from_rotor_frequency_in_Hz(new_plan, increment,
                                               rotor_frequency_in_Hz);
    event->plan = new_plan;
    return;
  }

  /* When only rotor angle is different, update the plan accordingly and return. */
  if (rotor_frequency_equal && !rotor_angle_equal) {
    MRS_plan *new_plan = MRS_copy_plan(plan);
    MRS_plan_update_from_rotor_angle_in_rad(new_plan, rotor_angle_in_rad,
                                            plan->allow_fourth_rank);
    event->plan = new_plan;
    return;
  }

  /* Otherwise, update plan for both rotor angle and freq and return. */
  MRS_plan *new_plan = MRS_copy_plan(plan);
  MRS_plan_update_from_rotor_frequency_in_Hz(new_plan, increment,
                                             rotor_frequency_in_Hz);
  MRS_plan_update_from_rotor_angle_in_rad(new_plan, rotor_angle_in_rad,
                                          plan->allow_fourth_rank);
  event->plan = new_plan;
}

/**
 * @brief Create a plan for every events within a dimension struct.
 *
 * @param dimension A pointer to the MRS_dimension.
 * @param scheme A pointer to the powder averaging scheme MRS_averaging_scheme.
 * @param count (int) The number of points along the dimension.
 * @param increment (double) The increment in Hz along the dimension.
 * @param coordinates_offset (double) The coordinates offset in Hz along the dimension.
 * @param n_events (int) The number of events within the dimension.
 * @param fraction A pointer to the fraction/weight of the events along the dimension.
 * @param rotor_frequency_in_Hz A pointer to rotor frequency in Hz per event.
 * @param rotor_angle_in_rad A pointer to rotor angle in rads per event.
 * @param magnetic_flux_density_in_T A pointer to field fluc density in T per event.
 * @param number_of_sidebands The total number of requested sidebands.
 */
static inline void create_plans_for_events_in_dimension(
    MRS_dimension *dim, MRS_averaging_scheme *scheme, int count, double increment,
    double coordinates_offset, int n_events, double *fraction,
    double *rotor_frequency_in_Hz, double *rotor_angle_in_rad,
    double *magnetic_flux_density_in_T, unsigned int number_of_sidebands) {
  int i;
  dim->count = count;
  dim->coordinates_offset = coordinates_offset;
  dim->increment = increment;
  dim->n_events = n_events;
  dim->events = (MRS_event *)malloc(n_events * sizeof(MRS_event));

  MRS_plan *the_plan =
      MRS_create_plan(scheme, number_of_sidebands, *rotor_frequency_in_Hz,
                      *rotor_angle_in_rad, increment, scheme->allow_fourth_rank);

  for (i = 0; i < n_events; i++) {
    dim->events[i].freq_amplitude = malloc_double(the_plan->size);
    vm_double_ones(the_plan->size, dim->events[i].freq_amplitude);
    MRS_set_event(&(dim->events[i]), *fraction++, *magnetic_flux_density_in_T++,
                  *rotor_frequency_in_Hz++, *rotor_angle_in_rad++, increment, the_plan);
  }

  dim->inverse_increment = 1.0 / increment;
  dim->normalize_offset = 0.5 - (coordinates_offset * dim->inverse_increment);
  dim->R0_offset = 0.0;

  /* buffer to hold the local frequencies and frequency offset. The buffer is useful
   * when the rotor angle is off magic angle (54.735 deg). */
  dim->local_frequency = malloc_double(scheme->total_orientations);
  dim->freq_offset = malloc_double(scheme->octant_orientations);
}

/**
 * Creates dimension and event structs and generates plans for every event with a
 * dimension.
 **/
MRS_dimension *MRS_create_dimensions(
    MRS_averaging_scheme *scheme, int *count, double *coordinates_offset,
    double *increment, double *fractions, double *magnetic_flux_density_in_T,
    double *rotor_frequency_in_Hz, double *rotor_angle_in_rad, int *n_events,
    unsigned int n_dim, unsigned int number_of_sidebands) {
  unsigned int i;
  MRS_dimension *dimension = (MRS_dimension *)malloc(n_dim * sizeof(MRS_dimension));

  for (i = 0; i < n_dim; i++) {
    create_plans_for_events_in_dimension(
        &(dimension[i]), scheme, count[i], increment[i], coordinates_offset[i],
        n_events[i], fractions, rotor_frequency_in_Hz, rotor_angle_in_rad,
        magnetic_flux_density_in_T, number_of_sidebands);

    fractions += n_events[i];
    rotor_frequency_in_Hz += n_events[i];
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
    //          dimension[i].events[j].rotor_frequency_in_Hz);
    //   printf("\t\trotor angle %f rad\n",
    //          dimension[i].events[j].rotor_angle_in_rad);
    //   printf("\t\tmagnetic flux density %f T\n",
    //          dimension[i].events[j].magnetic_flux_density_in_T);
    //   printf("\n");
    // }
  }
  return dimension;
}
