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
  if (the_event->plan != NULL) {
    if (DEBUG) printf("inside event->plan\n");
    MRS_free_plan(the_event->plan);
    the_event->plan = NULL;
  }
  free(the_event->event_freq_amplitude);
  if (DEBUG) printf("freed event->event_freq_amplitude\n");
}

/** Free the memory/buffer allocation for the MRS dimensions and events within. **/
void MRS_free_dimension(MRS_dimension *dimensions, unsigned int n) {
  unsigned int dim, evt;
  for (dim = 0; dim < n; dim++) {
    if (DEBUG) printf("post execution dimension %d cleanup\n", dim);
    for (evt = 0; evt < dimensions[dim].n_events; evt++) {
      if (DEBUG) printf("%d event\n", evt);
      MRS_free_event(&dimensions[dim].events[evt]);
    }
    free(dimensions[dim].events);
    free(dimensions[dim].local_frequency);
    free(dimensions[dim].local_phase);
    free(dimensions[dim].freq_offset);
    free(dimensions[dim].freq_amplitude);
    if (DEBUG) printf("freed events, local_frequency, freq_offset, freq_amplitude\n");
  }
  free(dimensions);
}

/**
 * @brief Populates event structs and creates/updates new plans based on the plan from
 * the first event along the dimension.
 *
 * @param event The pointer to the event.
 * @param fraction The fraction/weight of the event.
 * @param duration The duration of a delay event in µs.
 * @param is_spectral True (1) if the event is a SpectralEvent, False (0) if a
 * DelayEvent
 * @param magnetic_flux_density_in_T The external field flux density at the event.
 * @param rotor_frequency_in_Hz The rotor frequency at the event.
 * @param rotor_angle_in_rad The rotor angle at the event.
 * @param inverse_increment The inverse increment of the corresponding dimension.
 * @param plan The MRS plan of the first event along the dimension.
 */
static inline void MRS_set_event(MRS_event *event, double fraction, double duration,
                                 unsigned char is_spectral,
                                 double magnetic_flux_density_in_T,
                                 double rotor_frequency_in_Hz,
                                 double rotor_angle_in_rad, double inverse_increment,
                                 MRS_plan *plan) {
  event->fraction = fraction;
  event->duration = duration;
  event->is_spectral = is_spectral;
  event->rotor_frequency_in_Hz = rotor_frequency_in_Hz;
  event->rotor_angle_in_rad = rotor_angle_in_rad;
  event->magnetic_flux_density_in_T = magnetic_flux_density_in_T;

  bool rotor_frequency_equal = rotor_frequency_in_Hz == plan->rotor_frequency_in_Hz;
  bool rotor_angle_equal = rotor_angle_in_rad == plan->rotor_angle_in_rad;

  /* When both rotor angle and rotor freq is the same as the plan, return plan. */
  if (rotor_frequency_equal && rotor_angle_equal) {
    MRS_plan *new_plan = MRS_copy_plan(plan);
    new_plan->copy = true;
    event->plan = new_plan;
    return;
  }

  /* When only rotor freq is different, update the plan accordingly and return. */
  if (!rotor_frequency_equal && rotor_angle_equal) {
    MRS_plan *new_plan = MRS_copy_plan(plan);
    new_plan->copy = true;
    new_plan->copy_for_rotor_freq = true;
    MRS_plan_update_from_rotor_frequency_in_Hz(new_plan, rotor_frequency_in_Hz);

    // normalize the sideband frequencies.
    cblas_dscal(new_plan->number_of_sidebands, inverse_increment, new_plan->vr_freq, 1);

    event->plan = new_plan;
    MRS_plan_release_temp_storage(new_plan);
    return;
  }

  /* When only rotor angle is different, update the plan accordingly and return. */
  if (rotor_frequency_equal && !rotor_angle_equal) {
    MRS_plan *new_plan = MRS_copy_plan(plan);
    new_plan->copy = true;
    new_plan->copy_for_rotor_angle = true;
    MRS_plan_update_from_rotor_angle_in_rad(new_plan, rotor_angle_in_rad,
                                            plan->allow_4th_rank);
    event->plan = new_plan;
    return;
  }

  /* Otherwise, update plan for both rotor angle and freq and return. */
  MRS_plan *new_plan = MRS_copy_plan(plan);
  new_plan->copy = true;
  new_plan->copy_for_rotor_freq = true;
  new_plan->copy_for_rotor_angle = true;
  MRS_plan_update_from_rotor_frequency_in_Hz(new_plan, rotor_frequency_in_Hz);
  MRS_plan_update_from_rotor_angle_in_rad(new_plan, rotor_angle_in_rad,
                                          plan->allow_4th_rank);

  // normalize the sideband frequencies.
  cblas_dscal(new_plan->number_of_sidebands, inverse_increment, new_plan->vr_freq, 1);

  event->plan = new_plan;
  MRS_plan_release_temp_storage(new_plan);
}

/**
 * @brief Create a plan for every events within a dimension struct.
 *
 * @param dim A pointer to the MRS_dimension.
 * @param scheme A pointer to the powder averaging scheme MRS_averaging_scheme.
 * @param count (int) The number of points along the dimension.
 * @param increment (double) The increment in Hz along the dimension.
 * @param coordinates_offset (double) The coordinates offset in Hz along the dimension.
 * @param n_events (int) The number of events within the dimension.
 * @param fraction A pointer to the fraction/weight of the events along the dimension.
 * @param duration A pointer to the duration of delay events along the dimension.
 * @param is_spectral A pointer to a unsigned char array describing which event comes
 * from
 * @param rotor_frequency_in_Hz A pointer to rotor frequency in Hz per event.
 * @param rotor_angle_in_rad A pointer to rotor angle in rads per event.
 * @param magnetic_flux_density_in_T A pointer to field flux density in T per event.
 * @param number_of_sidebands The total number of requested sidebands.
 */
static inline void create_plans_for_events_in_dimension(
    MRS_dimension *dim, MRS_averaging_scheme *scheme, int count, double increment,
    double coordinates_offset, int n_events, double *fraction, double *duration,
    unsigned char *is_spectral, double *rotor_frequency_in_Hz,
    double *rotor_angle_in_rad, double *magnetic_flux_density_in_T,
    unsigned int number_of_sidebands) {
  int i;
  dim->count = count;
  dim->coordinates_offset = coordinates_offset;
  dim->increment = increment;
  dim->n_events = n_events;
  dim->events = (MRS_event *)malloc(n_events * sizeof(MRS_event));

  dim->inverse_increment = 1.0 / increment;
  dim->normalize_offset = 0.5 - (coordinates_offset * dim->inverse_increment);
  dim->R0_offset = 0.0;

  MRS_plan *plan =
      MRS_create_plan(scheme, number_of_sidebands, *rotor_frequency_in_Hz,
                      *rotor_angle_in_rad, increment, scheme->allow_4th_rank);
  // normalize the sideband frequencies.
  cblas_dscal(number_of_sidebands, dim->inverse_increment, plan->vr_freq, 1);

  for (i = 0; i < n_events; i++) {
    dim->events[i].event_freq_amplitude = malloc_complex128(plan->size);
    vm_double_ones(plan->size * 2, (double *)dim->events[i].event_freq_amplitude);
    cblas_dscal(plan->size, 0.0, (double *)dim->events[i].event_freq_amplitude + 1, 2);

    // if (*rotor_frequency_in_Hz != 0.0 && *rotor_frequency_in_Hz != 1.0e12) {
    //   dim->events[i].freq_amplitude = malloc_double(plan->size);
    //   vm_double_ones(plan->size, dim->events[i].freq_amplitude);
    // }
    MRS_set_event(&(dim->events[i]), *fraction++, *duration++, *is_spectral++,
                  *magnetic_flux_density_in_T++, *rotor_frequency_in_Hz++,
                  *rotor_angle_in_rad++, dim->inverse_increment, plan);
    if (i == 0) dim->events->plan->copy = false;
  }

  if (DEBUG) printf("Early memory release\n");
  MRS_plan_release_temp_storage(plan);

  /* buffer to hold the local frequencies and frequency offset. The buffer is useful
   * when the rotor angle is off magic angle (54.735 deg). */
  dim->local_frequency = malloc_double(scheme->n_gamma * scheme->total_orientations);
  dim->local_phase = malloc_double(scheme->n_gamma * scheme->total_orientations);
  dim->freq_offset = malloc_double(scheme->octant_orientations);
  dim->freq_amplitude = malloc_double(plan->size);
}

/**
 * Creates dimension and event structs and generates plans for every event with a
 * dimension.
 **/
MRS_dimension *MRS_create_dimensions(
    MRS_averaging_scheme *scheme, int *count, double *coordinates_offset,
    double *increment, double *fractions, double *durations, unsigned char *is_spectral,
    double *magnetic_flux_density_in_T, double *rotor_frequency_in_Hz,
    double *rotor_angle_in_rad, int *n_events, unsigned int n_dim,
    unsigned int *number_of_sidebands) {
  unsigned int i;
  MRS_dimension *dimension = (MRS_dimension *)malloc(n_dim * sizeof(MRS_dimension));

  for (i = 0; i < n_dim; i++) {
    create_plans_for_events_in_dimension(
        &(dimension[i]), scheme, count[i], increment[i], coordinates_offset[i],
        n_events[i], fractions, durations, is_spectral, rotor_frequency_in_Hz,
        rotor_angle_in_rad, magnetic_flux_density_in_T, number_of_sidebands[i]);

    fractions += n_events[i];
    durations += n_events[i];
    is_spectral += n_events[i];
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
