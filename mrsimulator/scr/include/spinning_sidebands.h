//
//  spinning_sidebands.h
//
//  Created by Deepansh J. Srivastava, Apr 11, 2019
//  Copyright © 2019 Deepansh J. Srivastava. All rights reserved.
//  Contact email = srivastava.89@osu.edu, deepansh2012@gmail.com
//

#include "mrsimulator.h"

// headerdoc

extern void __get_pre_phase_components(int number_of_sidebands,
                                       double spin_frequency,
                                       double complex *pre_phase);

extern void spinning_sideband_core(
    // spectrum information and related amplitude
    double *spec,              // The amplitude of the spectrum.
    double *cpu_time_,         // Execution time
    double spectral_start,     // The start of the frequency spectrum.
    double spectral_increment, // The increment of the frequency spectrum.
    int number_of_points,      // Number of points on the frequency spectrum.

    isotopomer_ravel *ravel_isotopomer, // isotopomer structure

    int quadSecondOrder,              // Quad theory for second order,
    int remove_second_order_quad_iso, // remove the isotropic contribution from
                                      // the second order quad Hamiltonian.

    // spin rate, spin angle and number spinning sidebands
    int number_of_sidebands,          // The number of sidebands
    double sample_rotation_frequency, // The rotor spin frequency
    double rotor_angle, // The rotor angle relative to lab-frame z-axis

    // Pointer to the transitions. transition[0] = mi and transition[1] = mf
    double *transition,

    // powder orientation average
    int geodesic_polyhedron_frequency // The number of triangle along the edge
                                      // of octahedron
);

// Return a vector ordered according to the fft output order.                //
// @params int n - The number of points                                      //
// @params double increment - The increment (sampling interval)              //
// @returns *double values = The pointer to the fft output order vector      //
static inline double *__get_frequency_in_FFT_order(int n, double increment) {
  double *vr_freq = createDouble1DArray(n);
  int i = 0, m, positive_limit, negative_limit;

  if (n % 2 == 0) {
    negative_limit = (int)(-n / 2);
    positive_limit = -negative_limit - 1;
  } else {
    negative_limit = (int)(-(n - 1) / 2);
    positive_limit = -negative_limit;
  }

  for (m = 0; m <= positive_limit; m++) {
    vr_freq[i] = (double)m * increment;
    i++;
  }
  for (m = negative_limit; m < 0; m++) {
    vr_freq[i] = (double)m * increment;
    i++;
  }
  return vr_freq;
}
