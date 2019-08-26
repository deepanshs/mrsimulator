//
//  frequencies.h
//
//  Created by Deepansh J. Srivastava, Apr 11, 2019
//  Copyright © 2019 Deepansh J. Srivastava. All rights reserved.
//  Contact email = srivastava.89@osu.edu, deepansh2012@gmail.com
//

#include "mrsimulator.h"
#include "transition_function.h"

/*
===============================================================================
            First order Nuclear shielding frequency in the PAS.
-------------------------------------------------------------------------------
The frequency includes the product of second rank tensor and the
spin transition functions.
*/
static inline void nuclear_shielding_frequency_to_1st_order(
    double *restrict R0, complex128 *restrict R2, const double iso,
    const double zeta, const double eta, const double *transition) {
  // Spin transition contribution
  double transition_fn = p(transition[1], transition[0]);

  // Scaled R00
  *R0 += iso * transition_fn;

  /*
  Scaled R2m containing the components of the quad second rank tensor in
  its principal axis frame.
  */
  double temp = -0.4082482905 * (zeta * eta) * transition_fn;

  vm_double_zeros(10, (double *)R2);
  double *R2_ = (double *)R2;

  R2_[0] = temp;                 // R2-2 real
  R2_[4] = zeta * transition_fn; // R2 0 real
  R2_[8] = temp;                 // R2 2 real
}

/*
===============================================================================
                First order Quadrupolar frequency in the PAS.
-------------------------------------------------------------------------------
The frequency includes the product of second rank tensor and the
spin transition functions.
*/
static inline void quadrupole_frequencies_to_first_order(
    double *restrict R0, complex128 *restrict R2, const double spin,
    const double Cq, const double eta, const double *transition) {
  // Spin transition contribution
  double transition_fn = d(transition[1], transition[0]);

  // Scaled R00
  *R0 += 0.0;

  /* vq is the Quadrupolar coupling constant given as vq = 3*Cq/(2I(2I-1)),
  where `I` is the spin quantum number. */
  double vq = 3.0 * Cq;
  double denominator = 2.0 * spin * (2.0 * spin - 1.0);
  vq /= denominator;

  /* Scaled R2m containing the components of the quad second rank tensor in
  its principal axis frame. */
  double temp = -0.1666666667 * (vq * eta) * transition_fn;

  vm_double_zeros(10, (double *)R2);
  double *R2_ = (double *)R2;

  R2_[0] = temp;                              // R2-2 real
  R2_[4] = 0.4082482905 * vq * transition_fn; // R2 0 real
  R2_[8] = temp;                              // R2 2 real
}

/*
===============================================================================
               Second order Quadrupolar frequency in the PAS.
-------------------------------------------------------------------------------
The frequency includes the product of second rank tensor and the
spin transition functions.
*/
static inline void quadrupole_frequencies_to_second_order(
    double *restrict R0, complex128 *restrict R2, complex128 *restrict R4,
    const double spin, const double Cq, const double eta,
    const double *transition, const double vo,
    const int remove_second_order_quad_isotropic) {
  // Spin transition contribution
  double *cl_value = malloc_double(3);
  cL(cl_value, transition[1], transition[0], spin);

  /* vq is the Quadrupolar coupling constant given as vq = 3*Cq/(2I(2I-1)),
  where `I` is the spin quantum number. */
  double vq = 3.0 * Cq;
  double denominator = 2.0 * spin * (2.0 * spin - 1.0);
  vq /= denominator;

  double scale = vq * vq / vo;
  double eta2 = eta * eta;

  // Scaled R00
  if (remove_second_order_quad_isotropic == 0) {
    *R0 += (eta2 * 0.33333333333 + 1.0) * 0.07453559925 * scale * *cl_value++;
  }

  /* Scaled R2m containing the components of the quad second rank tensor in
  its principal axis frame. */
  double temp = -eta * 0.07273929675 * scale * *cl_value;
  double temp2 = 0.08908708064 * scale * *cl_value;
  temp2 *= (eta2 * 0.33333333333 - 1.0);

  vm_double_zeros(10, (double *)R2);
  double *R2_ = (double *)R2;

  R2_[0] = temp;  // R2-2 real
  R2_[4] = temp2; // R2 0 real
  R2_[8] = temp;  // R2 2 real

  cl_value++;
  /* Scaled R4m containing the components of the quad second rank tensor in
  its principal axis frame. */

  temp = eta2 * 0.02777777778 * scale * *cl_value;
  temp2 = -0.06299407883 * eta * scale * *cl_value;
  double temp4 = 0.1195228609 * scale * *cl_value;
  temp4 *= (eta2 * 0.05555555556 + 1.0);

  vm_double_zeros(18, (double *)R4);
  double *R4_ = (double *)R4;

  R4_[0] = temp;   // R4-4 real
  R4_[4] = temp2;  // R4-2 real
  R4_[8] = temp4;  // R4 0 real
  R4_[12] = temp2; // R4 2 real
  R4_[16] = temp;  // R4 4 real
}

/*
===============================================================================
       First order Weakly coupled Magnetic Dipole frequency in the PAS.
-------------------------------------------------------------------------------
The frequency includes the product of second rank tensor and the
spin transition functions in the weak coupling limit.
*/
static inline void weakly_coupled_direct_dipole_frequencies_to_first_order(
    double *restrict R0, complex128 *restrict R2, const double D,
    const double *transition) {
  // Spin transition contribution
  double transition_fn = dIS(transition[0], transition[1], 0.5, 0.5);

  // Scaled R00
  *R0 += 0.0;

  /* Scaled R2m containing the components of the magnetic dipole second rank
  tensor in its principal axis frame. */
  vm_double_zeros(10, (double *)R2);
  double *R2_ = (double *)R2;
  R2_[4] = 2.0 * D * transition_fn; // R2 0 real
}
