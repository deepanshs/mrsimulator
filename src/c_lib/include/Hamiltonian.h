//
//  Hamiltonian.h
//
//  Created by Deepansh J. Srivastava, Apr 11, 2019
//  Copyright © 2019 Deepansh J. Srivastava. All rights reserved.
//  Contact email = srivastava.89@osu.edu, deepansh2012@gmail.com
//

#include "transition_function.h"
#include "mrsimulator.h"

/*
===============================================================================
            First order Nuclear shielding Hamiltonian in the PAS.
-------------------------------------------------------------------------------
The Hamiltonian includes the product of second rank tensor and the
spin transition functions.
*/
static inline void get_nuclear_shielding_hamiltonian_to_first_order(
    double *R0, complex128 *R2, double iso, double zeta, double eta,
    double *transition)
{
  // Spin transition contribution
  double transition_fn = p(transition[1], transition[0]);

  // Scaled R00
  R0[0] = iso * transition_fn;

  /*
  Scaled R2m containing the components of the quad second rank tensor in
  its principal axis frame.
  */
  double temp = -0.4082482905 * (zeta * eta) * transition_fn;
  self_cdadd(R2[0], temp);                 // R2-2
  self_cdadd(R2[1], 0.0);                  // R2-1
  self_cdadd(R2[2], zeta * transition_fn); // R2 0
  self_cdadd(R2[3], 0.0);                  // R2 1
  self_cdadd(R2[4], temp);                 // R2 2
}

/*
===============================================================================
                First order Quadrupolar Hamiltonian in the PAS.
-------------------------------------------------------------------------------
The Hamiltonian includes the product of second rank tensor and the
spin transition functions.
*/
static inline void get_quadrupole_hamiltonian_to_first_order(
    double *R0, complex128 *R2, double spin, double Cq, double eta,
    double *transition)
{
  // Spin transition contribution
  double transition_fn = d(transition[1], transition[0]);

  // Scaled R00
  R0[0] += 0.0;

  /* vq is the Quadrupolar coupling constant given as vq = 3*Cq/(2I(2I-1)),
  where `I` is the spin quantum number. */
  double vq = 3.0 * Cq;
  double denominator = 2.0 * spin * (2.0 * spin - 1.0);
  vq /= denominator;

  /* Scaled R2m containing the components of the quad second rank tensor in
  its principal axis frame. */
  double temp = -0.1666666667 * (vq * eta) * transition_fn;
  self_cdadd(R2[0], temp);                              // R2-2
  self_cdadd(R2[1], 0.0);                               // R2-1
  self_cdadd(R2[2], 0.4082482905 * vq * transition_fn); // R2 0
  self_cdadd(R2[3], 0.0);                               // R2 1
  self_cdadd(R2[4], temp);                              // R2 2
}

/*
===============================================================================
               Second order Quadrupolar Hamiltonian in the PAS.
-------------------------------------------------------------------------------
The Hamiltonian includes the product of second rank tensor and the
spin transition functions.
*/
static inline void get_quadrupole_hamiltonian_to_second_order(
    double *R0, complex128 *R2, complex128 *R4, double spin,
    double Cq, double eta, double *transition, double vo,
    int remove_second_order_quad_iso)
{
  // Spin transition contribution
  double c0, c2, c4;
  quad_ci(&c0, &c2, &c4, transition[1], transition[0], spin);

  /* vq is the Quadrupolar coupling constant given as vq = 3*Cq/(2I(2I-1)),
  where `I` is the spin quantum number. */
  double vq = 3.0 * Cq;
  double denominator = 2.0 * spin * (2.0 * spin - 1.0);
  vq /= denominator;

  double scale = vq * vq / vo;
  double eta2 = eta * eta;

  // Scaled R00
  if (remove_second_order_quad_iso == 0)
  {
    R0[0] += (eta2 * 0.33333333333 + 1.0) * 0.07453559925 * scale * c0;
  }

  /* Scaled R2m containing the components of the quad second rank tensor in
  its principal axis frame. */
  double temp = -eta * 0.07273929675 * scale * c2;
  double temp2 = 0.08908708064 * (eta2 * 0.33333333333 - 1.0) * scale * c2;
  self_cdadd(R2[0], temp);                                                      // R2-2
  self_cdadd(R2[1], 0.0);                                                       // R2-1
  self_cdadd(R2[2], temp2); // R2 0
  self_cdadd(R2[3], 0.0);                                                       // R2 1
  self_cdadd(R2[4], temp);                                                      // R2 2

  /* Scaled R4m containing the components of the quad second rank tensor in
  its principal axis frame. */
  temp = eta2 * 0.02777777778 * scale * c4;
  temp2 = -0.06299407883 * eta * scale * c4;
  double temp4 = 0.1195228609 * (eta2 * 0.05555555556 + 1.0) * scale * c4;
  self_cdadd(R4[0], temp);                                                     // R4-4
  self_cdadd(R4[2], temp2);                                                    // R4-2
  self_cdadd(R4[4], temp4); // R4 0
  self_cdadd(R4[6], temp2);                                                    // R4 2
  self_cdadd(R4[8], temp);                                                     // R4 4
}

/*
===============================================================================
       First order Weakly coupled Magnetic Dipole Hamiltonian in the PAS.
-------------------------------------------------------------------------------
The Hamiltonian includes the product of second rank tensor and the
spin transition functions in the weak coupling limit.
*/
static inline void get_weakly_coupled_direct_dipole_hamiltonian_to_first_order(
    double *R0, complex128 *R2, double D, double *transition)
{
  // Spin transition contribution
  double transition_fn = dIS(transition[0], transition[1], 0.5, 0.5);

  // Scaled R00
  R0[0] += 0.0;

  /* Scaled R2m containing the components of the magnetic dipole second rank
  tensor in its principal axis frame. */
  self_cdadd(R2[0], 0.0);                     // R2-2
  self_cdadd(R2[1], 0.0);                     // R2-1
  self_cdadd(R2[2], 2.0 * D * transition_fn); // R2 0
  self_cdadd(R2[3], 0.0);                     // R2 1
  self_cdadd(R2[4], 0.0);                     // R2 2
}
