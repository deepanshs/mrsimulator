//
//  frequencies.h
//
//  Created by Deepansh J. Srivastava, Apr 11, 2019
//  Copyright © 2019 Deepansh J. Srivastava. All rights reserved.
//  Contact email = srivastava.89@osu.edu, deepansh2012@gmail.com
//

#include "mrsimulator.h"
#include "spatial_function.h"
#include "transition_function.h"
/*
===============================================================================
            First order Nuclear shielding frequency in the PAS.
-------------------------------------------------------------------------------
*/

/**
 * @brief Frequency components from the first order nuclear shielding
 *      Hamiltonian is given as the product of the spatial and spin-transition
 *      functions,
 * @f[
 *
 * @f]
 *
 * spatial_part_nuclear_shielding_to_1st_order(double *restrict, complex128
 * restrict, const double, const double, const double),
 * p(const double mf, const double mi)
 * @param R0
 * The frequency includes the product of second rank tensor and the
 * spin transition functions.
 */
static inline void
frequency_components_from_1st_order_nuclear_shielding_Hamiltonian(
    double *restrict R0, complex128 *restrict R2, const double iso,
    const double zeta, const double eta, const double *transition) {
  // Spin transition contribution
  double transition_fn = p(transition[1], transition[0]);

  spatial_tensors_from_1st_order_nuclear_shielding_Hamiltonian(R0, R2, iso,
                                                               zeta, eta);

  // frequency contribution from zeroth rank tensor
  *R0 *= transition_fn;

  // frequency contribution from symmetric second rank tensor
  double *R2_ = (double *)R2;
  R2_[0] *= transition_fn; // R2-2 real
  R2_[4] *= transition_fn; // R2 0 real
  R2_[8] *= transition_fn; // R2 2 real
}

/*
===============================================================================
                First order Quadrupolar frequency in the PAS.
-------------------------------------------------------------------------------
The frequency includes the product of second rank tensor and the
spin transition functions.
*/
static inline void
frequency_components_from_1st_order_electric_quadrupole_Hamiltonian(
    complex128 *restrict R2, const double spin, const double Cq,
    const double eta, const double *transition) {
  // Spin transition contribution
  double transition_fn = d(transition[1], transition[0]);

  spatial_tensors_from_1st_order_electric_quadrupole_Hamiltonian(R2, spin, Cq,
                                                                 eta);

  // frequency contribution from symmetric second rank tensor
  double *R2_ = (double *)R2;
  R2_[0] *= transition_fn; // R2-2 real
  R2_[4] *= transition_fn; // R2 0 real
  R2_[8] *= transition_fn; // R2 2 real
}

/*
===============================================================================
               Second order Quadrupolar frequency in the PAS.
-------------------------------------------------------------------------------
The frequency includes the product of second rank tensor and the
spin transition functions.
*/
static inline void
frequency_components_from_2nd_order_electric_quadrupole_Hamiltonian(
    double *restrict R0, complex128 *restrict R2, complex128 *restrict R4,
    const double spin, const double Cq, const double eta,
    const double *transition, const double vo) {
  // Spin transition contribution
  double *cl_value = malloc_double(3);
  cL(cl_value, transition[1], transition[0], spin);

  spatial_tensors_from_2nd_order_electric_quadrupole_Hamiltonian(
      R0, R2, R4, spin, Cq, eta, vo);

  // frequency contribution from zeroth rank tensor
  *R0 *= *cl_value++;

  // frequency contribution from second rank tensor
  double *R2_ = (double *)R2;
  R2_[0] *= *cl_value; // R2-2 real
  R2_[4] *= *cl_value; // R2 0 real
  R2_[8] *= *cl_value; // R2 2 real

  cl_value++;

  // frequency contribution from fourth rank tensor
  double *R4_ = (double *)R4;
  R4_[0] *= *cl_value;  // R4-4 real
  R4_[4] *= *cl_value;  // R4-2 real
  R4_[8] *= *cl_value;  // R4 0 real
  R4_[12] *= *cl_value; // R4 2 real
  R4_[16] *= *cl_value; // R4 4 real
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
