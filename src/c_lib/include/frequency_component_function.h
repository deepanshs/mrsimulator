//
//  frequency_component_function.h
//
//  Created by Deepansh J. Srivastava, Apr 11, 2019
//  Copyright © 2019 Deepansh J. Srivastava. All rights reserved.
//  Contact email = srivastava.89@osu.edu, deepansh2012@gmail.com
//

#include "spatial_orientation_function.h"
#include "spin_transition_function.h"

/**
 * The frequency component functions, in the principal axis system, from the
 * first-order nuclear shielding Hamiltonian includes
 * @f[
 *    \Lambda_{0,0}^{(\sigma)}(i,j) &= \mathcal{R}_{0,0}^{(\sigma)}
 *                                    ~~  \mathbb{p}(i, j),~\text{and} \\
 *    \Lambda_{2,n}^{(\sigma)}(i,j) &= \mathcal{R}_{2,n}^{(\sigma)}
 *                                    ~~  \mathbb{p}(i, j).
 * @f]
 * where @f$\mathcal{R}_{L,n}^{(\sigma)}@f$ are the spatial orientation
 * functions, and @f$\mathbb{p}(i, j)@f$ is the spin transition
 * function for @f$i \rightarrow j@f$ transition.
 *
 * @param Lambda_0 A pointer to an array of length 1 where the frequency
 *      components from @f$\mathcal{\Lambda}_{0,0}^{(\sigma)}@f$ will be
 *      stored.
 * @param Lambda_2 A pointer to an array of length 5 where the frequency
 *      components from @f$\mathcal{\Lambda}_{2,n}^{(\sigma)}@f$ will be
 *      stored ordered according to
 *      @f$\left[\mathcal{\Lambda}_{2,n}^{(\sigma)}\right]_{n=-2}^2@f$.
 * @param iso_in_Hz The isotropic nuclear shielding given in Hz.
 * @param zeta_in_Hz The strength of the nuclear shielding anisotropy, in Hz,
 *      defined using Haeberlen convention.
 * @param eta The nuclear shielding asymmetry parameter,
 *      @f$\eta_\sigma \in [-1,1]@f$, defined using Haeberlen convention.
 * @param transition A pointer to an array of length 2 where the quantum
 *      numbers describing the two energy states involved in the transition
 *      @f$i \rightarrow j@f$ are stored ordered as @f$[i, j]@f$.
 */
static inline void
frequency_component_functions_from_1st_order_nuclear_shielding_Hamiltonian(
    double *restrict Lambda_0, void *restrict Lambda_2, const double iso_in_Hz,
    const double zeta_in_Hz, const double eta, const double *transition) {
  // Spin transition function
  double transition_fn = p(transition[1], transition[0]);

  // Spatial orientation function
  spatial_orientation_functions_from_1st_order_nuclear_shielding_Hamiltonian(
      Lambda_0, Lambda_2, iso_in_Hz, zeta_in_Hz, eta);

  // frequency component function from zeroth-rank irreducible tensor
  *Lambda_0 *= transition_fn;

  // frequency component function from second-rank irreducible tensor
  double *Lambda_2_ = (double *)Lambda_2;
  Lambda_2_[0] *= transition_fn; // Lambda_2-2 real
  Lambda_2_[4] *= transition_fn; // Lambda_2 0 real
  Lambda_2_[8] *= transition_fn; // Lambda_2 2 real
}

/**
 * The frequency component functions, in the principal axis system, from the
 * first-order electric quadrupole Hamiltonian include
 * @f[
 *    \Lambda_{2,n}^{(q)}(i,j) = \mathcal{R}_{2,n}^{(q)} ~~ \mathbb{d}(i, j).
 * @f]
 * where @f$\mathcal{R}_{L,n}^{(q)}@f$ are the spatial orientation
 * functions, and @f$\mathbb{d}(i, j)@f$ is the spin transition
 * function for @f$i \rightarrow j@f$ transition.
 *
 * @param Lambda_2 A pointer to an array of length 5 where the frequency
 *      components from @f$\mathcal{\Lambda}_{2,n}^{(q)}@f$ will be stored
 *      ordered according to
 *      @f$\left[\mathcal{\Lambda}_{2,n}^{(q)}\right]_{n=-2}^2@f$.
 * @param spin The spin quantum number, @f$I@f$.
 * @param Cq_in_Hz The quadrupole coupling constant, @f$C_q@f$, in Hz.
 * @param eta The quadrupole asymmetry parameter, @f$\eta_q \in [0, 1]@f$.
 * @param transition A pointer to an array of length 2 where the quantum
 *      numbers describing the two energy states involved in the transition
 *      @f$i \rightarrow j@f$ are stored ordered as @f$[i, j]@f$.
 */
static inline void
frequency_component_functions_from_1st_order_electric_quadrupole_Hamiltonian(
    void *restrict Lambda_2, const double spin, const double Cq_in_Hz,
    const double eta, const double *transition) {
  // Spin transition function
  double transition_fn = d(transition[1], transition[0]);

  // Spatial orientation function
  spatial_orientation_functions_from_1st_order_electric_quadrupole_Hamiltonian(
      Lambda_2, spin, Cq_in_Hz, eta);

  // frequency component function from second-rank irreducible tensor
  double *Lambda_2_ = (double *)Lambda_2;
  Lambda_2_[0] *= transition_fn; // Lambda_2-2 real
  Lambda_2_[4] *= transition_fn; // Lambda_2 0 real
  Lambda_2_[8] *= transition_fn; // Lambda_2 2 real
}

/**
 * The frequency component functions, in the principal axis system, from the
 * second-order electric quadrupole Hamiltonian includes
 * @f[
 *    \Lambda_{0,0}^{(qq)}(i,j) &= \mathcal{R}_{0,0}^{(qq)}
 *                                      ~~ \mathbb{c}_0(i, j), \\
 *    \Lambda_{2,n}^{(qq)}(i,j) &= \mathcal{R}_{2,n}^{(qq)}
 *                                      ~~ \mathbb{c}_2(i, j),~\text{and} \\
 *    \Lambda_{4,n}^{(qq)}(i,j) &= \mathcal{R}_{4,n}^{(qq)}
 *                                      ~~ \mathbb{c}_4(i, j),
 * @f]
 * where @f$\mathcal{R}_{L,n}^{(qq)}@f$ are the spatial orientation
 * functions, and @f$\mathbb{c}_i(i, j)@f$ are the composite spin transition
 * functions for @f$i \rightarrow j@f$ transition.
 *
 * @param Lambda_0 A pointer to an array of length 1 where the frequency
 *      components from @f$\mathcal{\Lambda}_{0,0}^{(qq)}@f$ will be stored.
 * @param Lambda_2 A pointer to an array of length 5 where the frequency
 *      components from @f$\mathcal{\Lambda}_{2,n}^{(qq)}@f$ will be stored
 *      ordered according to
 *      @f$\left[\mathcal{\Lambda}_{2,n}^{(qq)}\right]_{n=-2}^2@f$.
 * @param Lambda_4 A pointer to an array of length 5 where the frequency
 *      components from @f$\mathcal{\Lambda}_{4,n}^{(qq)}@f$ will be stored
 *      ordered according to
 *      @f$\left[\mathcal{\Lambda}_{4,n}^{(qq)}\right]_{n=-4}^4@f$.
 * @param spin The spin quantum number, @f$I@f$.
 * @param Cq_in_Hz The quadrupole coupling constant, @f$C_q@f$, in Hz.
 * @param eta The quadrupole asymmetry parameter, @f$\eta_q \in [0, 1]@f$.
 * @param v0_in_Hz The Larmor frequency, @f$\nu_0@f$, in Hz.
 * @param transition A pointer to an array of length 2 where the quantum
 *      numbers describing the two energy states involved in the transition
 *      @f$i \rightarrow j@f$ are stored ordered as @f$[i, j]@f$.
 */
static inline void
frequency_component_functions_from_2nd_order_electric_quadrupole_Hamiltonian(
    double *restrict Lambda_0, void *restrict Lambda_2, void *restrict Lambda_4,
    const double spin, const double Cq_in_Hz, const double eta,
    const double *transition, const double v0_in_Hz) {
  // Composite spin transition functions
  double *cl_value = malloc_double(3);
  cL(cl_value, transition[1], transition[0], spin);

  // Spatial orientation function
  spatial_orientation_functions_from_2nd_order_electric_quadrupole_Hamiltonian(
      Lambda_0, Lambda_2, Lambda_4, spin, Cq_in_Hz, eta, v0_in_Hz);

  // frequency component function from zeroth-rank irreducible tensor
  *Lambda_0 *= *cl_value++;

  // frequency component function from second-rank irreducible tensor
  double *Lambda_2_ = (double *)Lambda_2;
  Lambda_2_[0] *= *cl_value; // Lambda_2-2 real
  Lambda_2_[4] *= *cl_value; // Lambda_2 0 real
  Lambda_2_[8] *= *cl_value; // Lambda_2 2 real

  cl_value++;

  // frequency component function from fourth-rank irreducible tensor
  double *Lambda_4_ = (double *)Lambda_4;
  Lambda_4_[0] *= *cl_value;  // Lambda_4-4 real
  Lambda_4_[4] *= *cl_value;  // Lambda_4-2 real
  Lambda_4_[8] *= *cl_value;  // Lambda_4 0 real
  Lambda_4_[12] *= *cl_value; // Lambda_4 2 real
  Lambda_4_[16] *= *cl_value; // Lambda_4 4 real
}

/*
===============================================================================
       First order Weakly coupled Magnetic Dipole frequency in the PAS.
-------------------------------------------------------------------------------
The frequency includes the product of second rank tensor and the
spin transition functions in the weak coupling limit.
*/
static inline void weakly_coupled_direct_dipole_frequencies_to_first_order(
    double *restrict Lambda_0, void *restrict Lambda_2, const double D,
    const double *transition) {
  // Spin transition contribution
  double transition_fn = dIS(transition[0], transition[1], 0.5, 0.5);

  // Scaled R00
  *Lambda_0 += 0.0;

  /* Scaled R2m containing the components of the magnetic dipole second rank
  tensor in its principal axis frame. */
  vm_double_zeros(10, (double *)Lambda_2);
  double *Lambda_2_ = (double *)Lambda_2;
  Lambda_2_[4] = 2.0 * D * transition_fn; // Lambda_2 0 real
}
