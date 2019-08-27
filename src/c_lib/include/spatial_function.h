//
//  spatial_function.h
//
//  Created by Deepansh J. Srivastava, Aug 26, 2019
//  Copyright © 2019 Deepansh J. Srivastava. All rights reserved.
//  Contact email = srivastava.89@osu.edu, deepansh2012@gmail.com
//

#include "mrsimulator.h"

/**
 * @brief Contribution from the spatial part of nuclear shielding Hamiltonian to
 *        first-order in the principal axis system, (PAS), includes zeroth and
 *        second rank irreducible tensors, given as,
 * @f[\left. \rho_{0,0} = \sigma_\text{iso} \right\} \text{Rank-0},@f]
 * @f[ \left.
 *      \begin{array}{r l}
 *      \rho_{2,0} &= \zeta_\sigma, \\
 *      \rho_{2,\pm1} &= 0, \\
 *      \rho_{2,\pm2} &= -\frac{1}{\sqrt{6}} \eta_\sigma \zeta_\sigma, \\
 *      \end{array} \right\} \text{Rank-2},
 * @f]
 * where @f$\sigma_\text{iso}@f$ is the isotropic nuclear shielding, and,
 * @f$\zeta_\sigma@f$, @f$\eta_\sigma@f$ are the shielding anisotropy and
 * asymmetry from the symmetric second-rank irreducible nuclear shielding tensor
 * defined using the Haeberlen convention.
 *
 * @param rho_0 A pointer to the array of size 1 where the spatial contribution
 *      from the zeroth-rank irreducible tensor, in PAS, is stored.
 * @param rho_2 A pointer to the array of size 5 where the spatial contribution
 *      from the second-rank irreducible tensor, in PAS, is stored in order
 *      @f$\left[\rho_{2,m}\right]_{m=-2}^2@f$.
 * @param iso The isotropic nuclear shielding in Hz.
 * @param zeta The second rank shielding anisotropy in Hz.
 * @param eta The second rank shielding asymmetry, @f$\eta \in [0, 1]@f$.
 */
static inline void
spatial_tensors_from_1st_order_nuclear_shielding_Hamiltonian(
    double *restrict rho_0, complex128 *restrict rho_2, const double iso,
    const double zeta, const double eta) {

  // contribution from the zeroth rank
  *rho_0 = iso;

  // contribution from the shielding symmetric second rank
  vm_double_zeros(10, (double *)rho_2);
  double *rho_2_ = (double *)rho_2;

  double temp = -0.4082482905 * (zeta * eta);
  rho_2_[0] = temp; // R2-2 real
  rho_2_[4] = zeta; // R2 0 real
  rho_2_[8] = temp; // R2 2 real
}

/**
 * @brief Contribution from the spatial part of the first-order electric
 *      quadrupole Hamiltonian in the principal axis system, (PAS), includes the
 *      second rank irreducible tensor, given as,
 * @f[ \left.
 *      \begin{array}{rl}
 *      \rho_{2,0} &= \frac{1}{\sqrt{6}} \nu_q, \\
 *      \rho_{2,\pm1} &= 0, \\
 *      \rho_{2,\pm2} &= -\frac{1}{6} \eta_q \nu_q,
 *      \end{array} \right\} \text{Rank-2},
 * @f]
 * where @f$\nu_q = \frac{3 C_q}{2I(2I-1)}@f$ is the quadrupolar coupling
 * constant, and @f$\eta_q@f$ is the quadrupole asymmetry parameter defined
 * using the Haeberlen convension. Here, @f$I@f$ is the spin quantum number.
 *
 * @param rho_2 A pointer to the array of size 5 where the spatial contribution
 *      from the second-rank irreducible tensor, in PAS, is stored in order
 *      @f$\left[\rho_{2,m}\right]_{m=-2}^2@f$.
 * @param spin The spin quantum number.
 * @param Cq The quadrupolar coupling constant in Hz.
 * @param eta The quadrupolar asymmetry parameter, @f$\eta \in [0, 1]@f$.
 */
static inline void
spatial_tensors_from_1st_order_electric_quadrupole_Hamiltonian(
    complex128 *restrict rho_2, const double spin, const double Cq,
    const double eta) {

  /* vq is the Quadrupolar coupling constant given as vq = 3*Cq/(2I(2I-1)),
  where `I` is the spin quantum number. */
  double vq = 3.0 * Cq;
  double denominator = 2.0 * spin * (2.0 * spin - 1.0);
  vq /= denominator;

  // contribution from the symmetric second rank quadrupole tensor
  vm_double_zeros(10, (double *)rho_2);
  double *rho_2_ = (double *)rho_2;

  double temp = -0.1666666667 * (vq * eta);

  rho_2_[0] = temp;              // R2-2 real
  rho_2_[4] = 0.4082482905 * vq; // R2 0 real
  rho_2_[8] = temp;              // R2 2 real
}

/**
 * @brief Contribution from the spatial part of the second-order electric
 *      quadrupole Hamiltonian in the principal axis system, (PAS), includes the
 *      zeroth, second and fourth rank irreducible tensor, given as,
 * @f[\left.
 *    \rho_{0,0} = \frac{\nu_q^2}{\nu_0} \frac{1}{6\sqrt{5}}
 *                    \left(\frac{\eta_q^2}{3} + 1 \right)
 *    \right\} \text{Rank-0},@f]
 * @f[ \left.
 *      \begin{align}
 *      \rho_{2,0} &= \frac{\nu_q^2}{\nu_0} \frac{\sqrt{2}}{6\sqrt{7}}
 *                    \left(\frac{\eta_q^2}{3} - 1 \right), \\
 *      \rho_{2,\pm1} &= 0, \\
 *      \rho_{2,\pm2} &= -\frac{\nu_q^2}{\nu_0} \frac{1}{3\sqrt{21}} \eta_q, \\
 *      \end{align} \right\} \text{Rank-2},
 * @f]
 * @f[ \left.
 *      \begin{align}
 *      \rho_{4,0} &= \frac{\nu_q^2}{\nu_0} \frac{1}{\sqrt{70}}
 *                    \left(\frac{\eta_q^2}{18} + 1 \right), \\
 *      \rho_{4,\pm1} &= 0, \\
 *      \rho_{4,\pm2} &= -\frac{\nu_q^2}{\nu_0} \frac{1}{6\sqrt{7}} \eta_q, \\
 *      \rho_{4,\pm3} &= 0, \\
 *      \rho_{4,\pm4} &= \frac{\nu_q^2}{\nu_0} \frac{1}{36} \eta_q^2,
 *      \end{align} \right\} \text{Rank-4},
 * @f]
 * where @f$\nu_q = \frac{3 C_q}{2I(2I-1)}@f$ is the quadrupolar coupling
 * constant, @f$\nu_0@f$ is the Larmor frequency, @f$\eta_q@f$ is the
 * quadrupole asymmetry parameter defined using the Haeberlen convension, and
 * @f$I@f$ is the spin quantum number.
 *
 * @param rho_0 A pointer to the array of size 1 where the spatial contribution
 *      from the zeroth-rank irreducible tensor, in PAS, is stored.
 * @param rho_2 A pointer to the array of size 5 where the spatial contributions
 *      from the second-rank irreducible tensor, in PAS, is stored in order
 *      @f$\left[\rho_{2,m}\right]_{m=-2}^2@f$.
 * @param rho_4 A pointer to the array of size 9 where the spatial contributions
 *      from the fourth-rank irreducible tensor, in PAS, is stored in order
 *      @f$\left[\rho_{4,m}\right]_{m=-4}^4@f$.
 * @param spin The spin quantum number.
 * @param Cq The quadrupolar coupling constant in Hz.
 * @param eta The quadrupolar asymmetry parameter, @f$\eta \in [0, 1]@f$.
 * @param vo The Larmor frequency in Hz.
 */
static inline void
spatial_tensors_from_2nd_order_electric_quadrupole_Hamiltonian(
    double *restrict rho_0, complex128 *restrict rho_2,
    complex128 *restrict rho_4, const double spin, const double Cq,
    const double eta, const double vo) {

  /* vq is the Quadrupolar coupling constant given as vq = 3*Cq/(2I(2I-1)),
  where `I` is the spin quantum number. */
  double vq = 3.0 * Cq;
  double denominator = 2.0 * spin * (2.0 * spin - 1.0);
  vq /= denominator;

  double scale = vq * vq / vo;
  double eta2 = eta * eta;

  // contribution from the zeroth rank
  *rho_0 = (eta2 * 0.33333333333 + 1.0) * 0.07453559925 * scale;

  // contribution from the second rank
  vm_double_zeros(10, (double *)rho_2);
  double *rho_2_ = (double *)rho_2;

  double temp = -eta * 0.07273929675 * scale;
  double temp2 = 0.08908708064 * scale * (eta2 * 0.33333333333 - 1.0);

  rho_2_[0] = temp;  // R2-2 real
  rho_2_[4] = temp2; // R2 0 real
  rho_2_[8] = temp;  // R2 2 real

  // contribution from the fourth rank
  vm_double_zeros(18, (double *)rho_4);
  double *rho_4_ = (double *)rho_4;

  temp = eta2 * 0.02777777778 * scale;
  temp2 = -0.06299407883 * eta * scale;
  double temp4 = 0.1195228609 * scale * (eta2 * 0.05555555556 + 1.0);

  rho_4_[0] = temp;   // R4-4 real
  rho_4_[4] = temp2;  // R4-2 real
  rho_4_[8] = temp4;  // R4 0 real
  rho_4_[12] = temp2; // R4 2 real
  rho_4_[16] = temp;  // R4 4 real
}
