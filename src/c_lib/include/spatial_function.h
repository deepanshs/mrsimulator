//
//  spatial_function.h
//
//  Created by Deepansh J. Srivastava, Aug 26, 2019
//  Copyright © 2019 Deepansh J. Srivastava. All rights reserved.
//  Contact email = srivastava.89@osu.edu, deepansh2012@gmail.com
//

#include "mrsimulator.h"

/**
 * Contribution from the spatial part of the first-order nuclear shielding
 * Hamiltonian in the principal axis system, (PAS), includes the zeroth and
 * second-rank irreducible nuclear shielding tensors, given as
 * @f[ \left.
 *        \mathcal{R}_{0,0} = \nu_0\sigma_\text{iso}
 *      \right\} \text{Rank-0},
 * @f]
 * @f[ \left.
 *      \begin{array}{r l}
 *      \mathcal{R}_{2,0} &= \nu_0\zeta_\sigma, \\
 *      \mathcal{R}_{2,\pm1} &= 0, \\
 *      \mathcal{R}_{2,\pm2} &= -\frac{1}{\sqrt{6}}\nu_0 \eta_\sigma
 *                              \zeta_\sigma,
 *      \end{array}
 *     \right\} \text{Rank-2},
 * @f]
 * where @f$\sigma_\text{iso}@f$ is the isotropic nuclear shielding, and,
 * @f$\zeta_\sigma@f$, @f$\eta_\sigma@f$ are the shielding anisotropy and
 * asymmetry parameters from the symmetric second-rank irreducible nuclear
 * shielding tensor defined using Haeberlen convention. Here,
 * @f$\nu_0 = \gamma_I B_0@f$ is the Larmor frequency where, @f$\gamma_I@f$
 * and @f$B_0@f$ are the gyromagnetic ratio of the nucleus and the magnetic
 * flux density of the external magnetic field, respectively.
 *
 * @param R_0 A pointer to an array of length 1 where the spatial contribution
 *      from the zeroth-rank irreducible nuclear shielding tensor is stored in
 *      the principal axis system.
 * @param R_2 A pointer to an array of length 5 where the spatial contribution
 *      from the second-rank irreducible nuclear shielding tensor is stored in
 *      the principal axis system ordered according to
 *      @f$\left[\mathcal{R}_{2,m}\right]_{m=-2}^2@f$.
 * @param iso_in_Hz The product, @f$\nu_0\sigma_\text{iso}@f$, given in Hz.
 * @param zeta_in_Hz The product, @f$\nu_0\zeta_\sigma@f$, given in Hz.
 * @param eta The second rank shielding asymmetry, @f$\eta \in [0, 1]@f$.
 */
static inline void spatial_tensors_from_1st_order_nuclear_shielding_Hamiltonian(
    double *restrict R_0, void *restrict R_2, const double iso_in_Hz,
    const double zeta_in_Hz, const double eta) {

  // contribution from the zeroth rank
  *R_0 = iso_in_Hz;

  // contribution from the shielding symmetric second rank
  vm_double_zeros(10, (double *)R_2);
  double *R_2_ = (double *)R_2;

  double temp = -0.4082482905 * (zeta_in_Hz * eta);
  R_2_[0] = temp; // R2-2 real
  R_2_[4] = zeta_in_Hz; // R2 0 real
  R_2_[8] = temp; // R2 2 real
}

/**
 * Contribution from the spatial part of the first-order electric quadrupole
 * Hamiltonian in the principal axis system, (PAS), includes the second-rank
 * irreducible electric quadrupole tensor, given as
 * @f[ \left.
 *      \begin{array}{rl}
 *      \mathcal{R}_{2,0}^{(q)} &= \frac{1}{\sqrt{6}} \nu_q, \\
 *      \mathcal{R}_{2,\pm1}^{(q)} &= 0, \\
 *      \mathcal{R}_{2,\pm2}^{(q)} &= -\frac{1}{6} \eta_q \nu_q,
 *      \end{array}
 *     \right\} \text{Rank-2},
 * @f]
 * where @f$\nu_q = \frac{3 C_q}{2I(2I-1)}@f$ is the quadrupolar splitting
 * frequency, and @f$\eta_q@f$ is the quadrupole asymmetry parameter. Here,
 * @f$I@f$ is the spin quantum number, and @f$C_q@f$ is the quadrupole coupling
 * constant.
 *
 * @param R_2 A pointer to an array of length 5 where the spatial contribution
 *      from the second-rank irreducible electric quadrupole tensor is stored
 *      in the principal axis system ordered according to
 *      @f$\left[\mathcal{R}_{2,m}\right]_{m=-2}^2@f$.
 * @param spin The spin quantum number, @f$I@f$.
 * @param Cq_in_Hz The quadrupolar coupling constant, @f$C_q@f$, in Hz.
 * @param eta The quadrupolar asymmetry parameter, @f$\eta \in [0, 1]@f$.
 */
static inline void
spatial_tensors_from_1st_order_electric_quadrupole_Hamiltonian(
    void *restrict R_2, const double spin, const double Cq_in_Hz, const double eta) {

  /* vq is the Quadrupolar coupling constant given as vq = 3*Cq/(2I(2I-1)),
  where `I` is the spin quantum number. */
  double vq = 3.0 * Cq_in_Hz;
  double denominator = 2.0 * spin * (2.0 * spin - 1.0);
  vq /= denominator;

  // contribution from the symmetric second rank quadrupole tensor
  vm_double_zeros(10, (double *)R_2);
  double *R_2_ = (double *)R_2;

  double temp = -0.1666666667 * (vq * eta);

  R_2_[0] = temp;              // R2-2 real
  R_2_[4] = 0.4082482905 * vq; // R2 0 real
  R_2_[8] = temp;              // R2 2 real
}

/**
 * Contribution from the spatial part of the second-order electric quadrupole
 * Hamiltonian in the principal axis system, (PAS), includes the zeroth,
 * second and fourth rank irreducible electric quadrupole tensors, given as
 * @f[\left.
 *      \mathcal{R}_{0,0}^{(q)} = \frac{\nu_q^2}{\nu_0} \frac{1}{6\sqrt{5}}
 *                                \left(\frac{\eta_q^2}{3} + 1 \right)
 *    \right\} \text{Rank-0},
 * @f]
 * @f[ \left.
 *      \begin{align}
 *        \mathcal{R}_{2,0}^{(q)} &= \frac{\nu_q^2}{\nu_0}
 *                                   \frac{\sqrt{2}}{6\sqrt{7}}
 *                                   \left(\frac{\eta_q^2}{3} - 1 \right), \\
 *        \mathcal{R}_{2,\pm1}^{(q)} &= 0, \\
 *        \mathcal{R}_{2,\pm2}^{(q)} &= -\frac{\nu_q^2}{\nu_0}
 *                                       \frac{1}{3\sqrt{21}} \eta_q,
 *      \end{align}
 *     \right\} \text{Rank-2},
 * @f]
 * @f[ \left.
 *      \begin{align}
 *        \mathcal{R}_{4,0}^{(q)} &= \frac{\nu_q^2}{\nu_0} \frac{1}{\sqrt{70}}
 *                                   \left(\frac{\eta_q^2}{18} + 1 \right), \\
 *        \mathcal{R}_{4,\pm1}^{(q)} &= 0, \\
 *        \mathcal{R}_{4,\pm2}^{(q)} &= -\frac{\nu_q^2}{\nu_0}
 *                                       \frac{1}{6\sqrt{7}} \eta_q, \\
 *        \mathcal{R}_{4,\pm3}^{(q)} &= 0, \\
 *        \mathcal{R}_{4,\pm4}^{(q)} &= \frac{\nu_q^2}{\nu_0} \frac{1}{36}
 *                                      \eta_q^2,
 *      \end{align}
 *     \right\} \text{Rank-4},
 * @f]
 * where @f$\nu_q = \frac{3 C_q}{2I(2I-1)}@f$ is the quadrupolar splitting
 * frequency, @f$\nu_0@f$ is the Larmor frequency, and @f$\eta_q@f$ is the
 * quadrupole asymmetry parameter. Here,
 * @f$I@f$ is the spin quantum number, and @f$C_q@f$ is the quadrupole coupling
 * constant.
 *
 * @param R_0 A pointer to an array of length 1 where the spatial contribution
 *      from the zeroth-rank irreducible electric quadrupole tensor is stored
 *      in the principal axis system.
 * @param R_2 A pointer to an array of length 5 where the spatial contribution
 *      from the second-rank irreducible electric quadrupole tensor is stored
 *      in the principal axis system ordered according to
 *      @f$\left[\mathcal{R}_{2,m}\right]_{m=-2}^2@f$.
 * @param R_4 A pointer to an array of length 9 where the spatial contribution
 *      from the fourth-rank irreducible electric quadrupole tensor is stored
 *      in the principal axis system ordered according to
 *      @f$\left[\mathcal{R}_{4,m}\right]_{m=-4}^4@f$.
 * @param spin The spin quantum number, @f$I@f$.
 * @param Cq_in_Hz The quadrupolar coupling constant, @f$C_q@f$, in Hz.
 * @param eta The quadrupolar asymmetry parameter, @f$\eta \in [0, 1]@f$.
 * @param v0_in_Hz The Larmor frequency, @f$\nu_0@f$, in Hz.
 */
static inline void
spatial_tensors_from_2nd_order_electric_quadrupole_Hamiltonian(
    double *restrict R_0, void *restrict R_2, void *restrict R_4,
    const double spin, const double Cq_in_Hz, const double eta, const double v0_in_Hz) {

  /* vq is the Quadrupolar coupling constant given as vq = 3*Cq/(2I(2I-1)),
  where `I` is the spin quantum number. */
  double vq = 3.0 * Cq_in_Hz;
  double denominator = 2.0 * spin * (2.0 * spin - 1.0);
  vq /= denominator;

  double scale = vq * vq / v0_in_Hz;
  double eta2 = eta * eta;

  // contribution from the zeroth rank
  *R_0 = (eta2 * 0.33333333333 + 1.0) * 0.07453559925 * scale;

  // contribution from the second rank
  vm_double_zeros(10, (double *)R_2);
  double *R_2_ = (double *)R_2;

  double temp = -eta * 0.07273929675 * scale;
  double temp2 = 0.08908708064 * scale * (eta2 * 0.33333333333 - 1.0);

  R_2_[0] = temp;  // R2-2 real
  R_2_[4] = temp2; // R2 0 real
  R_2_[8] = temp;  // R2 2 real

  // contribution from the fourth rank
  vm_double_zeros(18, (double *)R_4);
  double *R_4_ = (double *)R_4;

  temp = eta2 * 0.02777777778 * scale;
  temp2 = -0.06299407883 * eta * scale;
  double temp4 = 0.1195228609 * scale * (eta2 * 0.05555555556 + 1.0);

  R_4_[0] = temp;   // R4-4 real
  R_4_[4] = temp2;  // R4-2 real
  R_4_[8] = temp4;  // R4 0 real
  R_4_[12] = temp2; // R4 2 real
  R_4_[16] = temp;  // R4 4 real
}
