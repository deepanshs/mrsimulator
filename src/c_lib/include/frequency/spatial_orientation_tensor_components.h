//
//  spatial_orientation_tensors.h
//
//  @copyright Deepansh J. Srivastava, 2019-2021.
//  Created by Deepansh J. Srivastava, Aug 26, 2019.
//  Contact email = srivastava.89@osu.edu
//

// The file contains the spatial orientation tensor components associated with
// lower order perturbation expansion of the nuclear spin interations
// Hamiltonian.

#ifndef __sptial_tensor_component__
#define __sptial_tensor_component__

#include "mrsimulator.h"

// =====================================================================================
//              First-order nuclear shielding spatial orientation tensor
// =====================================================================================

/**
 * The scaled spatial orientation tensors (sSOT) from the first-order perturbation
 * expansion of the nuclear shielding Hamiltonian, in the principal axis system (PAS),
 * include contributions from the zeroth and second-rank irreducible tensors which
 * follow,
 * @f[ \left.
 *        \varsigma_{0,0}^{(\sigma)} = \omega_0\delta_\text{iso}
 *     \right\} \text{Rank-0},
 * @f]
 * @f[ \left.
 *      \begin{aligned}
 *      \varsigma_{2,0}^{(\sigma)} &= -\omega_0\zeta_\sigma, \\
 *      \varsigma_{2,\pm1}^{(\sigma)} &= 0, \\
 *      \varsigma_{2,\pm2}^{(\sigma)} &= \frac{1}{\sqrt{6}}\omega_0
 *                                        \eta_\sigma \zeta_\sigma,
 *      \end{aligned}
 *     \right\} \text{Rank-2},
 * @f]
 * where @f$\sigma_\text{iso}@f$ is the isotropic nuclear shielding, and
 * @f$\zeta_\sigma@f$, @f$\eta_\sigma@f$ are the shielding anisotropy and asymmetry
 * parameters from the symmetric second-rank irreducible nuclear shielding tensor
 * defined using Haeberlen convention. Here, @f$\omega_0 = -\gamma_I B_0@f$ is the
 * Larmor frequency where, @f$\gamma_I@f$ and @f$B_0@f$ are the gyromagnetic ratio of
 * the nucleus and the magnetic flux density of the external magnetic field,
 * respectively.
 *
 * For non-zero Euler angles, @f$\Theta = [\alpha, \beta, \gamma]@f$, Wigner rotation of
 * @f$\varsigma_{2,n}^{(\sigma)}@f$ is applied following,
 * @f[ \mathcal{R'}_{2,n}^{(\sigma)}(\Theta) =
 *                                \sum_{m = -2}^2 D^2_{m, n}(\Theta)
 *                                \varsigma_{2,n}^{(\sigma)},
 * @f]
 * where @f$\mathcal{R'}_{2,n}^{(\sigma)}(\Theta)@f$ are the tensor components in the
 * frame defined by the Euler angles @f$\Theta@f$.
 *
 * @note
 *  - The method accepts frequency physical quantities, that is,
 *    @f$\omega_0\delta_\text{iso}/2\pi@f$ and @f$\omega_0\zeta_\sigma/2\pi@f$, as the
 *    isotropic chemical shift and nuclear shielding anisotropy, respectively.
 *  - When @f$\Theta = [0,0,0]@f$, @f$\mathcal{R'}_{2,n}^{(\sigma)}(\Theta) =
 *    \varsigma_{2,n}^{(\sigma)}@f$ where @f$ n \in [-2,2]@f$.
 *  - @f$\mathcal{R'}_{0,0}^{(\sigma)}(\Theta) = \varsigma_{0,0}^{(\sigma)} ~~~ \forall
 *    ~ \Theta@f$.
 *  - The method returns @f$\mathcal{R'}_{0,0}^{(\sigma)}(\Theta)/2\pi@f$ and
 *    @f$\mathcal{R'}_{2,n}^{(\sigma)}(\Theta)/2\pi@f$, that is, in **units of
 * frequency**.
 *
 * @param R_0 A pointer to an array of length 1, where the components of the zeroth-rank
 *    irreducible tensor, @f$\mathcal{R'}_{0,0}^{(\sigma)}(\Theta)/2\pi@f$, is stored.
 *
 * @param R_2 A pointer to a complex array of length 5, where the components of the
 *    second-rank irreducible tensor, @f$\mathcal{R'}_{2,n}^{(\sigma)}(\Theta)/2\pi@f$,
 *    is stored ordered as
 *    @f$\left[\mathcal{R'}_{2,n}^{(\sigma)}(\Theta)/2\pi\right]_{n=-2}^2@f$.
 *
 * @param omega_0_delta_iso_in_Hz The isotropic chemical shift in Hz,
 *    @f$\omega_0\sigma_\text{iso}/2\pi@f$.
 *
 * @param omega_0_zeta_sigma_in_Hz The shielding anisotropy in Hz,
 *    @f$\omega_0\zeta_\sigma/2\pi@f$.
 *
 * @param eta The shielding asymmetry, @f$\eta_\sigma \in [0, 1]@f$.
 *
 * @param Theta A pointer to an array of Euler angles, in radians, of length 3, ordered
 *    as @f$[\alpha, \beta, \gamma]@f$.
 */
static inline void sSOT_1st_order_nuclear_shielding_tensor_components(
    double *restrict R_0, void *restrict R_2, const double omega_0_delta_iso_in_Hz,
    const double omega_0_zeta_sigma_in_Hz, const double eta, const double *Theta) {
  // contribution from the zeroth-rank.
  *R_0 = omega_0_delta_iso_in_Hz;  // R0 0 real

  // contribution from the shielding symmetric second-rank.
  vm_double_zeros(10, (double *)R_2);
  double *R_2_ = (double *)R_2;

  double temp = 0.4082482905 * (omega_0_zeta_sigma_in_Hz * eta);
  R_2_[0] = temp;                       // R2-2 real
  R_2_[4] = -omega_0_zeta_sigma_in_Hz;  // R2 0 real
  R_2_[8] = temp;                       // R2 2 real

  // wigner rotations from PAS to common frame
  if (Theta[0] == 0.0 && Theta[1] == 0.0 && Theta[2] == 0.0) {
    return;
  }
  single_wigner_rotation(2, Theta, R_2, R_2);
}

// =====================================================================================
//             First-order electric quadrupolar spatial orientation tensor
// =====================================================================================

/**
 * The scaled spatial orientation tensors (sSOT) from the first-order perturbation
 * expansion of the electric quadrupole Hamiltonian, in the principal axis system (PAS),
 * include contributions from the second-rank irreducible tensor which follow,
 * @f[ \left.
 *      \begin{aligned}
 *      \varsigma_{2,0}^{(q)} &= \frac{1}{\sqrt{6}} \omega_q, \\
 *      \varsigma_{2,\pm1}^{(q)} &= 0, \\
 *      \varsigma_{2,\pm2}^{(q)} &= -\frac{1}{6} \eta_q \omega_q,
 *      \end{aligned}
 *     \right\} \text{Rank-2},
 * @f]
 * where @f$\omega_q = \frac{6\pi C_q}{2I(2I-1)}@f$ is the quadrupole splitting
 * frequency, and @f$\eta_q@f$ is the quadrupole asymmetry parameter. Here, @f$I@f$ is
 * the spin quantum number of the quadrupole nucleus, and @f$C_q@f$ is the quadrupole
 * coupling constant.
 *
 * As before, for non-zero Euler angles, @f$\Theta = [\alpha,\beta,\gamma]@f$, a Wigner
 * rotation of @f$\varsigma_{2,n}^{(q)}@f$ is applied following,
 * @f[ \mathcal{R'}_{2,n}^{(q)}(\Theta) = \sum_{m = -2}^2 D^2_{m, n}(\Theta)
 *                                            \varsigma_{2,n}^{(q)}.
 * @f]
 * where @f$\mathcal{R'}_{2,n}^{(q)}(\Theta)@f$ are the tensor components in the frame
 * defined by the Euler angles @f$\Theta@f$.
 *
 * @note
 *  - When @f$\Theta = [0,0,0]@f$, @f$\mathcal{R'}_{2,n}^{(q)}(\Theta) =
 *    \varsigma_{2,n}^{(q)}@f$ where @f$ n \in [-2,2]@f$.
 *  - The method returns @f$\mathcal{R'}_{2,0}^{(q)}(\Theta)/2\pi@f$, that is, in
 *    **units of frequency**.
 *
 * @param R_2 A pointer to a complex array of length 5, where the components of the
 *    second-rank irreducible tensor, @f$\mathcal{R'}_{2,n}^{(q)}(\Theta)/2\pi@f$, is
 *    stored ordered as
 *    @f$\left[\mathcal{R'}_{2,n}^{(q)}(\Theta)/2\pi\right]_{n=-2}^2@f$.
 *
 * @param spin The spin quantum number, @f$I@f$.
 *
 * @param Cq_in_Hz The quadrupole coupling constant, @f$C_q@f$, in Hz.
 *
 * @param eta The quadrupole asymmetry parameter, @f$\eta_q \in [0, 1]@f$.
 *
 * @param Theta A pointer to an array of Euler angles, in radians, of length 3, ordered
 *    as @f$[\alpha, \beta, \gamma]@f$.
 */
static inline void sSOT_1st_order_electric_quadrupole_tensor_components(
    void *restrict R_2, const double spin, const double Cq_in_Hz, const double eta,
    const double *Theta) {
  // R_2q = [-1/6 v_q η, 0, 1/√6 v_q, 0, -1/6 v_q η]
  //      = [-1/2 ζ_q η, 0, √3/2 ζ_q, 0, -1/2 ζ_q η] * (v_q/3ζ_q)

  /* vq is the Quadrupole coupling constant given as vq = 3*Cq/(2I(2I-1)), where `I` is
   * the spin quantum number. */
  double vq = 3.0 * Cq_in_Hz;
  double denominator = 2.0 * spin * (2.0 * spin - 1.0);
  vq /= denominator;

  // contribution from the symmetric second-rank quadrupole tensor.
  vm_double_zeros(10, (double *)R_2);
  double *R_2_ = (double *)R_2;

  double temp = -0.1666666667 * (vq * eta);

  R_2_[0] = temp;               // R2-2 real
  R_2_[4] = 0.4082482905 * vq;  // R2 0 real
  R_2_[8] = temp;               // R2 2 real

  // wigner rotations from PAS to common frame
  if (Theta[0] == 0.0 && Theta[1] == 0.0 && Theta[2] == 0.0) {
    return;
  }
  single_wigner_rotation(2, Theta, R_2, R_2);
}

// =====================================================================================
//              Second-order electric quadrupolar spatial orientation tensor
// =====================================================================================

/**
 * The scaled spatial orientation tensors (sSOT) from the second-order perturbation
 * expansion of the electric quadrupole Hamiltonian, in the principal axis system (PAS),
 * include contributions from the zeroth, second, and fourth-rank irreducible tensors
 * which follow,
 * @f[\left.
 *      \varsigma_{0,0}^{(qq)} = \frac{\omega_q^2}{\omega_0} \frac{1}{6\sqrt{5}}
 *                                \left(\frac{\eta_q^2}{3} + 1 \right)
 *    \right\} \text{Rank-0},
 * @f]
 * @f[ \left.
 *      \begin{aligned}
 *        \varsigma_{2,0}^{(qq)} &= \frac{\omega_q^2}{\omega_0}
 *                                   \frac{\sqrt{2}}{6\sqrt{7}}
 *                                   \left(\frac{\eta_q^2}{3} - 1 \right), \\
 *        \varsigma_{2,\pm1}^{(qq)} &= 0, \\
 *        \varsigma_{2,\pm2}^{(qq)} &= -\frac{\omega_q^2}{\omega_0}
 *                                      \frac{1}{3\sqrt{21}} \eta_q,
 *      \end{aligned}
 *     \right\} \text{Rank-2},
 * @f]
 * @f[ \left.
 *      \begin{aligned}
 *        \varsigma_{4,0}^{(qq)} &= \frac{\omega_q^2}{\omega_0}
 *              \frac{1}{\sqrt{70}} \left(\frac{\eta_q^2}{18} + 1 \right), \\
 *        \varsigma_{4,\pm1}^{(qq)} &= 0, \\
 *        \varsigma_{4,\pm2}^{(qq)} &= -\frac{\omega_q^2}{\omega_0}
 *                                       \frac{1}{6\sqrt{7}} \eta_q, \\
 *        \varsigma_{4,\pm3}^{(qq)} &= 0, \\
 *        \varsigma_{4,\pm4}^{(qq)} &= \frac{\omega_q^2}{\omega_0} \frac{1}{36}
 *                                      \eta_q^2,
 *      \end{aligned}
 *     \right\} \text{Rank-4},
 * @f]
 * where @f$\omega_q = \frac{6\pi C_q}{2I(2I-1)}@f$ is the quadrupole splitting
 * frequency, @f$\omega_0@f$ is the Larmor angular frequency, and @f$\eta_q@f$ is the
 * quadrupole asymmetry parameter. Here, @f$I@f$ is the spin quantum number, and
 * @f$C_q@f$ is the quadrupole coupling constant.
 *
 * For non-zero Euler angles, @f$\Theta = [\alpha,\beta,\gamma]@f$, Wigner rotation of
 * @f$\varsigma_{2,n}^{(qq)}@f$ and @f$\varsigma_{4,n}^{(qq)}@f$ are applied following,
 * @f[
 *    \mathcal{R'}_{2,n}^{(qq)}(\Theta) &= \sum_{m = -2}^2 D^2_{m, n}(\Theta)
 *                                     \varsigma_{2,n}^{(qq)}, \\
 *    \mathcal{R'}_{4,n}^{(qq)}(\Theta) &= \sum_{m = -4}^4 D^4_{m, n}(\Theta)
 *                                     \varsigma_{4,n}^{(qq)},
 * @f]
 * where @f$\mathcal{R'}_{2,n}^{(qq)}(\Theta)@f$ and
 * @f$\mathcal{R'}_{4,n}^{(qq)}(\Theta)@f$ are the second and fourth-rank tensor
 * components in the frame defined by the Euler angles @f$\Theta@f$, respectively.
 *
 * @note
 *  - When @f$\Theta = [0,0,0]@f$, @f$\mathcal{R'}_{2,n}^{(qq)}(\Theta) =
 *    \varsigma_{2,n}^{(qq)}@f$ where @f$ n \in [-2,2]@f$.
 *  - When @f$\Theta = [0,0,0]@f$, @f$\mathcal{R'}_{4,n}^{(qq)}(\Theta) =
 *    \varsigma_{4,n}^{(qq)}@f$ where @f$ n \in [-4,4]@f$.
 *  - @f$\mathcal{R'}_{0,0}^{(qq)}(\Theta) = \varsigma_{0,0}^{(qq)} ~~~ \forall ~
 *    \Theta@f$.
 *  - The method returns @f$\mathcal{R'}_{0,0}^{(qq)}(\Theta)/2\pi@f$,
 *    @f$\mathcal{R'}_{2,n}^{(qq)}(\Theta)/2\pi@f$, and
 *    @f$\mathcal{R'}_{4,n}^{(qq)}(\Theta)/2\pi@f$, that is, in **units of frequency**.
 *
 * @param R_0 A pointer to an array of length 1 where the zeroth-rank irreducible
 *    tensor, @f$\mathcal{R'}_{0,0}^{(qq)}(\Theta)/2\pi@f$, will be stored.
 *
 * @param R_2 A pointer to a complex array of length 5 where the second-rank irreducible
 *    tensor, @f$\mathcal{R'}_{2,n}^{(qq)}(\Theta)/2\pi@f$, will be stored ordered
 *    according to @f$\left[\mathcal{R'}_{2,n}^{(qq)}(\Theta)/2\pi\right]_{n=-2}^2@f$.
 *
 * @param R_4 A pointer to a complex array of length 9 where the fourth-rank irreducible
 *    tensor, @f$\mathcal{R'}_{4,n}^{(qq)}(\Theta)/2\pi@f$, will be stored ordered
 *    according to @f$\left[\mathcal{R'}_{4,n}^{(qq)}(\Theta)/2\pi\right]_{n=-4}^4@f$.
 *
 * @param spin The spin quantum number, @f$I@f$.
 *
 * @param v0_in_Hz The Larmor frequency, @f$\omega_0/2\pi@f$, in Hz.
 *
 * @param Cq_in_Hz The quadrupole coupling constant, @f$C_q@f$, in Hz.
 *
 * @param eta The quadrupole asymmetry parameter, @f$\eta_q \in [0, 1]@f$.
 *
 * @param Theta A pointer to an array of Euler angles, in radians, of length 3, ordered
 *    as @f$[\alpha, \beta, \gamma]@f$.
 */
static inline void sSOT_2nd_order_electric_quadrupole_tensor_components(
    double *restrict R_0, void *restrict R_2, void *restrict R_4, const double spin,
    const double v0_in_Hz, const double Cq_in_Hz, const double eta,
    const double *Theta) {
  /* vq is the Quadrupole coupling constant given as vq = 3*Cq/(2I(2I-1)), where `I` is
   * the spin quantum number. */
  double vq = 3.0 * Cq_in_Hz;
  double denominator = 2.0 * spin * (2.0 * spin - 1.0);
  vq /= denominator;

  double scale = vq * vq / v0_in_Hz;
  double eta2 = eta * eta;

  // contribution from the zeroth-rank.
  *R_0 = (eta2 * 0.33333333333 + 1.0) * 0.07453559925 * scale;  // R0 0 real

  // contribution from the second-rank.
  vm_double_zeros(10, (double *)R_2);
  double *R_2_ = (double *)R_2;

  double temp = -eta * 0.07273929675 * scale;
  double temp2 = 0.08908708064 * scale * (eta2 * 0.33333333333 - 1.0);

  R_2_[0] = temp;   // R2-2 real
  R_2_[4] = temp2;  // R2 0 real
  R_2_[8] = temp;   // R2 2 real

  // contribution from the fourth-rank.
  vm_double_zeros(18, (double *)R_4);
  double *R_4_ = (double *)R_4;

  temp = eta2 * 0.02777777778 * scale;
  temp2 = -0.06299407883 * eta * scale;
  double temp4 = 0.1195228609 * scale * (eta2 * 0.05555555556 + 1.0);

  R_4_[0] = temp;    // R4-4 real
  R_4_[4] = temp2;   // R4-2 real
  R_4_[8] = temp4;   // R4 0 real
  R_4_[12] = temp2;  // R4 2 real
  R_4_[16] = temp;   // R4 4 real

  // wigner rotations from PAS to common frame
  if (Theta[0] == 0.0 && Theta[1] == 0.0 && Theta[2] == 0.0) {
    return;
  }
  single_wigner_rotation(2, Theta, R_2, R_2);
  single_wigner_rotation(4, Theta, R_4, R_4);
}

// =====================================================================================
//      First-order J-coupling spatial orientation tensor for weakly coupled sites
// =====================================================================================

/**
 * The scaled spatial orientation tensors (sSOT) from the first-order perturbation
 * expansion of the @f$J@f$-coupling Hamiltonian under weak-coupling limit, in the
 * principal axis system (PAS), include contributions from the zeroth and second-rank
 * irreducible tensors which follow,
 * @f[ \left.
 *        \varsigma_{0,0}^{(J)} = 2\pi J_\text{iso}
 *      \right\} \text{Rank-0},
 * @f]
 * @f[ \left.
 *      \begin{aligned}
 *      \varsigma_{2,0}^{(J)} &= 2 \pi \zeta_J, \\
 *      \varsigma_{2,\pm1}^{(J)} &= 0, \\
 *      \varsigma_{2,\pm2}^{(J)} &= -2\pi \frac{1}{\sqrt{6}} \eta_J \zeta_J,
 *      \end{aligned}
 *     \right\} \text{Rank-2},
 * @f]
 * where @f$J_\text{iso}@f$ is the isotropic @f$J@f$-coupling, and @f$\zeta_J@f$,
 * @f$\eta_J@f$ are the @f$J@f$-coupling tensor anisotropy and asymmetry parameters from
 * the symmetric second-rank irreducible @f$J@f$ tensor, defined using Haeberlen
 * convention.
 *
 * For non-zero Euler angles, @f$\Theta = [\alpha, \beta, \gamma]@f$, Wigner rotation of
 * @f$\varsigma_{2,n}^{(J)}@f$ is applied following,
 * @f[ \mathcal{R'}_{2,n}^{(J)}(\Theta) =
 *                                \sum_{m = -2}^2 D^2_{m, n}(\Theta)
 *                                \varsigma_{2,n}^{(J)},
 * @f]
 * where @f$\mathcal{R'}_{2,n}^{(J)}(\Theta)@f$ are the tensor components in the frame
 * defined by the Euler angles, @f$\Theta@f$.
 *
 * @note
 *  - When @f$\Theta = [0,0,0]@f$, @f$\mathcal{R'}_{2,n}^{(J)}(\Theta) =
 *    \varsigma_{2,n}^{(J)}@f$ where @f$ n \in [-2,2]@f$.
 *  - @f$\mathcal{R'}_{0,0}^{(J)}(\Theta) = \varsigma_{0,0}^{(J)} ~~~ \forall ~
 *    \Theta@f$.
 *  - The method returns @f$\mathcal{R'}_{0,0}^{(J)}(\Theta)/2\pi@f$ and
 *    @f$\mathcal{R'}_{2,n}^{(J)}(\Theta)/2\pi@f$, that is, in **units of frequency**.
 *
 * @param R_0 A pointer to an array of length 1, where the zeroth-rank irreducible
 *    tensor, @f$\mathcal{R'}_{0,0}^{(J)}(\Theta)/2\pi@f$, is stored.
 *
 * @param R_2 A pointer to a complex array of length 5, where the second-rank
 *    irreducible tensor, @f$\mathcal{R'}_{2,n}^{(J)}(\Theta)/2\pi@f$, is stored ordered
 *    as @f$\left[\mathcal{R'}_{2,n}^{(J)}(\Theta)/2\pi\right]_{n=-2}^2@f$.
 *
 * @param J_iso_in_Hz The isotropic @f$J@f$-coupling, @f$J_\text{iso}@f$, in Hz.
 *
 * @param J_aniso_in_Hz The @f$J@f$-coupling anisotropy, @f$\zeta_J@f$, in Hz.
 *
 * @param J_eta The @f$J@f$-coupling asymmetry, @f$\eta_J \in [0, 1]@f$.
 *
 * @param Theta A pointer to an array of Euler angles, in radians, of length 3, ordered
 *    as @f$[\alpha, \beta, \gamma]@f$.
 */
static inline void sSOT_1st_order_weakly_coupled_J_tensor_components(
    double *restrict R_0, void *restrict R_2, const double J_iso_in_Hz,
    const double J_aniso_in_Hz, const double J_eta, const double *Theta) {
  // contribution from the zeroth-rank.
  *R_0 = J_iso_in_Hz;  // R0 0 real

  // contribution from the shielding symmetric second-rank.
  vm_double_zeros(10, (double *)R_2);
  double *R_2_ = (double *)R_2;

  double temp = -0.4082482905 * (J_aniso_in_Hz * J_eta);
  R_2_[0] = temp;           // R2-2 real
  R_2_[4] = J_aniso_in_Hz;  // R2 0 real
  R_2_[8] = temp;           // R2 2 real

  // wigner rotations from PAS to common frame
  if (Theta[0] == 0.0 && Theta[1] == 0.0 && Theta[2] == 0.0) {
    return;
  }
  single_wigner_rotation(2, Theta, R_2, R_2);
}

// =====================================================================================
//       First-order dipolar spatial orientation tensor for weakly coupled sites
// =====================================================================================

/**
 * The scaled spatial orientation tensors (sSOT) from the first-order perturbation
 * expansion of the dipolar-coupling Hamiltonian under weak-coupling limit, in the
 * principal axis system (PAS), include contributions from the second-rank irreducible
 * tensors which follow,
 * @f[ \left.
 *      \begin{aligned}
 *      \varsigma_{2,0}^{(d)} &= 4\pi D, \\
 *      \varsigma_{2,\pm1}^{(d)} &= 0, \\
 *      \varsigma_{2,\pm2}^{(d)} &= 0,
 *      \end{aligned}
 *     \right\} \text{Rank-2},
 * @f]
 * where @f$D@f$ is the dipolar-coupling.
 *
 * For non-zero Euler angles, @f$\Theta = [\alpha, \beta, \gamma]@f$, Wigner rotation of
 * @f$\varsigma_{2,n}^{(d)}@f$ is applied following,
 * @f[
 *    \mathcal{R'}_{2,n}^{(d)}(\Theta) = \sum_{m = -2}^2 D^2_{m, n}(\Theta)
 *                                \varsigma_{2,n}^{(d)},
 * @f]
 * where @f$\mathcal{R'}_{2,n}^{(d)}(\Theta)@f$ are the tensor components in the frame
 * defined by the Euler angles, @f$\Theta@f$.
 *
 * @note
 *  - When @f$\Theta = [0,0,0]@f$, @f$\mathcal{R'}_{2,n}^{(d)}(\Theta) =
 *    \varsigma_{2,n}^{(d)}@f$ where @f$ n \in [-2,2]@f$.
 *  - The method returns @f$\mathcal{R'}_{2,n}^{(d)}(\Theta)/2\pi@f$, that is, in
 *    **units of frequency**.
 *
 * @param R_2 A pointer to a complex array of length 5, where the second-rank
 *    irreducible tensor, @f$\mathcal{R'}_{2,n}^{(d)}(\Theta)/2\pi@f$, is stored ordered
 *    as @f$\left[\mathcal{R'}_{2,n}^{(d)}(\Theta)/2\pi\right]_{n=-2}^2@f$.
 *
 * @param D_in_Hz The dipolar coupling, @f$D@f$, in Hz.
 *
 * @param Theta A pointer to an array of Euler angles, in radians, of length 3, ordered
 *    as @f$[\alpha, \beta, \gamma]@f$.
 */
static inline void sSOT_1st_order_weakly_coupled_dipolar_tensor_components(
    void *restrict R_2, const double D_in_Hz, const double *Theta) {
  // contribution from the second-rank dipolar tensor.
  vm_double_zeros(10, (double *)R_2);
  double *R_2_ = (double *)R_2;

  R_2_[4] = 2 * D_in_Hz;  // R2 0 real

  // wigner rotations from PAS to common frame
  if (Theta[0] == 0.0 && Theta[1] == 0.0 && Theta[2] == 0.0) {
    return;
  }
  single_wigner_rotation(2, Theta, R_2, R_2);
}

// =====================================================================================
//                            Cross-term tensor components
// =====================================================================================

static inline void tensor_prod_element(const int rank, const int start_a,
                                       const double *cg_factors, const double *vec_a,
                                       const double *vec_b, double *res) {
  int l = (2 * rank) + 1;
  int size = l - (int)(start_a / 2);
  double vec_temp[10];
  cblas_dcopy(10, vec_a, 1, vec_temp, 1);

  vm_double_multiply_inplace(l, cg_factors, 1, vec_temp, 2);        // scale real
  vm_double_multiply_inplace(l, cg_factors, 1, vec_temp + 1, 2);    // scale imag
  cblas_zdotu_sub(size, &vec_temp[start_a], 1, &vec_b[0], 1, res);  // complex res
}

static inline void rank_2_tensor_products(const double *R_2a, const double *R_2b,
                                          double *R_0, double *R_2, double *R_4) {
  double R_2a_reverse[10];

  // reverse R_2a
  for (int i = 0; i < 10; i += 2) {
    R_2a_reverse[8 - i] = R_2a[i];      // real
    R_2a_reverse[9 - i] = R_2a[i + 1];  // imag
  }

  // zeroth-rank tensor components
  double CG_00[5] = {1.0, -1.0, 1.0, -1.0, 1.0};
  cblas_dscal(5, 0.4472135955, CG_00, 1);

  tensor_prod_element(2, 0, CG_00, R_2a_reverse, R_2b, R_0);  // R0 0

  // second-rank tensor components
  double CG_2m2[5] = {0.0, 0.0, 0.5345224838, -0.6546536707, 0.5345224838};
  double CG_2m1[5] = {0.0, 0.6546536707, -0.2672612419, -0.2672612419, 0.6546536707};
  double CG_20[5] = {0.5345224838, 0.2672612419, -0.5345224838, 0.2672612419,
                     0.5345224838};

  tensor_prod_element(2, 4, CG_2m2, R_2a_reverse, R_2b, &R_2[0]);  // R2 -2
  R_2[8] = R_2[0];                                                 // R2 2 real
  R_2[9] = -R_2[1];                                                // R2 2 imag

  tensor_prod_element(2, 2, CG_2m1, R_2a_reverse, R_2b, &R_2[2]);  // R2 -1
  R_2[6] = -R_2[2];                                                // R2 1 real
  R_2[7] = R_2[3];                                                 // R2 1 imag

  tensor_prod_element(2, 0, CG_20, R_2a_reverse, R_2b, &R_2[4]);  // R2 0

  // fourth-rank tensor components
  double CG_4m4[5] = {0.0, 0.0, 0.0, 0.0, 1.0};
  double CG_4m3[5] = {0.0, 0.0, 0.0, 0.7071067812, 0.7071067812};
  double CG_4m2[5] = {0.0, 0.0, 0.4629100499, 0.755928946, 0.4629100499};
  double CG_4m1[5] = {0.0, 0.2672612419, 0.6546536707, 0.6546536707, 0.2672612419};
  double CG_40[5] = {0.1195228609, 0.4780914437, 0.7171371656, 0.4780914437,
                     0.1195228609};

  tensor_prod_element(2, 8, CG_4m4, R_2a_reverse, R_2b, &R_4[0]);  // R4 -4
  R_4[16] = R_4[0];                                                // R4 4 real
  R_4[17] = -R_4[1];                                               // R4 4 imag

  tensor_prod_element(2, 6, CG_4m3, R_2a_reverse, R_2b, &R_4[2]);  // R4 -3
  R_4[14] = -R_4[2];                                               // R4 3 real
  R_4[15] = R_4[3];                                                // R4 3 imag

  tensor_prod_element(2, 4, CG_4m2, R_2a_reverse, R_2b, &R_4[4]);  // R4 -2
  R_4[12] = R_4[4];                                                // R4 2 real
  R_4[13] = -R_4[5];                                               // R4 2 imag

  tensor_prod_element(2, 2, CG_4m1, R_2a_reverse, R_2b, &R_4[6]);  // R4 -1
  R_4[10] = -R_4[6];                                               // R4 1 real
  R_4[11] = R_4[7];                                                // R4 1 imag

  tensor_prod_element(2, 0, CG_40, R_2a_reverse, R_2b, &R_4[8]);  // R4 0
}

// =====================================================================================
//                      Generic cross spatial orientation tensor
// =====================================================================================

static inline void sSOT_cross_tensor_components(double *restrict R_0,
                                                void *restrict R_2, void *restrict R_4,
                                                const double *R_2a,
                                                const double *R_2x) {
  // create a temp storage for R0, R0_res, because rank_2_tensor_products requires a
  // complex R0.
  double R0_res[2] = {0.0, 0.0};
  vm_double_zeros(10, (double *)R_2);
  vm_double_zeros(18, (double *)R_4);

  rank_2_tensor_products(R_2a, R_2x, R0_res, (double *)R_2, (double *)R_4);
  R_0[0] = R0_res[0];  // extract the real part and store it in R_0,
}

#endif /* __sptial_tensor_component__ */
