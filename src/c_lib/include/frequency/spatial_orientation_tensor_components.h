//
//  spatial_orientation_tensors.h
//
//  @copyright Deepansh J. Srivastava, 2019-2021.
//  Created by Deepansh J. Srivastava, Aug 26, 2019.
//  Contact email = srivastava.89@osu.edu
//

// The file defines function for the frequency-scaled spatial spherical tensor (fsSST)
// components associated with first and second-order perturbation expansion of the
// nuclear spin interations Hamiltonian.

#ifndef __sptial_tensor_component__
#define __sptial_tensor_component__

#include "mrsimulator.h"

// =====================================================================================
//                      First-order nuclear shielding fsSST
// =====================================================================================

/**
 * The fsSST from the first-order perturbation expansion of the nuclear shielding
 * Hamiltonian, in the principal axis system (PAS), include contributions from the
 * zeroth and second-rank irreducible tensors which follow,
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
 * @param Delta_0 A pointer to an array of length 1, where the components of the
 * zeroth-rank irreducible tensor, @f$\mathcal{R'}_{0,0}^{(\sigma)}(\Theta)/2\pi@f$, is
 * stored.
 *
 * @param Delta_2 A pointer to a complex array of length 5, where the components of the
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
static inline void fsSST_1st_order_nuclear_shielding_tensor_components(
    double *restrict Delta_0, void *restrict Delta_2,
    const double omega_0_delta_iso_in_Hz, const double omega_0_zeta_sigma_in_Hz,
    const double eta, const double *Theta) {
  // contribution from the zeroth-rank.
  *Delta_0 = omega_0_delta_iso_in_Hz;  // ς0 0 real = -ν_0 σ = ν_0 𝛿

  // contribution from the shielding symmetric second-rank.
  vm_double_zeros(10, (double *)Delta_2);
  double *varsigma_2 = (double *)Delta_2;

  double temp = 0.4082482905 * (omega_0_zeta_sigma_in_Hz * eta);
  varsigma_2[0] = temp;                       // ς2-2 real = (η ζ ν_0) / sqrt(6)
  varsigma_2[4] = -omega_0_zeta_sigma_in_Hz;  // ς2 0 real = (- ζ ν_0)
  varsigma_2[8] = temp;                       // ς2 2 real = (η ζ ν_0) / sqrt(6)

  // wigner rotations from PAS to common frame
  if (Theta[0] == 0.0 && Theta[1] == 0.0 && Theta[2] == 0.0) {
    return;
  }
  single_wigner_rotation(2, Theta, Delta_2, Delta_2);
}

// =====================================================================================
//                     First-order electric quadrupolar fsSST
// =====================================================================================

/**
 * The fsSST from the first-order perturbation expansion of the electric quadrupole
 * Hamiltonian, in the principal axis system (PAS), include contributions from the
 * second-rank irreducible tensor which follow,
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
 * @param Delta_2 A pointer to a complex array of length 5, where the components of the
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
static inline void fsSST_1st_order_electric_quadrupole_tensor_components(
    void *restrict Delta_2, const double spin, const double Cq_in_Hz, const double eta,
    const double *Theta) {
  // ς_2q = [-1/6 v_q η, 0, 1/√6 v_q, 0, -1/6 v_q η]
  //      = [-1/2 ζ_q η, 0, √3/2 ζ_q, 0, -1/2 ζ_q η] * (v_q/3ζ_q)

  /* vq is the Quadrupole coupling constant given as vq = 3*Cq/(2I(2I-1)), where `I` is
   * the spin quantum number. */
  double vq = 3.0 * Cq_in_Hz;
  double denominator = 2.0 * spin * (2.0 * spin - 1.0);
  vq /= denominator;

  // contribution from the symmetric second-rank quadrupole tensor.
  vm_double_zeros(10, (double *)Delta_2);
  double *varsigma_2 = (double *)Delta_2;

  double temp = -0.1666666667 * (vq * eta);

  varsigma_2[0] = temp;               // ς2-2 real = (-η ν_q) / 6
  varsigma_2[4] = 0.4082482905 * vq;  // ς2 0 real = ν_q / sqrt(6)
  varsigma_2[8] = temp;               // ς2 2 real = (-η ν_q) / 6

  // wigner rotations from PAS to common frame
  if (Theta[0] == 0.0 && Theta[1] == 0.0 && Theta[2] == 0.0) {
    return;
  }
  single_wigner_rotation(2, Theta, Delta_2, Delta_2);
}

// =====================================================================================
//                      Second-order electric quadrupolar fsSST
// =====================================================================================

/**
 * The fsSST from the second-order perturbation expansion of the electric quadrupole
 * Hamiltonian, in the principal axis system (PAS), include contributions from the
 * zeroth, second, and fourth-rank irreducible tensors which follow,
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
 * @param Delta_0 A pointer to an array of length 1 where the zeroth-rank irreducible
 *    tensor, @f$\mathcal{R'}_{0,0}^{(qq)}(\Theta)/2\pi@f$, will be stored.
 *
 * @param Delta_2 A pointer to a complex array of length 5 where the second-rank
 * irreducible tensor, @f$\mathcal{R'}_{2,n}^{(qq)}(\Theta)/2\pi@f$, will be stored
 * ordered according to
 * @f$\left[\mathcal{R'}_{2,n}^{(qq)}(\Theta)/2\pi\right]_{n=-2}^2@f$.
 *
 * @param Delta_4 A pointer to a complex array of length 9 where the fourth-rank
 * irreducible tensor, @f$\mathcal{R'}_{4,n}^{(qq)}(\Theta)/2\pi@f$, will be stored
 * ordered according to
 * @f$\left[\mathcal{R'}_{4,n}^{(qq)}(\Theta)/2\pi\right]_{n=-4}^4@f$.
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
static inline void fsSST_2nd_order_electric_quadrupole_tensor_components(
    double *restrict Delta_0, void *restrict Delta_2, void *restrict Delta_4,
    const double spin, const double v0_in_Hz, const double Cq_in_Hz, const double eta,
    const double *Theta) {
  /* vq is the Quadrupole coupling constant given as vq = 3*Cq/(2I(2I-1)), where `I` is
   * the spin quantum number. */
  double vq = 3.0 * Cq_in_Hz;
  double denominator = 2.0 * spin * (2.0 * spin - 1.0);
  vq /= denominator;

  double scale = vq * vq / v0_in_Hz;
  double eta2 = eta * eta;

  // contribution from the zeroth-rank.
  // ς0 0 = (η^2 / 3 + 1) (ν_q^2/ν_0) / (6 sqrt(5))
  *Delta_0 = (eta2 * 0.33333333333 + 1.0) * 0.07453559925 * scale;  // ς0 0 real

  // contribution from the second-rank.
  vm_double_zeros(10, (double *)Delta_2);
  double *varsigma_2 = (double *)Delta_2;

  double temp = -eta * 0.07273929675 * scale;
  double temp2 = 0.08908708064 * scale * (eta2 * 0.33333333333 - 1.0);

  varsigma_2[0] = temp;   // ς2-2 real = -(ν_q^2/ν_0) η / (3 sqrt(21))
  varsigma_2[4] = temp2;  // ς2 0 real = (η^2 / 3 - 1) (ν_q^2/ν_0) / (6 sqrt(7/2))
  varsigma_2[8] = temp;   // ς2 2 real = -(ν_q^2/ν_0) η / (3 sqrt(21))

  // contribution from the fourth-rank.
  vm_double_zeros(18, (double *)Delta_4);
  double *varsigma_4 = (double *)Delta_4;

  temp = eta2 * 0.02777777778 * scale;
  temp2 = -0.06299407883 * eta * scale;
  double temp4 = 0.1195228609 * scale * (eta2 * 0.05555555556 + 1.0);

  varsigma_4[0] = temp;    // ς4-4 real = (ν_q^2/ν_0) η^2 / 36
  varsigma_4[4] = temp2;   // ς4-2 real = -(ν_q^2/ν_0) η / (6 sqrt(7))
  varsigma_4[8] = temp4;   // ς4 0 real = (η^2 / 18 + 1) (ν_q^2/ν_0) / (sqrt(70))
  varsigma_4[12] = temp2;  // ς4 2 real = -(ν_q^2/ν_0) η / (6 sqrt(7))
  varsigma_4[16] = temp;   // ς4 4 real = (ν_q^2/ν_0) η^2 / 36

  // wigner rotations from PAS to common frame
  if (Theta[0] == 0.0 && Theta[1] == 0.0 && Theta[2] == 0.0) {
    return;
  }
  single_wigner_rotation(2, Theta, Delta_2, Delta_2);
  single_wigner_rotation(4, Theta, Delta_4, Delta_4);
}

// =====================================================================================
//                       Weakly coupled First-order J-coupling fsSST
// =====================================================================================

/**
 * The fsSST from the first-order perturbation expansion of the @f$J@f$-coupling
 * Hamiltonian under weak-coupling limit, in the principal axis system (PAS), include
 * contributions from the zeroth and second-rank irreducible tensors which follow,
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
 * @param Delta_0 A pointer to an array of length 1, where the zeroth-rank irreducible
 *    tensor, @f$\mathcal{R'}_{0,0}^{(J)}(\Theta)/2\pi@f$, is stored.
 *
 * @param Delta_2 A pointer to a complex array of length 5, where the second-rank
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
static inline void fsSST_1st_order_weakly_coupled_J_tensor_components(
    double *restrict Delta_0, void *restrict Delta_2, const double J_iso_in_Hz,
    const double J_aniso_in_Hz, const double J_eta, const double *Theta) {
  // contribution from the zeroth-rank.
  *Delta_0 = J_iso_in_Hz;  // ς0 0 real = J_iso

  // contribution from the shielding symmetric second-rank.
  vm_double_zeros(10, (double *)Delta_2);
  double *varsigma_2 = (double *)Delta_2;

  double temp = -0.4082482905 * (J_aniso_in_Hz * J_eta);
  varsigma_2[0] = temp;           // ς2-2 real = -ζ_J * η / sqrt(6)
  varsigma_2[4] = J_aniso_in_Hz;  // ς2 0 real = ζ_J
  varsigma_2[8] = temp;           // ς2 2 real = -ζ_J * η / sqrt(6)

  // wigner rotations from PAS to common frame
  if (Theta[0] == 0.0 && Theta[1] == 0.0 && Theta[2] == 0.0) {
    return;
  }
  single_wigner_rotation(2, Theta, Delta_2, Delta_2);
}

// =====================================================================================
//                    Weakly coupled First-order dipolar fsSST
// =====================================================================================

/**
 * The fsSST from the first-order perturbation expansion of the dipolar-coupling
 * Hamiltonian under weak-coupling limit, in the principal axis system (PAS), include
 * contributions from the second-rank irreducible tensors which follow,
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
 * @param Delta_2 A pointer to a complex array of length 5, where the second-rank
 *    irreducible tensor, @f$\mathcal{R'}_{2,n}^{(d)}(\Theta)/2\pi@f$, is stored ordered
 *    as @f$\left[\mathcal{R'}_{2,n}^{(d)}(\Theta)/2\pi\right]_{n=-2}^2@f$.
 *
 * @param D_in_Hz The dipolar coupling, @f$D@f$, in Hz.
 *
 * @param Theta A pointer to an array of Euler angles, in radians, of length 3, ordered
 *    as @f$[\alpha, \beta, \gamma]@f$.
 */
static inline void fsSST_1st_order_weakly_coupled_dipolar_tensor_components(
    void *restrict Delta_2, const double D_in_Hz, const double *Theta) {
  // contribution from the second-rank dipolar tensor.
  vm_double_zeros(10, (double *)Delta_2);
  double *varsigma_2 = (double *)Delta_2;

  varsigma_2[4] = 2 * D_in_Hz;  // ς2 0 real = 2 ν_d

  // wigner rotations from PAS to common frame
  if (Theta[0] == 0.0 && Theta[1] == 0.0 && Theta[2] == 0.0) {
    return;
  }
  single_wigner_rotation(2, Theta, Delta_2, Delta_2);
}

// =====================================================================================
//                          Cross-term tensor fsSST components
// =====================================================================================

static inline void tensor_prod_element(const int rank, const int start_a,
                                       const double *cg_factors, const double *vec_a,
                                       const double *vec_x, double *res) {
  int l = (2 * rank) + 1;
  int size = l - (int)(start_a / 2);
  double vec_temp[10];
  cblas_dcopy(10, vec_a, 1, vec_temp, 1);

  vm_double_multiply_inplace(l, cg_factors, 1, vec_temp, 2);        // scale real
  vm_double_multiply_inplace(l, cg_factors, 1, vec_temp + 1, 2);    // scale imag
  cblas_zdotu_sub(size, &vec_temp[start_a], 1, &vec_x[0], 1, res);  // complex res
}

static inline void rank_2_tensor_products(const double *Delta_2a,
                                          const double *Delta_2x, double *Delta_0,
                                          double *Delta_2, double *Delta_4) {
  double Delta_2a_rev[10];  // reverse order of elements in Delta_2a

  // reverse Delta_2a
  for (int i = 0; i < 10; i += 2) {
    Delta_2a_rev[8 - i] = Delta_2a[i];      // real
    Delta_2a_rev[9 - i] = Delta_2a[i + 1];  // imag
  }

  // zeroth-rank tensor components
  double CG_00[5] = {1.0, -1.0, 1.0, -1.0, 1.0};
  cblas_dscal(5, 0.4472135955, CG_00, 1);

  tensor_prod_element(2, 0, CG_00, Delta_2a_rev, Delta_2x, Delta_0);  // R0 0

  // second-rank tensor components
  double CG_2m2[5] = {0.0, 0.0, 0.5345224838, -0.6546536707, 0.5345224838};
  double CG_2m1[5] = {0.0, 0.6546536707, -0.2672612419, -0.2672612419, 0.6546536707};
  double CG_20[5] = {0.5345224838, 0.2672612419, -0.5345224838, 0.2672612419,
                     0.5345224838};

  tensor_prod_element(2, 4, CG_2m2, Delta_2a_rev, Delta_2x, &Delta_2[0]);  // R2 -2
  Delta_2[8] = Delta_2[0];                                                 // R2 2 real
  Delta_2[9] = -Delta_2[1];                                                // R2 2 imag

  tensor_prod_element(2, 2, CG_2m1, Delta_2a_rev, Delta_2x, &Delta_2[2]);  // R2 -1
  Delta_2[6] = -Delta_2[2];                                                // R2 1 real
  Delta_2[7] = Delta_2[3];                                                 // R2 1 imag

  tensor_prod_element(2, 0, CG_20, Delta_2a_rev, Delta_2x, &Delta_2[4]);  // R2 0

  // fourth-rank tensor components
  double CG_4m4[5] = {0.0, 0.0, 0.0, 0.0, 1.0};
  double CG_4m3[5] = {0.0, 0.0, 0.0, 0.7071067812, 0.7071067812};
  double CG_4m2[5] = {0.0, 0.0, 0.4629100499, 0.755928946, 0.4629100499};
  double CG_4m1[5] = {0.0, 0.2672612419, 0.6546536707, 0.6546536707, 0.2672612419};
  double CG_40[5] = {0.1195228609, 0.4780914437, 0.7171371656, 0.4780914437,
                     0.1195228609};

  tensor_prod_element(2, 8, CG_4m4, Delta_2a_rev, Delta_2x, &Delta_4[0]);  // R4 -4
  Delta_4[16] = Delta_4[0];                                                // R4 4 real
  Delta_4[17] = -Delta_4[1];                                               // R4 4 imag

  tensor_prod_element(2, 6, CG_4m3, Delta_2a_rev, Delta_2x, &Delta_4[2]);  // R4 -3
  Delta_4[14] = -Delta_4[2];                                               // R4 3 real
  Delta_4[15] = Delta_4[3];                                                // R4 3 imag

  tensor_prod_element(2, 4, CG_4m2, Delta_2a_rev, Delta_2x, &Delta_4[4]);  // R4 -2
  Delta_4[12] = Delta_4[4];                                                // R4 2 real
  Delta_4[13] = -Delta_4[5];                                               // R4 2 imag

  tensor_prod_element(2, 2, CG_4m1, Delta_2a_rev, Delta_2x, &Delta_4[6]);  // R4 -1
  Delta_4[10] = -Delta_4[6];                                               // R4 1 real
  Delta_4[11] = Delta_4[7];                                                // R4 1 imag

  tensor_prod_element(2, 0, CG_40, Delta_2a_rev, Delta_2x, &Delta_4[8]);  // R4 0
}

// =====================================================================================
//                             Generic cross-term fsSST
// =====================================================================================

static inline void fsSST_cross_tensor_components(double *restrict Delta_0,
                                                 void *restrict Delta_2,
                                                 void *restrict Delta_4,
                                                 const double *Delta_2a,
                                                 const double *Delta_2x) {
  // temp storage for Delta_0, Delta_0_tem, because rank_2_tensor_products requires a
  // complex Delta_0.
  double Delta_0_tem[2] = {0.0, 0.0};
  vm_double_zeros(10, (double *)Delta_2);
  vm_double_zeros(18, (double *)Delta_4);

  rank_2_tensor_products(Delta_2a, Delta_2x, Delta_0_tem, (double *)Delta_2,
                         (double *)Delta_4);
  Delta_0[0] = Delta_0_tem[0];  // extract the real part and store it in Delta_0,
}

#endif /* __sptial_tensor_component__ */
