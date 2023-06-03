// -*- coding: utf-8 -*-
//
//  frequency_tensor.h
//
//  @copyright Deepansh J. Srivastava, 2019-2021.
//  Created by Deepansh J. Srivastava, Apr 11, 2019.
//  Contact email = srivastava.89@osu.edu
//

#include "frequency/spatial_orientation_tensor_components.h"
#include "frequency/spin_transition_functions.h"

// =====================================================================================
//               First-order nuclear shielding frequency tensor components
// =====================================================================================

/**
 * The frequency tensors (FT) components from the first-order perturbation expansion of
 * the nuclear shielding Hamiltonian, in a given frame, @f$\mathcal{F}@f$, described by
 * the Euler angles @f$\Theta = [\alpha, \beta, \gamma]@f$ are
 * @f[
 *    {\Lambda'}_{0,0}^{(\sigma)} &= \mathcal{R'}_{0,0}^{(\sigma)}(\Theta)
 *                                    ~~  \mathbb{p}(i, j),~\text{and} \\
 *    {\Lambda'}_{2,n}^{(\sigma)} &= \mathcal{R'}_{2,n}^{(\sigma)}(\Theta)
 *                                    ~~  \mathbb{p}(i, j),
 * @f]
 * where @f$\mathcal{R'}_{0,0}^{(\sigma)}(\Theta)@f$ and
 * @f$\mathcal{R'}_{2,n}^{(\sigma)}(\Theta)@f$ are the spatial orientation functions in
 * frame @f$\mathcal{F}@f$, and @f$\mathbb{p}(i, j)@f$ is the spin transition function
 * for @f$\left|i\right> \rightarrow \left|j\right>@f$ transition.
 *
 * @param Lambda_0 A pointer to an array of length 1, where the frequency component from
 *      @f${\Lambda'}_{0,0}^{(\sigma)}@f$ is stored.
 * @param Lambda_2 A pointer to a complex array of length 5, where the frequency
 *      components from @f${\Lambda'}_{2,n}^{(\sigma)}@f$ is stored ordered as
 *      @f$\left[{\Lambda'}_{2,n}^{(\sigma)}\right]_{n=-2}^2@f$.
 * @param omega_0_delta_iso_in_Hz The isotropic chemical shift in Hz
 *      (@f$\omega_0\delta_\text{iso}/2\pi@f$).
 * @param omega_0_zeta_sigma_in_Hz The shielding anisotropy quantity in Hz
 *      (@f$\omega_0\zeta_\sigma/2\pi@f$) defined using Haeberlen convention.
 * @param eta The shielding asymmetry, @f$\eta_\sigma \in [-1,1]@f$, defined using
 *      Haeberlen convention.
 * @param Theta A pointer to an array of Euler angles, in radians, of length 3, ordered
 *      as @f$[\alpha, \beta, \gamma]@f$.
 * @param mf The spin quantum number of the final energy state.
 * @param mi The spin quantum number of the initial energy state.
 */
static inline void FCF_1st_order_nuclear_shielding_tensor_components(
    double *restrict Lambda_0, void *restrict Lambda_2,
    const double omega_0_delta_iso_in_Hz, const double omega_0_zeta_sigma_in_Hz,
    const double eta, const double *Theta, const float mf, const float mi) {
  // Spin transition function
  double transition_fn = STF_p(mf, mi);

  // Return if the transition is zero
  if (transition_fn == 0.0) {
    // zero the R0 and R2 components before populating with shielding components
    *Lambda_0 = 0.0;
    vm_double_zeros(10, (double *)Lambda_2);
    return;
  }

  // Spatial orientation function
  sSOT_1st_order_nuclear_shielding_tensor_components(
      Lambda_0, Lambda_2, omega_0_delta_iso_in_Hz, omega_0_zeta_sigma_in_Hz, eta,
      Theta);

  // frequency component function from the zeroth-rank irreducible tensor.
  *Lambda_0 *= transition_fn;

  // frequency component function from the second-rank irreducible tensor.
  cblas_dscal(10, transition_fn, (double *)Lambda_2, 1);
}

// =====================================================================================
//             First-order electric quadrupolar frequency tensor components
// =====================================================================================

/**
 * The frequency tensor (FT) components from the first-order perturbation expansion of
 * electric quadrupole Hamiltonian, in a given frame, @f$\mathcal{F}@f$, described by
 * the Euler angles @f$\Theta = [\alpha, \beta, \gamma]@f$ are
 * @f[
 *    {\Lambda'}_{2,n}^{(q)} = \mathcal{R'}_{2,n}^{(q)}(\Theta) ~~ \mathbb{d}(i, j),
 * @f]
 * where @f$\mathcal{R}_{2,n}^{(q)}(\Theta)@f$ are the spatial orientation functions in
 * frame @f$\mathcal{F}@f$, and @f$\mathbb{d}(i, j)@f$ is the spin transition function
 * for @f$\left|i\right> \rightarrow \left|j\right>@f$ transition.
 *
 * @param Lambda_2 A pointer to a complex array of length 5, where the frequency
 *      components from @f${\Lambda'}_{2,n}^{(q)}@f$ is stored ordered as
 *      @f$\left[{\Lambda'}_{2,n}^{(q)}\right]_{n=-2}^2@f$.
 * @param spin The spin quantum number, @f$I@f$.
 * @param Cq_in_Hz The quadrupole coupling constant, @f$C_q@f$, in Hz.
 * @param eta The quadrupole asymmetry parameter, @f$\eta_q \in [0, 1]@f$.
 * @param Theta A pointer to an array of Euler angles, in radians, of length 3, ordered
 *      as @f$[\alpha, \beta, \gamma]@f$.
 * @param mf The spin quantum number of the final energy state.
 * @param mi The spin quantum number of the initial energy state.
 */
static inline void FCF_1st_order_electric_quadrupole_tensor_components(
    void *restrict Lambda_2, const double spin, const double Cq_in_Hz, const double eta,
    const double *Theta, const float mf, const float mi) {
  // Spin transition function
  double transition_fn = STF_d(mf, mi);

  // Return if the transition is zero
  if (transition_fn == 0.0) {
    // zero the R2 components before populating with quad components
    vm_double_zeros(10, (double *)Lambda_2);
    return;
  }

  // Spatial orientation function
  sSOT_1st_order_electric_quadrupole_tensor_components(Lambda_2, spin, Cq_in_Hz, eta,
                                                       Theta);

  // frequency component function from the second-rank irreducible tensor.
  cblas_dscal(10, transition_fn, (double *)Lambda_2, 1);
}

// =====================================================================================
//             Second-order electric quadrupolar frequency tensor components
// =====================================================================================

/**
 * The frequency tensor (FCF) components from the second-order perturbation expansion of
 * electric quadrupole Hamiltonian, in a given frame, @f$\mathcal{F}@f$, described by
 * the Euler angles @f$\Theta = [\alpha, \beta, \gamma]@f$, are
 * @f[
 *    {\Lambda'}_{0,0}^{(qq)} &= \mathcal{R'}_{0,0}^{(qq)}(\Theta)
 *                                      ~~ \mathbb{c}_0(i, j), \\
 *    {\Lambda'}_{2,n}^{(qq)} &= \mathcal{R'}_{2,n}^{(qq)}(\Theta)
 *                                      ~~ \mathbb{c}_2(i, j),~\text{and} \\
 *    {\Lambda'}_{4,n}^{(qq)} &= \mathcal{R'}_{4,n}^{(qq)}(\Theta)
 *                                      ~~ \mathbb{c}_4(i, j),
 * @f]
 * where @f$\mathcal{R'}_{0,0}^{(qq)}(\Theta)@f$,
 * @f$\mathcal{R'}_{2,n}^{(qq)}(\Theta)@f$, and, @f$\mathcal{R'}_{4,n}^{(qq)}(\Theta)@f$
 * are the spatial orientation functions in frame @f$\mathcal{F}@f$, and
 * @f$\mathbb{c}_k(i, j)@f$ are the composite spin transition functions for
 * @f$\left|i\right> \rightarrow \left|j\right>@f$ transition.
 *
 * @param Lambda_0 A pointer to an array of length 1, where the frequency
 *      component from @f${\Lambda'}_{0,0}^{(qq)}@f$ is stored.
 * @param Lambda_2 A pointer to a complex array of length 5, where the frequency
 *      components from @f$\Lambda_{2,n}^{(qq)}@f$ are stored ordered as
 *      @f$\left[{\Lambda'}_{2,n}^{(qq)}\right]_{n=-2}^2@f$.
 * @param Lambda_4 A pointer to a complex array of length 9, where the frequency
 *      components from @f${\Lambda'}_{4,n}^{(qq)}@f$ are stored ordered as
 *      @f$\left[{\Lambda'}_{4,n}^{(qq)}\right]_{n=-4}^4@f$.
 * @param spin The spin quantum number, @f$I@f$.
 * @param Cq_in_Hz The quadrupole coupling constant, @f$C_q@f$, in Hz.
 * @param eta The quadrupole asymmetry parameter, @f$\eta_q \in [0, 1]@f$.
 * @param larmor_freq_in_Hz The Larmor frequency, @f$\nu_0@f$, in Hz.
 * @param Theta A pointer to an array of Euler angles, in radians, of length 3,
 *      ordered as @f$[\alpha, \beta, \gamma]@f$.
 * @param mf The spin quantum number of the final energy state.
 * @param mi The spin quantum number of the initial energy state.
 */
static inline void FCF_2nd_order_electric_quadrupole_tensor_components(
    double *restrict Lambda_0, void *restrict Lambda_2, void *restrict Lambda_4,
    const double spin, const double larmor_freq_in_Hz, const double Cq_in_Hz,
    const double eta, const double *Theta, const float mf, const float mi) {
  // Composite spin transition functions
  double cl_value[3];
  STF_cL(cl_value, mf, mi, spin);

  // Spatial orientation function
  sSOT_2nd_order_electric_quadrupole_tensor_components(
      Lambda_0, Lambda_2, Lambda_4, spin, larmor_freq_in_Hz, Cq_in_Hz, eta, Theta);

  // frequency component function from the zeroth-rank irreducible tensor.
  *Lambda_0 *= cl_value[0];

  // frequency component function from the second-rank irreducible tensor.
  cblas_dscal(10, cl_value[1], (double *)Lambda_2, 1);

  // frequency component function from the fourth-rank irreducible tensor.
  cblas_dscal(18, cl_value[2], (double *)Lambda_4, 1);
}

// =====================================================================================
//                Shielding-Quad cross frequency tensor components
// =====================================================================================

/**
 * The frequency tensor (FCF) components from the nuclear shielding and electric
 * quadrupolar cross term, in the crystallite frame, @f$\mathcal{F}@f$, are
 * @f[
 *    {\Lambda'}_{0,0}^{(\sigma q)} &= \mathcal{R'}_{0,0}^{(\sigma q)}(\Theta)
 *                                      ~~ \mathbb{d}(i, j), \\
 *    {\Lambda'}_{2,n}^{(\sigma q)} &= \mathcal{R'}_{2,n}^{(\sigma q)}(\Theta)
 *                                      ~~ \mathbb{d}(i, j),~\text{and} \\
 *    {\Lambda'}_{4,n}^{(\sigma q)} &= \mathcal{R'}_{4,n}^{(\sigma q)}(\Theta)
 *                                      ~~ \mathbb{d}(i, j),
 * @f]
 * where @f$\mathcal{R'}_{0,0}^{(\sigma q)}(\Theta)@f$,
 * @f$\mathcal{R'}_{2,n}^{(\sigma q)}(\Theta)@f$, and,
 * @f$\mathcal{R'}_{4,n}^{(\sigma q)}(\Theta)@f$
 * are the spatial orientation functions in frame @f$\mathcal{F}@f$, and
 * @f$\mathbb{d}(i, j)@f$ is spin transition functions for
 * @f$\left|i\right> \rightarrow \left|j\right>@f$ transition.
 *
 * @param Lambda_0 A pointer to an array of length 1, where the frequency
 *      component from @f${\Lambda'}_{0,0}^{(\sigma q)}@f$ is stored.
 * @param Lambda_2 A pointer to a complex array of length 5, where the frequency
 *      components from @f$\Lambda_{2,n}^{(\sigma q)}@f$ are stored ordered as
 *      @f$\left[{\Lambda'}_{2,n}^{(\sigma  q)}\right]_{n=-2}^2@f$.
 * @param Lambda_4 A pointer to a complex array of length 9, where the frequency
 *      components from @f${\Lambda'}_{4,n}^{(\sigma q)}@f$ are stored ordered as
 *      @f$\left[{\Lambda'}_{4,n}^{(\sigma q)}\right]_{n=-4}^4@f$.
 * @param R_2q A pointer to a complex array of length 5 holding the quadupolar
 *      tensor components.
 * @param R_2s A pointer to a complex array of length 5 holding the nuclear shileding
 *      tensor components.
 * @param mf The spin quantum number of the final energy state.
 * @param mi The spin quantum number of the initial energy state.
 */
static inline void FCF_NS_EQ_cross_tensor_components(
    double *restrict Lambda_0, void *restrict Lambda_2, void *restrict Lambda_4,
    const double *R_2q, const double *R_2s, const double larmor_freq_in_Hz,
    const float mf, const float mi) {
  // Spin transition function scalar

  // R_2q = [-1/6 v_q η, 0, 1/√6 v_q, 0 , -1/6 v_q η] * d(mf, mi)
  //      = [-1/2 ζ_q η, 0, √3/2 ζ_q, 0 , -1/2 ζ_q η] * (v_q/3ζ_q) * d(mf, mi)
  //
  // R_2s = [1/√6 v_0ζη, 0, -v_0ζ,  0 , 1/√6 v_0ζη] * p(mf, mi)
  //      = [-1/2 ζη,    0, √3/2 ζ, 0 ,  -1/2 ζη  ] * -√2/3 * v_0 * p(mf, mi)

  // Divide R_2s by (-√2/3 v_0 p(mf, mi))
  double transition_fn_scalar = -1.2247448714 / (STF_p(mf, mi) * larmor_freq_in_Hz);

  // Spatial orientation function.
  sSOT_cross_tensor_components(Lambda_0, Lambda_2, Lambda_4, R_2q, R_2s);

  // The following multiplicative factor are scaled by sqrt(3/2) to compensate the
  // sqrt(2/3) factor from R_2s
  // frequency component function from the zeroth-rank irreducible tensor.
  *Lambda_0 *= -0.4140393356054125 * transition_fn_scalar;

  // frequency component function from the second-rank irreducible tensor.
  cblas_dscal(10, 0.24743582965269675 * transition_fn_scalar, (double *)Lambda_2, 1);

  // frequency component function from the fourth-rank irreducible tensor.
  cblas_dscal(18, -0.4426266681379905 * transition_fn_scalar, (double *)Lambda_4, 1);
}

// =====================================================================================
//      First-order J-coupling frequency tensor components (weakly coupling limit)
// =====================================================================================

/**
 * The frequency tensor (FT) components from the first-order perturbation expansion of
 * the J-coupling Hamiltonian (weak coupling limit), in a given frame,
 * @f$\mathcal{F}@f$, described by the Euler angles @f$\Theta = [\alpha, \beta,
 * \gamma]@f$ are
 * @f[
 *    {\Lambda'}_{0,0}^{(J)} &= \mathcal{R'}_{0,0}^{(J)}(\Theta)
 *                     ~~  \mathbb{d}_{IS}(m_{i_I}, m_{i_S}, m_{f_I}, m_{f_S}),
 *                     ~\text{and} \\
 *    {\Lambda'}_{2,n}^{(J)} &= \mathcal{R'}_{2,n}^{(J)}(\Theta)
 *                     ~~  \mathbb{d}_{IS}(m_{i_I}, m_{i_S}, m_{f_I}, m_{f_S}),
 * @f]
 * where @f$\mathcal{R'}_{0,0}^{(J)}(\Theta)@f$ and
 * @f$\mathcal{R'}_{2,n}^{(J)}(\Theta)@f$ are the spatial orientation functions in frame
 * @f$\mathcal{F}@f$, and @f$\mathbb{d}_{IS}(m_{i_I}, m_{i_S}, m_{f_I}, m_{f_S})@f$ is
 * the spin transition function for @f$\left|m_{i_I}, m_{i_S}\right> \rightarrow
 * \left|m_{f_I}, m_{f_S}\right>@f$ transition.
 *
 * @param Lambda_0 A pointer to an array of length 1, where the frequency
 *      component from @f${\Lambda'}_{0,0}^{(J)}@f$ is stored.
 * @param Lambda_2 A pointer to a complex array of length 5, where the frequency
 *      components from @f${\Lambda'}_{2,n}^{(J)}@f$ is stored ordered as
 *      @f$\left[{\Lambda'}_{2,n}^{(J)}\right]_{n=-2}^2@f$.
 * @param J_iso_in_Hz The isotropic J-coupling, @f$J_\text{iso}@f$, in Hz.
 * @param J_aniso_in_Hz The J-coupling anisotropy, @f$\zeta_J@f$, in Hz and
 *      defined using Haeberlen convention.
 * @param J_eta The J-coupling anisotropy asymmetry parameter, @f$\eta_J \in
 *      [-1,1]@f$, defined using Haeberlen convention.
 * @param Theta A pointer to an array of Euler angles of length 3 ordered as
 *      @f$[\alpha, \beta, \gamma]@f$.
 * @param mIf The spin quantum number of the final energy state of site @f$I@f$.
 * @param mIi The spin quantum number of the initial energy state of site
 *      @f$I@f$.
 * @param mSf The spin quantum number of the final energy state of site @f$S@f$.
 * @param mSi The spin quantum number of the initial energy state of site
 *      @f$S@f$.
 */
static inline void FCF_1st_order_weak_J_coupling_tensor_components(
    double *restrict Lambda_0, void *restrict Lambda_2, const double J_iso_in_Hz,
    const double J_aniso_in_Hz, const double J_eta, const double *Theta,
    const float mIf, const float mIi, const float mSf, const float mSi) {
  // Spin transition function
  double transition_fn = STF_dIS(mIf, mIi, mSf, mSi);

  // Return if the transition is zero
  if (transition_fn == 0.0) {
    // zero the R0 and R2 components before populating with shielding components
    *Lambda_0 = 0.0;
    vm_double_zeros(10, (double *)Lambda_2);
    return;
  }

  // Spatial orientation function
  sSOT_1st_order_weakly_coupled_J_tensor_components(Lambda_0, Lambda_2, J_iso_in_Hz,
                                                    J_aniso_in_Hz, J_eta, Theta);

  // frequency component function from the zeroth-rank irreducible tensor.
  *Lambda_0 *= transition_fn;

  // frequency component function from the second-rank irreducible tensor.
  cblas_dscal(10, transition_fn, (double *)Lambda_2, 1);
}

// =====================================================================================
//           First-order Weakly coupled Magnetic Dipole frequency in the PAS.
// =====================================================================================

/**
 * The frequency tensor (FT) components from the first-order perturbation expansion of
 * the direct dipolar coupling Hamiltonian (weak coupling limit), in a given frame,
 * @f$\mathcal{F}@f$, described by the Euler angles @f$\Theta = [\alpha, \beta,
 * \gamma]@f$ are
 * @f[
 *    {\Lambda'}_{2,n}^{(d)} = \mathcal{R'}_{2,n}^{(d)}(\Theta)
 *                     ~~  \mathbb{d}_{IS}(m_{i_I}, m_{i_S}, m_{f_I}, m_{f_S}),
 * @f]
 * where @f$\mathcal{R'}_{2,n}^{(d)}(\Theta)@f$ are the spatial orientation functions in
 * frame @f$\mathcal{F}@f$, and @f$\mathbb{d}_{IS}(m_{i_I}, m_{i_S}, m_{f_I},
 * m_{f_S})@f$ is the spin transition function for @f$\left|m_{i_I}, m_{i_S}\right>
 * \rightarrow \left|m_{f_I}, m_{f_S}\right>@f$ transition.
 *
 * @param Lambda_2 A pointer to a complex array of length 5, where the frequency
 *      components from @f${\Lambda'}_{2,n}^{(d)}@f$ is stored ordered as
 *      @f$\left[{\Lambda'}_{2,n}^{(d)}\right]_{n=-2}^2@f$.
 * @param D_in_Hz The dipolar coupling, @f$D@f$, in Hz.
 * @param Theta A pointer to an array of Euler angles of length 3 ordered as
 *      @f$[\alpha, \beta, \gamma]@f$.
 * @param mIf The spin quantum number of the final energy state of site @f$I@f$.
 * @param mIi The spin quantum number of the initial energy state of site @f$I@f$.
 * @param mSf The spin quantum number of the final energy state of site @f$S@f$.
 * @param mSi The spin quantum number of the initial energy state of site @f$S@f$.
 */
static inline void FCF_1st_order_weak_dipolar_coupling_tensor_components(
    void *restrict Lambda_2, const double D_in_Hz, const double *Theta, const float mIf,
    const float mIi, const float mSf, const float mSi) {
  // Spin transition function
  double transition_fn = STF_dIS(mIf, mIi, mSf, mSi);

  // Return if the transition is zero
  if (transition_fn == 0.0) {
    // zero the R0 and R2 components before populating with shielding components
    vm_double_zeros(10, (double *)Lambda_2);
    return;
  }

  // Spatial orientation function
  sSOT_1st_order_weakly_coupled_dipolar_tensor_components(Lambda_2, D_in_Hz, Theta);
  // frequency component function from the second-rank irreducible tensor.
  cblas_dscal(10, transition_fn, (double *)Lambda_2, 1);
}

// =====================================================================================
//                Quad-2nd rank tensor cross frequency tensor components
// =====================================================================================

/**
 *
 * @param Lambda_0 A pointer to an array of length 1, where the frequency
 *      component from @f${\Lambda'}_{0,0}^{(\sigma q)}@f$ is stored.
 * @param Lambda_2 A pointer to a complex array of length 5, where the frequency
 *      components from @f$\Lambda_{2,n}^{(\sigma q)}@f$ are stored ordered as
 *      @f$\left[{\Lambda'}_{2,n}^{(\sigma  q)}\right]_{n=-2}^2@f$.
 * @param Lambda_4 A pointer to a complex array of length 9, where the frequency
 *      components from @f${\Lambda'}_{4,n}^{(\sigma q)}@f$ are stored ordered as
 *      @f$\left[{\Lambda'}_{4,n}^{(\sigma q)}\right]_{n=-4}^4@f$.
 * @param R_2q A pointer to a complex array of length 5 holding the quadupolar
 *      tensor components.
 * @param symm_aniso_in_Hz symmetric tensor anisotropy in Hz
 * @param symm_eta symmetric tensor asymmetry in Hz
 * @param symm_theta pointer to euler angle for symmetric tensor
 * @param mIf The spin quantum number of the final energy state of site @f$I@f$.
 * @param mIi The spin quantum number of the initial energy state of site @f$I@f$.
 * @param mSf The spin quantum number of the final energy state of site @f$S@f$.
 * @param mSi The spin quantum number of the initial energy state of site @f$S@f$.
 */
static inline void FCF_Quad_coupling_cross_tensor_components(
    double *restrict Lambda_0, void *restrict Lambda_2, void *restrict Lambda_4,
    const double *R_2q, const double symm_aniso_in_Hz, const double symm_eta,
    const double *symm_theta, const double larmor_freq_in_Hz, const float mIf,
    const float mIi, const float mSf, const float mSi) {
  // Spin transition function scalar

  double R_2tensor[10], temp;
  double transition_fn_scalar = STF_pdIS(mIf, mIi, mSf, mSi);

  // Return if the transition is zero
  if (transition_fn_scalar == 0.0 || symm_aniso_in_Hz == 0.0) {
    // zero the R0 and R2 components before populating with shielding components
    *Lambda_0 = 0.0;
    vm_double_zeros(10, (double *)Lambda_2);
    vm_double_zeros(18, (double *)Lambda_4);
    return;
  }

  vm_double_zeros(10, (double *)R_2tensor);
  // R_2tensor = [-1/2 ζη, 0, √3/2 ζ, 0, -1/2 ζη]
  temp = -0.5 * (symm_aniso_in_Hz * symm_eta);
  R_2tensor[0] = temp;                                  // R2-2 real
  R_2tensor[4] = 1.224744871391589 * symm_aniso_in_Hz;  // R2 0 real
  R_2tensor[8] = temp;                                  // R2 2 real
  single_wigner_rotation(2, symm_theta, R_2tensor, R_2tensor);

  sSOT_cross_tensor_components(Lambda_0, Lambda_2, Lambda_4, R_2q, R_2tensor);

  transition_fn_scalar /= larmor_freq_in_Hz;

  // frequency component function from the zeroth-rank irreducible tensor.
  *Lambda_0 *= 0.4140393356054125 * transition_fn_scalar;

  // frequency component function from the second-rank irreducible tensor.
  cblas_dscal(10, -0.24743582965269675 * transition_fn_scalar, (double *)Lambda_2, 1);

  // frequency component function from the fourth-rank irreducible tensor.
  cblas_dscal(18, 0.4426266681379905 * transition_fn_scalar, (double *)Lambda_4, 1);
}
