//
//  Hamiltonian.h
//
//  Created by Deepansh J. Srivastava, Apr 11, 2019
//  Copyright © 2019 Deepansh J. Srivastava. All rights reserved.
//  Contact email = srivastava.89@osu.edu, deepansh2012@gmail.com
//

#include "mrsimulator.h"

// Definning pi^{2,2}_{L,J} as piLJ //
#define pi01 = 0.3577708764
#define pi21 = 0.1069044968
#define pi41 = -0.1434274331

#define pi03 = 0.8485281374
#define pi23 = -1.0141851057
#define pi43 = -1.2850792082
// --------------------------------- //

/*
Return the p(mi, mf) transition element.
The expression follows,
      p(mf, mi) = < mf | T10 | mf > - < mi | T10 | mi >
                = mf - mi

@params double mi = quantum number of the initial state,
@params double mf = quantum number of the final state.
@returns double p = The value.
*/
static inline double __p__(double mf, double mi) { return (mf - mi); }


/*
Return the d(mi, mf) transition element.
The expression follows,
      d(mf, mi) = < mf | T20 | mf > - < mi | T20 | mi >
                = sqrt(3/2)*(mf^2 - mi^2)

@params double mi = quantum number of the initial state.
@params double mf = quantum number of the final state.
@returns double d = The value..
 */
static inline double __d__(double mf, double mi)
{
  return 1.2247448714 * (mf * mf - mi * mi);
}

/*
Return the dIS(mIi, mIf, mSi, mSf) transition element.
The expression follows,
      d(mf, mi) = < mIf mSf | T10(I) T10(S) | mIf mSf > -
                  < mIi mSi | T10(I) T10(S) | mTi mSi >
                = (mIf * mSf - mIi * mSi)

@params double mIi = quantum number of the initial state of spin I.
@params double mIf = quantum number of the final state of spin I.
@params double mSi = quantum number of the initial state of spin S.
@params double mSf = quantum number of the final state of spin S.
@returns double dIS = The value.
*/
static inline double __dIS__(double mIf, double mIi, double mSf, double mSi)
{
  return mIf * mSf - mIi * mSi;
}

/*
Return the f(mi, mf) transition element.
The expression follows
      f(mf, mi) = < mf | T30 | mf > - < mi | T30 | mi >
                = sqrt(1/10)*[5(mf^3 - mi^3)+(1-3I(I+1))(mf-mi)

@params double mi = quantum number of the initial state.
@params double mf = quantum number of the final state.
@returns double f = The value.
*/
static inline double __f__(double mf, double mi, double spin)
{
  double f = 1.0 - 3.0 * spin * (spin + 1.0);
  f *= (mf - mi);
  f += 5.0 * (mf * mf * mf - mi * mi * mi);
  f *= 0.316227766;
  return f;
}

static inline void __get_quad_ci__(double *c0, double *c2, double *c4,
                                   double mf, double mi, double spin)
{
  double f = __f__(mf, mi, spin);
  double p = __p__(mf, mi);

  double temp = spin * (spin + 1.0) - 0.75;
  c0[0] = 0.3577708764 * temp * p + 0.8485281374 * f;
  c2[0] = 0.1069044968 * temp * p + -1.0141851057 * f;
  c4[0] = -0.1434274331 * temp * p + -1.2850792082 * f;
}

/*
===============================================================================
            First order Nuclear shielding Hamiltonian in the PAS.
-------------------------------------------------------------------------------
The Hamiltonian includes the product of second rank tensor and the
spin transition functions.
*/
static inline void get_nuclear_shielding_hamiltonian_to_first_order(
    double *R0, double complex *R2, double iso, double zeta, double eta,
    double *transition)
{
  // Spin transition contribution
  double scale = __p__(transition[1], transition[0]);

  // Scaled R00
  R0[0] = iso * scale;

  /*
  Scaled R2m containing the components of the quad second rank tensor in
  its principal axis frame.
  */
  double temp = -0.4082482905 * (zeta * eta) * scale;
  R2[0] += temp;         // R2-2
  R2[1] += 0.0;          // R2-1
  R2[2] += zeta * scale; // R2 0
  R2[3] += 0.0;          // R2 1
  R2[4] += temp;         // R2 2
}

/*
===============================================================================
                First order Quadrupolar Hamiltonian in the PAS.
-------------------------------------------------------------------------------
The Hamiltonian includes the product of second rank tensor and the
spin transition functions.
*/
static inline void get_quadrupole_hamiltonian_to_first_order(
    double *R0, double complex *R2, double spin, double Cq, double eta,
    double *transition)
{
  // Spin transition contribution
  double transition_d_ = __d__(transition[1], transition[0]);

  // Scaled R00
  R0[0] += 0.0;

  /* vq is the Quadrupolar coupling constant given as vq = 3*Cq/(2I(2I-1)),
  where `I` is the spin quantum number. */
  double vq = 3.0 * Cq;
  double denominator = 2.0 * spin * (2.0 * spin - 1.0);
  vq /= denominator;

  /* Scaled R2m containing the components of the quad second rank tensor in
  its principal axis frame. */
  double temp = -0.1666666667 * (vq * eta) * transition_d_;
  R2[0] += temp;                              // R2-2
  R2[1] += 0.0;                               // R2-1
  R2[2] += 0.4082482905 * vq * transition_d_; // R2 0
  R2[3] += 0.0;                               // R2 1
  R2[4] += temp;                              // R2 2
}

/*
===============================================================================
               Second order Quadrupolar Hamiltonian in the PAS.
-------------------------------------------------------------------------------
The Hamiltonian includes the product of second rank tensor and the
spin transition functions.
*/
static inline void get_quadrupole_hamiltonian_to_second_order(
    double *R0, double complex *R2, double complex *R4, double spin,
    double Cq, double eta, double *transition, double vo,
    int remove_second_order_quad_iso)
{
  // Spin transition contribution
  double c0, c2, c4;
  __get_quad_ci__(&c0, &c2, &c4, transition[1], transition[0], spin);

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
  R2[0] += temp;                                                      // R2-2
  R2[1] += 0.0;                                                       // R2-1
  R2[2] += 0.08908708064 * (eta2 * 0.33333333333 - 1.0) * scale * c2; // R2 0
  R2[3] += 0.0;                                                       // R2 1
  R2[4] += temp;                                                      // R2 2

  /* Scaled R4m containing the components of the quad second rank tensor in
  its principal axis frame. */
  temp = eta2 * 0.02777777778 * scale * c4;
  double temp2 = -0.06299407883 * eta * scale * c4;
  R4[0] += temp;                                                     // R4-4
  R4[2] += temp2;                                                    // R4-2
  R4[4] += 0.1195228609 * (eta2 * 0.05555555556 + 1.0) * scale * c4; // R4 0
  R4[6] += temp2;                                                    // R4 2
  R4[8] += temp;                                                     // R4 4
}

/*
===============================================================================
       First order Weakly coupled Magnetic Dipole Hamiltonian in the PAS.
-------------------------------------------------------------------------------
The Hamiltonian includes the product of second rank tensor and the
spin transition functions in the weak coupling limit.
*/
static inline void get_weakly_coupled_direct_dipole_hamiltonian_to_first_order(
    double *R0, double complex *R2, double D, double *transition)
{
  // Spin transition contribution
  double transition_dIS_ = __dIS__(transition[0], transition[1], 0.5, 0.5);

  // Scaled R00
  R0[0] += 0.0;

  /* Scaled R2m containing the components of the magnetic dipole second rank
  tensor in its principal axis frame. */
  R2[0] += 0.0;                       // R2-2
  R2[1] += 0.0;                       // R2-1
  R2[2] += 2.0 * D * transition_dIS_; // R2 0
  R2[3] += 0.0;                       // R2 1
  R2[4] += 0.0;                       // R2 2
}
