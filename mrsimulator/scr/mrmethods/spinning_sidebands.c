
//
//  sideband_simulator.c
//
//  Created by Deepansh J. Srivastava, Apr 11, 2019
//  Copyright © 2019 Deepansh J. Srivastava. All rights reserved.
//  Contact email = srivastava.89@osu.edu, deepansh2012@gmail.com
//

// Definning pi^{2,2}_{L,J} as piLJ //
#define pi01 = 0.3577708764
#define pi21 = 0.1069044968
#define pi41 = -0.1434274331

#define pi03 = 0.8485281374
#define pi23 = -1.0141851057
#define pi43 = -1.2850792082
// --------------------------------- //

#include "spinning_sidebands.h"

// Return a vector ordered according to the fft output order.                //
// @params int n - The number of points                                      //
// @params double increment - The increment (sampling interval)              //
// @returns *double values = The pointer to the fft output order vector      //
inline double *__get_frequency_in_FFT_order(int n, double increment) {
  double *vr_freq = createDouble1DArray(n);
  int i = 0, m, positive_limit, negative_limit;

  if (n % 2 == 0) {
    negative_limit = (int)(-n / 2);
    positive_limit = -negative_limit - 1;
  } else {
    negative_limit = (int)(-(n - 1) / 2);
    positive_limit = -negative_limit;
  }

  for (m = 0; m <= positive_limit; m++) {
    vr_freq[i] = (double)m * increment;
    i++;
  }
  for (m = negative_limit; m < 0; m++) {
    vr_freq[i] = (double)m * increment;
    i++;
  }
  return vr_freq;
}

// Return the p(mi, mf) transition element                                   //
// The expression follows                                                    //
//        p(mf, mi) = < mf | T10 | mf > - < mi | T10 | mi >                  //
//                  = mf - mi                                                //
//                                                                           //
// @params double mi = quantum number of the initial state                   //
// @params double mf = quantum number of the final state                     //
// @returns double p = The value                                             //
static inline double __p__(double mf, double mi) { return (mf - mi); }

// Return the d(mi, mf) transition element                                   //
// The expression follows                                                    //
//        d(mf, mi) = < mf | T20 | mf > - < mi | T20 | mi >                  //
//                  = sqrt(3/2)*(mf^2 - mi^2)                                //
//                                                                           //
// @params double mi = quantum number of the initial state                   //
// @params double mf = quantum number of the final state                     //
// @returns double d = The value                                             //
static inline double __d__(double mf, double mi) {
  return 1.2247448714 * (mf * mf - mi * mi);
}

/*
  Return the dIS(mIi, mIf, mSi, mSf) transition element
  The expression follows
        d(mf, mi) = < mIf mSf | T10(I) T10(S) | mIf mSf > -
                    < mIi mSi | T10(I) T10(S) | mTi mSi >
                  = (mIf * mSf - mIi * mSi)

  @params double mIi = quantum number of the initial state of spin I
  @params double mIf = quantum number of the final state of spin I
  @params double mSi = quantum number of the initial state of spin S
  @params double mSf = quantum number of the final state of spin S
  @returns double dIS = The value
*/
static inline double __dIS__(double mIf, double mIi, double mSf, double mSi) {
  return mIf * mSf - mIi * mSi;
}

// Return the f(mi, mf) transition element                                   //
// The expression follows                                                    //
//        f(mf, mi) = < mf | T30 | mf > - < mi | T30 | mi >                  //
//                  = sqrt(1/10)*[5(mf^3 - mi^3)+(1-3I(I+1))(mf-mi)          //
//                                                                           //
// @params double mi = quantum number of the initial state                   //
// @params double mf = quantum number of the final state                     //
// @returns double f = The value                                             //
static inline double __f__(double mf, double mi, double spin) {
  double f = 1.0 - 3.0 * spin * (spin + 1.0);
  f *= (mf - mi);
  f += 5.0 * (mf * mf * mf - mi * mi * mi);
  f *= 0.316227766;
  return f;
}

static inline void __get_quad_ci__(double *c0, double *c2, double *c4,
                                   double mf, double mi, double spin) {
  double f = __f__(mf, mi, spin);
  double p = __p__(mf, mi);

  double temp = spin * (spin + 1.0) - 0.75;
  c0[0] = 0.3577708764 * temp * p + 0.8485281374 * f;
  c2[0] = 0.1069044968 * temp * p + -1.0141851057 * f;
  c4[0] = -0.1434274331 * temp * p + -1.2850792082 * f;
}

// ========================================================================= //
//          First order Nuclear shielding Hamiltonian in the PAS.            //
// ------------------------------------------------------------------------- //
// The Hamiltonian includes the product of second rank tensor and the        //
// spin transition functions.                                                //
// ------------------------------------------------------------------------- //
static inline void getNuclearShieldingHamiltonianUptoFirstOrder(
    double complex *R0, double complex *R2, double iso, double aniso,
    double eta, double *transition) {
  // Spin transition contribution
  double scale = __p__(transition[1], transition[0]);
  // printf("\n Entering CSA");
  // printf("\n Transition element %f \n", scale);
  // Scaled R00
  R0[0] = iso * scale;

  // Scaled R2m containing the components of the CSA second rank tensor in   //
  // its principal axis frame.                                               //
  double temp = -0.4082482905 * (aniso * eta) * scale;
  R2[0] += temp;          // R2-2
  R2[1] += 0.0;           // R2-1
  R2[2] += aniso * scale; // R2 0
  R2[3] += 0.0;           // R2 1
  R2[4] += temp;          // R2 2

  // printf("R00 %f \n", creal(R0[0]) );
  // printf("R2-2 %f \n", creal(R2[0]) );
  // printf("R2-1 %f \n", creal(R2[1]) );
  // printf("R20 %f \n", creal(R2[2]) );
  // printf("R21 %f \n", creal(R2[3]) );
  // printf("R22 %f \n", creal(R2[4]) );
}
// ========================================================================= //

// ========================================================================= //
//             First order Quadrupolar  Hamiltonian in the PAS.              //
// ------------------------------------------------------------------------- //
// The Hamiltonian includes the product of second rank tensor and the        //
// spin transition functions.                                                //
// ------------------------------------------------------------------------- //
static inline void getQuadrupoleHamiltonianUptoFirstOrder(double complex *R0,
                                                          double complex *R2,
                                                          double spin,
                                                          double Cq, double eta,
                                                          double *transition) {
  // Spin transition contribution
  double transition_d_ = __d__(transition[1], transition[0]);

  // printf("\n Entering Quad");
  // printf("\n Transition element %f \n", transition_d_);

  // Scaled R00
  R0[0] += 0.0;

  // vq is the Quadrupolar coupling constant                                 //
  // vq = 3*Cq/(2I(2I-1)), where I is the spin quantum number                //
  double vq = 3.0 * Cq;
  double denominator = 2.0 * spin * (2.0 * spin - 1.0);
  vq /= denominator;
  // printf("quad coupling constant, vq %f \n", vq);

  // Scaled R2m containing the components of the quad second rank tensor in  //
  // its principal axis frame. //
  double temp = -0.1666666667 * (vq * eta) * transition_d_;
  R2[0] += temp;                              // R2-2
  R2[1] += 0.0;                               // R2-1
  R2[2] += 0.4082482905 * vq * transition_d_; // R20
  R2[3] += 0.0;                               // R2 1
  R2[4] += temp;                              // R2 2

  // printf("R00 %f \n", creal(R0[0]) );
  // printf("R2-2 %f \n", creal(R2[0]) );
  // printf("R2-1 %f \n", creal(R2[1]) );
  // printf("R20 %f \n", creal(R2[2]) );
  // printf("R21 %f \n", creal(R2[3]) );
  // printf("R22 %f \n", creal(R2[4]) );
}
// ========================================================================= //

// ========================================================================= //
//            Second order Quadrupolar Hamiltonian in the PAS.               //
// ------------------------------------------------------------------------- //
// The Hamiltonian includes the product of second rank tensor and the        //
// spin transition functions.                                                //
// ------------------------------------------------------------------------- //
static inline void getQuadrupoleHamiltonianUptoSecondOrder(
    double complex *R0, double complex *R2, double complex *R4, double spin,
    double Cq, double eta, double *transition, double vo) {
  // Spin transition contribution
  double c0, c2, c4;
  __get_quad_ci__(&c0, &c2, &c4, transition[1], transition[0], spin);

  // printf("\n Entering Quad");
  // printf("\n Transition element c0=%f, c2=%f, c4=%f \n", c0, c2, c4);

  // vq is the Quadrupolar coupling constant                                 //
  // vq = 3*Cq/(2I(2I-1)), where I is the spin quantum number                //
  double vq = 3.0 * Cq;
  double denominator = 2.0 * spin * (2.0 * spin - 1.0);
  vq /= denominator;
  // printf("quad coupling constant, vq %f \n", vq);

  double scale = vq * vq / vo;
  double eta2 = eta * eta;

  // Scaled R00
  R0[0] += (eta2 * 0.33333333333 + 1.0) * 0.07453559925 * scale * c0;

  // Scaled R2m containing the components of the quad second rank tensor in  //
  // its principal axis frame. //
  double temp = -eta * 0.07273929675 * scale * c2;
  R2[0] += temp;                                                      // R2-2
  R2[1] += 0.0;                                                       // R2-1
  R2[2] += 0.08908708064 * (eta2 * 0.33333333333 - 1.0) * scale * c2; // R2 0
  R2[3] += 0.0;                                                       // R2 1
  R2[4] += temp;                                                      // R2 2

  // Scaled R2m containing the components of the quad fourth rank tensor in  //
  // its principal axis frame. //
  temp = eta2 * 0.02777777778 * scale * c4;
  double temp2 = -0.06299407883 * eta * scale * c4;
  R4[0] += temp;                                                     // R4-4
  R4[2] += temp2;                                                    // R4-2
  R4[4] += 0.1195228609 * (eta2 * 0.05555555556 + 1.0) * scale * c4; // R4 0
  R4[6] += temp2;                                                    // R4 2
  R4[8] += temp;                                                     // R4 4

  // printf("R00 %f \n", creal(R0[0]) );

  // printf("R2-2 %f \n", creal(R2[0]) );
  // printf("R2-1 %f \n", creal(R2[1]) );
  // printf("R20 %f \n", creal(R2[2]) );
  // printf("R21 %f \n", creal(R2[3]) );
  // printf("R22 %f \n", creal(R2[4]) );

  // printf("R4-4 %f \n", creal(R4[0]) );
  // printf("R4-2 %f \n", creal(R4[2]) );
  // printf("R40 %f \n", creal(R4[4]) );
  // printf("R42 %f \n", creal(R4[6]) );
  // printf("R44 %f \n", creal(R4[8]) );
}
// ========================================================================= //

// ========================================================================= //
//     First order Weakly coupled Magnetic Dipole Hamiltonian in the PAS.    //
// ------------------------------------------------------------------------- //
// The Hamiltonian includes the product of second rank tensor and the        //
// spin transition functions in the weak coupling limit.                     //
// ------------------------------------------------------------------------- //
static inline void getWeaklyCoupledMagneticDipoleHamiltonianUptoFirstOrder(
    double complex *R0, double complex *R2, double D, double *transition) {
  // Spin transition contribution
  double transition_dIS_ = __dIS__(transition[0], transition[1], 0.5, 0.5);

  // printf("\n Entering Quad");
  // printf("\n Transition element %f \n", transition_dIS_);

  // Scaled R00
  R0[0] += 0.0;

  // Scaled R2m containing the components of the magnetic dipole second rank //
  // tensor in its principal axis frame. //
  R2[0] += 0.0;                       // R2-2
  R2[1] += 0.0;                       // R2-1
  R2[2] += 2.0 * D * transition_dIS_; // R20
  R2[3] += 0.0;                       // R2 1
  R2[4] += 0.0;                       // R2 2

  // printf("R00 %f \n", creal(R0[0]) );
  // printf("R2-2 %f \n", creal(R2[0]) );
  // printf("R2-1 %f \n", creal(R2[1]) );
  // printf("R20 %f \n", creal(R2[2]) );
  // printf("R21 %f \n", creal(R2[3]) );
  // printf("R22 %f \n", creal(R2[4]) );
}
// ========================================================================= //

static inline void wigner_rotation(int l, double *wigner, double *cos_alpha,
                                   double *cos_gamma, double complex *scalex,
                                   double complex *initial_vector,
                                   double complex *final_vector) {
  int i = 0, m, mp, ll = 2 * l;
  double complex pha =
      cos_alpha[0] - I * sqrt(1.0 - cos_alpha[0] * cos_alpha[0]);
  double complex ph2 = pha;
  double complex temp_inital_vector[ll + 1];

  // copy the initial vector
  for (m = 0; m <= ll; m++) {
    temp_inital_vector[m] = initial_vector[m];
  }

  // scale the temp initial vector with exp[-I m alpha]
  for (m = 1; m <= l; m++) {
    temp_inital_vector[l + m] *= ph2;
    temp_inital_vector[l - m] *= conj(ph2);
    ph2 *= pha;
  }

  // Apply wigner rotation to the temp inital vector
  for (m = 0; m <= ll; m++) {
    final_vector[m] *= scalex[0];
    for (mp = 0; mp <= ll; mp++) {
      final_vector[m] += wigner[i++] * temp_inital_vector[mp];
    }
    // final_vector[ll-m] = creal(final_vector[m]) - I*cimag(final_vector[m]);
    // if(m%2!=0) final_vector[ll-m]*=-1.0;
  }

  // for(mp=0; mp<=ll; mp++){
  //     final_vector[l] += wigner[i++]*initial_vector[mp];
  // }
  // cblas_zgemv(CblasRowMajor, CblasNoTrans, ll, ll, &one, wigner, l,
  //               &temp_inital_vector[0], 1, scalex, &final_vector[0], 1);
}

static inline void __zero_components(double complex *R0, double complex *R2,
                                     double complex *R4) {
  int i;
  R0[0] = 0.0;
  for (i = 0; i <= 4; i++) {
    R2[i] = 0.0;
  }
  for (i = 0; i <= 8; i++) {
    R4[i] = 0.0;
  }
}

void __powder_averaging_setup(
    int nt, double *cosAlpha, double *cosBeta, double *amp,
    int space // 1 for octant, 2 for hemisphere and 4 for sphere
) {
  // unsigned int n_orientations = 1;
  if (space == 1) { // octant
    getPolarAngleTrigOverAnOctant(nt, cosAlpha, cosBeta, amp);
  }
}

void spinning_sideband_core(
    // spectrum information and related amplitude
    double *spec,              // The amplitude of the spectrum.
    double *cpu_time_,         // Execution time
    double spectral_start,     // The start of the frequency spectrum.
    double spectral_increment, // The increment of the frequency spectrum.
    int number_of_points,      // Number of points on the frequency spectrum.

    double spin_quantum_number, // Spin quantum numbers
    double larmor_frequency,    // Larmor frequency

    // Pointer to the array of CSA tensor information in the PAS.
    double *iso_n,   // The isotropic chemical shift.
    double *aniso_n, // The chemical shielding anisotropic.
    double *eta_n,   // Chemical shielding asymmetry

    // Pointer to the array of quadrupole tensor information in the PAS.
    double *Cq_e,        // The Cq of the quadrupole center.
    double *eta_e,       // The asymmetry term of the tensor.
    int quadSecondOrder, // Quad theory for second order,

    // Pointer to the array of dipolar tensor information in the PAS.
    double *D, // The dipolar coupling constant.

    // spin rate, spin angle and number spinning sidebands
    int number_of_sidebands, // The number of spinning sidebands to evaluate
    double spin_frequency,   // The rotor spin frequency
    double rotor_angle,      // The rotor angle relative to lab-frame z-axis

    // Pointer to the transitions. transition[0] = mi and transition[1] = mf
    double *transition,

    // The principal to molecular frame transformation euler angles.
    // double * omega_PM,

    // powder orientation average
    unsigned int n_orientations, // number of orientations
    double *cosAlpha,            // array of cosAlpha of orientations
    double *cosBeta,             // array of cosBeta of orientations
    double *amp,                 // array of amplitude of orientations
    int nt, // number of triangles along the edge of the octahedral face
    unsigned int number_of_sites // number of sites in the isotopomer
) {

  // The computation of the spinning sidebands is based on the method
  // described by Eden and Levitt et. al.
  // `Computation of Orientational Averages in Solid-State NMR by
  //  Gaussian Spherical Quadrature`
  //  JMR, 132, 1998. https://doi.org/10.1006/jmre.1998.1427

  // Time it
  clock_t start, end;
  start = clock();
  // mkl_set_num_threads( 4 );
  // fftw_init_threads();

  // Sampled over an octant
  // unsigned int n_orientations = (nt+1) * (nt+2)/2, orientation;
  unsigned int ji, ii, orientation;
  unsigned int size = n_orientations * number_of_sidebands;
  unsigned int site;
  int min_bound;

  // double* cosBeta = createDouble1DArray( n_orientations );
  // double* sinBeta = createDouble1DArray( n_orientations );
  // double* cosAlpha = createDouble1DArray( n_orientations );
  // double* sinAlpha = createDouble1DArray( n_orientations );

  // double* amp = createDouble1DArray(n_orientations);
  double *amp_temp = createDouble1DArray(n_orientations);

  int m, mp, step, i, allow_second_order_quad = 0;
  double tau, wrt, pht, spin_angular_freq, scale;
  // double number_of_sidebands_inverse = 1.0/((double) (number_of_sidebands));
  double spectral_increment_inverse = 1.0 / spectral_increment;
  double iso_n_, aniso_n_, eta_n_, Cq_e_, eta_e_, d_;

  double complex *pre_phase =
      createDoubleComplex1DArray(number_of_sidebands * 9);
  double complex *amp1;
  double *MR_full_DLM_2 = createDouble1DArray(25);
  double *MR_full_DLM_4 = createDouble1DArray(81);

  // set the rotor angle to 0 degree if the sample spin frequency is
  // less than 1 mHz
  if (spin_frequency < 1e-3) {
    spin_frequency = 1e6;
    rotor_angle = 0.0;
  }

  // Angle setup
  // getPolarAngleTrigOverAnOctant(nt, cosAlpha, cosBeta, amp);

  // __powder_averaging_setup(nt, &cosAlpha[0], &cosBeta[0], &amp[0], 1);   // 1
  // for octant, 2 for hemisphere and 4 for sphere Normalize the amplitudes
  // for(ii=0; ii<n_orientations; ii++){
  //   amp[ii]*=number_of_sidebands_inverse;
  // }
  // cblas_dscal(n_orientations, number_of_sidebands_inverse, amp, 1);

  // printf("\n Initializing system..........................");
  double complex R0[1] = {0.0};
  double complex R2[5] = {0.0, 0.0, 0.0, 0.0, 0.0};
  double complex R4[9] = {0., 0., 0., 0., 0., 0., 0., 0., 0.};

  double complex w_cs_2[5], rotor_lab_2[5];
  double complex w_cs_4[9], rotor_lab_4[9];
  double complex one = 1.0, zero = 0.0;
  double zero_f = 0.0;
  double shift = 0.0;
  if (number_of_points % 2 == 0) {
    shift = 0.5;
  }

  int spec_site;
  double *spec_site_ptr;

  // ----------------------------------------------------------------------- //
  // setup the fftw routine
  fftw_plan plan;
  fftw_complex *phi, *side_band;
  phi = (fftw_complex *)fftw_malloc(sizeof(fftw_complex) * number_of_sidebands);
  side_band =
      (fftw_complex *)fftw_malloc(sizeof(fftw_complex) * number_of_sidebands);
  plan = fftw_plan_dft_1d(number_of_sidebands, phi, side_band, FFTW_FORWARD,
                          FFTW_ESTIMATE);
  // fftw routine end
  // ----------------------------------------------------------------------- //

  // Calculate the spinning angular frequency
  spin_angular_freq = spin_frequency * PI2;

  // Generate the sideband order frequency relative to fft output order
  double *vr_freq =
      __get_frequency_in_FFT_order(number_of_sidebands, spin_frequency);
  // for(ii=0; ii<number_of_sidebands; ii++){
  //   vr_freq[ii]*=spectral_increment_inverse;
  // }
  cblas_dscal(number_of_sidebands, spectral_increment_inverse, vr_freq, 1);

  // create an empty array to hold the local spinning sideband frequencies.
  // This is useful when rotor angle is off the magic angle.
  double *local_frequency = createDouble1DArray(n_orientations);

  // Set up array for sideband_amplitudes
  double *sideband_amplitude = createDouble1DArray(size);

  // Calculate tau increments, where tau = (rotor period / number of phase
  // steps)
  tau = 1.0 / ((double)number_of_sidebands * spin_frequency);

  // pre-calculating the rotor to lab frame wigner terms
  for (mp = -4; mp <= 4; mp++) {
    rotor_lab_4[mp + 4] = wigner_d(4, mp, 0, rotor_angle);
  }
  for (mp = -2; mp <= 2; mp++) {
    rotor_lab_2[mp + 2] = wigner_d(2, mp, 0, rotor_angle);
  }

  // pre-calculate the m omega spinning frequencies
  double complex m_wr[9] = {-4., -3., -2., -1., 0., 1., 2., 3., 4.};
  // for(ii=0; ii<9; ii++){
  //   m_wr[ii]*=spin_angular_freq;
  // }
  cblas_zdscal(9, spin_angular_freq, &m_wr[0], 1);

  // pre-calculating the phase step exponents. --------------------------- //
  // --------------------------------------------------------------------- //
  //   phi = exp(sum_m w_cs_m * I 2pi [(exp(I m wr t) - 1)/(I m wr)])      //
  //   pre_phase(m, t) = I 2pi [(exp(I m wr t) - 1)/(I m wr)]              //
  //                   = (2 pi / m wr) (exp(I m wr t) - 1)                 //
  //                     ----scale----                                     //
  //                   = scale * (exp(I m wr t) - 1)                       //
  // --------------------------------------------------------------------- //
  i = 0;
  for (m = 0; m <= 8; m++) {
    if (m != 4) {
      wrt = m_wr[m] * tau;
      pht = 0.0;
      scale = PI2 / m_wr[m];
      // vzExp(number_of_sidebands, phi, pre_phase[m*number_of_sidebands]);
      for (step = 0; step < number_of_sidebands; step++) {
        // t = step * tau;
        // pht = m_wr[m] * t;
        pre_phase[i++] = scale * (cexp(I * pht) - 1.0);
        pht += wrt;
      }
    } else {
      i += number_of_sidebands;
    }
  }
  //  ---------------------------------------------------------------------- //

  // Per site base calculation
  for (site = 0; site < number_of_sites; site++) {
    spec_site = site * number_of_points;
    spec_site_ptr = &spec[spec_site];

    // Nuclear shielding terms
    iso_n_ = iso_n[site];
    aniso_n_ = aniso_n[site];
    eta_n_ = eta_n[site];

    // Electric quadrupolar terms
    Cq_e_ = Cq_e[site];
    eta_e_ = eta_e[site];

    // Magnetic dipole
    d_ = D[site];

    __zero_components(&R0[0], &R2[0], &R4[0]);

    getNuclearShieldingHamiltonianUptoFirstOrder(&R0[0], &R2[0], iso_n_,
                                                 aniso_n_, eta_n_, transition);

    getWeaklyCoupledMagneticDipoleHamiltonianUptoFirstOrder(&R0[0], &R2[0], d_,
                                                            transition);

    if (spin_quantum_number > 0.5) {
      getQuadrupoleHamiltonianUptoFirstOrder(
          &R0[0], &R2[0], spin_quantum_number, Cq_e_, eta_e_, transition);
      if (quadSecondOrder == 1) {
        allow_second_order_quad = 1;
        getQuadrupoleHamiltonianUptoSecondOrder(
            &R0[0], &R2[0], &R4[0], spin_quantum_number, Cq_e_, eta_e_,
            transition, larmor_frequency);
      }
    }

    // Equation [39] in the refernce above.
    //
    // w_cs^{m}(O_MR) = iso delta(m,0) + sum_{m', m" =-2}^{2} A[m"]
    // D^2_{m"m'}(O_PM) D^2_{m'm}(O_MR) d^2_{m'm}(b_RL)
    //

    for (orientation = 0; orientation < n_orientations; orientation++) {
      double *sideband_amplitude_f =
          &sideband_amplitude[number_of_sidebands * orientation];

      wigner_d_matrix(MR_full_DLM_2, 2, &cosBeta[orientation], 1);

      // ------------------------------------------------------------------- //
      //         Computing wigner rotation upto lab frame                    //

      // Second rank wigner rotation to rotor frame
      wigner_rotation(2, MR_full_DLM_2, &cosAlpha[orientation], &zero_f, &zero,
                      &R2[0], &w_cs_2[0]);

      // Second rank wigner rotation to lab frame cosidering alpha=gamma=0
      // vzMul( 5, &w_cs_2[0], &rotor_lab_2[0], &w_cs_2[0] );
      for (i = 0; i < 5; i++) {
        w_cs_2[i] *= rotor_lab_2[i];
      }

      // Fourth rank Wigner Rotation
      if (allow_second_order_quad) {

        wigner_d_matrix(MR_full_DLM_4, 4, &cosBeta[orientation], 1);
        wigner_rotation(4, MR_full_DLM_4, &cosAlpha[orientation], &zero_f,
                        &zero, R4, w_cs_4);

        // vzMul( 9, &w_cs_4[0], &rotor_lab_4[0], &w_cs_4[0] );
        for (i = 0; i < 9; i++) {
          w_cs_4[i] *= rotor_lab_4[i];
        }
      }

      // ------------------------------------------------------------------- //
      //    Computing phi = w_cs * I 2pi [(exp(i m wr t) - 1)/(i m wr)]      //
      //                           -------------- pre_phase------------      //
      // The pre_phase is calculated before.                                 //

      cblas_zgemv(CblasRowMajor, CblasTrans, 5, number_of_sidebands, &one,
                  &pre_phase[2 * number_of_sidebands], number_of_sidebands,
                  &w_cs_2[0], 1, &zero, phi, 1);

      if (allow_second_order_quad) {
        cblas_zgemv(CblasRowMajor, CblasTrans, 9, number_of_sidebands, &one,
                    pre_phase, number_of_sidebands, &w_cs_4[0], 1, &one, phi,
                    1);
      }

      // Computing exp(phi) ------------------------------------------------ //
      vzExp(number_of_sidebands, phi, phi);

      // Compute the fft --------------------------------------------------- //
      fftw_execute(plan);

      // ------------------------------------------------------------------- //
      // Taking the square of the the fft ampitudes
      for (m = 0; m < number_of_sidebands; m++) {
        amp1 = &side_band[m];
        sideband_amplitude_f[m] =
            creal(amp1[0]) * creal(amp1[0]) + cimag(amp1[0]) * cimag(amp1[0]);
      }

      // Multiplying the square amplitudes with the power scheme weights. -- //
      // And Normalizing with the number of phase steps squared ------------ //
      cblas_dscal(number_of_sidebands, amp[orientation], sideband_amplitude_f,
                  1);

      // adding the w_cs[0] term to the sideband frequencies before binning the
      // spectrum.
      local_frequency[orientation] =
          shift +
          (creal(w_cs_2[2]) + creal(w_cs_4[4]) + R0[0] - spectral_start) /
              spectral_increment;
    }

    // --------------------------------------------------------------------- //
    //              Calculating the tent for every sideband                  //
    // Allowing only sidebands that are within the spectral bandwidth ------ //
    for (i = 0; i < number_of_sidebands; i++) {
      min_bound = (int)(vr_freq[i] + iso_n[site] / spectral_increment);
      if (min_bound >= -number_of_points / 2 - 1 &&
          min_bound < number_of_points / 2 + 1) {
        ii = 0;
        ji = i;
        while (ii < n_orientations) {
          // ptr_ptr[ii] = &sideband_amplitude[ji];
          amp_temp[ii] = sideband_amplitude[ji];
          ii++;
          ji += number_of_sidebands;
        };

        // cblas_dcopy(n_orientations, &sideband_amplitude[0],
        // number_of_sidebands, &amp[0][0], 1);
        powderAverageWithTentingSchemeOverOctant(spec_site_ptr, local_frequency,
                                                 nt, amp_temp, &vr_freq[i],
                                                 number_of_points);
      }
    }
  }

  // clean up ------------------------------- //
  fftw_destroy_plan(plan);
  fftw_free(phi);
  fftw_free(side_band);
  destroyDoubleComplex1DArray(pre_phase);
  destroyDouble1DArray(MR_full_DLM_2);
  destroyDouble1DArray(MR_full_DLM_4);
  destroyDouble1DArray(vr_freq);
  destroyDouble1DArray(local_frequency);
  destroyDouble1DArray(sideband_amplitude);
  destroyDouble1DArray(amp_temp);

  end = clock();
  cpu_time_[0] += ((double)(end - start)) / (double)CLOCKS_PER_SEC;
}
