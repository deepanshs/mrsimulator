
//
//  sideband_simulator.c
//
//  Created by Deepansh J. Srivastava, Apr 11, 2019
//  Copyright © 2019 Deepansh J. Srivastava. All rights reserved.
//  Contact email = srivastava.89@osu.edu, deepansh2012@gmail.com
//

#include "spinning_sidebands.h"
// #include "mkl_dfti.h"

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

void __get_pre_phase_components(int number_of_sidebands,
                                double sample_rotation_frequency,
                                double complex *pre_phase) {
  double spin_angular_freq, tau, wrt, pht, scale;
  int step, i, m;
  // Calculate the spinning angular frequency
  spin_angular_freq = sample_rotation_frequency * PI2;

  // Calculate tau increments, where tau = (rotor period / number of phase
  // steps)
  tau = 1.0 / ((double)number_of_sidebands * sample_rotation_frequency);

  // pre-calculate the m omega spinning frequencies
  double m_wr[9] = {-4., -3., -2., -1., 0., 1., 2., 3., 4.};
  cblas_dscal(9, spin_angular_freq, &m_wr[0], 1);

  // pre-calculating the phase step exponents. --------------------------- //
  // --------------------------------------------------------------------- //
  //   phi = exp(sum_m v_cs_m * I 2pi [(exp(I m wr t) - 1)/(I m wr)])      //
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
      for (step = 0; step < number_of_sidebands; step++) {
        pre_phase[i++] = scale * (cexp(I * pht) - 1.0);
        pht += wrt;
      }
    } else {
      i += number_of_sidebands;
    }
  }
}
//  ---------------------------------------------------------------------- //

static inline void __spinning_sideband_core(
    // spectrum information and related amplitude
    double *spec,              // The amplitude of the spectrum.
    double spectral_start,     // The start of the frequency spectrum.
    double spectral_increment, // The increment of the frequency spectrum.
    int number_of_points,      // Number of points on the frequency spectrum.

    isotopomer_ravel *ravel_isotopomer, // isotopomer structure

    int quadSecondOrder,              // Quad theory for second order,
    int remove_second_order_quad_iso, // remove the isotropic contribution from
                                      // the second order quad Hamiltonian.

    // spin rate, spin angle and number spinning sidebands
    int number_of_sidebands, // The number of spinning sidebands to evaluate
    double sample_rotation_frequency, // The rotor spin frequency
    double rotor_angle, // The rotor angle relative to lab-frame z-axis

    // Pointer to the transitions. transition[0] = mi and transition[1] = mf
    double *transition,

    // fftw plan
    fftw_plan plan, fftw_complex *vector,

    // powder orientation average
    unsigned int n_orientations, // number of orientations
    int nt,

    // supplement
    double *local_frequency, // buffer for local frequencies
    double *freq_offset,     // buffer for frequencies + sideband frequencies
    double complex *w2,      // buffer for 2nd rank frequency calculation
    double complex *w4,      // buffer for 4nd rank frequency calculation
    double *vr_freq,         // sideband order frequencies in fft output order
    double *cos_alpha,       // array of cos_alpha of orientations
    double *amp,             // array of amplitude of orientations
    double *wigner_2j_matrices, // wigner-d 2j matrix for all orientations
    double *wigner_4j_matrices, // wigner-d 4j matrix for all orientations
    double complex *pre_phase_2, double complex *pre_phase_4,
    double *rotor_lab_2, double *rotor_lab_4) {

  clock_t start_site_time;
  double clock_time;
  // start = clock();
  /*
  The computation of the spinning sidebands is based on the method described by
  Eden and Levitt et. al.
    `Computation of Orientational Averages in Solid-State NMR by Gaussian
     Spherical Quadrature`
      JMR, 132, 1998. https://doi.org/10.1006/jmre.1998.1427
  */

  // Sampled over an octant
  unsigned int orientation, i, allow_second_order_quad = 0;
  unsigned int size = n_orientations * number_of_sidebands;
  unsigned int site;

  double spectral_increment_inverse = 1.0 / spectral_increment;
  double iso_n_, zeta_n_, eta_n_, Cq_e_, eta_e_, d_, offset, scale;

  double complex R0[1] = {0.0};
  double complex R2[5] = {0.0, 0.0, 0.0, 0.0, 0.0};
  double complex R4[9] = {0., 0., 0., 0., 0., 0., 0., 0., 0.};

  double complex one = 1.0, zero = 0.0;
  double shift_half_bin = 0.5;

  int spec_site;
  double *spec_site_ptr;
  double local_frequency_offset;

  // Per site base calculation
  for (site = 0; site < ravel_isotopomer[0].number_of_sites; site++) {
    start_site_time = clock();
    // start = clock();

    spec_site = site * number_of_points;
    spec_site_ptr = &spec[spec_site];

    // Nuclear shielding terms
    iso_n_ = ravel_isotopomer[0].isotropic_chemical_shift_in_Hz[site];
    zeta_n_ = ravel_isotopomer[0].shielding_anisotropy_in_Hz[site];
    eta_n_ = ravel_isotopomer[0].shielding_asymmetry[site];

    // Electric quadrupolar terms
    Cq_e_ = ravel_isotopomer[0].quadrupolar_constant_in_Hz[site];
    eta_e_ = ravel_isotopomer[0].quadrupolar_asymmetry[site];

    // Magnetic dipole
    d_ = ravel_isotopomer[0].dipolar_couplings[site];

    __zero_components(&R0[0], &R2[0], &R4[0]);

    get_nuclear_shielding_hamiltonian_to_first_order(
        &R0[0], &R2[0], iso_n_, zeta_n_, eta_n_, transition);

    get_weakly_coupled_direct_dipole_hamiltonian_to_first_order(&R0[0], &R2[0],
                                                                d_, transition);

    if (ravel_isotopomer[0].spin > 0.5) {
      get_quadrupole_hamiltonian_to_first_order(
          &R0[0], &R2[0], ravel_isotopomer[0].spin, Cq_e_, eta_e_, transition);
      if (quadSecondOrder == 1) {
        allow_second_order_quad = 1;
        get_quadrupole_hamiltonian_to_second_order(
            &R0[0], &R2[0], &R4[0], ravel_isotopomer[0].spin, Cq_e_, eta_e_,
            transition, ravel_isotopomer[0].larmor_frequency,
            remove_second_order_quad_iso);
      }
    }
    // ------------------------------------------------------------------- //

    // Equation [39] in the refernce above.
    //
    // w_cs^{m}(O_MR) = iso delta(m,0) + sum_{m', m" =-2}^{2} A[m"]
    // D^2_{m"m'}(O_PM) D^2_{m'm}(O_MR) d^2_{m'm}(b_RL)
    //

    local_frequency_offset = shift_half_bin + (creal(R0[0]) - spectral_start) *
                                                  spectral_increment_inverse;

    // -------------------------------------------------------------------
    //              Computing wigner rotation upto lab frame
    // Second rank wigner rotation
    __wigner_rotation(2, n_orientations, wigner_2j_matrices, cos_alpha, &R2[0],
                      &w2[0]);

    // Fourth rank wigner rotation
    if (allow_second_order_quad) {
      __wigner_rotation(4, n_orientations, wigner_4j_matrices, cos_alpha,
                        &R4[0], &w4[0]);
    }

    // -------------------------------------------------------------------
    //    Computing phi = w_cs * I 2pi [(exp(i m wr t) - 1)/(i m wr)]
    //                           -------------- pre_phase------------
    // second rank
    cblas_zgemm(CblasRowMajor, CblasTrans, CblasTrans, number_of_sidebands,
                n_orientations, 5, &one, pre_phase_2, number_of_sidebands, w2,
                5, &zero, vector, n_orientations);
    // fourth rank
    if (allow_second_order_quad) {
      cblas_zgemm(CblasRowMajor, CblasTrans, CblasTrans, number_of_sidebands,
                  n_orientations, 9, &one, pre_phase_4, number_of_sidebands, w4,
                  9, &one, vector, n_orientations);
    }

    // Compute exp(phi) ------------------------------------------------
    vmzExp(size, vector, vector, VML_EP);

    // Compute the fft ---------------------------------------------------
    fftw_execute(plan);
    // DftiComputeForward(plan, vector);

    // Taking the square of the the fft ampitudes
    vmdSqr(2 * size, (double *)&vector[0], (double *)&vector[0], VML_EP);
    cblas_daxpby(size, 1.0, (double *)&vector[0] + 1, 2, 1.0,
                 (double *)&vector[0], 2);

    // Multiplying the amplitude square with the power scheme weights.
    for (orientation = 0; orientation < n_orientations; orientation++) {
      cblas_dscal(number_of_sidebands, amp[orientation],
                  (double *)&vector[orientation], 2 * n_orientations);
    }

    // calculating local frequencies
    scale = spectral_increment_inverse * rotor_lab_2[2];
    cblas_daxpby(n_orientations, scale, (double *)&w2[2], 10, 0.0,
                 local_frequency, 1);
    if (allow_second_order_quad) {
      scale = spectral_increment_inverse * rotor_lab_4[4];
      cblas_daxpby(n_orientations, scale, (double *)&w4[4], 18, 1.0,
                   local_frequency, 1);
    }
    // ---------------------------------------------------------------------
    //              Calculating the tent for every sideband
    // Allowing only sidebands that are within the spectral bandwidth
    //
    for (i = 0; i < number_of_sidebands; i++) {
      offset = vr_freq[i] + local_frequency_offset;
      if ((int)offset >= 0 && (int)offset <= number_of_points) {

        vdLinearFrac(n_orientations, local_frequency, local_frequency, 1.0,
                     offset, 0.0, 1.0, freq_offset);
        octahedronInterpolation(spec_site_ptr, freq_offset, nt,
                                (double *)&vector[i * n_orientations], 2,
                                number_of_points);
      }
    }
    clock_time = ((double)(clock() - start_site_time)) / (double)CLOCKS_PER_SEC;
    printf("Total time per site %f \n", clock_time);
  }

  destroyDouble1DArray(local_frequency);
  destroyDouble1DArray(freq_offset);
}

void spinning_sideband_core(
    // spectrum information and related amplitude
    double *spec,              // The amplitude of the spectrum.
    double *cpu_time_,         // Execution time
    double spectral_start,     // The start of the frequency spectrum.
    double spectral_increment, // The increment of the frequency spectrum.
    int number_of_points,      // Number of points on the frequency spectrum.

    isotopomer_ravel *ravel_isotopomer, // isotopomer structure

    int quadSecondOrder,              // Quad theory for second order,
    int remove_second_order_quad_iso, // remove the isotropic contribution
                                      // from the second order quad
                                      // Hamiltonian.

    // spin rate, spin angle and number spinning sidebands
    int number_of_sidebands,          // The number of sidebands
    double sample_rotation_frequency, // The rotor spin frequency
    double rotor_angle, // The rotor angle relative to lab-frame z-axis

    // Pointer to the transitions. transition[0] = mi and transition[1] = mf
    double *transition,

    // powder orientation average
    int geodesic_polyhedron_frequency // The number of triangle along the edge
                                      // of octahedron
) {
  // Time it
  printf("start time %f s\n", cpu_time_[0]);
  clock_t begin, start1, all_site_time, all_c_time, fft_setup_time;
  all_c_time = clock();
  begin = clock();

  int nt = geodesic_polyhedron_frequency;
  unsigned int n_orientations = (nt + 1) * (nt + 2) / 2, i, size_2, size_4;
  int size = n_orientations * number_of_sidebands;

  // check for spinning speed
  if (sample_rotation_frequency < 1.0e-3) {
    sample_rotation_frequency = 1.0e9;
    rotor_angle = 0.0;
    number_of_sidebands = 1;
  }

  // setup the fftw routine
  fftw_plan plan;
  fftw_complex *vector;
  vector = (fftw_complex *)fftw_malloc(sizeof(fftw_complex) * size);
  fft_setup_time = clock();
  plan = fftw_plan_many_dft(1, &number_of_sidebands, n_orientations, vector,
                            NULL, n_orientations, 1, vector, NULL,
                            n_orientations, 1, FFTW_FORWARD, FFTW_ESTIMATE);
  // char *filename = "128_sidebands.wisdom";
  // int status = fftw_export_wisdom_to_filename(filename);
  // printf("file save status %i \n", status);

  // double complex *vector = createDoubleComplex1DArray(size);
  // DFTI_DESCRIPTOR_HANDLE plan;
  // MKL_LONG status;
  // MKL_LONG stride[2] = {0, n_orientations};
  // status = DftiCreateDescriptor(&plan, DFTI_DOUBLE, DFTI_COMPLEX, 1,
  //                               number_of_sidebands);
  // status = DftiSetValue(plan, DFTI_NUMBER_OF_TRANSFORMS, n_orientations);
  // status = DftiSetValue(plan, DFTI_INPUT_DISTANCE, 1);
  // status = DftiSetValue(plan, DFTI_OUTPUT_DISTANCE, 1);
  // status = DftiSetValue(plan, DFTI_INPUT_STRIDES, stride);
  // status = DftiSetValue(plan, DFTI_OUTPUT_STRIDES, stride);
  // DftiCommitDescriptor(plan);
  printf("fft time %f \n",
         ((double)(clock() - fft_setup_time)) / (double)CLOCKS_PER_SEC);

  double cost = fftw_cost(plan);
  printf("plan %f \n", cost);
  // fftw routine end

  start1 = clock();

  // create buffer for frequency calculations
  double complex *w2 = createDoubleComplex1DArray(5 * n_orientations);
  double complex *w4 = createDoubleComplex1DArray(9 * n_orientations);

  // Generate the sideband order frequency relative to fft output order
  double spectral_increment_inverse = 1.0 / spectral_increment;
  double *vr_freq = __get_frequency_in_FFT_order(number_of_sidebands,
                                                 sample_rotation_frequency);
  cblas_dscal(number_of_sidebands, spectral_increment_inverse, vr_freq, 1);

  // orientation set-up
  double *cos_alpha = createDouble1DArray(n_orientations);
  double *cos_beta = createDouble1DArray(n_orientations);
  double *amp_orientation = createDouble1DArray(n_orientations);
  __powder_averaging_setup(geodesic_polyhedron_frequency, cos_alpha, cos_beta,
                           amp_orientation, 1);
  // normalizing amplitudes
  double number_of_sideband_inverse = (1.0 / (double)number_of_sidebands);
  cblas_dscal(n_orientations, number_of_sideband_inverse, amp_orientation, 1);

  // setting up wigner 2j matrices for every orientation
  double *wigner_2j_matrices = createDouble1DArray(25 * n_orientations);
  __wigner_d_matrix_cosine(2, n_orientations, cos_beta, wigner_2j_matrices);
  double *wigner_4j_matrices = NULL;
  if (ravel_isotopomer[0].spin > 0.5) {
    double *wigner_4j_matrices = createDouble1DArray(81 * n_orientations);
    __wigner_d_matrix_cosine(4, n_orientations, cos_beta, wigner_4j_matrices);
  }
  destroyDouble1DArray(cos_beta);

  // rotor to lab frame tranformation setup
  double *rotor_lab_2 = createDouble1DArray(5);
  __wigner_dm0_vector(2, rotor_angle, &rotor_lab_2[0]);
  double *rotor_lab_4 = NULL;
  if (ravel_isotopomer[0].spin > 0.5) {
    double *rotor_lab_4 = createDouble1DArray(9);
    __wigner_dm0_vector(4, rotor_angle, rotor_lab_4);
  }

  // pre-calculating phase for side-band amplitudes
  size_4 = 9 * number_of_sidebands;
  double complex *pre_phase = createDoubleComplex1DArray(size_4);
  __get_pre_phase_components(number_of_sidebands, sample_rotation_frequency,
                             pre_phase);

  size_2 = 5 * number_of_sidebands;
  double complex *pre_phase_2 = createDoubleComplex1DArray(size_2);

  cblas_zcopy(size_2, &pre_phase[2 * number_of_sidebands], 1, pre_phase_2, 1);
  for (i = 0; i < 5; i++) {
    cblas_zdscal(number_of_sidebands, rotor_lab_2[i],
                 &pre_phase_2[i * number_of_sidebands], 1);
  }
  double complex *pre_phase_4 = NULL;
  if (ravel_isotopomer[0].spin > 0.5) {
    pre_phase_4 = &pre_phase[0];
    for (i = 0; i < 9; i++) {
      cblas_zdscal(number_of_sidebands, rotor_lab_4[i],
                   &pre_phase_4[i * number_of_sidebands], 1);
    }
  } else {
    destroyDoubleComplex1DArray(pre_phase);
  }

  // create an empty array to hold the local spinning sideband frequencies.
  // This is useful when rotor angle is off the magic angle.
  double *local_frequency = createDouble1DArray(n_orientations);
  double *freq_offset = createDouble1DArray(n_orientations);

  printf("setup time %f \n",
         ((double)(clock() - start1)) / (double)CLOCKS_PER_SEC);

  all_site_time = clock();
  __spinning_sideband_core(
      // spectrum information and related amplitude
      spec,               // The amplitude of the spectrum.
      spectral_start,     // The start of the frequency spectrum.
      spectral_increment, // The increment of the frequency spectrum.
      number_of_points,   // Number of points on the frequency spectrum.

      ravel_isotopomer, // isotopomer structure

      quadSecondOrder,              // Quad theory for second order,
      remove_second_order_quad_iso, // remove the isotropic contribution from
                                    // the second order quad Hamiltonian.

      // spin rate, spin angle and number spinning sidebands
      number_of_sidebands,       // The number of spinning sidebands to evaluate
      sample_rotation_frequency, // The rotor spin frequency
      rotor_angle,               // The rotor angle relative to lab-frame z-axis

      // Pointer to the transitions. transition[0] = mi and transition[1] = mf
      transition,

      // powder orientation average
      plan, vector,
      n_orientations,                // number of orientations
      geodesic_polyhedron_frequency, // number of triangles along the edge of
                                     // octahedron

      // supplement
      local_frequency,    // buffer for local frequencies
      freq_offset,        // buffer for frequencies + sideband frequencies
      w2,                 // buffer for 2nd rank frequency calculation
      w4,                 // buffer for 4nd rank frequency calculation
      vr_freq,            // sideband order frequencies in fft output order
      cos_alpha,          // array of cos_alpha of orientations
      amp_orientation,    // array of amplitude of orientations
      wigner_2j_matrices, // wigner-d 2j matrix for all orientations
      wigner_4j_matrices, // wigner-d 4j matrix for all orientations
      pre_phase_2, pre_phase_4, rotor_lab_2, rotor_lab_4);

  printf("all site time %f \n",
         ((double)(clock() - all_site_time)) / (double)CLOCKS_PER_SEC);

  printf("all c time %f \n",
         ((double)(clock() - all_c_time)) / (double)CLOCKS_PER_SEC);

  cpu_time_[0] += ((double)(clock() - begin)) / (double)CLOCKS_PER_SEC;
  printf("time %f s\n", cpu_time_[0]);

  // clean up ------------------------------- //
  fftw_destroy_plan(plan);
  fftw_free(vector);
  // DftiFreeDescriptor(&plan);
  // destroyDoubleComplex1DArray(vector);
  destroyDouble1DArray(rotor_lab_2);
  destroyDouble1DArray(rotor_lab_4);
  destroyDouble1DArray(wigner_2j_matrices);
  destroyDouble1DArray(wigner_4j_matrices);
  destroyDoubleComplex1DArray(pre_phase_2);
  destroyDoubleComplex1DArray(pre_phase_4);
  destroyDouble1DArray(amp_orientation);
  destroyDouble1DArray(cos_alpha);
  destroyDouble1DArray(vr_freq);
  destroyDoubleComplex1DArray(w2);
  destroyDoubleComplex1DArray(w4);
}
