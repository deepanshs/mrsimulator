

// #include "spinning_sidebands.h"


// void get_frequencies(
//           // spectrum information and related amplitude
// 					// List of frequencies. Total size = n_orientations*ph_step
// 					double *freq,                     
//           double * cpu_time_,               // Execution time

//           double *qunatum_number,           // Spin quantum numbers
//           double *wo,                       // omega_o

//           // Pointer to the array of CSA tensor information in the PAS. 
//           double *iso_n,                      // The isotropic chemical shift.
//           double *aniso_n,                    // The chemical shielding anisotropic.
//           double *eta_n,                      // Chemical shielding asymmetry

//           // Pointer to the array of quadrupole tensor information in the PAS. 
//           double *Cq_e,                       // The Cq of the quadrupole center.
//           double *eta_e,                      // The asymmetry term of the tensor.
//           int quadSecondOrder,                // Quad theory for second order, 

//           // spin rate, spin angle and number spinning sidebands
//           int ph_step,                      // The number of spinning sidebands to evaluate
//           double spin_frequency,            // The rotor spin frequency
//           double rotor_angle,               // The rotor angle relative to lab-frame z-axis

//           // Pointer to the transitions. transition[0] = mi and transition[1] = mf
//           double *transition,               

//           // The principal to molecular frame transformation euler angles.
//           // double * omega_PM,

//           directionCosines *cosines,                // Orientations
//           int n_orientations,               // number_of_orientations
//           unsigned int number_of_sites)
// {

//   // The following code is an adaption of Eden and Levitt et. al.
//   // `Computation of Orientational Averages in Solid-State NMR by
//   //  Gaussian Spherical Quadrature`
//   //  JMR, 132, 1998. https://doi.org/10.1006/jmre.1998.1427

//     // Time it
//     clock_t start, end;
//     start = clock();
//     // mkl_set_num_threads( 4 );
//     // fftw_init_threads();

//     // Sampled over an octant
//     int orientation;
//     // unsigned int n_orientations = (nt +1) * (2*nt+1), orientation;
//     // unsigned int increment = (2*nt+1) * ph_step;
//     unsigned int ji, ii;
//     int size = n_orientations * ph_step;
//     unsigned int site;

//     // double* cosBeta = createDouble1DArray( n_orientations );
//     // double* sinBeta = createDouble1DArray( n_orientations );
//     // double* cosAlpha = createDouble1DArray( n_orientations );
//     // double* sinAlpha = createDouble1DArray( n_orientations );

//     //
//     // double* amp = createDouble1DArray(n_orientations);
//     // double* amp_temp = createDouble1DArray(n_orientations);


//     // double** amp = createDouble2DMatrix(nt+1, 2*nt+1);
//     // double** amp_temp = createDouble2DMatrix(nt+1, 2*nt+1);
//     // double * amp_address_temp = &amp_temp[0];

//     // double *ptr[nt+1][2*nt+1];
//     // double **ptr_ptr = &ptr[0][0];

//     int m, mp, step, i, allow_second_order_quad=0;
//     double tau, t, pht, spin_angular_freq;
//     double ph_step_inverse = 1.0/((double) (ph_step));

//     // temporary interaction terms 
//     double iso_n_, aniso_n_, eta_n_, Cq_e_, eta_e_;

//     // double complex molecular_rotor;
//     double complex *pre_phase = createDoubleComplex1DArray(ph_step*9);
//     double complex *amp1;
//     double complex* MR_full_DLM_2 = createDoubleComplex1DArray(25);
//     double complex* MR_full_DLM_4 = createDoubleComplex1DArray(81);
//     // double complex alpha[9], phase_alpha;
//     // double complex* PM_full_DLM = createDoubleComplex1DArray(25);

//     // Angle setup
//     // getPolarAngleTrigOverAnOctant(nt, cosAlpha, sinAlpha, cosBeta, sinBeta, amp);

//     // Normalize the amplitudes
//     // cblas_dscal(n_orientations, ph_step_inverse, amp, 1);
//     // temp = ph_step_inverse_square*amp[orientation];

    
//     printf("\n Initializing system..........................");
//     double complex R0[1] = {0.0};
//     double complex R2[5] = {0.0, 0.0, 0.0, 0.0, 0.0};
//     double complex R4[9] = {0., 0., 0., 0., 0., 0., 0., 0., 0.};

//     double complex w_cs_2[5], rotor_lab_2[5];
//     double complex w_cs_4[9], rotor_lab_4[9];
//     double complex one=1.0, zero=0.0;

//     int freq_site;
//     double * freq_site_ptr;
// 		double sinBeta, sinAlpha;

//     // ----------------------------------------------------------------------- //
//     // setup the fftw routine
//     fftw_plan plan;
//     fftw_complex *phi, *side_band;
//     phi = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * ph_step);
//     side_band = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * ph_step);
//     plan = fftw_plan_dft_1d(ph_step, phi, side_band, FFTW_FORWARD, FFTW_ESTIMATE);
//     // fftw routine end
//     // ----------------------------------------------------------------------- //



//     // ----------------------------------------------------------------------- //
//     // prep for tenting --> histogram
//     // spin_frequency/=spectral_increment;
//     // spectral_start/=spectral_increment;

//     // Calcuate the spinning angular frequency
//     spin_angular_freq = spin_frequency * PI2;

//     // Generate the sideband order frequency relative to fft output order
//     double *vr_freq = __get_frequency_in_FFT_order(ph_step, spin_frequency);
//     // cblas_dscal(ph_step, 1./spectral_increment, vr_freq, 1);


//     // create an empty array to hold the local spinning sideband frequencies.
//     // This is useful when rotor angle is off the magic angle.
//     double *local_frequency = createDouble1DArray( n_orientations );


//     // Set up array for sideband_amplitudes
//     double *sideband_amplitude = createDouble1DArray( number_of_sites*size );


//     // Calculate tau increments, where tau = (rotor period / number of phase steps)
//     tau = 1.0/((double)ph_step*spin_frequency);

//     // Also pre-calculating the rotor to lab frame wigner terms
//     for(mp = -4; mp <= 4; mp++){
//       rotor_lab_4[mp+4] = wigner_d(4, mp, 0, rotor_angle);
//     }
//     for(mp = -2; mp <= 2; mp++){
//       rotor_lab_2[mp+2] = wigner_d(2, mp, 0, rotor_angle);
//     }

    

//     // pre-calculate the m omega spinning frequencies
//     double complex m_wr[9] = {-4., -3., -2., -1., 0., 1., 2., 3., 4.};
//     cblas_zdscal(9, spin_angular_freq, &m_wr[0], 1);

//     // pre-calculating the phase step exponents. ----------------------------- //
//     i = 0;
//     for(m =0; m<=8; m++){
//       if (m!=4){
//         // vzExp(ph_step, phi, pre_phase[m*ph_step]);
//         for(step=0; step<ph_step; step++){
//           t = step * tau;
//           pht = m_wr[m] * t;
//           pre_phase[i++] = PI2I * (cexp(I*pht) - 1.0)/ (I * m_wr[m]);
//         }
//       }
//       else{
//         i+=ph_step;
//       }
//     }
//     //  ---------------------------------------------------------------------- //


//     // Per site base calculation
//     for(site = 0; site < number_of_sites; site++){
//       freq_site = site*size;
//       freq_site_ptr = &freq[freq_site];

//       // Scaling all the interactions with the frequency interval.
//       // Nuclear shielding terms
//       iso_n_ = iso_n[site]; //spectral_increment;
//       aniso_n_ = aniso_n[site]; //spectral_increment;
//       eta_n_ = eta_n[site];

//       // Electric quadrupolar terms
//       Cq_e_ = Cq_e[site]; //spectral_increment;
//       eta_e_ = eta_e[site];


//       getNuclearShieldingUptoFirstRs(
//                     &R0[0],
//                     &R2[0],
//                     iso_n_,
//                     aniso_n_,
//                     eta_n_,
//                     transition
//       );

//       if (qunatum_number[site] > 0.5){
//         getQuandropuleShieldingUptoFirstOrderRs(
//                     &R0[0],
//                     &R2[0],
//                     qunatum_number[site],
//                     Cq_e_,
//                     eta_e_,
//                     transition
//           );
//         if (quadSecondOrder == 1){
//           allow_second_order_quad = 1;
//           getQuandropuleShieldingUptoSecondOrderRs(
//                     &R0[0],
//                     &R2[0],
//                     &R4[0],
//                     qunatum_number[site],
//                     Cq_e_,
//                     eta_e_,
//                     transition,
//                     wo[0]);
//         }
//       }

//       // Equation [39] in the refernce above.
//       //
//       // w_cs^{m}(O_MR) = iso delta(m,0) + sum_{m', m" =-2}^{2} A[m"] D^2_{m"m'}(O_PM) D^2_{m'm}(O_MR) d^2_{m'm}(b_RL)
//       //

//       for(orientation=0; orientation < n_orientations; orientation++){
//         double *sideband_amplitude_f = &sideband_amplitude[ph_step*orientation];

//         // phase_alpha = cosAlpha[orientation]-I*sinAlpha[orientation];
//         // alpha[4]=1.0;
//         // for(i=1; i<=4; i++){
//         //   alpha[i+4] = alpha[i+3] * phase_alpha;
//         //   alpha[-i+4] = conj(alpha[i+4]);
//         // }

//         // wigner_d_matrix(MR_full_DLM_2, 2, &cosBeta[orientation], 1);
        
//         // for(i=0;i<5;i++){
//         //   cblas_zscal(5, &alpha[i+2], &MR_full_DLM_2[i], 5);
//         // }
// 				sinAlpha = sqrt(1. - pow(cosines[orientation].cosAlpha, 2));
// 				sinBeta = sqrt(1. - pow(cosines[orientation].cosBeta, 2));
//         full_DLM_trig(MR_full_DLM_2, 2,
//               cosines[orientation].cosAlpha,
//               sinAlpha,
// 							cosines[orientation].cosBeta,
//               sinBeta);

//         // ------------------------------------------------------------------- //
//         //         Computing wigner rotation w_cs_PM to w_cs_MR                //
//         //         w_cs_PM are the A2m components in the Molecular frame.      //
//         //         w_cs are the A2m components in the rotor frame.             //

//         // cblas_zgemv(CblasRowMajor, CblasNoTrans, 5, 5, &one, MR_full_DLM, 5,
//         //               &w_cs_PM[0], 1, &zero, &w_cs[0], 1);

//         // Second rank wigner rotation
//         cblas_zgemv(CblasRowMajor, CblasNoTrans, 5, 5, &one, MR_full_DLM_2, 5,
//                       &R2[0], 1, &zero, &w_cs_2[0], 1);
        
        
        


//         // Fourth rank Wigner Rotation
//         if (allow_second_order_quad){
//           // get_even_DLM_4_from_2(MR_full_DLM_2, cosBeta[orientation]);
//           full_DLM_trig(MR_full_DLM_4, 4,
//               cosines[orientation].cosAlpha,
//               sinAlpha,
//               cosines[orientation].cosBeta,
//               sinBeta);

//           cblas_zgemv(CblasRowMajor, CblasNoTrans, 9, 9, &one, MR_full_DLM_4, 9,
//                       &R4[0], 1, &zero, &w_cs_4[0], 1);
//         }


//         // ------------------------------------------------------------------- //
//         //         Rotor to lab frame transformation.													 //
//         //         multiplying w_cs by lab frame wigner d elements. 					 //
//         //    Here, the alpha_RL = gamma_RL = 0 and beta is the rotor angle.	 //

//         vzMul( 5, &w_cs_2[0], &rotor_lab_2[0], &w_cs_2[0] );

//         if (allow_second_order_quad){
//           vzMul( 9, &w_cs_4[0], &rotor_lab_4[0], &w_cs_4[0] );
//         }


//         // ------------------------------------------------------------------- //
//         //    Computing phi = w_cs * I 2pi [(exp(i m wr t) - 1)/(i m wr)]      //
//         //                           -------------- pre_phase------------      //
//         //      The pre_phase is calculated before.                            //
//         cblas_zgemv(CblasRowMajor, CblasTrans, 5, ph_step, &one,
//                       &pre_phase[2*ph_step], ph_step, &w_cs_2[0], 
//                       1, &zero, phi, 1);

//         if (allow_second_order_quad){
//           cblas_zgemv(CblasRowMajor, CblasTrans, 9, ph_step, &one,
//                       pre_phase, ph_step, &w_cs_4[0], 
//                       1, &one, phi, 1);
//         }

//         // Computing exp(phi) ---------------------------------------------- //
//         vzExp(ph_step, phi, phi);

//         // Compute the fft --------------------------------------------------- //
//         fftw_execute(plan);


//         // ------------------------------------------------------------------- //
//         // Taking the square of the the fft ampitudes
//         // vzMulByConj(ph_step, side_band, side_band, sideband_amplitude_f);

//         // vzAbs( ph_step, side_band, sideband_amplitude_f);
//               // vdPowx( ph_step, sideband_amplitude_f, 2, sideband_amplitude_f);

//         // Taking the square of the the fft ampitudes
//         for(m=0; m<ph_step; m++){
//           amp1 = &side_band[m];
//           sideband_amplitude_f[m] = creal(amp1[0])*creal(amp1[0]) + cimag(amp1[0])*cimag(amp1[0]);
//                   // sideband_amplitude_f[m] = pow(cblas_dznrm2(1, &side_band[m], 1),2);
//         }


//         // Multiplying the square amplitudes with the power scheme weights. -- //
//         // And Normalizing with the number of phase steps squared ------------ //
//         // cblas_dscal(ph_step, amp[orientation], sideband_amplitude_f, 1);

//         // adding the w_cs[0] term to the sideband frequencies before binning the spectrum.
//         local_frequency[orientation] = creal(w_cs_2[2]) + creal(w_cs_4[4]) + R0[0];
//       }
//     }

//   // clean up ------------------------------- //
//   fftw_destroy_plan(plan);
//   fftw_free(phi); fftw_free(side_band);
//   destroyDoubleComplex1DArray(pre_phase);
//   // destroyDouble1DArray(cosAlpha);
//   // destroyDouble1DArray(cosBeta);
//   // destroyDouble1DArray(sinAlpha);
//   // destroyDouble1DArray(sinBeta);
//   destroyDoubleComplex1DArray(MR_full_DLM_2);
//   destroyDouble1DArray(vr_freq);
//   destroyDouble1DArray(local_frequency);
//   destroyDouble1DArray(sideband_amplitude);
//   // destroyDouble1DArray(amp);
//   // destroyDouble1DArray(amp_temp);

//   end = clock();
//   cpu_time_[0] = ((double) (end - start)) / (double) CLOCKS_PER_SEC;
// }


// void gettent(
// 		double *spec,
//     double spectral_start,
//     double spectral_increment,
//     double number_of_points,
//     double *freq,
// 		double *amp,
//     double number_of_sites,
// 		int n_orientations,
// 		int nt,
//     int ph_step)

// {
//     int site, i, ii, ji;
// 		double min_bound, zero=0.0;
// 		double *local_frequency = createDouble1DArray(n_orientations*ph_step);
// 		double *amp_temp = createDouble1DArray(n_orientations*ph_step);
// 		double *spec_site_ptr;
//     for(site = 0; site < number_of_sites; site++){
// 			spec_site_ptr = &spec[site*n_orientations*ph_step];
//     // --------------------------------------------------------------------- //
//     //              Calculating the tent for every sideband                  //
//     // Allowing only sidebands that are within the spectral bandwidth ------ //
//       freq[0] = (freq[0] - spectral_start)/spectral_increment;
//       for(i =0; i<ph_step; i++){
//         // min_bound = (int) (vr_freq[i] + iso_n[site]/spectral_increment);
//         // if (min_bound >= -number_of_points/2 -1 && min_bound < number_of_points/2 +1){
//           ii = 0;
//           ji = i;
//           while(ii<n_orientations){
//             // ptr_ptr[ii] = &sideband_amplitude[ji];
//             amp_temp[ii] = amp[ji];
// 						local_frequency[ii] = freq[ji];
//             ii++;
//             ji+=ph_step;
//           };

//           // cblas_dcopy(n_orientations, &sideband_amplitude[0], ph_step, &amp[0][0], 1);
//           powderAverageWithTentingSchemeOverOctant2(
//                       spec_site_ptr,
//                       local_frequency,
//                       nt,
//                       amp_temp,
//                       &zero,
//                       number_of_points);
//         }
//       }
// 		}
// }