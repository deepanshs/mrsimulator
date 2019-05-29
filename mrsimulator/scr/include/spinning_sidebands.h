//
//  spinning_sidebands.h
//
//  Created by Deepansh J. Srivastava, Apr 11, 2019
//  Copyright © 2019 Deepansh J. Srivastava. All rights reserved.
//  Contact email = srivastava.89@osu.edu, deepansh2012@gmail.com
//

#define PI2 6.2831853072
#define PI2I PI2*I

#include "MRAngularMomentum.h"
#include "powder_setup.h"
#include <complex.h>
#include <time.h>
#include "math.h"
#include "fftw/fftw3.h"
#include "fftw/fftw3_mkl.h"
#include "c_array.h"
#include "mkl.h"

#define MKL_Complex16 double complex

struct directionCosines
{
	double cosAlpha;
	double cosBeta;
}; 

extern void __powder_averaging_setup(
          int nt,
          double *cosAlpha, 
          double *cosBeta,
          double *amp,
          int space   // 1 for octant, 2 for hemisphere and 4 for sphere
);



// Return a vector ordered according to the fft output order.                //
// @params int n - The number of points                                      //
// @params double increment - The increment (sampling interval)              //
// @returns *double values = The pointer to the fft output order vector      //
extern inline double* __get_frequency_in_FFT_order(
                                      int n, 
                                      double increment
                                  );

extern void spinning_sideband_core(
          // spectrum information and related amplitude
          double * spec,                    // The amplitude of the spectrum.
          double * cpu_time_,               // Execution time
          double spectral_start,            // The start of the frequency spectrum.
          double spectral_increment,        // The bandwidth of the frequency spectrum.
          int number_of_points,             // Number of points on the frequency spectrum.

          double spin_quantum_number,       // Spin quantum numbers
          double larmor_frequency,          // Larmor frequency

          // Pointer to the array of CSA tensor information in the PAS. 
          double *iso_n,                      // The isotropic chemical shift.
          double *aniso_n,                    // The chemical shielding anisotropic.
          double *eta_n,                      // Chemical shielding asymmetry

          // Pointer to the array of quadrupole tensor information in the PAS. 
          double *Cq_e,                       // The Cq of the quadrupole center.
          double *eta_e,                      // The asymmetry term of the tensor.
          int quadSecondOrder,                // Quad theory for second order, 

          // Pointer to the array of dipolar tensor information in the PAS. 
          double *D,                          // The dipolar coupling constant.

          // spin rate, spin angle and number spinning sidebands
          int ph_step,                      // The number of spinning sidebands to evaluate
          double spin_frequency,            // The rotor spin frequency
          double rotor_angle,               // The rotor angle relative to lab-frame z-axis

          double *transition,               // The transition as transition[0] = mi and transition[1] = mf

          // The principal to molecular frame transformation euler angles.
          //   double * omega_PM,

          // powder orientation averager
          unsigned int n_orientations,
          double *cosAlpha, 
          double *cosBeta,
          double *amp,
          int nt,

          unsigned int number_of_sites);