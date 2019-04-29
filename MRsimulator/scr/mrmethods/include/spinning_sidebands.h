

#define PI2 6.2831853072
#define PI2I PI2*I

#include "MRAngularMomentum.h"
#include "OCEulerAngle.h"
#include "powder_setup.h"
// #include "csa_static_lineshape.h"
#include <complex.h>
#include <time.h>
#include "math.h"
#include "fftw3.h"
#include "fftw/fftw3_mkl.h"
#include "OCPowderScheme.h"
#include "c_array.h"
#include "histogram.h"
#include "mkl.h"

#define MKL_Complex16 double complex

struct directionCosines
{
	double cosAlpha;
	double cosBeta;
}; 


// Return a vector ordered according to the fft output order.                //
// @params int n - The number of points                                      //
// @params double increment - The increment (sampling interval)              //
// @returns *double values = The pointer to the fft output order vector      //
extern inline double* __get_frequency_in_FFT_order(
                                      int n, 
                                      double increment
                                  );


// static inline void get_frequency(int index,
// 				double cosAlpha, double sinAlpha,
// 				double cosBeta, double sinBeta,
// 				double complex * MR_full_DLM,
// 				double complex * w_cs_PM,
// 				double complex * w_cs,
// 				double complex * rotor_lab,
// 				double * iso,
// 				int * ph_step,
// 				double complex * pre_phase,
// 				fftw_complex * phi,
// 				fftw_complex * side_band,
// 				double complex * zero,
// 				double complex *one,
// 				double complex *i2pi,
// 				fftw_plan *plan,
// 				double **sideband_amplitude,
// 				double *vr_freq,
// 				double **local_vr_freq,
// 				double *amp
// 			);

extern void lineshape_cas_spinning_sideband_cython_angles(
						double * spec,
						double * cpu_time_,
						double frequency_start,
						double frequency_bandwidth,
						double number_of_points,

						double iso,
						double aniso,
						double eta,

						int ph_step,
						double spin_frequency,
						double rotor_angle,
						double * omega_PM_c,

						int averaging_scheme,
						int averaging_size);

extern void lineshape_cas_spinning_sideband_angles(
            // spectrum information and related amplitude
            double * spec, 								// The amplitude of the spectrum.
            double frequency_start, 			// The start of the frequency spectrum.
            double frequency_bandwidth, 	// The bondwidth of the frequency spectrum.
            int number_of_points, 				// Number of points on the frequency spectrum.

            // CSA tensor information
            double iso, 									// The isotropic chemical shift.
            double aniso, 								// The chemical shielding anisotropic.
            double eta, 									// Chemical shielding asymmetry

            // spin rate, spin angle and number spinning sidebands
            int ph_step, 									// The number of spinning sidebands to evaluate
            double spin_frequency, 				// The rotor spin frequency
            double rotor_angle, 					// The rotor angle relative to lab-frame z-axis

            // Euler angle -> principal to molecular frame
            OCEulerAngle omega_PM, 				// The principal to molecular frame transformation euler angles.

            // Euler angles for powder averaging scheme
            OCPowderScheme Omega); 				// Set of euler angles in powder averaging.


extern void spinning_sideband_core(
          // spectrum information and related amplitude
          double * spec,                    // The amplitude of the spectrum.
          double * cpu_time_,               // Execution time
          double spectral_start,            // The start of the frequency spectrum.
          double spectral_increment,        // The bandwidth of the frequency spectrum.
          int number_of_points,             // Number of points on the frequency spectrum.

          double *qunatum_number,           // Spin quantum numbers
          double *wo,                       // omega_o

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

          int nt,
          unsigned int number_of_sites);