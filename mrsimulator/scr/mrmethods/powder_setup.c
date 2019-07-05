
//
//  powder_setup.c
//
//  Created by Deepansh J. Srivastava, Apr 11, 2019
//  Copyright © 2019 Deepansh J. Srivastava. All rights reserved.
//  Contact email = srivastava.89@osu.edu, deepansh2012@gmail.com
//

#include "powder_setup.h"
#include "octahedron.h"

void __powder_averaging_setup(
    int nt, double *cos_alpha, double *cos_beta, double *amp,
    int space // 1 for octant, 2 for hemisphere and 4 for sphere
) {
  // unsigned int n_orientations = 1;
  if (space == 1) { // octant
    octahedronGetPolarAngleTrigOverAnOctant(nt, cos_alpha, cos_beta, amp);
  }
}

// void getDirectionCosineSquareOverHemishpereAndWeights(int nt, double **xr,
//                                                       double **yr, double
//                                                       **zr, double **amp) {

//   int i, j;
//   double x2, y2, z2, r2;

//   /* Do the (x + y + z = nt) face of the octahedron
//   !z -> 0 to nt-1
//   !y -> 0 to nt-z
//   !x -> nt - y - z
//   !*/

//   for (j = 0; j <= nt - 1; j++) {
//     for (i = 0; i <= nt - j; i++) {
//       // x = nt-i-j;
//       // y = i;
//       // z = j;
//       x2 = pow(nt - i - j, 2);
//       y2 = pow(i, 2);
//       z2 = pow(j, 2);
//       r2 = x2 + y2 + z2;
//       xr[i][j] = x2 / r2;
//       yr[i][j] = y2 / r2;
//       zr[i][j] = z2 / r2;
//       amp[i][j] = 1.0 / (r2 * sqrt(r2));
//     }

//     /* Do the (-x + y + z = nt) face of the octahedron
//     !z -> 0 to nt-1
//     !y -> 0 to nt-z
//     !x -> -nt + y + z
//     !*/
//     for (i = nt - j + 1; i <= nt; i++) {
//       // x = nt-i-j;
//       // y = nt-j;
//       // z = nt-i;
//       x2 = pow(nt - i - j, 2);
//       y2 = pow(nt - j, 2);
//       z2 = pow(nt - i, 2);
//       r2 = x2 + y2 + z2;
//       xr[i][j] = x2 / r2;
//       yr[i][j] = y2 / r2;
//       zr[i][j] = z2 / r2;
//       amp[i][j] = 1.0 / (r2 * sqrt(r2));
//     }
//   }

//   /* Do the (-x - y + z = nt) face of the octahedron
//   !*/

//   for (j = nt; j < 2 * nt; j++) {
//     for (i = j - nt + 1; i < nt; i++) {
//       // x = -nt-i+j;
//       // y = nt-j;
//       // z = nt-i;
//       x2 = pow(-nt - i + j, 2); // x*x;
//       y2 = pow(nt - j, 2);      // y*y;
//       z2 = pow(nt - i, 2);      // z*z;
//       r2 = x2 + y2 + z2;
//       xr[i][j] = x2 / r2;
//       yr[i][j] = y2 / r2;
//       zr[i][j] = z2 / r2;
//       amp[i][j] = 1.0 / (r2 * sqrt(r2));
//     }

//     /* Do the (x - y + z = nt) face of the octahedron
//     !*/

//     for (i = 1; i <= j - nt; i++) {
//       // x = -nt-i+j;
//       // y = -i;
//       // z = 2*nt-j;
//       x2 = pow(-nt - i + j, 2); // x*x;
//       y2 = pow(-i, 2);          // y*y;
//       z2 = pow(2 * nt - j, 2);  // z*z;
//       r2 = x2 + y2 + z2;
//       xr[i][j] = x2 / r2;
//       yr[i][j] = y2 / r2;
//       zr[i][j] = z2 / r2;
//       amp[i][j] = 1.0 / (r2 * sqrt(r2));
//     }
//   }

//   xr[0][nt] = 0.0;
//   yr[0][nt] = 0.0;
//   zr[0][nt] = 1.0;
//   r2 = nt;
//   amp[0][nt] = 1.0 / (r2 * r2 * r2);

//   for (j = 0; j < nt; j++) {
//     i = 2 * nt - j;
//     xr[0][i] = xr[0][j];
//     yr[0][i] = yr[0][j];
//     zr[0][i] = zr[0][j];
//     amp[0][i] = amp[0][j];
//   }

//   for (i = 0; i <= nt; i++) {
//     xr[nt][nt + i] = xr[i][0];
//     yr[nt][nt + i] = yr[i][0];
//     zr[nt][nt + i] = zr[i][0];
//     amp[nt][nt + i] = amp[i][0];
//   }

//   i = 2 * nt;
//   for (j = 1; j < nt; j++) {
//     xr[nt - j][i] = xr[nt][j];
//     yr[nt - j][i] = yr[nt][j];
//     zr[nt - j][i] = zr[nt][j];
//     amp[nt - j][i] = amp[nt][j];
//   }
// }

// void getPolarAngleTrigOverHemisphere(int nt, double *cos_alpha, double
// *sinAlpha,
//                                      double *cos_beta, double *sinBeta,
//                                      double **amp) {

//   int i, j;

//   double points = (nt + 1) * (2 * nt + 1);
//   double **xr = createDouble2DMatrix((nt + 1), (2 * nt + 1));
//   double **yr = createDouble2DMatrix((nt + 1), (2 * nt + 1));
//   double **zr = createDouble2DMatrix((nt + 1), (2 * nt + 1));

//   getDirectionCosineSquareOverHemishpereAndWeights(nt, xr, yr, zr, amp);

//   // Evaluate sqrt of zr to get cos(beta)
//   vdSqrt(points, &zr[0][0], &cos_beta[0]);

//   // evaluate A = x + y
//   vdAdd(points, &xr[0][0], &yr[0][0], &sinBeta[0]);
//   // Take sqrt of A to get sin(beta)
//   vdSqrt(points, &sinBeta[0], &sinBeta[0]);

//   // Evaluate sqrt of xr
//   vdSqrt(points, &xr[0][0], &xr[0][0]);

//   // Evaluate sqrt of xr
//   vdSqrt(points, &yr[0][0], &yr[0][0]);

//   int ii = 0;
//   for (i = 0; i <= nt; i++) {
//     for (j = 0; j <= 2 * nt; j++) {
//       // cos_beta[i][j] = sqrt(zr[i][j]);
//       // sinBeta[i][j] = sqrt(xr[i][j] + yr[i][j]);
//       if (sinBeta[ii] != 0.) {
//         cos_alpha[ii] = xr[i][j] / sinBeta[ii];
//         sinAlpha[ii] = yr[i][j] / sinBeta[ii];
//         ii++;
//       } else {
//         cos_alpha[ii] = 1.0;
//         sinAlpha[ii] = 0.0;
//         ii++;
//       }
//     }
//   }
//   destroyDouble2DMatrix(xr);
//   destroyDouble2DMatrix(yr);
//   destroyDouble2DMatrix(zr);
// }

// void powderAverageWithTentingSchemeOverHemisphere2(double *spec, double
// *freq,
//                                                    int nt, double *amp,
//                                                    double *offset, int m) {

//   int i = 0, j = 0, local_index, n_pts = (nt + 1) * (nt + 2) / 2;
//   double amp1 = 0.0, temp;
//   double *amp_address, *freq_address;

//   /* Interpolate between frequencies by setting up tents */
//   int l = nt;
//   local_index = nt - 1;
//   amp_address = &amp[4 * nt];
//   freq_address = &freq[4 * nt];

//   while (i < n_pts - 1) {
//     temp = amp[i + 1] + amp_address[j];
//     amp1 = temp;
//     amp1 += amp[i];

//     triangle_interpolation(&freq[i], &freq[i + 1], &freq_address[j], offset,
//                            &amp1, spec, &m);

//     if (i < local_index) {
//       amp1 = temp;
//       if (j == 4 * l) {
//         amp1 += amp_address[j - 4 * l];
//         triangle_interpolation(&freq[i + 1], &freq_address[j],
//         &freq_address[0],
//                                offset, &amp1, spec, &m);
//       } else
//         amp1 += amp_address[j + 1];
//       triangle_interpolation(&freq[i + 1], &freq_address[j],
//                              &freq_address[j + 1], offset, &amp1, spec, &m);
//       i++;
//       j++;
//     } else {
//       local_index = i + nt;
//       i++;
//     }
//   }
// }

// void powderAverageWithTentingSchemeOverHemisphere(double *spec,
//                                                   double **powfreq, int nt,
//                                                   double **amp, double
//                                                   *offset, int m) {

//   int i, j;
//   double amp1, amp2, temp;

//   /* Interpolate between frequencies by setting up tents */

//   for (i = 0; i <= nt - 1; i++) {
//     for (j = 0; j <= nt - 1; j++) {
//       temp = amp[i][j + 1] + amp[i + 1][j];
//       amp1 = amp[i][j] + temp;
//       amp2 = temp + amp[i + 1][j + 1];

//       triangle_interpolation(&powfreq[i + 1][j], &powfreq[i][j + 1],
//                              &powfreq[i][j], offset, &amp1, spec, &m);
//       triangle_interpolation(&powfreq[i + 1][j], &powfreq[i][j + 1],
//                              &powfreq[i + 1][j + 1], offset, &amp2, spec,
//                              &m);

//       // tent_amp(&powfreq[i+1][j], &powfreq[i][j+1], &powfreq[i][j], offset, \
//       //          &amp[i+1][j], &amp[i][j+1], &amp[i][j], spec, m);
//       // tent_amp(&powfreq[i+1][j], &powfreq[i][j+1], &powfreq[i+1][j+1], offset, \
//       //          &amp[i+1][j], &amp[i][j+1], &amp[i+1][j+1], spec, m);
//     }
//   }

//   // if (octa == 0){
//   for (i = 0; i <= nt - 1; i++) {
//     for (j = nt; j <= 2 * nt - 1; j++) {
//       temp = amp[i][j] + amp[i + 1][j + 1];
//       amp1 = temp + amp[i + 1][j];
//       amp2 = temp + amp[i][j + 1];

//       triangle_interpolation(&powfreq[i][j], &powfreq[i + 1][j + 1],
//                              &powfreq[i + 1][j], offset, &amp1, spec, &m);
//       triangle_interpolation(&powfreq[i][j], &powfreq[i + 1][j + 1],
//                              &powfreq[i][j + 1], offset, &amp2, spec, &m);

//       // tent_amp(&powfreq[i][j], &powfreq[i+1][j+1], &powfreq[i][j+1], offset, \
//         //          &amp[i][j], &amp[i+1][j+1], &amp[i][j+1], spec, m);
//       // tent_amp(&powfreq[i][j], &powfreq[i+1][j+1], &powfreq[i][j+1], offset, \
//         //          &amp[i][j], &amp[i+1][j+1], &amp[i][j+1], spec, m);
//     }
//   }
//   // }
// };

// static inline double two_times_triangle_area(
//                         double *a,
//                         double *b,
//                         double *c)
// {
//     return ((b[0]-a[0])*(c[1]-a[1]) - (b[1]-a[1])*(c[0]-a[0]));
// }
