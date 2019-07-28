
//
//  octahedron.c
//
//  Created by Deepansh J. Srivastava, Apr 11, 2019
//  Copyright © 2019 Deepansh J. Srivastava. All rights reserved.
//  Contact email = srivastava.89@osu.edu, deepansh2012@gmail.com
//

#include "octahedron.h"

void octahedronGetDirectionCosineSquareOverOctantAndWeights(int nt, double *xr,
                                                            double *yr,
                                                            double *zr,
                                                            double *amp) {

  int i, j, k = 0;
  double x2, y2, z2, r2, scale = 1.0;

  /* Do the (x + y + z = nt) face of the octahedron
  !z -> 0 to nt-1
  !y -> 0 to nt-z
  !x -> nt - y - z
  !*/

  for (j = 0; j <= nt - 1; j++) {
    for (i = 0; i <= nt - j; i++) {
      // x = nt-i-j;
      // y = i;
      // z = j;
      x2 = pow(nt - i - j, 2);
      y2 = pow(i, 2);
      z2 = pow(j, 2);
      r2 = x2 + y2 + z2;
      xr[k] = x2 / r2;
      yr[k] = y2 / r2;
      zr[k] = z2 / r2;
      amp[k] = scale / (r2 * sqrt(r2));
      k++;
    }
  }

  xr[k] = 0.0;
  yr[k] = 0.0;
  zr[k] = 1.0;
  r2 = (double)nt;
  amp[k] = scale / (r2 * r2 * r2);
}

void octahedronGetPolarAngleTrigOverAnOctant(int nt, double *cos_alpha,
                                             double *cos_beta, double *amp) {

  int points = (nt + 1) * (nt + 2) / 2;
  double *xr = malloc_double(points);
  double *yr = malloc_double(points);
  double *zr = malloc_double(points);
  double *sin_beta = malloc_double(points);

  octahedronGetDirectionCosineSquareOverOctantAndWeights(nt, xr, yr, zr, amp);

  // Evaluate sqrt of zr to get cos(beta)
  vmd_sqrt(points, zr, cos_beta);

  // Evaluate A = x + y
  vmd_add(points, xr, yr, sin_beta);

  // Take sqrt of A to get sin(beta)
  vmd_sqrt(points, sin_beta, sin_beta);

  // Evaluate sqrt of xr
  vmd_sqrt(points, xr, xr);

  // Evaluate sqrt of xr
  // vmd_sqrt(points, &yr[0], &yr[0]);

  vmd_div(points - 1, xr, sin_beta, cos_alpha);

  // vmd_div(points-1, yr, sinBeta, sinAlpha );

  cos_alpha[points - 1] = 1.0;
  // sinAlpha[points-1] = 0.0;

  free_double(xr);
  free_double(yr);
  free_double(zr);
  free_double(sin_beta);
}

void octahedronGetPolarAngleCosineAzimuthalAnglePhaseOverOctant(
    int nt, complex128 *exp_I_alpha, double *cos_beta, double *amp) {

  int points = (nt + 1) * (nt + 2) / 2;
  double *xr = malloc_double(points);
  double *yr = malloc_double(points);
  double *zr = malloc_double(points);
  double *sin_beta = malloc_double(points);

  // The values xr = x^2, yr = y^2, zr = z^2, where x, y, and z are the
  // direction cosines.
  octahedronGetDirectionCosineSquareOverOctantAndWeights(nt, xr, yr, zr, amp);

  // Evaluate sqrt of zr to get cos(beta)
  // cos(beta) = sqrt(z^2)
  vmd_sqrt(points, zr, cos_beta);

  // Evaluate A = x + y
  // sin^2(beta) = x^2 + y^2
  vmd_add(points, xr, yr, sin_beta);
  // Take sqrt of A to get sin(beta)
  // sin(beta) = sqrt(x^2 + y^2)
  vmd_sqrt(points, sin_beta, sin_beta);

  // Evaluate sqrt of xr
  // value of xr is updated to sqrt(x^2)
  vmd_sqrt(points, xr, xr);

  // Evaluate sqrt of yr
  // value of yr is updated to sqrt(y^2)
  vmd_sqrt(points, yr, yr);

  // Evaluate cos_alpha = x/sqrt(x^2 + y^2) = xr/sin_beta
  // value of xr is updated to hold cos_alpha
  vmd_div(points - 1, xr, sin_beta, xr);

  // Evaluate sin_alpha = y/sqrt(x^2 + y^2) = yr/sin_beta
  // value of yr is updated to hold sin_alpha
  vmd_div(points - 1, yr, sin_beta, yr);

  xr[points - 1] = 1.0;
  yr[points - 1] = 0.0;

  // copy cos_alpha and sin_alpha to alternating address of exp_I_alpha
  // to emulate a complex array.
  cblas_dcopy(points, xr, 1, (double *)&exp_I_alpha[0], 2);
  cblas_dcopy(points, yr, 1, (double *)&exp_I_alpha[0] + 1, 2);

  free_double(xr);
  free_double(yr);
  free_double(zr);
  free_double(sin_beta);
}

void octahedronInterpolation(double *spec, double *freq, int nt, double *amp,
                             int stride, int m) {

  int i = 0, j = 0, local_index, n_pts = (nt + 1) * (nt + 2) / 2;
  unsigned int int_i_stride, int_j_stride;
  double amp1 = 0.0, temp;
  double *amp_address, *freq_address;

  /* Interpolate between 1d points by setting up triangles of unit area */

  local_index = nt - 1;
  amp_address = &amp[(nt + 1) * stride];
  freq_address = &freq[nt + 1];

  while (i < n_pts - 1) {
    int_i_stride = i * stride;
    int_j_stride = j * stride;
    temp = amp[int_i_stride + stride] + amp_address[int_j_stride];
    amp1 = temp;
    amp1 += amp[int_i_stride];

    triangle_interpolation(&freq[i], &freq[i + 1], &freq_address[j], &amp1,
                           spec, &m);

    if (i < local_index) {
      amp1 = temp;
      amp1 += amp_address[int_j_stride + stride];
      triangle_interpolation(&freq[i + 1], &freq_address[j],
                             &freq_address[j + 1], &amp1, spec, &m);
    } else {
      local_index = j + nt;
      i++;
    }
    i++;
    j++;
  }
}
