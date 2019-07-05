
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
  double *xr = createDouble1DArray(points);
  double *yr = createDouble1DArray(points);
  double *zr = createDouble1DArray(points);
  double *sin_beta = createDouble1DArray(points);

  octahedronGetDirectionCosineSquareOverOctantAndWeights(nt, xr, yr, zr, amp);

  // Evaluate sqrt of zr to get cos(beta)
  vdSqrt(points, &zr[0], &cos_beta[0]);

  // Evaluate A = x + y
  vdAdd(points, &xr[0], &yr[0], &sin_beta[0]);
  // Take sqrt of A to get sin(beta)
  vdSqrt(points, &sin_beta[0], &sin_beta[0]);

  // Evaluate sqrt of xr
  vdSqrt(points, &xr[0], &xr[0]);

  // Evaluate sqrt of xr
  // vdSqrt(points, &yr[0], &yr[0]);

  vdDiv(points - 1, xr, sin_beta, cos_alpha);
  // vdDiv(points-1, yr, sinBeta, sinAlpha );

  cos_alpha[points - 1] = 1.0;
  // sinAlpha[points-1] = 0.0;

  destroyDouble1DArray(xr);
  destroyDouble1DArray(yr);
  destroyDouble1DArray(zr);
  destroyDouble1DArray(sin_beta);
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
