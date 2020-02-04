// -*- coding: utf-8 -*-
//
//  octahedron.c
//
//  @copyright Deepansh J. Srivastava, 2019-2020.
//  Created by Deepansh J. Srivastava, Apr 11, 2019
//  Contact email = deepansh2012@gmail.com
//

#include "octahedron.h"

void octahedronGetDirectionCosineSquareOverOctantAndWeights(int nt, double *xr,
                                                            double *yr,
                                                            double *zr,
                                                            double *amp) {
  int i, j, k = 0;
  double x2, y2, z2, r2, scale = (double)nt;

  /* Do the (x + y + z = nt) face of the octahedron
  !z -> 0 to nt-1
  !y -> 0 to nt-z
  !x -> nt - y - z
  !*/

  for (j = 0; j <= nt - 1; j++) {
    for (i = 0; i <= nt - j; i++) {
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

  // The values xr = x^2, yr = y^2, zr = z^2, where x, y, and z are the
  // direction cosines.
  octahedronGetDirectionCosineSquareOverOctantAndWeights(nt, xr, yr, zr, amp);

  // Evaluate sqrt of zr to get cos(beta) ---> cos(beta) = sqrt(z^2)
  vm_double_square_root(points, zr, cos_beta);

  // Evaluate zr = xr + yr ---> sin^2(beta) = x^2 + y^2. Here zr is sin^2(beta).
  vm_double_add(points, xr, yr, zr);

  // Take sqrt of zr to get sin(beta) ---> sin(beta) = sqrt(x^2 + y^2)
  vm_double_square_root(points, zr, sin_beta);

  // Evaluate sqrt of xr and set the value to zr, zr = sqrt(x^2) = sqrt(xr)
  vm_double_square_root(points, xr, zr);

  // Evaluate sqrt of xr
  // vdSqrt(points, &yr[0], &yr[0]);

  // Evaluate cos_alpha = x/sqrt(x^2 + y^2) = zr/sin_beta
  vm_double_divide(points - 1, zr, sin_beta, cos_alpha);
  // vdDiv(points-1, yr, sinBeta, sinAlpha );

  cos_alpha[points - 1] = 1.0;
  // sinAlpha[points-1] = 0.0;

  free(xr);
  free(yr);
  free(zr);
  free(sin_beta);
}

void octahedronGetPolarAngleCosineAzimuthalAnglePhaseOverOctant(
    int nt, void *exp_I_alpha, void *exp_I_beta, double *amp) {
  int points = (nt + 1) * (nt + 2) / 2;
  double *xr = malloc_double(points);
  double *yr = malloc_double(points);
  double *zr = malloc_double(points);

  octahedronGetDirectionCosineSquareOverOctantAndWeights(nt, xr, yr, zr, amp);
  // At this point the variables
  //    xr = x^2,
  //    yr = y^2,
  //    zr = z^2,
  // where x, y, and z are the direction cosines.

  // Cos beta .............................................................. //
  // cos(beta) = sqrt(z^2) ... (1)
  // In terms of the variables, Eq (1) is given as sqrt(zr)
  // Evaluate zr = sqrt(zr)         ==>> sqrt(z^2)
  vm_double_square_root(points, zr, zr);
  // Copy cos(beta), aka zr, to the even addresses of exp_I_beta
  cblas_dcopy(points, zr, 1, (double *)exp_I_beta, 2);

  // Sin beta .............................................................. //
  // sin^2(beta) = x^2 + y^2 ... (2)
  // In terms of the variables, Eq (2) is given as xr + yr
  // Evaluate zr = xr + yr.
  vm_double_add(points, xr, yr, zr);
  // Take the square root of zr to get sin(beta).
  // Evaluate zr = sqrt(zr)         ==>> sqrt(x^2 + y^2)
  vm_double_square_root(points, zr, zr);
  // Copy sin(beta), aka zr, to the odd addresses of exp_I_beta
  cblas_dcopy(points, zr, 1, (double *)exp_I_beta + 1, 2);

  // Evaluate squate root of xr
  // xr = sqrt(xr)        ==>> sqrt(x^2)
  vm_double_square_root(points, xr, xr);

  // Evaluate squate root of yr
  // yr = sqrt(yr)        ==>> sqrt(y^2)
  vm_double_square_root(points, yr, yr);

  // .. note
  // At this point the variables,
  //    xr = x,
  //    yr = y,
  //    zr = sin(beta) = sqrt(x^2 + y^2)

  // Cos alpha ............................................................. //
  // cos(alpha) = x/sqrt(x^2 + y^2) = x/sin(beta) ... (3)
  // In terms of the variables, Eq (3) is given as xr/zr
  // Evaluate xr = xr/zr
  vm_double_divide(points - 1, xr, zr, xr);
  xr[points - 1] = 0.0;
  // Copy cos(alpha), aka xr, to the even addresses of exp_I_alpha
  cblas_dcopy(points, xr, 1, (double *)exp_I_alpha, 2);

  // Sin alpha ............................................................. //
  // sin(alpha) = y/sqrt(x^2 + y^2) = y/sin(beta) ... (4)
  // In terms of the variables, Eq (4) is given as yr/zr
  // Evaluate yr = yr/zr
  vm_double_divide(points - 1, yr, zr, yr);
  yr[points - 1] = 0.0;
  // Copy sin(alpha), aka yr, to the odd addresses of exp_I_alpha
  cblas_dcopy(points, yr, 1, (double *)exp_I_alpha + 1, 2);

  // free temporary memory
  free(xr);
  free(yr);
  free(zr);
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
