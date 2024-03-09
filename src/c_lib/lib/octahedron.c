// -*- coding: utf-8 -*-
//
//  octahedron.c
//
//  @copyright Deepansh J. Srivastava, 2019-2021.
//  Created by Deepansh J. Srivastava, Apr 11, 2019.
//  Contact email = srivastava.89@osu.edu
//

#include "octahedron.h"

void octahedronGetDirectionCosineSquareAndWeightsOverOctant(const unsigned int nt,
                                                            double *restrict xr,
                                                            double *restrict yr,
                                                            double *restrict zr,
                                                            double *restrict amp) {
  unsigned int i, j;
  // scale is divided by pi to remove angular dependence.
  // The factor 6 is because each point in the triangle interpolation is shared by
  // 6 neighboring triangles.
  // Additional fudge factor of 2 (overall factor is now 3) to match zcw, step surface
  // area.
  double x2, y2, z2, r2, factor = 3.0 * CONST_PI, scale = (double)nt / factor;

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
      *xr++ = x2 / r2;
      *yr++ = y2 / r2;
      *zr++ = z2 / r2;
      *amp++ = scale / (r2 * sqrt(r2));
    }
  }

  *xr = 0.0;
  *yr = 0.0;
  *zr = 1.0;
  r2 = (double)nt;
  *amp = 1.0 / (r2 * r2 * factor);
}

void octahedronGetPolarAngleTrigOverOctant(const unsigned int nt,
                                           double *restrict cos_alpha,
                                           double *restrict cos_beta,
                                           double *restrict amp) {
  unsigned int points = (nt + 1) * (nt + 2) / 2;
  double *xr = malloc_double(points);
  double *yr = malloc_double(points);
  double *zr = malloc_double(points);
  double *sin_beta = malloc_double(points);

  // The values xr = x^2, yr = y^2, zr = z^2, where x, y, and z are the
  // direction cosines.
  octahedronGetDirectionCosineSquareAndWeightsOverOctant(nt, xr, yr, zr, amp);

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

void octahedronGetComplexExpOfPolarAngleOverOctant(const unsigned int nt,
                                                   void *restrict exp_I_alpha,
                                                   void *restrict exp_I_beta,
                                                   double *restrict amp) {
  unsigned int points = (nt + 1) * (nt + 2) / 2;
  double *xr = malloc_double(points);
  double *yr = malloc_double(points);
  double *zr = malloc_double(points);

  octahedronGetDirectionCosineSquareAndWeightsOverOctant(nt, xr, yr, zr, amp);
  // At this point the variables
  //    xr = x^2,
  //    yr = y^2,
  //    zr = z^2,
  // where x, y, and z are the direction cosines.

  // Cos beta .............................................................. //
  // cos(beta) = sqrt(z^2) ... (1)
  // In terms of the variables, Eq (1) is given as sqrt(zr).
  // Evaluate zr = sqrt(zr)         ==>> sqrt(z^2).
  vm_double_square_root_inplace(points, zr);
  // Copy cos(beta), aka zr, to the even addresses of exp_I_beta.
  cblas_dcopy(points, zr, 1, (double *)exp_I_beta, 2);

  // Sin beta .............................................................. //
  // sin^2(beta) = x^2 + y^2 ... (2)
  // In terms of the variables, Eq (2) is given as xr + yr.
  // Evaluate zr = xr + yr.
  vm_double_add(points, xr, yr, zr);
  // Take the square root of zr to get sin(beta).
  // Evaluate zr = sqrt(zr)         ==>> sqrt(x^2 + y^2).
  vm_double_square_root_inplace(points, zr);
  // Copy sin(beta), aka zr, to the odd addresses of exp_I_beta.
  cblas_dcopy(points, zr, 1, (double *)exp_I_beta + 1, 2);

  // Evaluate squate root of xr
  // xr = sqrt(xr)        ==>> sqrt(x^2).
  vm_double_square_root_inplace(points, xr);

  // Evaluate squate root of yr
  // yr = sqrt(yr)        ==>> sqrt(y^2).
  vm_double_square_root_inplace(points, yr);

  // .. note
  // At this point the variables,
  //    xr = x,
  //    yr = y,
  //    zr = sin(beta) = sqrt(x^2 + y^2)

  // Cos alpha ............................................................. //
  // cos(alpha) = x/sqrt(x^2 + y^2) = x/sin(beta) ... (3)
  // In terms of the variables, Eq (3) is given as xr/zr.
  // Evaluate xr = xr/zr
  vm_double_divide_inplace(points - 1, zr, xr);
  xr[points - 1] = 0.0;
  // Copy cos(alpha), aka xr, to the even addresses of exp_I_alpha
  cblas_dcopy(points, xr, 1, (double *)exp_I_alpha, 2);

  // Sin alpha ............................................................. //
  // sin(alpha) = y/sqrt(x^2 + y^2) = y/sin(beta) ... (4)
  // In terms of the variables, Eq (4) is given as yr/zr.
  // Evaluate yr = yr/zr
  vm_double_divide_inplace(points - 1, zr, yr);
  yr[points - 1] = 0.0;
  // Copy sin(alpha), aka yr, to the odd addresses of exp_I_alpha.
  cblas_dcopy(points, yr, 1, (double *)exp_I_alpha + 1, 2);

  // free temporary memory
  free(xr);
  free(yr);
  free(zr);
}

void get_total_amplitude(const unsigned int nt, double *amp, double *amp_sum) {
  unsigned int i = 0, j = 0, local_index = nt - 1, n_pts = (nt + 1) * (nt + 2) / 2;
  double *amp_address, temp;

  amp_address = &amp[nt + 1];

  *amp_sum = 0.0;
  while (i < n_pts - 1) {
    temp = amp[i + 1] + amp_address[j];
    *amp_sum += temp + amp[i];

    if (i < local_index) {
      temp += amp_address[j + 1];
      *amp_sum += temp;
    } else {
      local_index = j + nt;
      i++;
    }
    i++;
    j++;
  }
}

// fix amplitude for binning case:
// octant face edge points are divided by 2
// octant vertexes are divided by 4
void fix_amplitude_for_binning(unsigned int nt, double *amp) {
  unsigned int i = 0, i_max = (nt + 1) * (nt + 2) / 2, tt = nt - 1;
  bool inc_one = false;
  while (i <= nt) amp[i++] /= 2.0;
  while (i < i_max) {
    amp[i] /= 2.0;
    if (inc_one) {
      i++;
      inc_one = false;
    } else {
      i += tt--;
      inc_one = true;
    }
  }
  amp[0] /= 2.0;
  amp[nt] /= 2.0;
  amp[i_max - 1] /= 2.0;

  // rescale all amps by 6 to remove the 1point shared by 6 triangle factor
  cblas_dscal(i_max, 6.0, amp, 1);
}

void averaging_setup(unsigned int nt, void *exp_I_alpha, void *exp_I_beta, double *amp,
                     bool interpolation) {
  // octahedronGetPolarAngleTrigOverOctant(nt, cos_alpha, cos_beta, amp);
  octahedronGetComplexExpOfPolarAngleOverOctant(nt, exp_I_alpha, exp_I_beta, amp);
  if (!interpolation) {
    fix_amplitude_for_binning(nt, amp);
  }
}
