// -*- coding: utf-8 -*-
//
//  spherical.c
//
//  @copyright Deepansh J. Srivastava, 2019-2021.
//  Created by Deepansh J. Srivastava, Feb 12, 2021
//  Contact email = srivastava.89@osu.edu
//

#include "spherical.h"

void sphericalGetCoordinates(int nt, double *restrict x, double *restrict y,
                             double *restrict z) {
  int i, j;
  double temp1, temp2, scale = (double)nt, x_o, y_o, z_o;

  for (j = 0; j <= nt - 1; j++) {
    for (i = 0; i <= nt - j; i++) {
      z_o = (double)j / scale;
      y_o = (double)i / scale;
      x_o = (double)(nt - i - j) / scale;
      temp1 = PIby2 * (y_o / (y_o + x_o));
      *z = 2.0 * z_o * (1.0 - 0.5 * z_o);
      temp2 = sqrt(1 - *z * *z);
      *x++ = temp2 * cos(temp1);
      *y++ = temp2 * sin(temp1);
      z++;
    }
  }
  *x = 0.0;
  *y = 0.0;
  *z = 1.0;
}

void sphericalGetComplexExpOfPolarAngleOverOctant(int nt, void *exp_I_alpha,
                                                  void *exp_I_beta,
                                                  double *amp) {
  /** Uniform area trioangle coordinates over the Octahedron.
   *
   *      z_o = (double)j / scale                                (1a)
   *      y_o = (double)i / scale                                (1b)
   *      x_o = (double)(nt-i-j)  / scale                        (1c)
   *
   * Area preservong mapping of the (x_o, y_o, z_o) coordinates over the
   * octahedron to (x, y, z) over the sphere.
   *
   *      z = 2.0 * z_o * (1.0 - 0.5 * z_o)                      (2a)
   *      x = sqrt(1 - z^2) cos( (π/2) (y_o / (y_o + x_o)) )     (2b)
   *      y = sqrt(1 - z^2) sin( (π/2) (y_o / (y_o + x_o)) )     (2c)
   *
   * Cosines of alpha and beta from the (x, y, z) coordinates
   *
   *      cos_beta = z                                           (3a)
   *      sin_beta = sqrt(1 - z^2)                               (3b)
   *      cos_alpha = sqrt(x^2 / (x^2 + y^2))                    (3c)
   *      sin_alpha = sqrt(y^2 / (x^2 + y^2))                    (3d)
   *
   * From Eqs. (2b), (2c), and (3b)
   *      x^2 + y^2 = 1 - z^2 = (sin_beta)^2                     (4)
   *
   * Substituting Eqs. (4) and (2b) into (3c) gives
   *      cos_alpha = x / sin_beta                               (5a)
   *      cos_alpha = sin_beta cos( arg ) / sin_beta             (5b)
   *      cos_alpha = cos( arg )                                 (5b)
   *
   * where,
   *      arg = (π/2) (y_o / (y_o + x_o))                        (6a)
   *          = (π/2) ( i / (nt - j) )                           (6b)
   *
   * Similarly,
   *      sin_alpha = y / sin_beta                               (7a)
   *      sin_alpha = sin_beta sin( arg ) / sin_beta             (7b)
   *      sin_alpha = sin( arg )                                 (7c)
   **/

  int i, j, points = (nt + 1) * (nt + 2) / 2, index;
  double temp, arg, z, scale = (double)nt;
  double *exp_I_alpha_ = (double *)exp_I_alpha;
  double *exp_I_beta_ = (double *)exp_I_beta;

  vm_double_ones(points, amp);

  for (j = 0; j <= nt - 1; j++) {
    for (i = 0; i <= nt - j; i++) {
      z = (double)j / scale;    // z-coordinate over octahedron
      temp = (double)(nt - j);  // x + y

      z = 2.0 * z * (1.0 - 0.5 * z);  // z-coordinate over sphere
      *exp_I_beta_++ = z;
      *exp_I_beta_++ = sqrt(1.0 - z * z);

      arg = PIby2 * (double)i / temp;
      *exp_I_alpha_++ = cos(arg);
      *exp_I_alpha_++ = sin(arg);
    }
  }
  *exp_I_beta_++ = 1.0;
  *exp_I_beta_ = 0.0;
  *exp_I_alpha_++ = 1.0;
  *exp_I_alpha_ = 0.0;
}
