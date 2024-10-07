// -*- coding: utf-8 -*-
//
//  interpolation.c
//
//  @copyright Deepansh J. Srivastava, 2019-2021.
//  Created by Deepansh J. Srivastava, Apr 11, 2019.
//  Contact email = srivastava.89@osu.edu
//

#include "interpolation.h"

/**
 * ================================================================================== *
 *                         One-dimensional interpolation                              *
 * ================================================================================== *
 */

/**
 * @brief Linear interpolation scheme for isotropic values.
 *
 * @param freq The frequency coordinate.
 * @param points Total number of points in spec.
 * @param amp The area corresponding to the frequency coordinate.
 * @param spec1 A pointer to the intensity vector.
 */
static void inline delta_fn_linear_interpolation(const double *freq, const int *points,
                                                 double *amp, double *spec1) {
  int p = (int)*freq;
  if (p >= *points || p < 0) return;

  double diff, delta;
  bool left;
  diff = *freq - (double)p;
  delta = diff - 0.5;
  double *spec = &spec1[2 * p];

  // Do not interpolate the intensity if the difference < 1.0e-6.
  // Ensures that the sideband dimension frequencies do not get interpolated.
  if (absd(delta) < TOL) {
    *spec += *amp;
    return;
  }

  // Linear interpolation
  left = delta < 0;
  if (left) {
    if (p != 0) *(spec - 2) -= *amp * delta;
    *spec += *amp * (1.0 + delta);
  } else {
    if (p + 1 != *points) *(spec + 2) += *amp * delta;
    *spec += *amp * (1.0 - delta);
  }
}

static void inline delta_fn_gauss_interpolation(const double *freq, const int *points,
                                                double *amp, double *spec) {
  double res, w, a0, a1, a2, a3, a4, sum, temp;
  int p = (int)(floor(*freq)), index, ip2, ip1, im1, im2, pad = 2, p2 = 2 * p;
  if (p >= *points + pad || p < -pad + 1) return;

  // for sideband delta freq. It avoids round-off errors.
  res = *freq - (double)p;
  if (absd(res - 0.5) < TOL && p >= 0 && p < *points) {
    spec[p2] += *amp;
    return;
  }

  p = (int)(floor(*freq - 0.5));
  p2 = 2 * p;
  res = *freq - (double)p - 0.5;
  res *= gauss_table_precision_inverse;
  index = (int)res;
  w = res - (double)index;

  ip2 = 2 * gauss_table_precision_inverse - index;
  ip1 = gauss_table_precision_inverse - index;
  im1 = gauss_table_precision_inverse + index;
  im2 = 2 * gauss_table_precision_inverse + index;

  a0 = lerp_plus(w, im2);
  a1 = lerp_plus(w, im1);
  a2 = lerp_plus(w, index);
  a3 = lerp_minus(w, ip1);
  a4 = lerp_minus(w, ip2);

  sum = a0;
  sum += a1;
  sum += a2;
  sum += a3;
  sum += a4;

  temp = *amp / sum;

  // double *spec = &spec1[2 * p];
  if (p - 2 >= 0 && p - 2 < *points) spec[p2 - 4] += temp * a0;
  if (p - 1 >= 0 && p - 1 < *points) spec[p2 - 2] += temp * a1;
  if (p >= 0 && p < *points) spec[p2] += temp * a2;
  if (p + 1 >= 0 && p + 1 < *points) spec[p2 + 2] += temp * a3;
  if (p + 2 >= 0 && p + 2 < *points) spec[p2 + 4] += temp * a4;
}

/**
 * @brief Get the clipping conditions object
 *
 * @param p Start index of the triangle.
 * @param pmid Apex index of the triangle.
 * @param pmax End index of the triangle.
 * @param points Max attainable index.
 * @param clips Pointer to bool array to store the clip conditions.
 */
static inline void get_clipping_conditions(int *p, int *pmid, int *pmax, int *points,
                                           bool *clips) {
  *clips = *p <= 0;  // left triangle left clip
  *p = *clips++ ? 0 : *p;

  *clips = *pmid > *points - 1;  // left triangle right clip
  *pmid = *clips++ ? *points - 1 : *pmid;

  *clips = *pmid <= 0;  // right triangle left clip
  *pmid = *clips++ ? 0 : *pmid;

  *clips = *pmax > *points - 1;  // right triangle right clip
  *pmax = *clips ? *points - 1 : *pmax;
}

/**
 * @brief Interpolation scheme of a right-angled triangle with clipping. (left face)
 *
 * @param p Start bin index for triangle interpolation.
 * @param pmid End bin index for triangle interpolation.
 * @param l_clip Left clip boolean.
 * @param r_clip Right clip booleaan.
 * @param top The height of the triangle.
 * @param f Pointer to the 1D coordinates of the triangle.
 * @param spec1 Pointer to a vector where the interpolated intensities are added.
 */
// A pictorial representation of the interpolation scheme.
//
//        /|                  /|                  / |                       / |      *
//       /*|                /**|                /|  |                    /    |      *
//      /**|              /****|              /**|  |                 /*|     |      *
//     /***|            /|*****|            /****|  |              /|***|     |      *
//    /****|          /__|*****|          /******|__|           /___|***|_____|      *
//   p    pmid           p    pmid       p     pmid                 p  pmid          *
//  no clipping        left clip           right clip         left and right clip    *
static inline void left_triangle_interpolate(int p, int pmid, bool l_clip, bool r_clip,
                                             bool r_clip_r, double top, double *f,
                                             double *spec1) {
  double f10 = f[1] - f[0], df1 = top / f10, diff;
  double *spec = &spec1[2 * p];

  if (p == pmid) {
    *spec += (!r_clip && !l_clip) ? f10 * top * 0.5 : 0.0;
    // *spec += (!r_clip_r && l_clip) ? f[1] * (f10 - f[0]) * 0.5 * df1 : 0.0;
    return;
  }

  diff = (double)p + 1.0 - f[0];
  *spec += (l_clip && f[0] <= 0.0) ? (diff - 0.5) * df1 : 0.5 * diff * diff * df1;

  diff -= 0.5;
  diff *= df1;
  while (++p < pmid) *(spec += 2) += diff += df1;

  double df = (double)p;
  *(spec += 2) += (!r_clip) ? (f[1] - df) * (f10 + df - f[0]) * 0.5 * df1 : diff + df1;
}

/**
 * @brief Interpolation scheme of a right-angled triangle with clipping. (right face)
 *
 * @param p Start bin index for triangle interpolation.
 * @param pmax End bin index for triangle interpolation.
 * @param l_clip Left clip boolean.
 * @param r_clip Right clip booleaan.
 * @param top The height of the triangle.
 * @param f Pointer to the 1D coordinates of the triangle.
 * @param spec1 Pointer to a vector where the interpolated intensities are added.
 */
// A pictorial representation of the interpolation scheme.
//
//   |\               |\                  | \                   | \                   *
//   |*\              |**\                |  |\                 |    \                *
//   |**\             |****\              |  |**\               |     |*\             *
//   |***\            |*****|\            |  |****\             |     |***|\          *
//   |****\           |*****|__\          |__|******\           |_____|***|___\       *
//   p   pmax         p       pmax        p        pmax         p            pmax     *
//  no clipping        right clip           left clip         left and right clip     *
//
static inline void right_triangle_interpolate(int p, int pmax, bool l_clip, bool r_clip,
                                              bool r_clip_l, double top, double *f,
                                              double *spec1) {
  double f21 = f[2] - f[1], df1 = top / f21;
  double *spec = &spec1[2 * p];
  double diff = f[2] - (double)p - 1.0;

  if (p == pmax) {
    *spec += (!r_clip) ? f21 * top * 0.5 : 0.0;
    *spec += (r_clip && !r_clip_l) ? (f21 - diff) * (diff + f21) * 0.5 * df1 : 0.0;
    return;
  }

  *spec += (!l_clip) ? (f21 - diff) * (diff + f21) * 0.5 * df1 : (diff + 0.5) * df1;

  diff += 0.5;
  diff *= df1;
  while (++p < pmax) *(spec += 2) += diff -= df1;

  double df = f[2] - (double)p;
  *(spec += 2) += (!r_clip) ? df * df * 0.5 * df1 : diff - df1;
}

static inline void __triangle_interpolation(double *freq1, double *freq2, double *freq3,
                                            double *amp, double *spec, int *points) {
  int p, pmid, pmax, i, j;
  double top, t;
  p = (int)(*freq1);

  // check if the three points lie within a bin interval.
  if (p == (int)(*freq2) && p == (int)(*freq3)) {
    if (p >= *points || p < 0) return;
    spec[2 * p] += *amp;
    return;
  }

  double f[3] = {freq1[0], freq2[0], freq3[0]};

  // arrange the numbers in ascending order (sort)

  for (j = 1; j <= 2; j++) {
    t = f[j];
    i = j - 1;
    while (i >= 0 && f[i] > t) {
      f[i + 1] = f[i];
      i--;
    }
    f[i + 1] = t;
  }

  // if min frequency is higher than the last bin, return
  p = (int)f[0];
  if (p >= *points) return;

  // if max frequency is lower than the first bin, return
  pmax = (int)f[2];
  if (pmax < 0) return;

  pmid = (int)f[1];
  top = *amp * 2.0 / (f[2] - f[0]);

  bool clips[4];
  get_clipping_conditions(&p, &pmid, &pmax, points, clips);

  // printf("\n%d %d %d %d", clips[0], clips[1], clips[2], clips[3]);
  if (f[1] >= 0.0) {
    left_triangle_interpolate(p, pmid, clips[0], clips[1], clips[2], top, f, spec);
  }
  if (f[2] >= 0.0) {
    right_triangle_interpolate(pmid, pmax, clips[2], clips[3], clips[1], top, f, spec);
  }
}

void triangle_interpolation1D(double *freq1, double *freq2, double *freq3, double *amp,
                              double *spec, int *points, unsigned int iso_intrp) {
  if (absd(*freq1 - *freq2) < TOL && absd(*freq1 - *freq3) < TOL) {
    if (iso_intrp == 0) return delta_fn_linear_interpolation(freq1, points, amp, spec);
    if (iso_intrp == 1) return delta_fn_gauss_interpolation(freq1, points, amp, spec);
  }
  __triangle_interpolation(freq1, freq2, freq3, amp, spec, points);
}

void triangle_interpolation1D_linear(double *freq1, double *freq2, double *freq3,
                                     double *amp, double *spec, int *points) {
  if (absd(*freq1 - *freq2) < TOL && absd(*freq1 - *freq3) < TOL)
    return delta_fn_linear_interpolation(freq1, points, amp, spec);

  __triangle_interpolation(freq1, freq2, freq3, amp, spec, points);
}

void triangle_interpolation1D_gaussian(double *freq1, double *freq2, double *freq3,
                                       double *amp, double *spec, int *points) {
  if (absd(*freq1 - *freq2) < TOL && absd(*freq1 - *freq3) < TOL)
    return delta_fn_gauss_interpolation(freq1, points, amp, spec);

  __triangle_interpolation(freq1, freq2, freq3, amp, spec, points);
}

/**
 * ================================================================================== *
 *                         Two-dimensional interpolation                              *
 * ================================================================================== *
 */

// static inline void rectangle_bin(double x0, double x1, double amp, double *spec,
//                                  int m1) {
//   double da, diff;
//   bool test = (x0 < x1);

//   diff = x0;
//   x0 = test ? x0 : x1;
//   x1 = test ? x1 : diff;

//   int p = (int)x0, pmax = (int)x1;

//   da = amp / (x1 - x0);

//   if (p < 0) {
//     p = 0;
//     x0 = 0.0;
//   }
//   if (pmax >= m1) {
//     pmax = m1;
//     x1 = (double)m1;
//   }
//   diff = (double)p + 1. - x0;
//   spec[2 * p++] += diff * da;

//   while (p < pmax) spec[2 * p++] += da;
//   spec[2 * p] += (x1 - (double)p) * da;
// }

// //    f10 ----------- f11  top-line
// //    /|             | \.
// //   / |             |    \.
// // f00 ------------------ f01  bottom line
// static inline void quadrilateral_bin_2(double *f00, double *f11, double *f10,
//                                        double *f01, double top, double bottom,
//                                        double total, double amp, double *spec, int
//                                        m1) {
//   double amp_left, amp_mid, amp_right, a, a1, a2, a3;

//   // if (*f00 == *f01) {
//   //   triangle_interpolation1D(f00, f11, f10, &amp, spec, &m1);
//   // } else if (*f10 == *f11) {
//   //   triangle_interpolation1D(f00, f11, f01, &amp, spec, &m1);
//   // } else {

//   a1 = 0.5 * absd(*f10 - *f00);
//   a2 = absd(*f11 - *f10);
//   a3 = 0.5 * absd(*f01 - *f11);
//   a = a1 + a2 + a3;
//   amp_left = amp * a1 / a;
//   amp_mid = amp * a2 / a;
//   amp_right = amp * a3 / a;
//   triangle_interpolation1D(f00, f10, f10, &amp_left, spec, &m1);
//   rectangle_bin(*f10, *f11, amp_mid, spec, m1);
//   triangle_interpolation1D(f11, f11, f01, &amp_right, spec, &m1);
//   // }
// }

//    f00 ----------- f01  top-line
//    /   \ _____        \.
//   /            \_____  \.
// f10 ------------------ f11  bottom line
static inline void quadrilateral_bin(double *f00, double *f11, double *f10, double *f01,
                                     double top, double bottom, double total,
                                     double amp, double *spec, int m1,
                                     unsigned int iso_intrp) {
  double amp_bottom, amp_top, norm;

  if (total != 0) {
    norm = amp / total;
    amp_bottom = bottom * norm;
    amp_top = top * norm;
    triangle_interpolation1D(f00, f11, f10, &amp_bottom, spec, &m1, iso_intrp);
    triangle_interpolation1D(f00, f11, f01, &amp_top, spec, &m1, iso_intrp);
  } else {
    triangle_interpolation1D(f00, f11, f10, &amp, spec, &m1, iso_intrp);
  }
}

static inline void lower_triangle_interpolation_2d(int p, int pmid, bool l_clip,
                                                   bool r_clip, bool r_clip_r,
                                                   double top, double *f1, double *f2,
                                                   double *x10, double *x11, int m1,
                                                   double *spec,
                                                   unsigned int iso_intrp) {
  double f10, slope_diff, abs_sdiff, abs_sdiff_2, temp, df1, diff;
  double amp, denom, line_up, line_down;
  double x00, x01, f01_slope, f02_slope;

  f10 = f1[1] - f1[0];
  df1 = top / f10;

  f01_slope = (f2[1] - f2[0]) / f10;
  f02_slope = (f2[2] - f2[0]) / (f1[2] - f1[0]);
  slope_diff = f02_slope - f01_slope;
  abs_sdiff = absd(slope_diff);
  abs_sdiff_2 = 2.0 * abs_sdiff;

  spec += 2 * p * m1;
  if (p == pmid) {
    *x10 = f2[1];
    *x11 = f02_slope * f10 + f2[0];
    if (!r_clip && !l_clip) {
      amp = f10 * top * 0.5;
      triangle_interpolation1D(&f2[0], x11, x10, &amp, spec, &m1, iso_intrp);
    }
    if (!r_clip_r && l_clip) {
      amp = f1[1] * (f10 - f1[0]) * 0.5 * df1;
      triangle_interpolation1D(&f2[0], x11, x10, &amp, spec, &m1, iso_intrp);
    }
    return;
  }

  // interpolate f2 dimensions here.  // clip_left1 = l_clip
  // Calculate the f2 coordinates from lines connecting f0-f1 and f0-f2

  // Part 1: At start bin.
  diff = (double)p + 1. - f1[0];

  *x10 = f01_slope * diff + f2[0];
  *x11 = f02_slope * diff + f2[0];
  if (!l_clip) {
    amp = 0.5 * diff * diff * df1;
    triangle_interpolation1D(&f2[0], x11, x10, &amp, spec, &m1, iso_intrp);
  } else {
    amp = (diff - 0.5) * df1;
    x00 = *x10 - f01_slope;
    x01 = *x11 - f02_slope;
    line_down = absd(*x11 - *x10);
    line_up = absd(x01 - x00);
    denom = line_down + line_up;
    quadrilateral_bin(&x00, x11, x10, &x01, line_up, line_down, denom, amp, spec, m1,
                      iso_intrp);
  }
  p++;
  spec += 2 * m1;

  // Part 2: After start to before mid bin.
  temp = diff * slope_diff;
  line_up = absd(temp);
  line_down = absd(temp + slope_diff);
  denom = line_up + line_down;

  diff -= 0.5;
  diff *= df1;
  while (p < pmid) {
    diff += df1;
    x00 = *x10;
    x01 = *x11;
    *x10 += f01_slope;
    *x11 += f02_slope;
    quadrilateral_bin(&x00, x11, x10, &x01, line_up, line_down, denom, diff, spec, m1,
                      iso_intrp);
    line_up += abs_sdiff;
    line_down += abs_sdiff;
    denom += abs_sdiff_2;
    p++;
    spec += 2 * m1;
  }

  // Part 3: population to the left of the mid bin.
  if (!r_clip) {
    amp = (f1[1] - (double)p) * (f10 + ((double)p - f1[0])) * 0.5 * df1;
    x00 = *x10;
    x01 = *x11;
    *x10 = f2[1];
    *x11 = f02_slope * f10 + f2[0];
    line_down = absd(*x11 - *x10);
    denom = line_up + line_down;
    quadrilateral_bin(&x00, x11, x10, &x01, line_up, line_down, denom, amp, spec, m1,
                      iso_intrp);
  } else {
    diff += df1;
    x00 = *x10;
    x01 = *x11;
    *x10 += f01_slope;
    *x11 += f02_slope;
    quadrilateral_bin(&x00, x11, x10, &x01, line_up, line_down, denom, diff, spec, m1,
                      iso_intrp);
  }
}

static inline void upper_triangle_interpolation_2d(int p, int pmax, bool l_clip,
                                                   bool r_clip, bool r_clip_l,
                                                   double top, double *f1, double *f2,
                                                   double *x10, double *x11, int m1,
                                                   double *spec,
                                                   unsigned int iso_intrp) {
  double f21, slope_diff, abs_sdiff, abs_sdiff_2, temp, df2, diff;
  double amp, denom, line_up, line_down;
  double x00, x01, f12_slope, f02_slope;

  diff = f1[2] - (double)p - 1.;
  f21 = f1[2] - f1[1];
  df2 = top / f21;

  f02_slope = (f2[2] - f2[0]) / (f1[2] - f1[0]);
  f12_slope = (f2[2] - f2[1]) / f21;
  slope_diff = f02_slope - f12_slope;
  abs_sdiff = absd(slope_diff);
  abs_sdiff_2 = 2 * abs_sdiff;

  spec += 2 * p * m1;
  if (p == pmax) {
    if (!r_clip) {
      amp = f21 * top * 0.5;
      triangle_interpolation1D(x11, &f2[1], &f2[2], &amp, spec, &m1, iso_intrp);
    }
    if (r_clip && !r_clip_l) {
      amp = (f21 - diff) * (diff + f21) * 0.5 * df2;
      triangle_interpolation1D(x11, &f2[1], &f2[2], &amp, spec, &m1, iso_intrp);
    }
    return;
  }

  // interpolate f2 dimensions here.
  // Calculate the f2 coordinates from lines connecting f1-f2 and f0-f2

  // Part 4: population to the right of the mid bin.
  if (!l_clip) {
    amp = (f21 - diff) * (diff + f21) * 0.5 * df2;
    x00 = *x10;
    x01 = *x11;
    *x10 = f12_slope * ((double)(p + 1) - f1[1]) + f2[1];
    *x11 = f02_slope * ((double)(p + 1) - f1[0]) + f2[0];
    line_up = absd(x01 - x00);
    line_down = absd(*x11 - *x10);
    denom = line_down + line_up;
    quadrilateral_bin(&x00, x11, x10, &x01, line_up, line_down, denom, amp, spec, m1,
                      iso_intrp);
  } else {
    amp = (diff + 0.5) * df2;
    x00 = f12_slope * ((double)p - f1[1]) + f2[1];
    x01 = f02_slope * ((double)p - f1[0]) + f2[0];
    *x10 = x00 + f12_slope;
    *x11 = x01 + f02_slope;
    line_up = absd(x01 - x00);
    line_down = absd(*x11 - *x10);
    denom = line_down + line_up;
    quadrilateral_bin(&x00, x11, x10, &x01, line_up, line_down, denom, amp, spec, m1,
                      iso_intrp);
  }
  p++;
  spec += 2 * m1;

  // Part 6: After mid to before end bin.
  diff += 0.5;
  diff *= df2;
  line_up = line_down;
  line_down -= abs_sdiff;
  denom = line_up + line_down;
  while (p < pmax) {
    diff -= df2;
    x00 = *x10;
    x01 = *x11;
    *x10 += f12_slope;
    *x11 += f02_slope;
    quadrilateral_bin(&x00, x11, x10, &x01, line_up, line_down, denom, diff, spec, m1,
                      iso_intrp);
    line_up -= abs_sdiff;
    line_down -= abs_sdiff;
    denom -= abs_sdiff_2;
    p++;
    spec += 2 * m1;
  }

  // Part 7: The end bin.
  if (!r_clip) {
    temp = (f1[2] - (double)p);
    amp = temp * temp * 0.5 * df2;
    x01 = *x11;
    *x11 = f2[2];
    triangle_interpolation1D(x10, x11, &x01, &amp, spec, &m1, iso_intrp);
  } else {
    diff -= df2;
    x00 = *x10;
    x01 = *x11;
    *x10 += f12_slope;
    *x11 += f02_slope;
    quadrilateral_bin(&x00, x11, x10, &x01, line_up, line_down, denom, diff, spec, m1,
                      iso_intrp);
  }
}

void triangle_interpolation2D(double *freq11, double *freq12, double *freq13,
                              double *freq21, double *freq22, double *freq23,
                              double *amp, double *spec, int m0, int m1,
                              unsigned int iso_intrp) {
  double top, t1, t2, diff, temp, n_i;
  int p, pmid, pmax, i, j;
  double freq10_01, freq11_02;

  p = (int)(*freq11);

  if (absd(freq11[0] - freq12[0]) < TOL && absd(freq11[0] - freq13[0]) < TOL) {
    if (p >= m0 || p < 0) return;

    diff = freq11[0] - (double)p;
    n_i = 0.5;
    if (absd(diff - n_i) < TOL) {
      triangle_interpolation1D(freq21, freq22, freq23, amp, &spec[2 * p * m1], &m1,
                               iso_intrp);
      return;
    }
    if (diff < n_i) {
      if (p != 0) {
        temp = *amp * (n_i - diff);
        triangle_interpolation1D(freq21, freq22, freq23, &temp, &spec[2 * (p - 1) * m1],
                                 &m1, iso_intrp);
      }
      temp = *amp * (n_i + diff);
      triangle_interpolation1D(freq21, freq22, freq23, &temp, &spec[2 * p * m1], &m1,
                               iso_intrp);
      return;
    }
    if (diff > n_i) {
      if (p + 1 != m0) {
        temp = *amp * (diff - n_i);
        triangle_interpolation1D(freq21, freq22, freq23, &temp, &spec[2 * (p + 1) * m1],
                                 &m1, iso_intrp);
      }
      temp = *amp * (1 + n_i - diff);
      triangle_interpolation1D(freq21, freq22, freq23, &temp, &spec[2 * p * m1], &m1,
                               iso_intrp);
      return;
    }
    return;
  }

  double f1[3] = {freq11[0], freq12[0], freq13[0]};
  double f2[3] = {freq21[0], freq22[0], freq23[0]};

  // arrange the numbers in ascending order
  for (j = 1; j <= 2; j++) {
    t1 = f1[j];
    t2 = f2[j];
    i = j - 1;
    while (i >= 0 && f1[i] > t1) {
      f1[i + 1] = f1[i];
      f2[i + 1] = f2[i];
      i--;
    }
    f1[i + 1] = t1;
    f2[i + 1] = t2;
  }

  // if min frequency is higher than the last bin, return
  p = (int)f1[0];
  if (p >= m0) return;

  // if max frequency is lower than the first bin, return
  pmax = (int)f1[2];
  if (pmax < 0) return;

  pmid = (int)f1[1];
  top = *amp * 2.0 / (f1[2] - f1[0]);

  bool clips[4];
  get_clipping_conditions(&p, &pmid, &pmax, &m0, clips);
  if (f1[1] >= 0.0) {
    lower_triangle_interpolation_2d(p, pmid, clips[0], clips[1], clips[2], top, f1, f2,
                                    &freq10_01, &freq11_02, m1, spec, iso_intrp);
  }
  if (f1[2] >= 0.0) {
    upper_triangle_interpolation_2d(pmid, pmax, clips[2], clips[3], clips[1], top, f1,
                                    f2, &freq10_01, &freq11_02, m1, spec, iso_intrp);
  }
}

void rasterization(double *grid, double *v0, double *v1, double *v2, int rows,
                   int columns) {
  double A12, B12, C12, A20, B20, C20, A01, B01, C01;
  double minX, minY, maxX, maxY, w0, w1, w2;
  double w0_row, w1_row, w2_row;
  int i, j, i_, minX_, minY_, maxX_, maxY_;

  minX = fmin(fmin(v0[0], v1[0]), v2[0]);
  minY = fmin(fmin(v0[1], v1[1]), v2[1]);
  maxX = fmax(fmax(v0[0], v1[0]), v2[0]);
  maxY = fmax(fmax(v0[1], v1[1]), v2[1]);

  // clip against screen bounds
  minX_ = (int)fmax(minX, 0.);
  minY_ = (int)fmax(minY, 0.);
  maxX_ = (int)fmin(maxX, (double)rows - 1.);
  maxY_ = (int)fmin(maxY, (double)columns - 1.);

  A12 = (v2[0] - v1[0]);
  B12 = (v2[1] - v1[1]);
  C12 = -A12 * v1[1] + B12 * v1[0];
  A20 = (v0[0] - v2[0]);
  B20 = (v0[1] - v2[1]);
  C20 = -A20 * v2[1] + B20 * v2[0];
  A01 = (v1[0] - v0[0]);
  B01 = (v1[1] - v0[1]);
  C01 = -A01 * v0[1] + B01 * v0[0];

  w0_row = A12 * minY - B12 * minX + C12;
  w1_row = A20 * minY - B20 * minX + C20;
  w2_row = A01 * minY - B01 * minX + C01;

  // Rasterizer
  for (i = minY_; i <= maxY_; i++) {
    // Determine barycentric coordinates
    w0 = w0_row;
    w1 = w1_row;
    w2 = w2_row;

    i_ = rows * i;
    for (j = minX_; j <= maxX_; j++) {
      // If p is on or inside all edges, render pixel.
      if ((int)w0 >= 0 && (int)w1 >= 0 && (int)w2 >= 0) {
        grid[i_ + j] += 1.;  //(w0+w1+w2);
      }
      if ((int)w0 <= 0 && (int)w1 <= 0 && (int)w2 <= 0) {
        grid[i_ + j] += -1.;  //(w0+w1+w2);
      }
      // i_++;

      w0 -= B12;
      w1 -= B20;
      w2 -= B01;
    }
    w0_row += A12;
    w1_row += A20;
    w2_row += A01;
  }
}

// static inline void get_KL(
//             double *v0,
//             double *v1,
//             double *x,
//             double *y)
// {
//   double K[2] =
// }
// def get_KL(self, x0, y0, x1, y1, eps = 1e-5):
//     self.x0, self.y0, self.x1, self.y1 = x0, y0, x1, y1
//     Kx, Ky, Lx, Ly = 0, 0, 0, 0
//     for sec in self.clip(0, 1, 1, 0):
//         v1 = Point(sec[0], sec[1])
//         v0 = Point(sec[2], sec[3])
//         if abs(v0.x - 1) < eps and abs(v1.x - 1) < eps \
//         or abs(v0.y - 1) < eps and abs(v1.y - 1) < eps:
//             continue

//         Kx += 1./4 * (v0.y - v1.y)
//         Ky += 1./4 * (v1.x - v0.x)
//         Lx += 1./8 * (v0.y-v1.y) * (v0.x+v1.x)
//         Ly += 1./8 * (v1.x-v0.x) * (v0.y+v1.y)
//     return Point(Kx, Ky), Point(Lx, Ly)

// static inline void clip(
//             double left,
//             double right,
//             double bottom,
//             double top,
//             double *p0,
//             double *p1,
//             double *clipped_lines
//             )
// {
//   int edge, i;
//   double t0=0, t1=1, r;
//   double delta_x = p1[0] - p0[0], delta_y=p1[1]-p0[1], p, q;

//   for(edge=0; edge<4; edge++){
//     if(edge ==0){
//       p = -delta_x;
//       q = -(left - p0[0]);
//     }
//     else if(edge ==1){
//       p = delta_x;
//       q = (right - p0[0]);
//     }
//     else if(edge ==2){
//       p = delta_y;
//       q = (bottom - p0[1]);
//     }
//     else if(edge ==3){
//       p = -delta_y;
//       q = -(top - p0[1]);
//     }
//     if(p == 0 && q < 0) return 0;
//     if (p < 0){
//       r = q / (double)p;
//       if(r > t1)  return 0;
//       if(r > t0) t0 = r;  // line is clipped!
//     }
//     else if(p > 0){
//       r = q / (double)p;
//       if(r < t0) return 0;
//       if(r < t1)  t1 = r;   // line is clipped!
//     }
//   }
// }

//     def clip(self, left, right, bottom, top):
//         t0, t1 = 0, 1
//         xdelta = self.x1 - self.x0
//         ydelta = self.y1 - self.y0
//         for edge in range(4): #traverse through left, right, bottom, top
//         edges.
//             if   edge == 0:   p, q = -xdelta, -(left-self.x0)
//             elif edge == 1:   p, q =  xdelta,  (right-self.x0)
//             elif edge == 2:   p, q =  ydelta,  (bottom-self.y0)
//             elif edge == 3:   p, q = -ydelta, -(top-self.y0)
//             if p == 0 and q < 0:    return []
//             if p < 0:
//                 r = q / float(p)
//                 if r > t1:          return []
//                 elif r > t0:        t0 = r   # line is clipped!
//             elif p > 0:
//                 r = q / float(p)
//                 if r < t0:          return []
//                 elif r < t1:        t1 = r   # line is clipped!
//         clipped_line = (self.x0 + t0*xdelta, self.y0 + t0*ydelta,
//                         self.x0 + t1*xdelta, self.y0 + t1*ydelta)
//         return [clipped_line]

void octahedronDeltaInterpolation(const unsigned int nt, double *freq, double *amp,
                                  int stride, int n_spec, double *spec,
                                  unsigned int iso_intrp) {
  int i = 0, j = 0, local_index, n_pts = (nt + 1) * (nt + 2) / 2;
  unsigned int int_i_stride = 0, int_j_stride = 0;
  double amp1, temp, *amp_address;

  local_index = nt - 1;
  amp_address = &amp[(nt + 1) * stride];
  amp1 = 0.0;
  while (i < n_pts - 1) {
    temp = amp[int_i_stride + stride] + amp_address[int_j_stride];
    amp1 += temp + amp[int_i_stride];

    if (i < local_index) {
      temp += amp_address[int_j_stride + stride];
      amp1 += temp;
    } else {
      local_index = j + nt;
      i++;
      int_i_stride += stride;
    }
    i++;
    j++;
    int_i_stride += stride;
    int_j_stride += stride;
  }
  if (iso_intrp == 0) return delta_fn_linear_interpolation(freq, &n_spec, &amp1, spec);
  if (iso_intrp == 1) return delta_fn_gauss_interpolation(freq, &n_spec, &amp1, spec);
}

void octahedronInterpolation(double *spec, double *freq, const unsigned int nt,
                             double *amp, int stride, int m) {
  int i = 0, j = 0, local_index, n_pts = (nt + 1) * (nt + 2) / 2;
  unsigned int int_i_stride = 0, int_j_stride = 0;
  double amp1, temp, *amp_address, *freq_address;

  /* Interpolate between 1d points by setting up triangles of unit area */
  local_index = nt - 1;
  amp_address = &amp[(nt + 1) * stride];
  freq_address = &freq[nt + 1];

  // Note, amp is the sum of amplitude at three vertexes because the factor 3 is already
  // applied to the amplitude vector.
  while (i < n_pts - 1) {
    temp = amp[int_i_stride + stride] + amp_address[int_j_stride];
    amp1 = temp + amp[int_i_stride];

    __triangle_interpolation(&freq[i], &freq[i + 1], &freq_address[j], &amp1, spec, &m);

    if (i < local_index) {
      temp += amp_address[int_j_stride + stride];
      __triangle_interpolation(&freq[i + 1], &freq_address[j], &freq_address[j + 1],
                               &temp, spec, &m);
    } else {
      local_index = j + nt;
      i++;
      int_i_stride += stride;
    }
    i++;
    j++;
    int_i_stride += stride;
    int_j_stride += stride;
  }
}

void generic_1d_triangle_interpolation(double *spec, const unsigned int freq_size,
                                       double *freq, double *amp, int m,
                                       const unsigned int position_size,
                                       int32_t *positions) {
  unsigned int pos_size = position_size;
  int32_t p1, p2, p3;
  double amp_sum;

  while (pos_size-- > 0) {
    p1 = *positions++;
    p2 = *positions++;
    p3 = *positions++;
    // we do amp_sum because amps are already scaled to account for the factor 3
    amp_sum = amp[p1] + amp[p2] + amp[p3];
    __triangle_interpolation(&freq[p1], &freq[p2], &freq[p3], &amp_sum, spec, &m);
  }
}

void hist1d(double *spec, const unsigned int freq_size, double *freq, double *amp,
            int m, const unsigned int nt) {
  unsigned int i = 0, ix;
  double temp_freq;

  for (i = 0; i < freq_size; i++) {
    temp_freq = freq[i];
    if (temp_freq >= 0 && temp_freq < m) {
      ix = (unsigned int)temp_freq;
      spec[2 * ix] += amp[i];
    }
  }
}

void one_d_averaging(double *spec, const unsigned int freq_size, double *freq,
                     double *amp_real, double *amp_imag, int dimension_count,
                     const unsigned int position_size, int32_t *positions,
                     const unsigned int nt, bool user_defined, bool interpolation,
                     bool is_complex) {
  if (!user_defined) {
    if (interpolation) {
      octahedronInterpolation(spec, freq, nt, amp_real, 1, dimension_count);
      if (is_complex) {
        octahedronInterpolation(spec + 1, freq, nt, amp_imag, 1, dimension_count);
      }
    } else {
      hist1d(spec, freq_size, freq, amp_real, dimension_count, nt);
      if (is_complex) {
        hist1d(spec + 1, freq_size, freq, amp_imag, dimension_count, nt);
      }
    }
  } else {
    generic_1d_triangle_average(spec, freq_size, freq, amp_real, dimension_count,
                                position_size, positions, nt);
    if (is_complex) {
      generic_1d_triangle_average(spec + 1, freq_size, freq, amp_imag, dimension_count,
                                  position_size, positions, nt);
    }
  }
}

void generic_1d_triangle_average(double *spec, const unsigned int freq_size,
                                 double *freq, double *amp, int m,
                                 const unsigned int position_size, int32_t *positions,
                                 const unsigned int nt) {
  if (positions == NULL) {
    hist1d(spec, freq_size, freq, amp, m, nt);
  } else {
    generic_1d_triangle_interpolation(spec, freq_size, freq, amp, m, position_size,
                                      positions);
  }
}

void two_d_averaging(double *spec, const unsigned int freq_size, double *freq1,
                     double *freq2, double *amp, int amp_stride,
                     const unsigned int position_size, int32_t *positions,
                     int dimension0_count, int dimension1_count, unsigned int iso_intrp,
                     const unsigned int nt, bool user_defined, bool interpolation) {
  if (!user_defined) {
    if (interpolation) {
      // Perform tenting on every sideband order over all orientations
      octahedronInterpolation2D(spec, freq1, freq2, nt, amp, amp_stride,
                                dimension0_count, dimension1_count, iso_intrp);
    } else {
      hist2d(spec, freq_size, freq1, freq2, amp, amp_stride, dimension0_count,
             dimension1_count, nt);
    }
  } else {
    generic_2d_triangle_average(spec, freq_size, freq1, freq2, amp, amp_stride,
                                dimension0_count, dimension1_count, position_size,
                                positions, nt, iso_intrp);
  }
}

void generic_2d_triangle_interpolation(double *spec, const unsigned int freq_size,
                                       double *freq1, double *freq2, double *amp,
                                       int amp_stride, const unsigned int position_size,
                                       int32_t *positions, int m0, int m1,
                                       unsigned int iso_intrp) {
  unsigned int pos_size = position_size;
  int32_t p1, p2, p3;
  double amp_sum;

  while (pos_size-- > 0) {
    p1 = *positions++;
    p2 = *positions++;
    p3 = *positions++;
    // we do amp_sum because amps are already scaled to account for the factor 3
    amp_sum = amp[amp_stride * p1] + amp[amp_stride * p2] + amp[amp_stride * p3];
    triangle_interpolation2D(&freq1[p1], &freq1[p2], &freq1[p3], &freq2[p1], &freq2[p2],
                             &freq2[p3], &amp_sum, spec, m0, m1, iso_intrp);
  }
}

void hist2d(double *spec, const unsigned int freq_size, double *freq1, double *freq2,
            double *amp, int amp_stride, int m0, int m1, const unsigned int nt) {
  unsigned int i = 0, ix, iy, hist_index;
  double temp_freq_1, temp_freq_2;

  for (i = 0; i < freq_size; i++) {
    temp_freq_1 = freq1[i];
    temp_freq_2 = freq2[i];
    if (temp_freq_1 >= 0 && temp_freq_1 < m0 && temp_freq_2 >= 0 && temp_freq_2 < m1) {
      ix = (unsigned int)temp_freq_1;
      iy = (unsigned int)temp_freq_2;
      hist_index = iy + m1 * ix;
      spec[2 * hist_index] += amp[i * amp_stride];
    }
  }
}

void generic_2d_triangle_average(double *spec, const unsigned int freq_size,
                                 double *freq1, double *freq2, double *amp,
                                 int amp_stride, int m0, int m1,
                                 const unsigned int position_size, int32_t *positions,
                                 const unsigned int nt, unsigned int iso_intrp) {
  if (positions == NULL) {
    hist2d(spec, freq_size, freq1, freq2, amp, amp_stride, m0, m1, nt);
  } else {
    generic_2d_triangle_interpolation(spec, freq_size, freq1, freq2, amp, amp_stride,
                                      position_size, positions, m0, m1, iso_intrp);
  }
}

void octahedronInterpolation2D(double *spec, double *freq1, double *freq2, int nt,
                               double *amp, int stride, int m0, int m1,
                               unsigned int iso_intrp) {
  int i = 0, j = 0, local_index, n_pts = (nt + 1) * (nt + 2) / 2;
  unsigned int int_i_stride = 0, int_j_stride = 0;
  double amp1, temp, *amp_address, *freq1_address, *freq2_address;

  /* Interpolate between 1d points by setting up triangles of unit area */

  local_index = nt - 1;
  amp_address = &amp[(nt + 1) * stride];
  freq1_address = &freq1[nt + 1];
  freq2_address = &freq2[nt + 1];

  // Note, amp is the sum of amplitude at three vertexes because the factor 3 is already
  // applied to the amplitude vector.
  while (i < n_pts - 1) {
    temp = amp[int_i_stride + stride] + amp_address[int_j_stride];
    amp1 = temp + amp[int_i_stride];

    triangle_interpolation2D(&freq1[i], &freq1[i + 1], &freq1_address[j], &freq2[i],
                             &freq2[i + 1], &freq2_address[j], &amp1, spec, m0, m1,
                             iso_intrp);

    if (i < local_index) {
      temp += amp_address[int_j_stride + stride];
      triangle_interpolation2D(&freq1[i + 1], &freq1_address[j], &freq1_address[j + 1],
                               &freq2[i + 1], &freq2_address[j], &freq2_address[j + 1],
                               &temp, spec, m0, m1, iso_intrp);
    } else {
      local_index = j + nt;
      i++;
      int_i_stride += stride;
    }
    i++;
    j++;
    int_i_stride += stride;
    int_j_stride += stride;
  }
}
