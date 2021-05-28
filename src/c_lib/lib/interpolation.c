// -*- coding: utf-8 -*-
//
//  interpolation.c
//
//  @copyright Deepansh J. Srivastava, 2019-2021.
//  Created by Deepansh J. Srivastava, Apr 11, 2019.
//  Contact email = srivastava.89@osu.edu
//

#include "interpolation.h"

// triangle_interpolation is an optimized version of tent. Still plenty of
// room for optimization.

static void inline delta_fn_interpolation(const double *freq1, const int *points,
                                          double *amp, double *spec) {
  double diff, n_i;
  int p = (int)*freq1;
  if (p >= *points || p < 0) return;

  diff = *freq1 - (double)p;
  n_i = 0.5;
  if (fabs(diff - n_i) < TOL) {
    spec[p] += *amp;
    return;
  }
  if (diff < n_i) {
    if (p != 0) spec[p - 1] += *amp * (n_i - diff);
    spec[p] += *amp * (n_i + diff);
    return;
  }
  if (diff > n_i) {
    if (p + 1 != *points) spec[p + 1] += *amp * (diff - n_i);
    spec[p] += *amp * (1 + n_i - diff);
    return;
  }
}

// interpolate first half of the triangle.
//        /|                /|                    / |
//       /*|              /**|                 /    |
//      /**|            /****|              /*|     |
//     /***|          /|*****|           /|***|     |
//    /____|        /  |_____|        /   |___|_____|
//   p    pmid         p    pmid          p        pmid
// clips are off     left clip on   both left and right clips on
static inline void left_triangle_interpolate(int p, int pmid, bool l_clip, bool r_clip,
                                             double top, double *f, double *spec) {
  double f10 = f[1] - f[0], dp;

  if (p == pmid) {
    spec[p] += (!r_clip && !l_clip) ? f10 * top * 0.5 : 0.0;
    return;
  }

  double df1 = top / f10, diff = (double)p + 1. - f[0];

  spec[p++] += (!l_clip) ? 0.5 * diff * diff * df1 : (diff - 0.5) * df1;

  diff -= 0.5;
  diff *= df1;
  while (p != pmid) spec[p++] += diff += df1;

  dp = (double)p;
  spec[p] += (!r_clip) ? (f[1] - dp) * (f10 + (dp - f[0])) * 0.5 * df1 : 0.0;
}

// interpolate second half of the triangle.
static inline void right_triangle_interpolate(int p, int pmax, bool l_clip, bool r_clip,
                                              double top, double *f, double *spec) {
  double f21 = f[2] - f[1];
  if (p == pmax) {
    spec[p] += (!r_clip) ? f21 * top * 0.5 : 0.0;
    return;
  }

  double df1 = top / f21, diff = f[2] - (double)p - 1.0;

  spec[p++] += (!l_clip) ? (f21 - diff) * (diff + f21) * 0.5 * df1 : (diff + 0.5) * df1;

  diff += 0.5;
  diff *= df1;
  while (p != pmax) spec[p++] += diff -= df1;

  if (!r_clip) {
    diff = (f[2] - (double)p);
    spec[p] += diff * diff * 0.5 * df1;
  }
}

static inline void __triangle_interpolation(double *freq1, double *freq2, double *freq3,
                                            double *amp, double *spec, int *points) {
  double top, t;

  int p, pmid, pmax, i, j;
  bool clip_right1 = false, clip_left1 = false, clip_right2 = false, clip_left2 = false;

  p = (int)(*freq1);

  // check if the three points lie within a bin interval.
  if ((int)freq1[0] == (int)freq2[0] && (int)freq1[0] == (int)freq3[0]) {
    if (p >= *points || p < 0) return;
    spec[p] += *amp;
    return;
  }

  // arrange the numbers in ascending order (sort)
  double f[3] = {freq1[0], freq2[0], freq3[0]};
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
  if (p > *points) return;

  // if max frequency is lower than the first bin, return
  pmax = (int)f[2];
  if (pmax < 0) return;

  pmid = (int)f[1];
  if (pmid >= *points) {
    pmid = *points;
    clip_right1 = true;
  }

  if (pmax >= *points) {
    pmax = *points;
    clip_right2 = true;
  }

  if (p < 0) {
    p = 0;
    clip_left1 = true;
  }

  if (pmid < 0) {
    pmid = 0;
    clip_left2 = true;
  }

  top = (*amp * 2.0 / (f[2] - f[0]));

  left_triangle_interpolate(p, pmid, clip_left1, clip_right1, top, f, spec);
  right_triangle_interpolate(pmid, pmax, clip_left2, clip_right2, top, f, spec);

  return;
}

int triangle_interpolation2D(double *freq11, double *freq12, double *freq13,
                             double *freq21, double *freq22, double *freq23,
                             double *amp, double *spec, int m0, int m1) {
  double df1, df2, top = 0.0, t1, t2, diff, f10 = 0.0, f21 = 0.0, temp, n_i;
  double slope_diff, abs_slope_diff, line_up, line_down;
  int p, pmid, pmax, i, j;
  int clip_right1 = 0, clip_left1 = 0, clip_right2 = 0, clip_left2 = 0;
  double f01_slope, f02_slope, f12_slope, amp_section, denom;
  double area_up_triangle, area_down_triangle;
  double freq00_01, freq01_02, freq10_01, freq11_02, freq00_12, freq10_12;

  p = (int)(freq11[0]);

  if (fabs(freq11[0] - freq12[0]) < TOL && fabs(freq11[0] - freq13[0]) < TOL) {
    if (p >= m0 || p < 0) {
      return 0;
    }
    diff = freq11[0] - (double)p;
    n_i = 0.5;
    if (fabs(diff - n_i) < TOL) {
      triangle_interpolation1D(freq21, freq22, freq23, amp, &spec[p * m1], &m1);
      return 0;
    }
    if (diff < n_i) {
      if (p != 0) {
        temp = amp[0] * (n_i - diff);
        triangle_interpolation1D(freq21, freq22, freq23, &temp, &spec[(p - 1) * m1],
                                 &m1);
      }
      temp = amp[0] * (n_i + diff);
      triangle_interpolation1D(freq21, freq22, freq23, &temp, &spec[p * m1], &m1);
      return 0;
    }
    if (diff > n_i) {
      if (p + 1 != m0) {
        temp = amp[0] * (diff - n_i);
        triangle_interpolation1D(freq21, freq22, freq23, &temp, &spec[(p + 1) * m1],
                                 &m1);
      }
      temp = amp[0] * (1 + n_i - diff);
      triangle_interpolation1D(freq21, freq22, freq23, &temp, &spec[p * m1], &m1);
      return 0;
    }
    return 0;
  }

  if ((int)freq11[0] == (int)freq12[0] && (int)freq11[0] == (int)freq13[0]) {
    if (p >= m0 || p < 0) {
      return 0;
    }
    triangle_interpolation1D(freq21, freq22, freq23, amp, &spec[p * m1], &m1);
    return 0;
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

  p = (int)f1[0];
  if (p > m0) {
    return 0;
  }

  pmax = (int)f1[2];
  if (pmax < 0) {
    return 0;
  }

  pmid = (int)f1[1];
  if (pmid >= m0) {
    pmid = m0;
    clip_right1 = 1;
  }

  if (pmax >= m0) {
    pmax = m0;
    clip_right2 = 1;
  }

  if (p < 0) {
    p = 0;
    clip_left1 = 1;
  }

  if (pmid < 0) {
    pmid = 0;
    clip_left2 = 1;
  }

  top += (amp[0] * 2.0 / (f1[2] - f1[0]));
  f10 += f1[1] - f1[0];
  f21 += f1[2] - f1[1];

  f01_slope = (f2[1] - f2[0]) / (f1[1] - f1[0]);
  f02_slope = (f2[2] - f2[0]) / (f1[2] - f1[0]);
  slope_diff = f02_slope - f01_slope;
  abs_slope_diff = fabs(slope_diff);

  if (clip_left2 == 0) {
    if (p != pmid) {
      // interpolate f2 dimensions here.
      // Calculate the f2 coordinates from lines connecting f0-f1 and f0-f2

      // Part 1: At start bin.
      df1 = top / f10;
      diff = (double)p + 1. - f1[0];
      if (clip_left1 == 0) {
        amp_section = 0.5 * diff * diff * df1;
        freq00_01 = f2[0];
        freq10_01 = f01_slope * diff + f2[0];
        freq11_02 = f02_slope * diff + f2[0];
        triangle_interpolation1D(&freq00_01, &freq11_02, &freq10_01, &amp_section,
                                 &spec[p * m1], &m1);
        p++;
      } else {
        amp_section = (diff - 0.5) * df1;
        freq10_01 = f01_slope * diff + f2[0];
        freq11_02 = f02_slope * diff + f2[0];
        freq00_01 = freq10_01 - f01_slope;
        freq01_02 = freq11_02 - f02_slope;

        line_down = fabs(freq11_02 - freq10_01);
        line_up = fabs(freq01_02 - freq00_01);
        denom = line_down + line_up;
        if (denom != 0) {
          area_down_triangle = line_down / denom * amp_section;
          area_up_triangle = line_up / denom * amp_section;
          triangle_interpolation1D(&freq00_01, &freq11_02, &freq10_01,
                                   &area_down_triangle, &spec[p * m1], &m1);
          triangle_interpolation1D(&freq00_01, &freq11_02, &freq01_02,
                                   &area_up_triangle, &spec[p * m1], &m1);
        } else {
          triangle_interpolation1D(&freq00_01, &freq11_02, &freq10_01, &amp_section,
                                   &spec[p * m1], &m1);
        }
        p++;
      }

      // Part 2: After start to before mid bin.
      temp = diff * slope_diff;
      line_up = fabs(temp);
      line_down = fabs(temp + slope_diff);
      denom = line_up + line_down;

      diff -= 0.5;
      diff *= df1;
      while (p != pmid) {
        diff += df1;
        amp_section = diff;
        freq00_01 = freq10_01;
        freq01_02 = freq11_02;
        freq10_01 += f01_slope;
        freq11_02 += f02_slope;

        if (denom != 0) {
          area_down_triangle = line_down / denom * amp_section;
          area_up_triangle = line_up / denom * amp_section;

          triangle_interpolation1D(&freq00_01, &freq11_02, &freq10_01,
                                   &area_down_triangle, &spec[p * m1], &m1);
          triangle_interpolation1D(&freq00_01, &freq11_02, &freq01_02,
                                   &area_up_triangle, &spec[p * m1], &m1);
        } else {
          triangle_interpolation1D(&freq00_01, &freq11_02, &freq10_01, &amp_section,
                                   &spec[p * m1], &m1);
        }
        line_up += abs_slope_diff;
        line_down += abs_slope_diff;
        denom += 2 * abs_slope_diff;
        p++;
      }

      // Part 3: population to the left of the mid bin.
      if (clip_right1 == 0) {
        amp_section = (f1[1] - (double)p) * (f10 + ((double)p - f1[0])) * 0.5 * df1;
        freq00_01 = freq10_01;
        freq01_02 = freq11_02;
        freq10_01 = f2[1];
        freq11_02 = f02_slope * (f1[1] - f1[0]) + f2[0];

        if (denom != 0) {
          denom = fabs(freq11_02 - freq10_01) + line_up;
          area_down_triangle = fabs(freq11_02 - freq10_01) / denom * amp_section;
          area_up_triangle = line_up / denom * amp_section;
          triangle_interpolation1D(&freq00_01, &freq11_02, &freq10_01,
                                   &area_down_triangle, &spec[p * m1], &m1);
          triangle_interpolation1D(&freq00_01, &freq11_02, &freq01_02,
                                   &area_up_triangle, &spec[p * m1], &m1);
        } else {
          triangle_interpolation1D(&freq00_01, &freq11_02, &freq10_01, &amp_section,
                                   &spec[p * m1], &m1);
        }
      }
    } else {
      if (clip_right1 == 0 && clip_left1 == 0) {
        amp_section = f10 * top * 0.5;
        freq10_01 = f2[1];
        freq11_02 = f02_slope * (f1[1] - f1[0]) + f2[0];
        triangle_interpolation1D(&f2[0], &freq10_01, &freq11_02, &amp_section,
                                 &spec[p * m1], &m1);
      }
    }
  }

  f12_slope = (f2[2] - f2[1]) / (f1[2] - f1[1]);
  slope_diff = f02_slope - f12_slope;
  abs_slope_diff = fabs(slope_diff);

  if (p != pmax) {
    // interpolate f2 dimensions here.
    // Calculate the f2 coordinates from lines connecting f1-f2 and f0-f2

    // Part 4: population to the right of the mid bin.
    df2 = top / f21;
    diff = f1[2] - (double)p - 1.;

    if (clip_left2 == 0) {
      amp_section = (f21 - diff) * (diff + f21) * 0.5 * df2;
      freq00_12 = freq10_01;
      freq01_02 = freq11_02;
      freq10_12 = f12_slope * ((double)(p + 1) - f1[1]) + f2[1];
      freq11_02 = f02_slope * ((double)(p + 1) - f1[0]) + f2[0];

      line_up = freq01_02 - freq00_12;
      line_down = freq11_02 - freq10_12;
      denom = fabs(line_down) + fabs(line_up);
      if (denom != 0) {
        area_down_triangle = fabs(line_down) / denom * amp_section;
        area_up_triangle = fabs(line_up) / denom * amp_section;
        triangle_interpolation1D(&freq00_12, &freq11_02, &freq10_12,
                                 &area_down_triangle, &spec[p * m1], &m1);
        triangle_interpolation1D(&freq00_12, &freq11_02, &freq01_02, &area_up_triangle,
                                 &spec[p * m1], &m1);
      } else {
        triangle_interpolation1D(&freq00_12, &freq11_02, &freq10_12, &amp_section,
                                 &spec[p * m1], &m1);
      }
      p++;
    } else {
      amp_section = (diff + 0.5) * df2;
      freq00_12 = f12_slope * ((double)p - f1[1]) + f2[1];
      freq01_02 = f02_slope * ((double)p - f1[0]) + f2[0];
      freq10_12 = freq00_12 + f12_slope;
      freq11_02 = freq01_02 + f02_slope;

      line_down = fabs(freq11_02 - freq10_12);
      line_up = fabs(freq01_02 - freq00_12);
      denom = line_down + line_up;
      if (denom != 0) {
        area_down_triangle = line_down / denom * amp_section;
        area_up_triangle = line_up / denom * amp_section;
        triangle_interpolation1D(&freq00_12, &freq11_02, &freq10_12,
                                 &area_down_triangle, &spec[p * m1], &m1);
        triangle_interpolation1D(&freq00_12, &freq11_02, &freq01_02, &area_up_triangle,
                                 &spec[p * m1], &m1);
      } else {
        triangle_interpolation1D(&freq00_12, &freq11_02, &freq10_12, &amp_section,
                                 &spec[p * m1], &m1);
      }
      p++;
    }

    // Part 6: After mid to before end bin.
    diff += 0.5;
    diff *= df2;
    line_up = fabs(line_down);
    line_down += slope_diff;
    line_down = fabs(line_down);
    denom = line_up + line_down;
    while (p != pmax) {
      diff -= df2;
      amp_section = diff;
      freq00_12 = freq10_12;
      freq01_02 = freq11_02;
      freq10_12 += f12_slope;
      freq11_02 += f02_slope;

      if (denom != 0) {
        area_down_triangle = line_down / denom * amp_section;
        area_up_triangle = line_up / denom * amp_section;
        triangle_interpolation1D(&freq00_12, &freq11_02, &freq10_12,
                                 &area_down_triangle, &spec[p * m1], &m1);
        triangle_interpolation1D(&freq00_12, &freq11_02, &freq01_02, &area_up_triangle,
                                 &spec[p * m1], &m1);
      } else {
        triangle_interpolation1D(&freq00_12, &freq11_02, &freq10_12, &amp_section,
                                 &spec[p * m1], &m1);
      }
      line_up -= abs_slope_diff;
      line_down -= abs_slope_diff;
      denom -= 2 * abs_slope_diff;
      p++;
    }

    // Part 7: The end bin.
    if (clip_right2 == 0) {
      temp = (f1[2] - (double)p);
      amp_section = temp * temp * 0.5 * df2;
      freq00_12 = freq10_12;
      freq01_02 = freq11_02;
      freq11_02 = f2[2];
      triangle_interpolation1D(&freq00_12, &freq11_02, &freq01_02, &amp_section,
                               &spec[p * m1], &m1);
    }
  } else {
    if (clip_right2 == 0) {
      amp_section = f21 * top * 0.5;
      triangle_interpolation1D(&freq11_02, &f2[1], &f2[2], &amp_section, &spec[p * m1],
                               &m1);
    }
  }
  return 0;
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

void triangle_interpolation1D(double *freq1, double *freq2, double *freq3, double *amp,
                              double *spec, int *points) {
  if (fabs(*freq1 - *freq2) < TOL && fabs(*freq1 - *freq3) < TOL)
    return delta_fn_interpolation(freq1, points, amp, spec);

  __triangle_interpolation(freq1, freq2, freq3, amp, spec, points);
}

void octahedronDeltaInterpolation(const unsigned int nt, double *freq, double *amp,
                                  int stride, int n_spec, double *spec) {
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
  return delta_fn_interpolation(freq, &n_spec, &amp1, spec);
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

void octahedronInterpolation2D(double *spec, double *freq1, double *freq2, int nt,
                               double *amp, int stride, int m0, int m1) {
  int i = 0, j = 0, local_index, n_pts = (nt + 1) * (nt + 2) / 2;
  unsigned int int_i_stride = 0, int_j_stride = 0;
  double amp1, temp, *amp_address, *freq1_address, *freq2_address;

  /* Interpolate between 1d points by setting up triangles of unit area */

  local_index = nt - 1;
  amp_address = &amp[(nt + 1) * stride];
  freq1_address = &freq1[nt + 1];
  freq2_address = &freq2[nt + 1];

  while (i < n_pts - 1) {
    temp = amp[int_i_stride + stride] + amp_address[int_j_stride];
    amp1 = temp + amp[int_i_stride];

    triangle_interpolation2D(&freq1[i], &freq1[i + 1], &freq1_address[j], &freq2[i],
                             &freq2[i + 1], &freq2_address[j], &amp1, spec, m0, m1);

    if (i < local_index) {
      temp += amp_address[int_j_stride + stride];
      triangle_interpolation2D(&freq1[i + 1], &freq1_address[j], &freq1_address[j + 1],
                               &freq2[i + 1], &freq2_address[j], &freq2_address[j + 1],
                               &temp, spec, m0, m1);
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
