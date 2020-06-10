// -*- coding: utf-8 -*-
//
//  interpolation.c
//
//  @copyright Deepansh J. Srivastava, 2019-2020.
//  Created by Deepansh J. Srivastava, Apr 11, 2019.
//  Contact email = deepansh2012@gmail.com
//

#include "interpolation.h"

// triangle_interpolation is an optimized version of tent. Still plenty of
// room for optimization.

int triangle_interpolation(double *freq1, double *freq2, double *freq3,
                           double *amp, double *spec, int *points) {
  double df1, df2, top = 0.0, t, diff, f10 = 0.0, f21 = 0.0, temp;
  int p, pmid, pmax, i, j;
  int clip_right1 = 0, clip_left1 = 0, clip_right2 = 0, clip_left2 = 0;

  p = (int)(freq1[0]);
  if ((int)freq1[0] == (int)freq2[0] && (int)freq1[0] == (int)freq3[0]) {
    if (p >= points[0] || p < 0) {
      return 0;
    }
    spec[p] += amp[0];
    return 0;
  }

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

  p = (int)f[0];
  if (p > points[0]) {
    return 0;
  }

  pmax = (int)f[2];
  if (pmax < 0) {
    return 0;
  }

  pmid = (int)f[1];
  if (pmid >= points[0]) {
    pmid = points[0];
    clip_right1 = 1;
  }

  if (pmax >= points[0]) {
    pmax = points[0];
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

  top += (amp[0] * 2.0 / (f[2] - f[0]));
  f10 += f[1] - f[0];
  f21 += f[2] - f[1];

  if (p != pmid) {
    df1 = top / f10;
    diff = (double)p + 1. - f[0];
    if (clip_left1 == 0) {
      spec[p++] += 0.5 * diff * diff * df1;
    } else {
      spec[p++] += (diff - 0.5) * df1;
    }
    diff -= 0.5;
    diff *= df1;
    while (p != pmid) {
      diff += df1;
      spec[p++] += diff;
    }
    if (clip_right1 == 0) {
      spec[p] += (f[1] - (double)p) * (f10 + ((double)p - f[0])) * 0.5 * df1;
    }
  } else {
    if (clip_right1 == 0 && clip_left1 == 0) {
      spec[p] += f10 * top * 0.5;
    }
  }

  if (p != pmax) {
    df2 = top / f21;
    diff = f[2] - (double)p - 1.;

    if (clip_left2 == 0) {
      spec[p++] += (f21 - diff) * (diff + f21) * 0.5 * df2;
    } else {
      spec[p++] += (diff + 0.5) * df2;
    }
    diff += 0.5;
    diff *= df2;
    while (p != pmax) {
      diff -= df2;
      spec[p++] += diff;
    }
    if (clip_right2 == 0) {
      temp = (f[2] - (double)p);
      spec[p] += temp * temp * 0.5 * df2;
    }
  } else {
    if (clip_right2 == 0) {
      spec[p] += f21 * top * 0.5;
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
