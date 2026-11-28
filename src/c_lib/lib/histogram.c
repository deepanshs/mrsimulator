// -*- coding: utf-8 -*-
//
//  histogram.h
//
//  @copyright Deepansh J. Srivastava, 2019-2024.
//  Created by Deepansh J. Srivastava, Mar 6, 2024.
//  Contact email = deepansh2012@gmail.com
//

#include "histogram.h"

#include <math.h>
#include <stdbool.h>
#include <stdio.h>

void histogram1d_c(int x_count, const double x_min, const double x_max,
                   double *restrict hist, int sample_count,
                   const double *restrict sample_x, const int stride_x,
                   const double *restrict weights, const int stride_w) {
  int ix, counter = sample_count;
  double temp_sample, normx;
  double *x_pt = sample_x, *w_pt = weights;

  normx = (double)x_count / (x_max - x_min);

  while (counter-- > 0) {
    temp_sample = *x_pt;
    if (temp_sample >= x_min && temp_sample < x_max) {
      ix = (int)((temp_sample - x_min) * normx);
      if (ix >= x_count) ix = x_count - 1;  // round-off just below x_max
      hist[ix] += *w_pt;
    }
    x_pt += stride_x;
    w_pt += stride_w;
  }
}

void histogram2d_c(int x_count, const double x_min, const double x_max, int y_count,
                   const double y_min, const double y_max, double *restrict hist,
                   int sample_count, const double *restrict sample_x,
                   const int stride_x, const double *restrict sample_y,
                   const int stride_y, const double *restrict weights,
                   const int stride_w) {
  int ix, iy, hist_index, counter = sample_count;
  double temp_sample_x, temp_sample_y, normx, normy;
  double *x_pt = sample_x, *y_pt = sample_y, *w_pt = weights;

  normx = (double)(x_count) / (x_max - x_min);
  normy = (double)(y_count) / (y_max - y_min);

  while (counter-- > 0) {
    temp_sample_x = *x_pt;
    if (temp_sample_x >= x_min && temp_sample_x < x_max) {
      temp_sample_y = *y_pt;
      if (temp_sample_y >= y_min && temp_sample_y < y_max) {
        ix = (int)((temp_sample_x - x_min) * normx);
        iy = (int)((temp_sample_y - y_min) * normy);
        if (ix >= x_count) ix = x_count - 1;  // round-off just below x_max
        if (iy >= y_count) iy = y_count - 1;  // round-off just below y_max
        hist_index = iy + y_count * ix;
        hist[hist_index] += *w_pt;
      }
    }
    x_pt += stride_x;
    y_pt += stride_y;
    w_pt += stride_w;
  }
}

/**
 * Bilinear deposit of the samples onto the bin centres. A sample between the outer bin
 * centre and the edge of the range is deposited onto the outer bin.
 */
void histogram2d_interp_c(int x_count, const double x_min, const double x_max,
                          int y_count, const double y_min, const double y_max,
                          double *restrict hist, int sample_count,
                          const double *restrict sample_x, const int stride_x,
                          const double *restrict sample_y, const int stride_y,
                          const double *restrict weights, const int stride_w) {
  int ix0, iy0, ix1, iy1, counter = sample_count;
  double temp_sample_x, temp_sample_y, normx, normy, fx, fy, dx, dy, w;
  const double *x_pt = sample_x, *y_pt = sample_y, *w_pt = weights;

  normx = (double)(x_count) / (x_max - x_min);
  normy = (double)(y_count) / (y_max - y_min);

  while (counter-- > 0) {
    temp_sample_x = *x_pt;
    temp_sample_y = *y_pt;
    if (temp_sample_x >= x_min && temp_sample_x < x_max && temp_sample_y >= y_min &&
        temp_sample_y < y_max) {
      // position in units of bins, relative to the centre of the first bin
      fx = (temp_sample_x - x_min) * normx - 0.5;
      fy = (temp_sample_y - y_min) * normy - 0.5;
      ix0 = (int)floor(fx);
      iy0 = (int)floor(fy);
      dx = fx - ix0;
      dy = fy - iy0;

      // clamp the neighbours to the range
      ix1 = ix0 + 1;
      iy1 = iy0 + 1;
      if (ix0 < 0) ix0 = 0;
      if (iy0 < 0) iy0 = 0;
      if (ix1 > x_count - 1) ix1 = x_count - 1;
      if (iy1 > y_count - 1) iy1 = y_count - 1;
      if (ix0 > x_count - 1) ix0 = x_count - 1;
      if (iy0 > y_count - 1) iy0 = y_count - 1;

      w = *w_pt;
      hist[iy0 + y_count * ix0] += w * (1.0 - dx) * (1.0 - dy);
      hist[iy1 + y_count * ix0] += w * (1.0 - dx) * dy;
      hist[iy0 + y_count * ix1] += w * dx * (1.0 - dy);
      hist[iy1 + y_count * ix1] += w * dx * dy;
    }
    x_pt += stride_x;
    y_pt += stride_y;
    w_pt += stride_w;
  }
}
