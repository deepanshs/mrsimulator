// -*- coding: utf-8 -*-
//
//  histogram.h
//
//  @copyright Deepansh J. Srivastava, 2019-2024.
//  Created by Deepansh J. Srivastava, Mar 6, 2024.
//  Contact email = deepansh2012@gmail.com
//

#include "histogram.h"

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
        hist_index = iy + y_count * ix;
        hist[hist_index] += *w_pt;
      }
    }
    x_pt += stride_x;
    y_pt += stride_y;
    w_pt += stride_w;
  }
}

void histogram2d_interp_c(int x_count, const double x_min, const double x_max,
                          int y_count, const double y_min, const double y_max,
                          double *restrict hist, int sample_count,
                          const double *restrict sample_x, const int stride_x,
                          const double *restrict sample_y, const int stride_y,
                          const double *restrict weights, const int stride_w) {
  int ix, iy, p, counter = sample_count;
  bool left;
  double temp_sample_x, temp_sample_y, normx, normy;
  double *x_pt = sample_x, *y_pt = sample_y, *w_pt = weights;
  double delta_x, delta_y, x_min_adj, y_min_adj;

  normx = (double)(x_count) / (x_max - x_min);
  normy = (double)(y_count) / (y_max - y_min);

  x_min_adj = x_min;  // - 1.0 / normx;
  y_min_adj = y_min;  // - 1.0 / normy;

  while (counter-- > 0) {
    temp_sample_x = *x_pt;
    if (temp_sample_x >= x_min_adj && temp_sample_x < x_max) {
      temp_sample_y = *y_pt;
      if (temp_sample_y >= y_min_adj && temp_sample_y < y_max) {
        ix = (int)((temp_sample_x - x_min) * normx);
        iy = (int)((temp_sample_y - y_min) * normy);
        p = iy + y_count * ix;

        delta_y = temp_sample_y - iy / normy;
        delta_x = temp_sample_x - ix / normx;

        left = (delta_y) < 0;
        if (left) {
          if (iy != 0) hist[p - 1] -= *w_pt * delta_y;
          hist[p] += *w_pt * (1.0 + delta_y);
        } else {
          if (iy != y_count) hist[p + 1] += *w_pt * delta_y;
          hist[p] += *w_pt * (1.0 - delta_y);
        }

        left = (delta_x) < 0;
        if (left) {
          if (iy != 0) hist[p - y_count] -= *w_pt * delta_x;
          hist[p] += *w_pt * (1.0 + delta_x);
        } else {
          if (iy != y_count) hist[p + y_count] += *w_pt * delta_x;
          hist[p] += *w_pt * (1.0 - delta_x);
        }
      }
    }
    x_pt += stride_x;
    y_pt += stride_y;
    w_pt += stride_w;
  }
}
