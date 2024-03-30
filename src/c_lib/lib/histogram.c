// -*- coding: utf-8 -*-
//
//  histogram.h
//
//  @copyright Deepansh J. Srivastava, 2019-2024.
//  Created by Deepansh J. Srivastava, Mar 6, 2024.
//  Contact email = deepansh2012@gmail.com
//

#include "histogram.h"

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

// void histogram2d_linterp(int x_count, const double hist_x_min, const double
// hist_x_max,
//                  int y_count, const double hist_y_min, const double hist_y_max,
//                  double *restrict hist, int sample_count,
//                  const double *restrict sample_x, const double *restrict sample_y,
//                  const double *restrict weights) {
//   int i, ix, iy, hist_index;
//   double temp_sample_x, temp_sample_y, normx, normy;
//   double interval_x, interval_y;

//   interval_x = (hist_x_max - hist_x_min) / (x_count - 1);
//   interval_y = (hist_y_max - hist_y_min) / (y_count - 1);

//   normx = 1.0 / interval_x;
//   normy = 1.0 / interval_y;

//   for (i = 0; i <= sample_count - 1; i++) {
//     temp_sample_x = sample_x[i];
//     temp_sample_y = sample_y[i];

//     if (temp_sample_x >= hist_x_min &&
//         temp_sample_x < hist_x_max &&
//         temp_sample_y >= hist_y_min &&
//         temp_sample_y < hist_y_max) {
//       ix = (int)((temp_sample_x - hist_x_min) * normx);
//       iy = (int)((temp_sample_y - hist_y_min) * normy);
//       hist_index = iy + y_count * ix;
//       hist[hist_index] += weights[i];
//     }
//   }
// }
