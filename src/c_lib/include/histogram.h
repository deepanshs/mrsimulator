// -*- coding: utf-8 -*-
//
//  histogram.h
//
//  @copyright Deepansh J. Srivastava, 2019-2024.
//  Created by Deepansh J. Srivastava, Mar 6, 2024.
//  Contact email = deepansh2012@gmail.com
//

extern void histogram1d_c(int x_count, const double x_min, const double x_max,
                          double *restrict hist, int sample_count,
                          const double *restrict sample_x, const int stride_x,
                          const double *restrict weights, const int stride_w);

extern void histogram2d_c(int x_count, const double x_min, const double x_max,
                          int y_count, const double y_min, const double y_max,
                          double *restrict hist, int sample_count,
                          const double *restrict sample_x, const int stride_x,
                          const double *restrict sample_y, const int stride_y,
                          const double *restrict weights, const int stride_w);

extern void histogram2d_interp_c(int x_count, const double x_min, const double x_max,
                                 int y_count, const double y_min, const double y_max,
                                 double *restrict hist, int sample_count,
                                 const double *restrict sample_x, const int stride_x,
                                 const double *restrict sample_y, const int stride_y,
                                 const double *restrict weights, const int stride_w);
