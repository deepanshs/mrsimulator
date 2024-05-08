# -*- coding: utf-8 -*-
#
#  nmr_method.pxd
#
#  @copyright Deepansh J. Srivastava, 2019-2024.
#  Created by Deepansh J. Srivastava.
#  Contact email = deepansh2012@gmail.com
#


cdef extern from "histogram.h":
    void histogram1d_c(
        int x_count,
        const double x_min,
        const double x_max,
        double *hist,
        int sample_count,
        const double *sample_x,
        const int stride_x,
        const double *weights,
        const int stride_w,
    )

    void histogram2d_c(
        int x_count,
        const double x_min,
        const double x_max,
        int y_count,
        const double y_min,
        const double y_max,
        double *hist,
        int sample_count,
        const double *sample_x,
        const int stride_x,
        const double *sample_y,
        const int stride_y,
        const double *weights,
        const int stride_w,
    )

    void histogram2d_interp_c(
        int x_count,
        const double x_min,
        const double x_max,
        int y_count,
        const double y_min,
        const double y_max,
        double *hist,
        int sample_count,
        const double *sample_x,
        const int stride_x,
        const double *sample_y,
        const int stride_y,
        const double *weights,
        const int stride_w,
    )
