cimport clib as clib
cimport numpy as cnp
from numpy cimport ndarray
import numpy as np
import cython
import ctypes

__author__ = "Deepansh J. Srivastava"
__email__ = "deepansh2012@gmail.com"


@cython.profile(False)
@cython.boundscheck(False)
@cython.wraparound(False)
def histogram1d(
        int x_count,
        double x_min,
        double x_max,
        ndarray[double, ndim=1] sample_x,
        ndarray[double, ndim=1] weights,
    ):
    cdef ndarray[double] hist = np.zeros(x_count, dtype=np.float64)
    cdef int sample_count = sample_x.size
    cdef int c_int64 = ctypes.sizeof(ctypes.c_int64)
    cdef int stride_x = sample_x.ctypes.strides[0] / c_int64
    cdef int stride_w = weights.ctypes.strides[0] / c_int64

    clib.histogram1d_c(x_count, x_min, x_max, &hist[0], sample_count,
                       &sample_x[0], stride_x, &weights[0], stride_w)

    dx = 0.5 * (x_max - x_min) / x_count
    bins = np.linspace(x_min, x_max, x_count, endpoint=False) + dx
    return bins, hist


@cython.profile(False)
@cython.boundscheck(False)
@cython.wraparound(False)
def histogram2d(
        int x_count,
        double x_min,
        double x_max,
        int y_count,
        double y_min,
        double y_max,
        ndarray[double, ndim=1] sample_x,
        ndarray[double, ndim=1] sample_y,
        ndarray[double, ndim=1] weights,
        int interp=0,
    ):
    cdef ndarray[double, ndim=2] hist = np.zeros((x_count, y_count), dtype=np.float64)
    cdef int sample_count = sample_x.size
    cdef int c_int64 = ctypes.sizeof(ctypes.c_int64)
    cdef int stride_x = sample_x.ctypes.strides[0] / c_int64
    cdef int stride_y = sample_y.ctypes.strides[0] / c_int64
    cdef int stride_w = weights.ctypes.strides[0] / c_int64

    if interp == 1:
        clib.histogram2d_interp_c(x_count, x_min, x_max, y_count, y_min, y_max,
                                  &hist[0, 0], sample_count, &sample_x[0], stride_x,
                                  &sample_y[0], stride_y, &weights[0], stride_w)
    else:
        clib.histogram2d_c(x_count, x_min, x_max, y_count, y_min, y_max,
                           &hist[0, 0], sample_count, &sample_x[0], stride_x,
                           &sample_y[0], stride_y, &weights[0], stride_w)

    dx = 0.5 * (x_max - x_min) / x_count
    dy = 0.5 * (y_max - y_min) / y_count

    bins_x = np.linspace(x_min, x_max, x_count, endpoint=False) + dx
    bins_y = np.linspace(y_min, y_max, y_count, endpoint=False) + dy
    return bins_x, bins_y, hist
