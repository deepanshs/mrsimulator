# -*- coding: utf-8 -*-
from scipy.ndimage import gaussian_filter


def line_broadening(datum, sigma):
    # normalize = datum.max()
    return gaussian_filter(datum, sigma, 0)
    # return (datum / datum.max()) * normalize


if __name__ == "__main__":
    import numpy

    line_broadening(numpy.asarray([1, 2]), 2)
