# -*- coding: utf-8 -*-
import numpy as np

# from . import transition_function as tf

__author__ = "Deepansh J. Srivastava"
__email__ = ["srivastava.89@osu.edu", "deepansh2012@gmail.com"]

__all__ = ["nuclear_shielding"]


def nuclear_shielding(iso, zeta, eta, order=1):
    # scale = tf.p(transition[1], transition[0])
    R0 = iso  # * scale

    R2 = np.zeros(5)
    temp = -0.4082482905 * (zeta * eta)  # * scale
    R2[0] += temp  # R2-2
    R2[1] += 0.0  # R2-1
    R2[2] += zeta  # * scale  # R2 0
    R2[3] += 0.0  # R2 1
    R2[4] += temp  # R2 2

    return R0, R2.astype(np.complex128)
