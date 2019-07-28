import numpy as np
from . import transition_function as tf

__author__ = "Deepansh J. Srivastava"
__email__ = ["srivastava.89@osu.edu", "deepansh2012@gmail.com"]

__all__ = ["nuclear_shielding"]


def nuclear_shielding(iso, zeta, eta, order=1):
    # scale = tf.p(transition[1], transition[0])
    R0 = iso  # * scale

    R2 = np.zeros(5)
    R4 = np.zeros(9)
    temp = -0.4082482905 * (zeta * eta)  # * scale
    R2[0] += temp  # R2-2
    R2[1] += 0.0  # R2-1
    R2[2] += zeta  # * scale  # R2 0
    R2[3] += 0.0  # R2 1
    R2[4] += temp  # R2 2

    return R0, R2.astype(np.complex128), R4.astype(np.complex128)


def quadrupolar(spin, Cq, eta, order=1, vo=None):
    vq = 3.0 * Cq
    denominator = 2.0 * spin * (2.0 * spin - 1.0)
    vq /= denominator

    R2 = np.zeros(5)
    R4 = np.zeros(9)

    if order == 1:
        R0 = 0

        temp = -0.1666666667 * (vq * eta)
        R2[0] = temp  # R2-2
        R2[1] = 0.0  # R2-1
        R2[2] = 0.4082482905 * vq  # R2 0
        R2[3] = 0.0  # R2 1
        R2[4] = temp  # R2 2

    if order == 2:
        scale = vq * vq / vo
        eta2 = eta * eta

        R0 = (eta2 * 0.33333333333 + 1.0) * 0.07453559925 * scale

        temp = -eta * 0.07273929675 * scale
        R2[0] = temp  # R2-2
        R2[1] = 0.0  # R2-1
        R2[2] = 0.08908708064 * (eta2 * 0.33333333333 - 1.0) * scale  # R2 0
        R2[3] = 0.0  # R2 1
        R2[4] = temp  # R2 2

        temp = eta2 * 0.02777777778 * scale
        temp2 = -0.06299407883 * eta * scale
        R4[0] = temp  # R4-4
        R4[2] = temp2  # R4-2
        R4[4] = 0.1195228609 * (eta2 * 0.05555555556 + 1.0) * scale  # R4 0
        R4[6] = temp2  # R4 2
        R4[8] = temp  # R4 4

    return R0, R2.astype(np.complex128), R4.astype(np.complex128)
