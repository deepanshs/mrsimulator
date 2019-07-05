import numpy as np


def _p_(mf, mi):
    return mf - mi


def _d_(mf, mi):
    return 1.2247448714 * (mf * mf - mi * mi)


def _dIS_(mIf, mIi, mSf, mSi):
    return mIf * mSf - mIi * mSi


def _f_(mf, mi, spin):
    f = 1.0 - 3.0 * spin * (spin + 1.0)
    f *= mf - mi
    f += 5.0 * (mf * mf * mf - mi * mi * mi)
    f *= 0.316227766
    return f


def _qqad_ci_(mf, mi, spin):
    f = _f_(mf, mi, spin)
    p = _p_(mf, mi)
    temp = spin * (spin + 1.0) - 0.75
    c0 = 0.3577708764 * temp * p + 0.8485281374 * f
    c2 = 0.1069044968 * temp * p + -1.0141851057 * f
    c4 = -0.1434274331 * temp * p + -1.2850792082 * f
    return c0, c2, c4


def nuclear_shielding(iso, zeta, eta, transition, order=1):
    scale = _p_(transition[1], transition[0])
    R0 = iso * scale

    R2 = np.zeros(5)
    temp = -0.4082482905 * (zeta * eta) * scale
    R2[0] += temp  # R2-2
    R2[1] += 0.0  # R2-1
    R2[2] += zeta * scale  # R2 0
    R2[3] += 0.0  # R2 1
    R2[4] += temp  # R2 2

    return R0, R2.astype(np.complex128)
