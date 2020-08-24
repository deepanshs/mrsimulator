# -*- coding: utf-8 -*-
import numpy as np


def get_principal_components(zeta, eta):
    """
    Return the principal components of a traceless second-rank symmetric
    Cartesian tensor.

    Args:
        zeta: The zeta parameter in PAS, according to the Haeberlen convention.
        eta: The eta parameter in PAS, according to the Haeberlen convention.
    """
    xx = -0.5 * zeta * (eta + 1.0)
    yy = 0.5 * zeta * (eta - 1.0)
    zz = zeta

    return [xx, yy, zz]


def get_Haeberlen_components(tensors):
    """Return zeta and eta parameters of the tensor using the Haeberlen convention.

    Args:
        ndarray tensors: A `N x 3 x 3` ndarray of `N` traceless symmetric second-rank
            Cartesian tensors.
    """
    n = tensors.shape[0]
    eig_val = np.linalg.eigvalsh(tensors)
    eig_val_sort_ = np.argsort(np.abs(eig_val), axis=1, kind="mergesort")
    eig_val_sort_ = (eig_val_sort_.T + 3 * np.arange(n)).T.ravel()
    eig_val_sorted = eig_val.ravel()[eig_val_sort_].reshape(n, 3)

    eig_val_sort_ = eig_val = None
    del eig_val_sort_, eig_val

    zeta = eig_val_sorted[:, -1]
    eta = (eig_val_sorted[:, 0] - eig_val_sorted[:, 1]) / zeta

    return zeta, eta


def x_y_from_zeta_eta(zeta, eta):
    """Convert the zeta, eta coordinates from the Haeberlen convention to the
    x-y notation."""
    xa = np.empty(zeta.size)
    ya = np.empty(zeta.size)

    index = np.where(zeta >= 0)
    temp = np.tan(0.7853981634 * eta[index])
    ya[index] = np.sqrt(zeta[index] * zeta[index] / (temp * temp + 1.0))
    xa[index] = temp * ya[index]

    index = np.where(zeta < 0)
    temp = np.tan(0.7853981634 * eta[index])
    xa[index] = np.sqrt(zeta[index] * zeta[index] / (temp * temp + 1.0))
    ya[index] = temp * xa[index]

    zeta = eta = None
    del zeta, eta

    return xa, ya
