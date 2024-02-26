import numpy as np


__author__ = "Deepansh J. Srivastava"
__email__ = "dsrivastava@hyperfine.io"


def czjzek(sigma, zeta, eta):
    """Analytical czjzek distribution

    Args:
        sigma: czjzek standard distribution
        zeta: zeta array
        eta: eta array
    """
    sigma_ = 2 * sigma
    denom = np.sqrt(2 * np.pi) * sigma_**5
    res = (zeta**4 * eta) * (1 - eta**2 / 9) / denom
    res *= np.exp(-(zeta**2 * (1 + (eta**2 / 3))) / (2 * sigma_**2))
    res /= res.sum()
    return res
