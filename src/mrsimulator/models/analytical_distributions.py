import numpy as np
from mrsimulator.clib import histogram2d
from mrsimulator.models.utils import zeta_eta_to_x_y

__author__ = "Deepansh J. Srivastava"
__email__ = "dsrivastava@hyperfine.io"


def czjzek(sigma: float, pos: list, polar: bool):
    """Analytical czjzek distribution on polar or non-polar gird

    Args:
        sigma: czjzek standard distribution
        pos: zeta-eta or x-y array based on polar=False ot True, respectively
        polar: grid type
    """
    if polar:
        return czjzek_polar(sigma, pos)
    return czjzek_zeta_eta(sigma, pos)


def czjzek_distribution(sigma: float, zeta: np.ndarray, eta: np.ndarray):
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


def czjzek_zeta_eta(sigma: float, pos: list):
    """Haeberlen czjzek distribution

    Args:
        sigma: standard deviation
        pos: x and y numpy arrays
    """
    sub_zero = np.where(pos[1] < 0)[0]
    super_one = np.where(pos[1] > 1)[0]
    zeta, eta = np.meshgrid(pos[0], pos[1])
    pdf_model = czjzek_distribution(sigma, zeta, eta)
    pdf_model[sub_zero] = 0
    pdf_model[super_one] = 0
    eta_idx = np.where(eta.ravel() == 1)
    pdf_model = pdf_model.ravel()
    pdf_model[eta_idx] /= 2.0
    pdf_model.shape = zeta.shape
    return pos[0], pos[1], pdf_model


def czjzek_polar(sigma: float, pos: list):
    """Polar czjzek distribution

    Args:
        sigma: standard deviation
        pos: x and y numpy arrays
    """
    bins = [pos[0].size, pos[1].size]

    max_val = np.sqrt(pos[0][-1] ** 2 + pos[1][-1] ** 2)
    max_size = int(np.max(bins) * np.sqrt(2)) * 8
    zeta_1d = 2 * max_val * (np.arange(max_size) / max_size)
    eta_1d = np.arange(max_size) / (max_size - 1)
    zeta, eta = np.meshgrid(zeta_1d, eta_1d)

    pdf_model = czjzek_distribution(sigma, zeta, eta)
    eta = eta.ravel()
    zeta = zeta.ravel()
    x, y = zeta_eta_to_x_y(zeta, eta)

    eta_idx = np.where(eta == 1)
    pdf_model = pdf_model.ravel()
    pdf_model[eta_idx] /= 2.0

    delta_z = (pos[0][1] - pos[0][0]) / 2
    delta_e = (pos[1][1] - pos[1][0]) / 2
    range_x = [pos[0][0] - delta_z, pos[0][-1] + delta_z]
    range_y = [pos[1][0] - delta_e, pos[1][-1] + delta_e]

    _, _, hist_x_y = histogram2d(
        sample_x=x,
        sample_y=y,
        weights=pdf_model,
        x_count=bins[0],
        y_count=bins[1],
        x_min=range_x[0],
        x_max=range_x[1],
        y_min=range_y[0],
        y_max=range_y[1],
    )
    hist_x_y += hist_x_y.T
    hist_x_y /= hist_x_y.sum()
    return pos[0], pos[1], hist_x_y
