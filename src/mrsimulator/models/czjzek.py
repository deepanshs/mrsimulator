# -*- coding: utf-8 -*-
import numpy as np

from .utils import get_Haeberlen_components
from .utils import get_principal_components


def _czjzek_random_distribution_tensors(sigma, n):
    r"""Czjzek random distribution model.

    Args:
        float sigma: The standard deviation of the five-dimensional multi-variate normal
            distribution.
        int n: Number of samples drawn from the Czjzek random distribution model.

    Description
    -----------

    U is an array of the coordinates randomly drawn from an uncorrelated five-dimensional
    multivariate normal distribution with standard deviation `sigma` and zero mean.

    The components of the traceless second-rank symmetric cartesian tensor, S_ij, follows

    Sxx = sqrt(3) * U5 - U1

    Sxy = Syx = sqrt(3) * U4

    Syy = -sqrt(3) * U5 - U1

    Sxz = Szx = sqrt(3) * U2

    Szz = 2 * U1

    Syz = Szy = sqrt(3) * U3
    """

    # The random sampling U1, U2, ... U5
    U1 = np.random.normal(0.0, sigma, n)

    sqrt_3_sigmal = np.sqrt(3) * sigma
    sqrt_3_U2 = np.random.normal(0.0, sqrt_3_sigmal, n)
    sqrt_3_U3 = np.random.normal(0.0, sqrt_3_sigmal, n)
    sqrt_3_U4 = np.random.normal(0.0, sqrt_3_sigmal, n)
    sqrt_3_U5 = np.random.normal(0.0, sqrt_3_sigmal, n)

    # Create N random tensors
    tensors = np.zeros((n, 3, 3))  # n x 3 x 3 tensors

    tensors[:, 0, 0] = sqrt_3_U5 - U1  # xx
    tensors[:, 0, 1] = sqrt_3_U4  # xy
    tensors[:, 0, 2] = sqrt_3_U2  # xz

    tensors[:, 1, 0] = sqrt_3_U4  # yx
    tensors[:, 1, 1] = -sqrt_3_U5 - U1  # yy
    tensors[:, 1, 2] = sqrt_3_U3  # yz

    tensors[:, 2, 0] = sqrt_3_U2  # zx
    tensors[:, 2, 1] = sqrt_3_U3  # zy
    tensors[:, 2, 2] = 2 * U1  # zz

    return tensors


def czjzek_distribution(sigma, n):
    r"""Draw `N` samples of zeta and eta from the Czjzek distribution model.

    Args:
        float sigma: The Gaussian standard deviation.
        int n:  Number of samples drawn from the Czjzek distribution model.

    Description
    -----------

    The Czjzek distribution model is a random sampling of second-rank traceless
    symmetric tensors whose explicit matrix form follows

    .. math::
        {\bf S} = \left[
        \begin{array}{l l l}
        \sqrt{3} U_5 - U_1   & \sqrt{3} U_4          & \sqrt{3} U_2 \\
        \sqrt{3} U_4         & -\sqrt{3} U_5 - U_1   & \sqrt{3} U_3 \\
        \sqrt{3} U_2         & \sqrt{3} U_3          & 2 U_1
        \end{array}
        \right],

    where the components, $U_i$, are randomly drawn from a five-dimensional
    multivariate normal distribution.
    Each component, $U_i$, is a dimension of the five-dimensional uncorrelated
    multivariate normal distribution with the mean of :math:`<U_i>=0` and the variance
    :math:`<U_iU_i>=\sigma^2`. Because the :math:`U_i`s are uncorrelated, the
    covariance matrix, :math:`\Lambda=\sigma^2 {\bf I}_5`, where :math:`{\bf I}_5` is a
    5 x 5 identity matrix.

    .. math::
        S_T = S_C(\sigma),

    Example
    -------

    >>> Vzz_dist, eta_dist = czjzek_distribution(0.5, n=1000000)
    """
    tensors = _czjzek_random_distribution_tensors(sigma, n)
    return get_Haeberlen_components(tensors)


def extended_czjzek_distribution(zeta, eta, eps, n):
    r"""Draw `N` samples of zeta and eta from the extended czjzek distribution model.

    Description
    -----------

    The extended Czjzek random distribution model is an extension of the Czjzek model,
    given as

    .. math::
        S_T = S(0) + \rho S_C(\sigma=1),

    where :math:`S_T` is the total tensor, :math:`S(0)` is the dominant tensor,
    :math:`S_C(\sigma=1)` is the Czjzek random model attributing to the random
    perturbation of the tensor about the dominant tensor, :math:`S(0)`, and
    :math:`\rho` is the size of the perturbation. Note, in the above equation, the
    :math:`\sigma` parameter from the Czjzek random model, :math:`S_C`, has no meaning
    and is set to one. The factor, :math:`\rho`, is defined as

    .. math::
        \rho = \frac{||S(0)|| \epsilon}{\sqrt{30}},

    where :math:`\|S(0)\|` is the 2-norm of the tensor, and :math:`\epsilon` is a
    fraction.

    Args:
        float zeta: The size of the anisotropy of the dominant tensor.
        float eta: The asymmetry parameter of the dominant tensor.
        float eps: A fraction determining the extent of perturbation.
        int n: Number of samples drawn from the extended Czjzek distribution model.

    Reference
    ---------

    Gérard Le Caër, Bruno Bureau, and Dominique Massiot,
    An extension of the Czjzek model for the distributions of electric field gradients
    in disordered solids and an application to NMR spectra of 71Ga in chalcogenide
    glasses. Journal of Physics: Condensed Matter, 2010, 22, 065402.
    DOI: 10.1088/0953-8984/22/6/065402

    Example
    -------

    >>> from mrsimulator.models import extended_czjzek_distribution
    >>> Vzz_dist, eta_dist = extended_czjzek_distribution(1, 0.3, eps=0.35, n=1000000)
    """
    # czjzek_random_distribution model
    tensors = _czjzek_random_distribution_tensors(1, n)

    # the traceless second-rank symmetric Cartesian tensor in PAS
    T0 = [0.0, 0.0, 0.0]
    if zeta != 0 and eta != 0:
        T0 = get_principal_components(zeta, eta)

    # 2-norm of the tensor
    norm_T0 = np.linalg.norm(T0)

    # the perturbation factor
    rho = eps * norm_T0 / np.sqrt(30)

    # total tensor
    total_tensors = np.diag(T0) + rho * tensors

    return get_Haeberlen_components(total_tensors)
