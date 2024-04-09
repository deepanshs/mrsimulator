import numpy as np
from mrsimulator.base_model import get_Haeberlen_components as haeberlen_c
from mrsimulator.models.czjzek import _czjzek_random_distribution_tensors
from mrsimulator.models.czjzek import get_p_q_basis
from mrsimulator.models.utils import get_Haeberlen_components


def test_roots():
    n = 400_000

    tensors = _czjzek_random_distribution_tensors(sigma=1.0, n=n)
    u1 = tensors[:, 2, 2] / 2.0
    u2 = tensors[:, 2, 0]
    u3 = tensors[:, 2, 1]
    u4 = tensors[:, 1, 0]
    u5 = (tensors[:, 0, 0] - tensors[:, 1, 1]) / 2.0

    p_q_basis = get_p_q_basis(u1, u2, u3, u4, u5)

    zeta = 1.4e6
    eta = 0.7
    eps = 0.3

    T0 = [-0.5 * zeta * (eta + 1.0), 0.5 * zeta * (eta - 1.0), zeta]
    norm_T0 = np.linalg.norm(T0)
    rho = eps * norm_T0 / 5.4772255751

    zeta_eig, eta_eig = get_Haeberlen_components(np.diag(T0) + tensors * rho)

    zeta_c, eta_c = haeberlen_c(*p_q_basis, zeta, eta, rho)

    assert np.allclose(zeta_eig, zeta_c)
    assert np.allclose(eta_eig, eta_c)
