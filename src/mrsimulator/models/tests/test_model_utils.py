import numpy as np
from mrsimulator.models.utils import get_Haeberlen_components
from mrsimulator.models.utils import get_principal_components
from mrsimulator.models.utils import x_y_from_zeta_eta
from mrsimulator.models.utils import x_y_to_zeta_eta


__author__ = "Deepansh J. Srivastava"
__email__ = "srivastava.89@osu.edu"


def test_get_principal_components():
    # 1
    zeta = 30
    eta = 0
    components = get_principal_components(zeta, eta)
    assert components == [-15, -15, 30]

    # 2
    zeta = 30
    eta = 1
    components = get_principal_components(zeta, eta)
    assert components == [-30, 0, 30]

    # 3
    zeta = 30
    eta = 0.5
    components = get_principal_components(zeta, eta)
    assert components == [-22.5, -7.5, 30]


def test_get_Haeberlen_components():
    # 1
    pas = [-15, -15, 30]
    tensors = np.diag(pas)[np.newaxis, :, :]
    components = get_Haeberlen_components(tensors)
    assert components == (30, 0)

    # 2
    pas = [-30, 0, 30]
    tensors = np.diag(pas)[np.newaxis, :, :]
    components = get_Haeberlen_components(tensors)
    assert components == (30, 1)

    # 3
    pas = [-22.5, -7.5, 30]
    tensors = np.diag(pas)[np.newaxis, :, :]
    components = get_Haeberlen_components(tensors)
    assert components == (30, 0.5)


def test_x_y_from_zeta_eta():
    zeta = np.random.rand(20) * 40
    eta = np.random.rand(20)
    x, y = x_y_from_zeta_eta(zeta, eta)

    theta = np.pi * eta / 4.0
    x_ = zeta * np.sin(theta)
    y_ = zeta * np.cos(theta)

    assert np.allclose(x, x_)
    assert np.allclose(y, y_)


def test_x_y_to_zeta_eta():
    x = np.random.rand(16) * 3000
    y = np.random.rand(16) * 3000
    x[-1] = y[-1] = 56.0
    x[-2] = y[-2] = 0.0
    factor_ = 4 / np.pi
    zeta_ = []
    eta_ = []
    for x_, y_ in zip(x, y):
        z = np.sqrt(x_**2 + y_**2)
        if x_ < y_:
            eta_.append(factor_ * np.arctan(x_ / y_))
            zeta_.append(z)
        elif x_ > y_:
            eta_.append(factor_ * np.arctan(y_ / x_))
            zeta_.append(-z)
        else:
            zeta_.append(z)
            eta_.append(1.0)

        z_temp, e_temp = x_y_to_zeta_eta([x_], [y_])
        assert np.allclose(zeta_[-1], z_temp)
        assert np.allclose(eta_[-1], e_temp)

    zeta, eta = x_y_to_zeta_eta(x, y)
    assert np.allclose(zeta, np.asarray(zeta_))
    assert np.allclose(eta, np.asarray(eta_))
