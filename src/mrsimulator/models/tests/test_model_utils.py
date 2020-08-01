# -*- coding: utf-8 -*-
import numpy as np
from mrsimulator.models.utils import get_Haeberlen_components
from mrsimulator.models.utils import get_principal_components
from mrsimulator.models.utils import x_y_from_zeta_eta


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
