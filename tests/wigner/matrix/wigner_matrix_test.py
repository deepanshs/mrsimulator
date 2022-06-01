# -*- coding: utf-8 -*-
import mrsimulator.tests.tests as clib
import numpy as np
from sympy.physics.quantum.spin import Rotation

from tests.python_test_for_c_code.angular_momentum import wigner_d_matrix_cosines
from tests.python_test_for_c_code.angular_momentum import wigner_dm0_vector


def wigner_dm0_vector_sympy(ang_momentum_l, angle):
    R_out = np.empty(2 * ang_momentum_l + 1, dtype=np.float64)
    for i in range(2 * ang_momentum_l + 1):
        rotate = Rotation.d(ang_momentum_l, -ang_momentum_l + i, 0, angle).doit()
        R_out[i] = complex(rotate).real
    return R_out


def wigner(ang_momentum_l, cos_beta):
    # python test
    cos_beta = np.asarray([cos_beta]).ravel()
    sin_beta = np.sqrt(1.0 - cos_beta**2)
    exp_I_beta = cos_beta + 1j * sin_beta
    wigner_py = wigner_d_matrix_cosines(ang_momentum_l, cos_beta)

    # c test
    wigner_c = clib.wigner_d_matrices_from_exp_I_beta(ang_momentum_l, False, exp_I_beta)
    wigner_c_half = clib.wigner_d_matrices_from_exp_I_beta(
        ang_momentum_l, True, exp_I_beta
    )
    return wigner_py, wigner_c, wigner_c_half


# All wigner matrix are tested against the wigner matrix computed using
# Sympy Rotation methods. See the main function in this file.


def test_wigner_2j_matrix_angle_00():
    ang_momentum_l = 2
    beta = np.arccos(0.5)
    array = np.load("tests/wigner/matrix/l=2_cx=0.5.npy")
    wigner = clib.wigner_d_matrices(ang_momentum_l, np.asarray([beta]))
    np.testing.assert_almost_equal(array, wigner, decimal=8)


def test_wigner_2j_matrix_cosine_00():
    ang_momentum_l = 2
    cos_beta = 0.5
    wigner_py, wigner_c, wigner_c_half = wigner(ang_momentum_l, cos_beta)
    array = np.load("tests/wigner/matrix/l=2_cx=0.5.npy")
    np.testing.assert_almost_equal(array, wigner_py.ravel(), decimal=8)
    np.testing.assert_almost_equal(array, wigner_c.ravel(), decimal=8)
    np.testing.assert_almost_equal(wigner_c[:, :3, :], wigner_c_half, decimal=8)


def test_wigner_2j_matrix_angle_01():
    ang_momentum_l = 2
    cos_beta = [-0.5498, 0.230]
    beta = np.arccos(cos_beta)
    array = np.load("tests/wigner/matrix/l=2_cx=[-0.5498, 0.230].npy")
    wigner = clib.wigner_d_matrices(ang_momentum_l, beta)
    np.testing.assert_almost_equal(array, wigner, decimal=8)


def test_wigner_2j_matrix_cosine_01():
    ang_momentum_l = 2
    cos_beta = [-0.5498, 0.230]
    wigner_py, wigner_c, wigner_c_half = wigner(ang_momentum_l, cos_beta)
    array = np.load("tests/wigner/matrix/l=2_cx=[-0.5498, 0.230].npy")
    np.testing.assert_almost_equal(array, wigner_py.ravel(), decimal=8)
    np.testing.assert_almost_equal(array, wigner_c.ravel(), decimal=8)
    np.testing.assert_almost_equal(wigner_c[:, :3, :], wigner_c_half, decimal=8)


def test_wigner_4j_matrix_angle_02():
    ang_momentum_l = 4
    cos_beta = [-0.8459]
    beta = np.arccos(cos_beta)
    array = np.load("tests/wigner/matrix/l=4_cx=-0.8459.npy")
    wigner = clib.wigner_d_matrices(ang_momentum_l, beta)
    np.testing.assert_almost_equal(array, wigner, decimal=8)


def test_wigner_4j_matrix_cosine_02():
    ang_momentum_l = 4
    cos_beta = [-0.8459]
    wigner_py, wigner_c, wigner_c_half = wigner(ang_momentum_l, cos_beta)
    array = np.load("tests/wigner/matrix/l=4_cx=-0.8459.npy")
    np.testing.assert_almost_equal(array, wigner_py.ravel(), decimal=8)
    np.testing.assert_almost_equal(array, wigner_c.ravel(), decimal=8)
    np.testing.assert_almost_equal(wigner_c[:, :5, :], wigner_c_half, decimal=8)


def test_wigner_4j_matrix_angle_03():
    ang_momentum_l = 4
    cos_beta = [-0.934, 0.4958]
    beta = np.arccos(cos_beta)
    array = np.load("tests/wigner/matrix/l=4_cx=[-0.934, 0.4958].npy")
    wigner = clib.wigner_d_matrices(ang_momentum_l, beta)
    np.testing.assert_almost_equal(array, wigner, decimal=8)


def test_wigner_4j_matrix_cosine_03():
    ang_momentum_l = 4
    cos_beta = [-0.934, 0.4958]
    wigner_py, wigner_c, wigner_c_half = wigner(ang_momentum_l, cos_beta)
    array = np.load("tests/wigner/matrix/l=4_cx=[-0.934, 0.4958].npy")
    np.testing.assert_almost_equal(array, wigner_py.ravel(), decimal=8)
    np.testing.assert_almost_equal(array, wigner_c.ravel(), decimal=8)
    np.testing.assert_almost_equal(wigner_c[:, :5, :], wigner_c_half, decimal=8)


def test_wigner_2j_dm0_vector():
    R_py = wigner_dm0_vector(2, 0.235)
    R_c = clib.wigner_dm0_vector(2, 0.235)
    R_sympy = wigner_dm0_vector_sympy(2, 0.235)
    np.testing.assert_almost_equal(R_py, R_sympy, decimal=8)
    np.testing.assert_almost_equal(R_c, R_sympy, decimal=8)


def test_wigner_4j_dm0_vector():
    R_py = wigner_dm0_vector(4, 0.235)
    R_c = clib.wigner_dm0_vector(4, 0.235)
    R_sympy = wigner_dm0_vector_sympy(4, 0.235)
    np.testing.assert_almost_equal(R_py, R_sympy, decimal=8)
    np.testing.assert_almost_equal(R_c, R_sympy, decimal=8)


#
# This code was used to generate and store wigner d-matrices from sympy rotation.
# All tests are compared with these stored values to speed up the test runs.
#
# if __name__ == "__main__":
#     value = []
#     l = 2
#     cos_beta = [-0.5498, 0.230]
#     for k, cos_b in enumerate(cos_beta):
#         for i in range(2 * l + 1):
#             for j in range(2 * l + 1):
#                 value.append(
#                     complex(
#                       Rotation.d(l, -l + j, -l + i, np.arccos(cos_b)).doit()
#                     ).real
#                 )
#     value = np.asarray(value)
#     np.save("tests/wigner/l=2_cx=[-0.5498, 0.230]", value)

#     value = []
#     l = 2
#     cos_beta = 0.5
#     for i in range(2 * l + 1):
#         for j in range(2 * l + 1):
#             value.append(
#                 complex(
#                   Rotation.d(l, -l + j, -l + i, np.arccos(cos_beta)).doit()).real
#             )
#     value = np.asarray(value)
#     np.save("tests/wigner/l=2_cx=0.5", value)

#     value = []
#     l = 4
#     cos_beta = -0.8459
#     for i in range(2 * l + 1):
#         for j in range(2 * l + 1):
#             value.append(
#                 complex(
#                   Rotation.d(l, -l + j, -l + i, np.arccos(cos_beta)).doit()).real
#             )
#     value = np.asarray(value)
#     np.save("tests/wigner/l=4_cx=-0.8459", value)

#     value = []
#     l = 4
#     cos_beta = [-0.934, 0.4958]
#     for k, cos_b in enumerate(cos_beta):
#         for i in range(2 * l + 1):
#             for j in range(2 * l + 1):
#                 value.append(
#                     complex(
#                       Rotation.d(l, -l + j, -l + i, np.arccos(cos_b)).doit()).real
#                 )
#     value = np.asarray(value)
#     np.save("tests/wigner/l=4_cx=[-0.934, 0.4958]", value)
