# -*- coding: utf-8 -*-
import mrsimulator.tests.tests as clib
import numpy as np

from tests.python_test_for_c_code.angular_momentum import wigner_rotation


def test__batch_wigner_rotation():
    n = 40
    n_octants = 1
    alpha = np.random.rand(n) * np.pi / 2.0
    beta = np.random.rand(n) * np.pi / 2.0

    cos_beta = np.cos(beta)
    sin_beta = np.sqrt(1.0 - cos_beta ** 2)
    exp_I_beta = cos_beta + 1j * sin_beta
    wigner_2j_matrices = clib.wigner_d_matrices_from_exp_I_beta(2, exp_I_beta).ravel()
    wigner_4j_matrices = clib.wigner_d_matrices_from_exp_I_beta(4, exp_I_beta).ravel()

    R4 = np.asarray(
        [0 - 0.2j, 0, 0 + 0.5j, 0, 0 + 0.1j, 0, 0 - 0.5j, 0, 0 + 0.2j],
        dtype=np.complex128,
    )
    R2 = np.asarray([0 + 0.5j, 0, 0 + 0.1j, 0, 0 - 0.5j], dtype=np.complex128)

    cos_alpha = np.cos(alpha)
    exp_im_alpha = clib.get_exp_Im_alpha(n, cos_alpha, True).ravel()
    exp_im_alpha_in = np.empty(4 * n, dtype=np.complex128)
    exp_im_alpha_in[:] = exp_im_alpha.copy()

    w2, w4 = clib.__batch_wigner_rotation(
        n, n_octants, wigner_2j_matrices, R2, wigner_4j_matrices, R4, exp_im_alpha
    )

    # assert np.allclose(exp_im_alpha_in, exp_im_alpha, atol=1e-15)

    alpha_octants = []
    for i in range(n_octants):
        alpha_octants.append(alpha + i * np.pi / 2.0)
    alpha_octants = np.asarray(alpha_octants).ravel()
    cos_alpha_octants = np.cos(alpha_octants)

    beta_octants = []
    for i in range(n_octants):
        beta_octants.append(beta)
    beta_octants = np.asarray(beta_octants).ravel()
    cos_beta_octants = np.cos(beta_octants)

    w2_1 = clib.__wigner_rotation_2(2, cos_alpha_octants, cos_beta_octants, R2).ravel()
    w4_1 = clib.__wigner_rotation_2(4, cos_alpha_octants, cos_beta_octants, R4).ravel()

    # print(np.argmax(w2 - w2_1))
    assert np.allclose(w2, w2_1, atol=1e-12)
    np.testing.assert_almost_equal(w4, w4_1, decimal=8)


def test_single_2j_rotation_00():
    ang_momentum_l = 2
    R_in = np.asarray([0 + 0.5j, 0, 0 + 0.1j, 0, 0 - 0.5j], dtype=np.complex128)
    indexes = np.arange(5) - 2.0

    for _ in range(10):
        euler_angle = np.random.rand(3) * 2.0 * np.pi

        # single rotation
        R_out = clib.single_wigner_rotation(ang_momentum_l, euler_angle, R_in)

        exp_im_alpha = np.exp(-1j * indexes * euler_angle[0])
        R_out_p = R_in * exp_im_alpha
        wigner = clib.wigner_d_matrices(ang_momentum_l, np.asarray([euler_angle[1]]))
        R_out_p = np.dot(wigner.reshape(5, 5), R_out_p)
        exp_im_gamma = np.exp(-1j * indexes * euler_angle[2])
        R_out_p *= exp_im_gamma

        np.testing.assert_almost_equal(R_out, R_out_p, decimal=8)


def test_single_4j_rotation_00():
    ang_momentum_l = 4
    R_in = np.asarray(
        [0 - 0.2j, 0, 0 + 0.5j, 0, 0 + 0.1j, 0, 0 - 0.5j, 0, 0 + 0.2j],
        dtype=np.complex128,
    )
    indexes = np.arange(9) - 4.0

    for _ in range(10):
        euler_angle = np.random.rand(3) * 2.0 * np.pi

        # single rotation
        R_out = clib.single_wigner_rotation(ang_momentum_l, euler_angle, R_in)

        exp_im_alpha = np.exp(-1j * indexes * euler_angle[0])
        R_out_p = R_in * exp_im_alpha
        wigner = clib.wigner_d_matrices(ang_momentum_l, np.asarray([euler_angle[1]]))
        R_out_p = np.dot(wigner.reshape(9, 9), R_out_p)
        exp_im_gamma = np.exp(-1j * indexes * euler_angle[2])
        R_out_p *= exp_im_gamma

        np.testing.assert_almost_equal(R_out, R_out_p, decimal=8)


def test_wigner_2j_rotation_00():
    cos_beta = 0.5
    cos_alpha = 0.5
    cos_beta = np.asarray([cos_beta], dtype=np.float64).ravel()
    cos_alpha = np.asarray([cos_alpha], dtype=np.float64).ravel()

    R_in = np.asarray([0 + 0.5j, 0, 0 + 0.1j, 0, 0 - 0.5j], dtype=np.complex128)

    # R_out_c = clib.wigner_rotation(2, R_in, cos_alpha, cos_beta)
    R_out_c2 = clib.__wigner_rotation_2(2, cos_alpha, cos_beta, R_in)
    R_out_py = wigner_rotation(2, R_in, cos_alpha, cos_beta)

    # np.testing.assert_almost_equal(R_out_py, R_out_c, decimal=8)
    np.testing.assert_almost_equal(R_out_py, R_out_c2, decimal=8)


def test_wigner_2j_rotation_01():
    cos_beta = [0.15, 0.51, 0.223]
    cos_alpha = [0.53, 0.95, 0.391]
    cos_beta = np.asarray([cos_beta], dtype=np.float64).ravel()
    cos_alpha = np.asarray([cos_alpha], dtype=np.float64).ravel()

    R_in = np.asarray([0 + 0.5j, 0, 0 + 0.1j, 0, 0 - 0.5j], dtype=np.complex128)

    # R_out_c = clib.wigner_rotation(2, R_in, cos_alpha, cos_beta)
    R_out_c2 = clib.__wigner_rotation_2(2, cos_alpha, cos_beta, R_in)
    R_out_py = wigner_rotation(2, R_in, cos_alpha, cos_beta)

    # np.testing.assert_almost_equal(R_out_py, R_out_c, decimal=8)
    np.testing.assert_almost_equal(R_out_py, R_out_c2, decimal=8)


def test_wigner_2j_rotation_02():
    cos_beta = 2.0 * np.random.rand(1024) - 1.0
    cos_alpha = 2.0 * np.random.rand(1024) - 1.0
    cos_beta = np.asarray([cos_beta], dtype=np.float64).ravel()
    cos_alpha = np.asarray([cos_alpha], dtype=np.float64).ravel()

    R_in = np.asarray([0 + 0.5j, 0, 0 + 0.1j, 0, 0 - 0.5j], dtype=np.complex128)

    # R_out_c = clib.wigner_rotation(2, R_in, cos_alpha, cos_beta)
    R_out_c2 = clib.__wigner_rotation_2(2, cos_alpha, cos_beta, R_in)
    R_out_py = wigner_rotation(2, R_in, cos_alpha, cos_beta)

    # np.testing.assert_almost_equal(R_out_py, R_out_c, decimal=8)
    np.testing.assert_almost_equal(R_out_py, R_out_c2, decimal=8)


def test_wigner_4j_rotation_03():
    cos_beta = 0.35
    cos_alpha = 0.598
    cos_beta = np.asarray([cos_beta], dtype=np.float64).ravel()
    cos_alpha = np.asarray([cos_alpha], dtype=np.float64).ravel()

    R_in = np.asarray(
        [0 - 0.2j, 0, 0 + 0.5j, 0, 0 + 0.1j, 0, 0 - 0.5j, 0, 0 + 0.2j],
        dtype=np.complex128,
    )

    # R_out_c = clib.wigner_rotation(4, R_in, cos_alpha, cos_beta)
    R_out_c2 = clib.__wigner_rotation_2(4, cos_alpha, cos_beta, R_in)
    R_out_py = wigner_rotation(4, R_in, cos_alpha, cos_beta)

    # np.testing.assert_almost_equal(R_out_py, R_out_c, decimal=8)
    np.testing.assert_almost_equal(R_out_py, R_out_c2, decimal=8)


def test_wigner_4j_rotation_04():
    cos_beta = 2.0 * np.random.rand(1024) - 1.0
    cos_alpha = 2.0 * np.random.rand(1024) - 1.0
    cos_beta = np.asarray([cos_beta], dtype=np.float64).ravel()
    cos_alpha = np.asarray([cos_alpha], dtype=np.float64).ravel()

    R_in = np.asarray(
        [0 - 0.2j, 0, 0 + 0.5j, 0, 0 + 0.1j, 0, 0 - 0.5j, 0, 0 + 0.2j],
        dtype=np.complex128,
    )

    # R_out_c = clib.wigner_rotation(4, R_in, cos_alpha, cos_beta)
    R_out_c2 = clib.__wigner_rotation_2(4, cos_alpha, cos_beta, R_in)
    R_out_py = wigner_rotation(4, R_in, cos_alpha, cos_beta)

    # np.testing.assert_almost_equal(R_out_py, R_out_c, decimal=8)
    np.testing.assert_almost_equal(R_out_py, R_out_c2, decimal=8)
