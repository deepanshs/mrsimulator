import mrsimulator.sandbox as clib
from mrsimulator.python.angular_momentum import wigner_rotation
import numpy as np


def test_wigner_2j_rotation_00():
    cos_beta = 0.5
    cos_alpha = 0.5
    cos_beta = np.asarray([cos_beta], dtype=np.float64).ravel()
    cos_alpha = np.asarray([cos_alpha], dtype=np.float64).ravel()

    R_in = np.asarray([0 + 0.5j, 0, 0 + 0.1j, 0, 0 - 0.5j], dtype=np.complex128)

    R_out_c = clib.wigner_rotation(2, R_in, cos_alpha, cos_beta)
    R_out_py = wigner_rotation(2, R_in, cos_alpha, cos_beta)

    assert np.allclose(R_out_py, R_out_c, atol=1e-15)


def test_wigner_2j_rotation_01():
    cos_beta = [0.15, 0.51, 0.223]
    cos_alpha = [0.53, 0.95, 0.391]
    cos_beta = np.asarray([cos_beta], dtype=np.float64).ravel()
    cos_alpha = np.asarray([cos_alpha], dtype=np.float64).ravel()

    R_in = np.asarray([0 + 0.5j, 0, 0 + 0.1j, 0, 0 - 0.5j], dtype=np.complex128)

    R_out_c = clib.wigner_rotation(2, R_in, cos_alpha, cos_beta)
    R_out_py = wigner_rotation(2, R_in, cos_alpha, cos_beta)

    assert np.allclose(R_out_py, R_out_c, atol=1e-15)


def test_wigner_2j_rotation_02():
    cos_beta = 2.0 * np.random.rand(1024) - 1.0
    cos_alpha = 2.0 * np.random.rand(1024) - 1.0
    cos_beta = np.asarray([cos_beta], dtype=np.float64).ravel()
    cos_alpha = np.asarray([cos_alpha], dtype=np.float64).ravel()

    R_in = np.asarray([0 + 0.5j, 0, 0 + 0.1j, 0, 0 - 0.5j], dtype=np.complex128)

    R_out_c = clib.wigner_rotation(2, R_in, cos_alpha, cos_beta)
    R_out_py = wigner_rotation(2, R_in, cos_alpha, cos_beta)

    assert np.allclose(R_out_py, R_out_c, atol=1e-15)


def test_wigner_4j_rotation_03():
    cos_beta = 0.35
    cos_alpha = 0.598
    cos_beta = np.asarray([cos_beta], dtype=np.float64).ravel()
    cos_alpha = np.asarray([cos_alpha], dtype=np.float64).ravel()

    R_in = np.asarray(
        [0 - 0.2j, 0, 0 + 0.5j, 0, 0 + 0.1j, 0, 0 - 0.5j, 0, 0 + 0.2j],
        dtype=np.complex128,
    )

    R_out_c = clib.wigner_rotation(4, R_in, cos_alpha, cos_beta)
    R_out_py = wigner_rotation(4, R_in, cos_alpha, cos_beta)

    assert np.allclose(R_out_py, R_out_c, atol=1e-15)


def test_wigner_4j_rotation_04():
    cos_beta = 2.0 * np.random.rand(1024) - 1.0
    cos_alpha = 2.0 * np.random.rand(1024) - 1.0
    cos_beta = np.asarray([cos_beta], dtype=np.float64).ravel()
    cos_alpha = np.asarray([cos_alpha], dtype=np.float64).ravel()

    R_in = np.asarray(
        [0 - 0.2j, 0, 0 + 0.5j, 0, 0 + 0.1j, 0, 0 - 0.5j, 0, 0 + 0.2j],
        dtype=np.complex128,
    )

    R_out_c = clib.wigner_rotation(4, R_in, cos_alpha, cos_beta)
    R_out_py = wigner_rotation(4, R_in, cos_alpha, cos_beta)

    assert np.allclose(R_out_py, R_out_c, atol=1e-15)
