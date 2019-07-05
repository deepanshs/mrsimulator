import mrsimulator.sandbox as clib
from mrsimulator.python.angular_momentum import (
    wigner_d_matrix_cosines,
    wigner_dm0_vector,
)
import numpy as np
from sympy.physics.quantum.spin import Rotation


def wigner_dm0_vector_sympy(l, angle):
    R_out = np.empty(2 * l + 1, dtype=np.float64)
    for i in range(2 * l + 1):
        R_out[i] = complex(Rotation.d(l, -l + i, 0, angle).doit()).real
    return R_out


def wigner(l, cos_beta):
    # python test
    cos_beta = np.asarray([cos_beta]).ravel()
    wigner_py = wigner_d_matrix_cosines(l, cos_beta)

    # c test
    wigner_c = clib.wigner_d_matrix_cosines(l, cos_beta)
    return wigner_py.ravel(), wigner_c.ravel()


# All wigner matrix are tested against the wigner matrix computed using
# Sympy Rotaion methods. See the main function in this file.


def test_wigner_2j_matrix_cosine_00():
    l = 2
    cos_beta = 0.5
    wigner_py, wigner_c = wigner(l, cos_beta)
    array = np.load("tests/wigner/l=2_cx=0.5.npy")
    assert np.allclose(array, wigner_py, atol=1e-15)
    assert np.allclose(array, wigner_c, atol=1e-15)


def test_wigner_2j_matrix_cosine_01():
    l = 2
    cos_beta = [-0.5498, 0.230]
    wigner_py, wigner_c = wigner(l, cos_beta)
    array = np.load("tests/wigner/l=2_cx=[-0.5498, 0.230].npy")
    assert np.allclose(array, wigner_py, atol=1e-15)
    assert np.allclose(array, wigner_c, atol=1e-15)


def test_wigner_4j_matrix_cosine_02():
    l = 4
    cos_beta = [-0.8459]
    wigner_py, wigner_c = wigner(l, cos_beta)
    array = np.load("tests/wigner/l=4_cx=-0.8459.npy")
    assert np.allclose(array, wigner_py, atol=1e-15)
    assert np.allclose(array, wigner_c, atol=1e-15)


def test_wigner_4j_matrix_cosine_03():
    l = 4
    cos_beta = [-0.934, 0.4958]
    wigner_py, wigner_c = wigner(l, cos_beta)
    array = np.load("tests/wigner/l=4_cx=[-0.934, 0.4958].npy")
    assert np.allclose(array, wigner_py, atol=1e-15)
    assert np.allclose(array, wigner_c, atol=1e-15)


def test_wigner_2j_dm0_vector():
    R_py = wigner_dm0_vector(2, 0.235)
    R_c = clib.wigner_dm0_vector(2, 0.235)
    R_sympy = wigner_dm0_vector_sympy(2, 0.235)
    assert np.allclose(R_py, R_sympy, atol=1e-15)
    assert np.allclose(R_c, R_sympy, atol=1e-15)


def test_wigner_4j_dm0_vector():
    R_py = wigner_dm0_vector(4, 0.235)
    R_c = clib.wigner_dm0_vector(4, 0.235)
    R_sympy = wigner_dm0_vector_sympy(4, 0.235)
    assert np.allclose(R_py, R_sympy, atol=1e-15)
    assert np.allclose(R_c, R_sympy, atol=1e-15)


# def triangle(f, n_points):
#     f = np.sort(f)
#     h = 2.0 / (f[2] - f[0])
#     x = f
#     y = [0, h, 0]

#     p = np.arange(int(x[1] - x[0]))
#     spec = np.arange(n_points)
#     y[int(x[0]):int(x[0]) + p] = y[0] + (y[1] - y[0]) / (x[1] - x[0]) * (p - x[0])


if __name__ == "__main__":
    value = []
    l = 2
    cos_beta = [-0.5498, 0.230]
    for k, cos_b in enumerate(cos_beta):
        for i in range(2 * l + 1):
            for j in range(2 * l + 1):
                value.append(
                    complex(Rotation.d(l, -l + j, -l + i, np.arccos(cos_b)).doit()).real
                )
    value = np.asarray(value)
    np.save("tests/wigner/l=2_cx=[-0.5498, 0.230]", value)

    value = []
    l = 2
    cos_beta = 0.5
    for i in range(2 * l + 1):
        for j in range(2 * l + 1):
            value.append(
                complex(Rotation.d(l, -l + j, -l + i, np.arccos(cos_beta)).doit()).real
            )
    value = np.asarray(value)
    np.save("tests/wigner/l=2_cx=0.5", value)

    value = []
    l = 4
    cos_beta = -0.8459
    for i in range(2 * l + 1):
        for j in range(2 * l + 1):
            value.append(
                complex(Rotation.d(l, -l + j, -l + i, np.arccos(cos_beta)).doit()).real
            )
    value = np.asarray(value)
    np.save("tests/wigner/l=4_cx=-0.8459", value)

    value = []
    l = 4
    cos_beta = [-0.934, 0.4958]
    for k, cos_b in enumerate(cos_beta):
        for i in range(2 * l + 1):
            for j in range(2 * l + 1):
                value.append(
                    complex(Rotation.d(l, -l + j, -l + i, np.arccos(cos_b)).doit()).real
                )
    value = np.asarray(value)
    np.save("tests/wigner/l=4_cx=[-0.934, 0.4958]", value)
