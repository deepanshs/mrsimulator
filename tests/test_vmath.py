import mrsimulator.tests.tests as clib
import numpy as np


def test_c_arithmetic():
    for _ in range(10):
        a = np.random.rand(500)
        b = np.random.rand(500)

        c = [clib.vm_absd(_a) for _a in a]
        np.testing.assert_allclose(c, np.abs(a))

        c = [clib.vm_absd(_b) for _b in b]
        np.testing.assert_allclose(c, np.abs(b))

        c = clib.vm_add(a, b)
        np.testing.assert_allclose(c, a + b)

        c = clib.vm_sub(a, b)
        np.testing.assert_allclose(c, a - b)

        c = clib.vm_mult(a, b)
        np.testing.assert_allclose(c, a * b)

        zero_index = np.where(b == 0)[0]
        b[zero_index] = 1
        c = clib.vm_div(a, b)
        np.testing.assert_allclose(c, a / b)

        b_copy = b.copy()
        clib.vm_add_inplace(a, b_copy)
        np.testing.assert_allclose(b_copy, a + b)

        b_copy = b.copy()
        clib.vm_sub_inplace(a, b_copy)
        np.testing.assert_allclose(b_copy, b - a)

        b_copy = b.copy()
        clib.vm_mult_inplace(a, b_copy)
        np.testing.assert_allclose(b_copy, a * b)

        b_copy = b.copy()
        zero_index = np.where(a == 0)[0]
        a[zero_index] = 1
        c = clib.vm_div_inplace(a, b_copy)
        np.testing.assert_allclose(b_copy, b / a)


def test_complex():
    for _ in range(10):
        a = np.random.rand(500) + 1j * np.random.rand(500)
        b = np.random.rand(500) + 1j * np.random.rand(500)
        c = clib.vm_cmult(a, b)
        np.testing.assert_allclose(c, a * b)

        c = clib.vm_cmult_conj(a, b)
        np.testing.assert_allclose(c, a * b.conj())

        c = clib.vm_I_exp(a)
        np.testing.assert_almost_equal(c, np.exp(a), decimal=3)


def test_c_power():
    for _ in range(10):
        a = np.abs(np.random.rand(500))

        c = clib.vm_sq(a)
        np.testing.assert_allclose(c, a**2)

        c = clib.vm_sqrt(a)
        np.testing.assert_allclose(c, np.sqrt(a))

        a_copy = a.copy()
        clib.vm_sq_inplace(a_copy)
        np.testing.assert_allclose(a_copy, a**2)

        a_copy = a.copy()
        clib.vm_sqrt_inplace(a_copy)
        np.testing.assert_allclose(a_copy, np.sqrt(a))


def test_c_trig():
    a = ((np.arange(12000) / 2000) - 3) * np.pi * 2
    print(a.min(), a.max())

    c = clib.vm_cos(a)
    np.testing.assert_almost_equal(c, np.cos(a), decimal=9)

    c = clib.vm_sin(a)
    np.testing.assert_almost_equal(c, np.sin(a), decimal=9)

    c = clib.vm_cos_I_sin(a)
    np.testing.assert_almost_equal(c, np.cos(a) + 1j * np.sin(a), decimal=9)
