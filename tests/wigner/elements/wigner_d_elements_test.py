# -*- coding: utf-8 -*-
import mrsimulator.tests.tests as clib
import numpy as np


def prep_assertion(ang_momentum_l, filename):
    beta = np.asarray([0, 30, 45, 60, 90]) * np.pi / 180.0
    n1 = int(2 * ang_momentum_l + 1)
    s1 = np.load(filename)
    s2 = [
        clib.wigner_d_element(
            ang_momentum_l, -ang_momentum_l + i, -ang_momentum_l + j, b
        )
        for b in beta
        for i in range(n1)
        for j in range(n1)
    ]
    np.testing.assert_almost_equal(s1, s2, decimal=8)


def test_wigner_l_half_elements():
    prep_assertion(1 / 2, "tests/wigner/elements/2l=1_beta=[0,30,45,60,90].npy")


def test_wigner_l_one_elements():
    prep_assertion(1, "tests/wigner/elements/2l=2_beta=[0,30,45,60,90].npy")


def test_wigner_l_three_half_elements():
    prep_assertion(3 / 2, "tests/wigner/elements/2l=3_beta=[0,30,45,60,90].npy")


def test_wigner_l_two_elements():
    prep_assertion(2, "tests/wigner/elements/2l=4_beta=[0,30,45,60,90].npy")


#
# This code was used to generate and store wigner d-elements from sympy rotation.
# All tests are compared with these stored values to speed up the test runs.
#
# if __name__ == "__main__":
#     from sympy import Integer
#     from sympy.physics.quantum.spin import Rotation

#     def prep_setup(l, filename):
#         beta = np.asarray([0, 30, 45, 60, 90]) * np.pi / 180.0
#         n1 = int(2 * l + 1)
#         s = [
#             complex(Rotation.d(l, -l + i, -l + j, b).doit()).real
#             for b in beta
#             for i in range(n1)
#             for j in range(n1)
#         ]
#         np.save(filename, s)

#     prep_setup(1 / Integer(2), "tests/wigner/elements/2l=1_beta=[0,30,45,60,90].npy")
#     prep_setup(1, "tests/wigner/elements/2l=2_beta=[0,30,45,60,90].npy")
#     prep_setup(3 / Integer(2), "tests/wigner/elements/2l=3_beta=[0,30,45,60,90].npy")
#     prep_setup(2, "tests/wigner/elements/2l=4_beta=[0,30,45,60,90].npy")
