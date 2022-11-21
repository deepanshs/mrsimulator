import numpy as np
from mrsimulator.tests.tests import rank_2_tensor_products
from mrsimulator.tests.tests import single_wigner_rotation


def test_rank_2_tensor_products():
    pi2 = np.pi / 2
    euler_angles = np.array(
        [[0, 0, 0], [0, pi2, 0], [1.81, pi2, 0.321], [pi2, pi2, 0], [pi2, pi2, pi2]]
    )

    for angles in euler_angles:
        for zeta_q, eta_q in zip([10, 1.5], [0.1, 0.8]):
            R2_q = [
                -eta_q * zeta_q / 2,
                0,
                np.sqrt(3 / 2) * zeta_q,
                0,
                -eta_q * zeta_q / 2,
            ]
            R2_q = np.array(R2_q, dtype=np.complex128)
            R2_q = single_wigner_rotation(l=2, euler_angles=angles, R_in=R2_q)
            R0, R2, R4 = rank_2_tensor_products(R2_q, R2_q)

            error = "tensor product error in"
            zq2 = zeta_q**2
            eq2 = eta_q**2
            q_0 = (9 * zq2) * (0.333333333 * eq2 + 1) / (6 * np.sqrt(5))
            np.testing.assert_almost_equal(
                complex(q_0), R0[0], decimal=5, err_msg=f"{error} R0"
            )

            q_2 = np.asarray(
                [
                    -zq2 * eta_q * 3 / np.sqrt(21),
                    0,
                    3 * zq2 * (0.333333333 * eq2 - 1) / (2 * np.sqrt(3.5)),
                    0,
                    -zq2 * eta_q * 3 / np.sqrt(21),
                ],
                dtype=complex,
            )
            q_2 = single_wigner_rotation(l=2, euler_angles=angles, R_in=q_2)
            np.testing.assert_almost_equal(q_2, R2, decimal=5, err_msg=f"{error} R2")

            q_4 = np.asarray(
                [
                    zq2 * eq2 / 4,
                    0,
                    -zq2 * eta_q * 3 / (2 * np.sqrt(7)),
                    0,
                    9 * zq2 * ((eq2 / 18) + 1) / np.sqrt(70),
                    0,
                    -zq2 * eta_q * 3 / (2 * np.sqrt(7)),
                    0,
                    zq2 * eq2 / 4,
                ],
                dtype=complex,
            )
            q_4 = single_wigner_rotation(l=4, euler_angles=angles, R_in=q_4)
            np.testing.assert_almost_equal(q_4, R4, decimal=5, err_msg=f"{error} R4")
