# -*- coding: utf-8 -*-
import numpy as np

__author__ = "Deepansh J. Srivastava"
__email__ = "deepansh2012@gmail.com"


def wigner_rotation(
    l, R_in, cos_alpha=None, cos_beta=None, wigner_matrix=None, phase_alpha=None
):
    n = 2 * l + 1
    if wigner_matrix is None:
        n_orientation = cos_beta.size
        wigner = wigner_d_matrix_cosines(l, cos_beta)
    else:
        wigner = wigner_matrix
        n_orientation = wigner.shape[0]

    pha = cos_alpha - 1j * np.sqrt(1.0 - cos_alpha ** 2)
    ph2 = np.copy(pha)

    R_vec = np.tile(R_in, n_orientation).reshape(n_orientation, n)

    for m in range(1, l + 1):
        R_vec[:, l + m] *= ph2
        R_vec[:, l - m] *= ph2.conj()
        ph2 *= pha

    R_out = wigner.ravel() * np.tile(R_vec, (1, n)).ravel()
    return R_out.reshape(n_orientation, n, n).sum(axis=-1)


def wigner_d_matrix(l, beta):
    """
    Evaluates a wigner matrix of rank `l` ev
    """
    return wigner_d_matrix_cosines(l, np.cos(beta))


def wigner_d_matrix_cosines(l, cos_beta):
    r"""
    Returns a $(2l+1) \times (2l+1)$ wigner-d(cos_beta) matrix for rank $l$ at
    a given `cos_beta`. Currently only rank l=2 and l=4 is supported.

    If `cos_beta` is a 1D-numpy array of size n, a
    `n x (2l+1) x (2l+1)` matrix is returned instead.

    :ivar l: The angular momentum quantum number.
    :ivar cos_beta: An 1D numpy array or a scalar representing the cosine of
                    $\beta$ angles.
    """
    if not isinstance(cos_beta, np.ndarray):
        cos_beta = np.asarray([cos_beta]).ravel()

    cx = cos_beta
    cx2 = cx * cx
    sx = np.sqrt(1.0 - cx2)

    if l == 2:
        wigner = np.empty((25, cx.size), dtype=np.float64)

        t1 = 1.0 + cx
        temp = -sx * t1 / 2.0
        wigner[19] = temp  # 2,  1 # 19
        wigner[5] = -temp  # -2, -1 #  5
        wigner[23] = -temp  # 1,  2 # 23
        wigner[1] = temp  # -1, -2 #  1

        temp = t1 * t1 / 4.0
        wigner[24] = temp  # 2,  2  # 24
        wigner[0] = temp  # -2, -2 #  0

        t1 = 1.0 - cx
        temp = -sx * t1 / 2.0
        wigner[9] = temp  # 2, -1 #  9
        wigner[15] = -temp  # -2,  1 # 15
        wigner[3] = temp  # 1, -2 #  3
        wigner[21] = -temp  # -1,  2 # 21

        temp = t1 * t1 / 4.0
        wigner[4] = temp  # 2, -2 #  4
        wigner[20] = temp  # -2,  2 # 20

        temp = 0.6123724355 * sx * sx
        wigner[14] = temp  # 2,  0 # 14
        wigner[10] = temp  # -2,  0 # 10
        wigner[22] = temp  # 0,  2 # 22
        wigner[2] = temp  # 0, -2 #  2

        temp = 1.224744871 * sx * cx
        wigner[13] = -temp  # 1,  0 # 13
        wigner[17] = temp  # 0,  1 # 17
        wigner[7] = -temp  # 0, -1 #  7
        wigner[11] = temp  # -1,  0 # 11

        temp = (2.0 * cx2 + cx - 1.0) / 2.0
        wigner[18] = temp  # 1,  1 # 18
        wigner[6] = temp  # -1, -1 #  6

        temp = -(2.0 * cx2 - cx - 1.0) / 2.0
        wigner[8] = temp  # 1, -1 #  8
        wigner[16] = temp  # -1,  1 # 16

        wigner[12] = 1.5 * cx2 - 0.5  # 0,  0 # 12

        wigner = wigner.reshape(5, 5, cx.size)
        return np.moveaxis(wigner, -1, 0)

    if l == 4:
        wigner = np.empty((81, cx.size), dtype=np.float64)

        sx2 = sx * sx
        sx3 = sx2 * sx

        cxp1 = 1.0 + cx
        cxm1 = 1.0 - cx
        cxp12 = cxp1 * cxp1
        cxm12 = cxm1 * cxm1
        cxm13 = cxm12 * cxm1
        cxp13 = cxp12 * cxp1

        temp = 0.0625 * cxp12 * cxp12
        wigner[0] = temp  # -4, -4 #  0
        wigner[80] = temp  # 4,  4 # 80

        temp = 0.0625 * cxm12 * cxm12
        wigner[72] = temp  # -4,  4 # 72
        wigner[8] = temp  # 4, -4 #  8

        temp = -0.1767766953 * cxp13 * sx
        wigner[1] = temp  # -3, -4 #  1
        wigner[9] = -temp  # -4, -3 #  9
        wigner[79] = -temp  # 3,  4 # 79
        wigner[71] = temp  # 4,  3 #  9

        temp = -0.1767766953 * cxm13 * sx
        wigner[7] = temp  # 3, -4 #  7
        wigner[63] = -temp  # -4,  3 # 63
        wigner[73] = -temp  # -3,  4 # 73
        wigner[17] = temp  # 4, -3 # 17

        temp = -0.4677071733 * cxp1 * sx3
        wigner[53] = temp  # 4,  1 # 53
        wigner[27] = -temp  # -4, -1 # 27
        wigner[77] = -temp  # 1,  4 # 77
        wigner[3] = temp  # -1, -4 #  3

        temp = -0.4677071733 * cxm1 * sx3
        wigner[35] = temp  # 4, -1 # 35
        wigner[45] = -temp  # -4,  1 # 45
        wigner[75] = -temp  # -1,  4 # 75
        wigner[5] = temp  # 1, -4 #  5

        temp = 0.5229125166 * sx3 * sx
        wigner[44] = temp  # 4,  0 # 44
        wigner[36] = temp  # -4,  0 # 36
        wigner[76] = temp  # 0,  4 # 76
        wigner[4] = temp  # 0, -4 #  4

        temp = -1.4790199458 * sx3 * cx
        wigner[43] = temp  # 3,  0 # 43
        wigner[37] = -temp  # -3,  0 # 37
        wigner[67] = -temp  # 0,  3 # 67
        wigner[13] = temp  # 0, -3 # 13

        temp = 0.3307189139 * sx2 * cxp12
        wigner[78] = temp  # 2,  4 # 78
        wigner[2] = temp  # -2, -4 #  2
        wigner[62] = temp  # 4,  2 # 62
        wigner[18] = temp  # -4, -2 # 18

        temp = 0.3307189139 * sx2 * cxm12
        wigner[6] = temp  # 2, -4 #  6
        wigner[74] = temp  # -2,  4 # 74
        wigner[54] = temp  # -4,  2 # 54
        wigner[26] = temp  # 4, -2 # 26

        temp = 0.4677071733 * cxp12 * sx * (2.0 * cx - 1.0)
        wigner[69] = temp  # 2,  3 # 69
        wigner[11] = -temp  # -2, -3 # 11
        wigner[61] = -temp  # 3,  2 # 61
        wigner[19] = temp  # -3, -2 # 19

        temp = 0.4677071733 * cxm12 * sx * (-2.0 * cx - 1.0)
        wigner[15] = temp  # 2, -3 # 15
        wigner[65] = -temp  # -2,  3 # 65
        wigner[55] = -temp  # -3,  2 # 55
        wigner[25] = temp  # 3, -2 # 25

        temp = 0.25 * cxp12 * (1.0 - 7.0 * cxm1 + 7.0 * cxm12)
        wigner[60] = temp  # 2,  2 # 60
        wigner[20] = temp  # -2, -2 # 20

        temp = 0.25 * cxm12 * (1.0 - 7.0 * cxp1 + 7.0 * cxp12)
        wigner[56] = temp  # -2,  2 # 56
        wigner[24] = temp  # 2, -2 # 24

        temp = 0.3952847075 * sx2 * (7.0 * cx2 - 1)
        wigner[42] = temp  # 2,  0 # 42
        wigner[38] = temp  # -2,  0 # 38
        wigner[58] = temp  # 0,  2 # 58
        wigner[22] = temp  # 0, -2 # 22

        temp = 0.125 * cxp13 * (-3.0 + 4.0 * cx)
        wigner[10] = temp  # -3, -3 # 10
        wigner[70] = temp  # 3,  3 # 70

        temp = 0.125 * cxm13 * (3.0 + 4.0 * cx)
        wigner[64] = temp  # -3,  3 # 64
        wigner[16] = temp  # 3, -3 # 16

        temp = 0.3307189139 * cxm1 * cxp12 * (-1.0 + 4.0 * cx)
        wigner[12] = temp  # -1, -3 # 12
        wigner[28] = temp  # -3, -1 # 28
        wigner[68] = temp  # 1,  3 # 68
        wigner[52] = temp  # 3,  1 # 52

        temp = 0.3307189139 * cxm12 * cxp1 * (1.0 + 4.0 * cx)
        wigner[14] = temp  # 1, -3 # 14
        wigner[46] = temp  # -3,  1 # 46
        wigner[66] = temp  # -1,  3 # 66
        wigner[34] = temp  # 3, -1 # 34

        temp = -0.5590169944 * (4.0 - 18.0 * cxm1 + 21.0 * cxm12 - 7.0 * cxm13) * sx
        wigner[41] = temp  # 1,  0 # 41
        wigner[39] = -temp  # -1,  0 # 39
        wigner[49] = -temp  # 0,  1 # 49
        wigner[31] = temp  # 0, -1 # 31

        temp = -0.3535533906 * (3.0 - 10.5 * cxm1 + 7.0 * cxm12) * sx * cxp1
        wigner[51] = temp  # 2,  1 # 51
        wigner[29] = -temp  # -2, -1 # 29
        wigner[59] = -temp  # 1,  2 # 59
        wigner[21] = temp  # -1, -2 # 21

        temp = -0.3535533906 * (10.0 - 17.5 * cxm1 + 7.0 * cxm12) * sx * cxm1
        wigner[23] = temp  # 1, -2 # 23
        wigner[57] = -temp  # -1,  2 # 57
        wigner[47] = -temp  # -2,  1 # 47
        wigner[33] = temp  # 2, -1 # 33

        temp = 0.5 * (1.0 - 9.0 * cxm1 + 15.75 * cxm12 - 7.0 * cxm13) * cxp1
        wigner[30] = temp  # -1, -1 # 30
        wigner[50] = temp  # 1,  1 # 50

        temp = 0.5 * (10.0 - 30.0 * cxm1 + 26.25 * cxm12 - 7.0 * cxm13) * cxm1
        wigner[32] = temp  # 1, -1 # 32
        wigner[48] = temp  # -1,  1 # 48

        temp = 0.125 * (3.0 - 30.0 * cx2 + 35.0 * cx2 * cx2)
        wigner[40] = temp  # 0,  0 # 40

        # factor = (np.arange(9) - 4.)*(cos_alpha -
        #                               1j * np.sqrt(1.0 - cos_alpha ** 2))
        wigner = wigner.reshape(9, 9, cx.size)
        # wigner *= factor[np.newaxis, :, np.newaxis]
        return np.moveaxis(wigner, -1, 0)


def wigner_dm0_vector(l: int, angle):
    R_out = np.empty(2 * l + 1, dtype=np.float64)
    sx = np.sin(angle)
    cx = np.cos(angle)
    if l == 2:
        R_out[0] = 0.6123724355 * sx * sx
        R_out[1] = 1.224744871 * sx * cx
        R_out[2] = 1.5 * cx * cx - 0.5
        R_out[3] = -R_out[1] * 1
        R_out[4] = R_out[0] * 1
        return R_out

    if l == 4:
        sx2 = sx * sx
        sx3 = sx2 * sx
        cx2 = 1.0 - sx2
        cxm1 = 1.0 - cx
        cxm12 = cxm1 * cxm1
        temp = 4.0 - 18.0 * cxm1 + 21.0 * cxm12 - 7.0 * cxm12 * cxm1
        R_out[0] = 0.5229125166 * sx3 * sx
        R_out[1] = 1.4790199458 * sx3 * cx
        R_out[2] = 0.3952847075 * sx2 * (7.0 * cx2 - 1)
        R_out[3] = 0.5590169944 * temp * sx
        R_out[4] = 0.125 * (3.0 - 30.0 * cx2 + 35 * cx2 * cx2)
        R_out[5] = -R_out[3]
        R_out[6] = R_out[2]
        R_out[7] = -R_out[1]
        R_out[8] = R_out[0]
        return R_out
