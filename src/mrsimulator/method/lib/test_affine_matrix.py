import numpy as np
from mrsimulator.method.lib import FiveQ_VAS
from mrsimulator.method.lib import SevenQ_VAS
from mrsimulator.method.lib import ST1_VAS
from mrsimulator.method.lib import ST2_VAS
from mrsimulator.method.lib import ThreeQ_VAS


def test_mqmas_affine_matrix():
    assert np.allclose(
        ThreeQ_VAS(channels=["23Na"]).affine_matrix, [9 / 16, 7 / 16, 0, 1]
    )
    assert np.allclose(
        ThreeQ_VAS(channels=["27Al"]).affine_matrix, [12 / 31, 19 / 31, 0, 1]
    )
    assert np.allclose(
        ThreeQ_VAS(channels=["51V"]).affine_matrix, [45 / 146, 101 / 146, 0, 1]
    )
    assert np.allclose(
        ThreeQ_VAS(channels=["93Nb"]).affine_matrix, [36 / 127, 91 / 127, 0, 1]
    )

    assert np.allclose(
        FiveQ_VAS(channels=["27Al"]).affine_matrix, [12 / 37, 25 / 37, 0, 1]
    )
    assert np.allclose(
        FiveQ_VAS(channels=["51V"]).affine_matrix, [9 / 20, 11 / 20, 0, 1]
    )
    assert np.allclose(
        FiveQ_VAS(channels=["93Nb"]).affine_matrix, [36 / 131, 95 / 131, 0, 1]
    )

    assert np.allclose(
        SevenQ_VAS(channels=["51V"]).affine_matrix, [45 / 206, 161 / 206, 0, 1]
    )
    assert np.allclose(
        SevenQ_VAS(channels=["93Nb"]).affine_matrix, [18 / 25, 7 / 25, 0, 1]
    )


def test_stvas_affine_matrix():
    assert np.allclose(ST1_VAS(channels=["23Na"]).affine_matrix, [9 / 17, 8 / 17, 0, 1])
    assert np.allclose(
        ST1_VAS(channels=["27Al"]).affine_matrix, [24 / 31, 7 / 31, 0, 1]
    )
    assert np.allclose(
        ST1_VAS(channels=["51V"]).affine_matrix, [45 / 73, 28 / 73, 0, 1]
    )
    assert np.allclose(
        ST1_VAS(channels=["93Nb"]).affine_matrix, [72 / 127, 55 / 127, 0, 1]
    )

    assert np.allclose(
        ST2_VAS(channels=["27Al"]).affine_matrix, [6 / 17, 11 / 17, 0, 1]
    )
    assert np.allclose(
        ST2_VAS(channels=["51V"]).affine_matrix, [45 / 68, 23 / 68, 0, 1]
    )
    assert np.allclose(
        ST2_VAS(channels=["93Nb"]).affine_matrix, [18 / 19, 1 / 19, 0, 1]
    )
