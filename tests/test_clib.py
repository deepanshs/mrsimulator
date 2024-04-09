import numpy as np
from mrsimulator.clib import histogram1d
from mrsimulator.clib import histogram2d


def test_hist1s():
    vec = np.random.normal(loc=0.5, scale=0.2, size=100_000)
    hist_np, x_np = np.histogram(vec, bins=10, range=[-0.1, 1.1])

    x_np = (x_np[1:] + x_np[:-1]) / 2

    x_c, hist_c = histogram1d(
        sample_x=vec, x_count=10, x_min=-0.1, x_max=1.1, weights=np.ones(vec.size)
    )

    assert np.allclose(x_c, x_np)
    assert np.allclose(hist_c, hist_np)


def test_hist2d():
    s1 = 0.2
    s2 = 0.4
    r = 0.5
    cov = [
        [s1**2, r * s1 * s2],
        [r * s1 * s2, s2**2],
    ]
    vec_2s = np.random.multivariate_normal(mean=[0.3, 0.5], cov=cov, size=100_000)
    vec_x = vec_2s[:, 0]
    vec_y = vec_2s[:, 1]

    hist_np, x_np, y_np = np.histogram2d(
        vec_x, vec_y, bins=[50, 50], range=[[-0.2, 2.1], [-0.2, 2.1]]
    )
    x_np = (x_np[1:] + x_np[:-1]) / 2
    y_np = (y_np[1:] + y_np[:-1]) / 2

    x_c, y_c, hist_c = histogram2d(
        sample_x=vec_x,
        sample_y=vec_y,
        x_count=50,
        y_count=50,
        x_min=-0.2,
        x_max=2.1,
        y_min=-0.2,
        y_max=2.1,
        weights=np.ones(vec_x.size),
    )

    assert np.allclose(x_c, x_np)
    assert np.allclose(y_c, y_np)
    assert np.allclose(hist_c, hist_np)
