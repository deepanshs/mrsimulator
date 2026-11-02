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


def test_histogram_rounding_at_upper_edge():
    """A sample just below x_max, where (x - x_min) * count / (x_max - x_min) rounds up
    to count, is counted in the last bin."""
    import numpy as np
    from mrsimulator import clib

    x_max, n = 6.5395468350513815, 34
    x = np.array([np.nextafter(x_max, 0)])
    hist = clib.histogram1d(n, 0.0, x_max, x, np.ones(1))[1]
    assert hist.sum() == 1 and hist[-1] == 1

    hist = clib.histogram2d(n, 0.0, x_max, n, 0.0, x_max, x, x, np.ones(1))[2]
    assert hist.sum() == 1 and hist[-1, -1] == 1


def test_histogram2d_interpolation():
    """Bilinear deposit between bin centres: a sample at a bin centre goes to that
    bin, a sample midway between two centres is split equally, and the total weight is
    conserved."""
    import numpy as np
    from mrsimulator import clib

    def hist2d(x, y):
        x, y = np.atleast_1d(x) * 1.0, np.atleast_1d(y) * 1.0
        return clib.histogram2d(4, 0.0, 4.0, 5, 0.0, 5.0, x, y, np.ones(x.size), 1)[2]

    h = hist2d(1.5, 2.5)  # bin centre
    expected = np.zeros((4, 5))
    expected[1, 2] = 1
    assert np.allclose(h, expected)

    h = hist2d(2.0, 2.5)  # midway between x centres 1.5 and 2.5
    expected = np.zeros((4, 5))
    expected[1, 2] = expected[2, 2] = 0.5
    assert np.allclose(h, expected)

    h = hist2d(2.0, 3.0)  # midway in both
    expected = np.zeros((4, 5))
    expected[1:3, 2:4] = 0.25
    assert np.allclose(h, expected)

    rng = np.random.default_rng(0)
    x, y = rng.uniform(0, 4, 1000), rng.uniform(0, 5, 1000)
    h = hist2d(x, y)
    assert np.allclose(h.sum(), 1000)
    assert np.all(h >= 0)


def test_octahedron_binning_vertex_weights():
    """In binning mode, the three octant vertices (x, y, and z axes) have equal weight
    by symmetry."""
    import numpy as np
    import mrsimulator.tests.tests as tst

    for nt in [1, 2, 3, 5, 10]:
        amp = tst.cosine_of_polar_angles_and_amplitudes(nt, interpolation=False)[2]
        assert np.allclose([amp[nt], amp[-1]], amp[0]), nt


def _triangles_per_vertex(nt):
    """Number of mesh triangles sharing each vertex of the octant, enumerated as in
    octahedronInterpolation and get_total_amplitude."""
    import numpy as np

    n_pts = (nt + 1) * (nt + 2) // 2
    count = np.zeros(n_pts, dtype=int)
    i = j = 0
    local_index = nt - 1
    base = nt + 1
    while i < n_pts - 1:
        count[[i, i + 1, base + j]] += 1
        if i < local_index:
            count[[i + 1, base + j, base + j + 1]] += 1
        else:
            local_index = j + nt
            i += 1
        i += 1
        j += 1
    return count


def test_octahedron_binning_weights_conserve_triangle_shares():
    """In binning mode, the weight of each orientation equals its total share over the
    triangles of the interpolation mode: six shares for an interior point, three for an
    edge point, and one for an octant vertex."""
    import numpy as np
    import mrsimulator.tests.tests as tst

    for nt in [1, 2, 3, 4, 6, 12, 24]:
        amp_i = tst.cosine_of_polar_angles_and_amplitudes(nt, interpolation=True)[2]
        amp_b = tst.cosine_of_polar_angles_and_amplitudes(nt, interpolation=False)[2]
        np.testing.assert_allclose(
            amp_b, amp_i * _triangles_per_vertex(nt), rtol=1e-12, err_msg=f"nt={nt}"
        )
