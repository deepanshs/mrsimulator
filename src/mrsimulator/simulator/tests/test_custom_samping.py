"""Test for the custum sampling class."""
from pathlib import Path

import numpy as np
from mrsimulator.simulator.sampling_scheme import generate_custom_sampling
from mrsimulator.simulator.sampling_scheme import zcw_averaging


def test_generate_custom_sampling():
    alpha = np.random.rand(10)
    beta = np.random.rand(10)
    weight = np.random.rand(10)

    res = generate_custom_sampling(alpha, beta, weight)

    np.testing.assert_allclose(res.alpha, alpha)
    np.testing.assert_allclose(res.beta, beta)
    np.testing.assert_allclose(res.weight, weight)


def test_custom_sampling():
    sampling = zcw_averaging(M=4, integration_volume="hemisphere", triangle_mesh=False)

    size = sampling.alpha.size
    sampling.save(f"zcw_h_{size}_deg.txt", units="deg")

    p = Path(f"zcw_h_{size}_deg.txt")
    assert p.is_file()
    p.unlink(missing_ok=False)

    sampling.save(f"zcw_h_{size}_rad.txt", units="rad")

    p = Path(f"zcw_h_{size}_rad.txt")
    assert p.is_file()
    p.unlink(missing_ok=False)
