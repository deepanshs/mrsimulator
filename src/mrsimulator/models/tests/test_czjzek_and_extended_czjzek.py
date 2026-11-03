from os import path

import numpy as np
from mrsimulator.models import CzjzekDistribution
from mrsimulator.models import ExtCzjzekDistribution
from mrsimulator.models.utils import x_y_to_zeta_eta
from mrsimulator.models.utils import zeta_eta_to_x_y

__author__ = "Deepansh J. Srivastava"
__email__ = "srivastava.89@osu.edu"

MODULE_DIR = path.dirname(path.abspath(__file__))
COUNT = int(1e6)


def test_extended_czjzek_eta_distribution_1():
    filename = path.join(MODULE_DIR, "test_data", "eps=0.05.npy")
    with open(filename, "rb") as f:
        data = np.load(f)

    S0 = {"zeta": 1, "eta": 0.1}
    _, eta1 = ExtCzjzekDistribution(S0, eps=0.05).rvs(size=COUNT)
    hist1, _ = np.histogram(eta1, bins=100, range=[0, 1])

    message = "failed to compare values with file eps=0.05.npy"
    np.testing.assert_almost_equal(hist1 / COUNT, data[0], decimal=2, err_msg=message)


def test_extended_czjzek_polar():
    S0 = {"zeta": 1, "eta": 0.1}
    x, y = ExtCzjzekDistribution(S0, eps=0.05, polar=True).rvs(size=COUNT)
    x1, y1 = zeta_eta_to_x_y(*x_y_to_zeta_eta(x, y))
    np.testing.assert_almost_equal(x, x1)
    np.testing.assert_almost_equal(y, y1)


def test_extended_czjzek_eta_distribution_2():
    filename = path.join(MODULE_DIR, "test_data", "eps=0.2.npy")
    with open(filename, "rb") as f:
        data = np.load(f)

    S0 = {"Cq": 1e6, "eta": 0.3}
    _, eta1 = ExtCzjzekDistribution(S0, eps=0.2).rvs(size=COUNT)
    hist1, _ = np.histogram(eta1, bins=100, range=[0, 1])

    message = "failed to compare values with file eps=0.2.npy"
    np.testing.assert_almost_equal(hist1 / COUNT, data[1], decimal=2, err_msg=message)


def test_extended_czjzek_eta_distribution_3():
    filename = path.join(MODULE_DIR, "test_data", "eps=0.65.npy")
    with open(filename, "rb") as f:
        data = np.load(f)

    S0 = {"Cq": 1e6, "eta": 0.7}
    _, eta1 = ExtCzjzekDistribution(S0, eps=0.65).rvs(size=COUNT)
    hist1, _ = np.histogram(eta1, bins=100, range=[0, 1])

    message = "failed to compare values with file eps=0.05.npy"
    np.testing.assert_almost_equal(hist1 / COUNT, data[3], decimal=2, err_msg=message)


def test_czjzek_distribution():
    sigma = 0.5

    # numerical Czjzek distribution
    count_ = COUNT
    zeta, eta = CzjzekDistribution(sigma).rvs(size=count_)

    # eta projection
    e_hist, ran_e = np.histogram(eta, bins=15, range=[0, 1])
    e_vector = e_hist / count_
    e_range = (ran_e[1:] + ran_e[:-1]) / 2

    # zeta projection
    z_hist, ran_z = np.histogram(zeta, bins=20, range=[-20, 20])
    z_vector = z_hist / count_
    z_range = (ran_z[1:] + ran_z[:-1]) / 2

    # czjzek distribution from analytical formula
    res = CzjzekDistribution(sigma).pdf(pos=[z_range, e_range], pack_as_csdm=True)

    eta_pro = res.sum(axis=0).y[0].components[0]
    zeta_pro = res.sum(axis=1).y[0].components[0]

    # eta test
    message = "failed to compare eta projection for Czjzek distribution"
    np.testing.assert_almost_equal(e_vector, eta_pro, decimal=2, err_msg=message)

    # zeta test
    message = "failed to compare zeta projection for Czjzek distribution"
    np.testing.assert_almost_equal(z_vector, zeta_pro, decimal=2, err_msg=message)


def test_czjzek_pdf_full_eta():
    sigma = 0.5
    z_range = np.arange(100) * 30 / 100 - 15
    e_range = np.arange(21) / 20

    # czjzek distribution from analytical formula
    _, _, res = CzjzekDistribution(sigma).pdf(pos=[z_range, e_range], analytical=True)

    _, _, amp = CzjzekDistribution(sigma).pdf([z_range, e_range], analytical=False)

    error = "Czjzek analytical is not equal to numerical"
    np.testing.assert_almost_equal(res, amp, decimal=2, err_msg=error)


def test_czjzek_polar_analytical_vs_numerical():
    x_range = np.linspace(0, 12, 100)
    y_range = np.linspace(0, 12, 100)

    czjzek_dist = CzjzekDistribution(sigma=1.25, polar=True)
    # The 2D amplitude grid of Cq and eta is sampled from the Czjzek distribution model.
    _, _, amp_ana = czjzek_dist.pdf(pos=[x_range, y_range], analytical=True)
    _, _, amp_num = czjzek_dist.pdf(
        pos=[x_range, y_range], analytical=False, size=1_000_000
    )

    error = "Czjzek analytical is not equal to numerical"
    np.testing.assert_almost_equal(amp_ana, amp_num, decimal=2, err_msg=error)


def test_czjzek_pdf_sub_eta():
    sigma = 0.5
    z_range = np.arange(100) * 30 / 100 - 15
    e_range = np.arange(21) / 21

    # czjzek distribution from analytical formula
    _, _, res = CzjzekDistribution(sigma).pdf(pos=[z_range, e_range], analytical=True)

    _, _, amp = CzjzekDistribution(sigma).pdf([z_range, e_range], analytical=False)

    error = "Czjzek analytical is not equal to numerical"
    np.testing.assert_almost_equal(res, amp, decimal=2, err_msg=error)


def test_czjzek_polar():
    x, y = CzjzekDistribution(sigma=0.5, polar=True).rvs(size=COUNT)
    x1, y1 = zeta_eta_to_x_y(*x_y_to_zeta_eta(x, y))
    np.testing.assert_almost_equal(x, x1)
    np.testing.assert_almost_equal(y, y1)


def test_extended_czjzek_axially_symmetric_dominant_tensor():
    """eta=0 dominant tensor must still be perturbed (rho != 0)."""
    n = 200000
    for S0 in [{"Cq": 1e6, "eta": 0.0}, {"zeta": 10, "eta": 0.0}]:
        zeta0, eta0 = ExtCzjzekDistribution(S0, eps=0.3).rvs(size=n)
        S1 = {**S0, "eta": 1e-9}
        zeta1, eta1 = ExtCzjzekDistribution(S1, eps=0.3).rvs(size=n)

        assert np.std(zeta0) > 0
        assert np.allclose(np.std(zeta0), np.std(zeta1), rtol=0.05)
        assert np.allclose(np.mean(eta0), np.mean(eta1), rtol=0.05)


def test_czjzek_sigma_update_with_cache():
    """With the cached basis, a change of sigma scales the sampled zeta, and a change of
    size returns the requested number of samples."""
    model = CzjzekDistribution(sigma=1.0)
    zeta1, eta1 = model.rvs(size=5000)

    model.sigma = 5.0
    zeta5, eta5 = model.rvs(size=5000)
    assert np.allclose(zeta5, 5 * zeta1)
    assert np.allclose(eta5, eta1)

    zeta, eta = model.rvs(size=100)
    assert zeta.size == 100 and eta.size == 100

    model = ExtCzjzekDistribution({"zeta": 10, "eta": 0.3}, eps=0.2)
    zeta, eta = model.rvs(size=1000)
    zeta, eta = model.rvs(size=50)
    assert zeta.size == 50 and eta.size == 50


def test_czjzek_sigma_without_cache():
    """Without the cache, sigma sets the scale of the distribution."""
    std = [
        np.std(CzjzekDistribution(sigma=s, cache=False).rvs(size=200000)[0])
        for s in [1.0, 3.0]
    ]
    assert np.allclose(std[1] / std[0], 3.0, rtol=0.02)


def test_haeberlen_components_axial_and_zero_tensor():
    """Axially symmetric and zero tensors (rho = 0) give finite Haeberlen components
    and do not crash. Run in a subprocess, since the failure mode is a segfault."""
    import subprocess
    import sys

    code = """
import numpy as np
from mrsimulator.base_model import get_Haeberlen_components
from mrsimulator.models.czjzek import czjzek_random_components, get_p_q_basis

np.random.seed(1)
basis = get_p_q_basis(*czjzek_random_components(1, 2000))
for zeta in np.linspace(-50, 50, 41):
    for rho in [0.0, 1e-14, 1e-12]:
        z, e = get_Haeberlen_components(*basis, zeta, 0.0, rho)
        assert np.all(np.isfinite(z)) and np.all(np.isfinite(e)), (zeta, rho)
        if rho == 0:
            assert np.allclose(z, zeta, atol=1e-9) and np.allclose(e, 0, atol=1e-6)
print("ok")
"""
    result = subprocess.run(
        [sys.executable, "-c", code], capture_output=True, text=True
    )
    assert result.returncode == 0, result.stderr[-2000:]
    assert result.stdout.strip() == "ok"
