from lmfit import Parameters
from mrsimulator.models import CzjzekDistribution
from mrsimulator.models import ExtCzjzekDistribution


def test_dist_1():
    dist = CzjzekDistribution(sigma=1.5, mean_isotropic_chemical_shift=23.4)
    prefix = "czjzek"

    params = Parameters()
    params = dist.add_lmfit_params(params=params, i=0)
    assert len(params.keys()) == 3

    assert params[f"{prefix}_0_sigma"] == 1.5
    assert params[f"{prefix}_0_mean_isotropic_chemical_shift"] == 23.4
    assert params[f"{prefix}_0_abundance"] == 100.0

    params[f"{prefix}_0_sigma"].value = 10
    params[f"{prefix}_0_mean_isotropic_chemical_shift"].value = -63.4
    params[f"{prefix}_0_abundance"].value = 0.15

    dist.update_lmfit_params(params=params, i=0)
    assert dist.sigma == 10.0
    assert dist.mean_isotropic_chemical_shift == -63.4
    assert dist.abundance == 0.15


def test_dist_2():
    dist = ExtCzjzekDistribution(
        symmetric_tensor={"zeta": 4.1, "eta": 0.2},
        eps=0.3,
        mean_isotropic_chemical_shift=0.0,
    )
    prefix = "ext_czjzek"

    params = Parameters()
    params = dist.add_lmfit_params(params=params, i=0)
    assert len(params.keys()) == 5

    assert params[f"{prefix}_0_symmetric_tensor_zeta"] == 4.1
    assert params[f"{prefix}_0_symmetric_tensor_eta"] == 0.2
    assert params[f"{prefix}_0_eps"] == 0.3
    assert params[f"{prefix}_0_mean_isotropic_chemical_shift"] == 0.0
    assert params[f"{prefix}_0_abundance"] == 100.0

    params[f"{prefix}_0_symmetric_tensor_zeta"].value = -30.2
    params[f"{prefix}_0_symmetric_tensor_eta"].value = 0.65
    params[f"{prefix}_0_eps"].value = 0.6
    params[f"{prefix}_0_mean_isotropic_chemical_shift"].value = -3.4
    params[f"{prefix}_0_abundance"].value = 0.9

    dist.update_lmfit_params(params=params, i=0)
    assert dist.symmetric_tensor.zeta == -30.2
    assert dist.symmetric_tensor.eta == 0.65
    assert dist.eps == 0.6
    assert dist.mean_isotropic_chemical_shift == -3.4
    assert dist.abundance == 0.9
