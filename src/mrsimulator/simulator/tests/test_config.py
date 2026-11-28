"""Test for the base ConfigSimulator class."""
import pytest
from mrsimulator import Simulator

__author__ = "Deepansh Srivastava"
__email__ = "srivastava.89@osu.edu"


def test_config():

    # set config
    b = Simulator()

    error = "value is not a valid dict"
    with pytest.raises(ValueError, match=f".*{error}.*"):
        b.config = 5

    a = Simulator()

    assert b == a

    # number of sidebands
    assert a.config.number_of_sidebands == 64

    a.config.number_of_sidebands = 1.4
    assert a.config.number_of_sidebands == 1

    a.config.number_of_sidebands = 10
    assert a.config.number_of_sidebands == 10

    error = "ensure this value is greater than 0"
    with pytest.raises(ValueError, match=f".*{error}.*"):
        a.config.number_of_sidebands = 0

    # integration density
    assert a.config.integration_density == 70
    a.config.integration_density = 20
    assert a.config.integration_density == 20

    error = "ensure this value is greater than 0"
    with pytest.raises(ValueError, match=f".*{error}.*"):
        a.config.integration_density = -12

    error = "value is not a valid integer"
    with pytest.raises(ValueError, match=f".*{error}.*"):
        a.config.integration_density = {}

    # integration volume
    assert a.config.integration_volume == "octant"
    a.config.integration_volume = "hemisphere"
    assert a.config.integration_volume == "hemisphere"

    error = "unexpected value; permitted: 'octant', 'hemisphere'"
    with pytest.raises(ValueError, match=f".*{error}.*"):
        a.config.integration_volume = "pentagon"

    # decompose spectrum
    assert a.config.decompose_spectrum == "none"
    a.config.decompose_spectrum = "spin_system"
    assert a.config.decompose_spectrum == "spin_system"

    error = "unexpected value; permitted: 'none', 'spin_system'"
    with pytest.raises(ValueError, match=f".*{error}.*"):
        a.config.decompose_spectrum = "haha"

    # isotropic interpolation
    assert a.config.isotropic_interpolation == "linear"
    a.config.isotropic_interpolation = "gaussian"
    assert a.config.isotropic_interpolation == "gaussian"

    error = "unexpected value; permitted: 'linear', 'gaussian'"
    with pytest.raises(ValueError, match=f".*{error}.*"):
        a.config.isotropic_interpolation = "haha"

    # number of gamma angles
    assert a.config.number_of_gamma_angles == 1
    a.config.number_of_gamma_angles = 14
    assert a.config.number_of_gamma_angles == 14

    error = "ensure this value is greater than 0"
    with pytest.raises(ValueError, match=f".*{error}.*"):
        a.config.integration_density = -1

    # overall hemisphere
    assert a.config.dict(exclude={"property_units"}) == {
        "custom_sampling": None,
        "decompose_spectrum": "spin_system",
        "number_of_sidebands": 10,
        "number_of_gamma_angles": 14,
        "integration_volume": "hemisphere",
        "integration_density": 20,
        "isotropic_interpolation": "gaussian",
        "name": None,
        "description": None,
        "label": None,
        "is_complex": True,
    }

    assert a.config.get_int_dict() == {
        "decompose_spectrum": 1,
        "number_of_sidebands": 10,
        "number_of_gamma_angles": 14,
        "integration_volume": 1,
        "integration_density": 20,
        "isotropic_interpolation": 1,
        "is_complex": True,
    }

    assert b != a

    # get orientation count
    assert a.config.get_orientations_count() == 4 * 21 * 22 * 14 / 2

    # overall sphere
    a.config.integration_volume = "sphere"
    assert a.config.integration_volume == "sphere"
    assert a.config.dict(exclude={"property_units"}) == {
        "custom_sampling": None,
        "decompose_spectrum": "spin_system",
        "number_of_sidebands": 10,
        "number_of_gamma_angles": 14,
        "integration_volume": "sphere",
        "integration_density": 20,
        "isotropic_interpolation": "gaussian",
        "name": None,
        "description": None,
        "label": None,
        "is_complex": True,
    }

    assert a.config.get_int_dict() == {
        "decompose_spectrum": 1,
        "number_of_sidebands": 10,
        "number_of_gamma_angles": 14,
        "integration_volume": 2,
        "integration_density": 20,
        "isotropic_interpolation": 1,
        "is_complex": True,
    }

    assert b != a

    # get orientation count
    assert a.config.get_orientations_count() == 4 * 21 * 22 * 14


def test_custom_sampling_validation():
    """The arrays of a custom sampling are of the same size, and the vertex indexes
    refer to the orientations. Lists and arrays of any numeric type are accepted."""
    import numpy as np
    from mrsimulator.simulator.config import CustomSampling

    sampling = CustomSampling(
        alpha=[0, 1, 2, 3],
        beta=np.arange(4, dtype=np.float32),
        weight=np.ones(4, dtype=int),
        vertex_indexes=[[0, 1, 2], [1, 2, 3]],
    )
    assert sampling.alpha.dtype == np.float64 and sampling.beta.dtype == np.float64
    assert sampling.weight.dtype == np.float64
    assert sampling.vertex_indexes.dtype == np.int32
    assert sampling.vertex_indexes.shape == (2, 3)

    with pytest.raises(ValueError, match=".*same size.*"):
        CustomSampling(alpha=np.zeros(200), beta=np.zeros(10), weight=np.ones(10))
    with pytest.raises(ValueError, match=".*vertex_indexes.*"):
        CustomSampling(
            alpha=np.zeros(4), beta=np.zeros(4), weight=np.ones(4),
            vertex_indexes=[[0, 1, 4]],
        )
    with pytest.raises(ValueError, match=".*vertex_indexes.*"):
        CustomSampling(
            alpha=np.zeros(4), beta=np.zeros(4), weight=np.ones(4),
            vertex_indexes=[[0, 1], [2, 3]],
        )


def test_custom_sampling_json():
    """A simulator with a custom sampling is serialized and parsed back."""
    import numpy as np
    from mrsimulator import Site
    from mrsimulator import SpinSystem
    from mrsimulator.method.lib import BlochDecaySpectrum
    from mrsimulator.simulator.sampling_scheme import zcw_averaging

    sim = Simulator(
        spin_systems=[SpinSystem(sites=[Site(isotope="1H")])],
        methods=[BlochDecaySpectrum(channels=["1H"])],
    )
    for mesh in [False, True]:
        sim.config.custom_sampling = zcw_averaging(M=3, triangle_mesh=mesh)
        py_dict = sim.json()
        sampling = Simulator.parse(py_dict).config.custom_sampling
        assert np.allclose(sampling.alpha, sim.config.custom_sampling.alpha)
        assert np.allclose(sampling.weight, sim.config.custom_sampling.weight)
        if mesh:
            assert np.array_equal(
                sampling.vertex_indexes, sim.config.custom_sampling.vertex_indexes
            )
