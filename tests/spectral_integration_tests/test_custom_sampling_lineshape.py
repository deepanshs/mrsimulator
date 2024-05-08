import numpy as np
import pytest
from mrsimulator import Simulator
from mrsimulator import Site
from mrsimulator import SpinSystem
from mrsimulator.method import SpectralDimension
from mrsimulator.method.lib import BlochDecaySpectrum
from mrsimulator.simulator.sampling_scheme import step_averaging
from mrsimulator.simulator.sampling_scheme import zcw_averaging
from mrsimulator.spin_system.tensors import SymmetricTensor


def setup_sim():
    site = Site(
        isotope="2H",
        isotropic_chemical_shift=-96,  # in ppm
        shielding_symmetric=SymmetricTensor(
            zeta=-566, eta=0.0, alpha=0.7, beta=2.0, gamma=3.14
        ),
        quadrupolar=SymmetricTensor(Cq=75.2e3, eta=0.97),
    )
    spin_systems = [SpinSystem(sites=[site])]
    method = BlochDecaySpectrum(
        channels=["2H"],
        magnetic_flux_density=9.395,  # in T
        rotor_frequency=0,  # in Hz
        rotor_angle=0,  # in rads
        spectral_dimensions=[
            SpectralDimension(count=1024, spectral_width=400000),
        ],
    )
    sim = Simulator(spin_systems=spin_systems, methods=[method])
    return sim


def test_one_d():
    sim = setup_sim()
    sim.config.integration_volume = "sphere"
    sim.config.integration_density = 24
    sim.run()
    spec_asg = sim.methods[0].simulation.y[0].components[0].real

    sim.config.custom_sampling = zcw_averaging(M=12)
    assert sim.config.get_orientations_count() == 2584
    sim.run()
    spec_zcw = sim.methods[0].simulation.y[0].components[0].real

    sim.config.custom_sampling = step_averaging(N_alpha=51, N_beta=51)
    assert sim.config.get_orientations_count() == 2601
    sim.run()
    spec_step = sim.methods[0].simulation.y[0].components[0].real

    np.testing.assert_almost_equal(spec_asg, spec_zcw, decimal=3)
    np.testing.assert_almost_equal(spec_asg, spec_step, decimal=3)


def test_internal_external_averaging_spectrum():
    site = Site(
        isotope="2H",
        isotropic_chemical_shift=-96,  # in ppm
        shielding_symmetric=SymmetricTensor(
            zeta=-566, eta=0.0, alpha=0.7, beta=2.0, gamma=3.14
        ),
        quadrupolar=SymmetricTensor(Cq=75.2e3, eta=0.97),  # Cq in Hz
    )

    spin_systems = [SpinSystem(sites=[site])]

    method = BlochDecaySpectrum(
        channels=["2H"],
        magnetic_flux_density=9.395,  # in T
        rotor_frequency=0,  # in Hz
        rotor_angle=0,  # in rads
        spectral_dimensions=[
            SpectralDimension(count=2046, spectral_width=400000),
        ],
    )

    sim = Simulator(spin_systems=spin_systems, methods=[method])
    sim.config.integration_volume = "hemisphere"

    # ASG interpolation
    sim.config.integration_density = 44
    sim.run()
    spec_asg_interp = sim.methods[0].simulation.y[0].components[0].real

    # ASG binning
    sim.config.integration_volume = "hemisphere"
    sim.config.integration_density = 400
    sim.run(interpolation=False)
    spec_asg_bin = sim.methods[0].simulation.y[0].components[0].real
    np.testing.assert_almost_equal(spec_asg_interp, spec_asg_bin, decimal=2)

    # ZCW binning
    sim.config.custom_sampling = zcw_averaging(M=15)
    sim.run()
    spec_zcw_interp = sim.methods[0].simulation.y[0].components[0].real

    # ZCW interpolation
    sim.config.custom_sampling = zcw_averaging(M=23, triangle_mesh=False)
    sim.run()
    spec_zcw_bin = sim.methods[0].simulation.y[0].components[0].real

    np.testing.assert_almost_equal(spec_zcw_interp, spec_zcw_bin, decimal=2)

    # STEP interpolate
    sim.config.custom_sampling = step_averaging(N_alpha=100, N_beta=100)
    sim.run()
    spec_step_interp = sim.methods[0].simulation.y[0].components[0].real

    # STEP binning
    sim.config.custom_sampling = step_averaging(
        N_alpha=1160, N_beta=1160, triangle_mesh=False
    )
    sim.run()
    spec_step_bin = sim.methods[0].simulation.y[0].components[0].real

    np.testing.assert_almost_equal(spec_step_interp, spec_step_bin, decimal=2)
    np.testing.assert_almost_equal(spec_asg_interp, spec_step_interp, decimal=3)
    np.testing.assert_almost_equal(spec_asg_interp, spec_zcw_interp, decimal=3)


def test_sampling_triangulation():
    error = "Triangulation of non sphere geometry is not implemented."
    for vol in ["octant", "hemisphere"]:
        with pytest.raises(NotImplementedError, match=error):
            _ = step_averaging(
                N_alpha=50, N_beta=50, triangle_mesh=True, integration_volume=vol
            )

        with pytest.raises(NotImplementedError, match=error):
            _ = zcw_averaging(M=5, triangle_mesh=True, integration_volume=vol)


def test_sampling_non_sphere():
    for vol in ["octant", "hemisphere"]:
        step_s = step_averaging(
            N_alpha=50, N_beta=50, triangle_mesh=False, integration_volume=vol
        )
        assert step_s.alpha.size == 2500
        zcw_s = zcw_averaging(M=5, triangle_mesh=False, integration_volume=vol)
        assert zcw_s.alpha.size == 89
