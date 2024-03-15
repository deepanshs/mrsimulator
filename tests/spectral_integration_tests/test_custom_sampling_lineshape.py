import numpy as np
from mrsimulator import Simulator
from mrsimulator import Site
from mrsimulator import SpinSystem
from mrsimulator.method import SpectralDimension
from mrsimulator.method.lib import BlochDecaySpectrum
from mrsimulator.simulator.sampling_scheme import STEP_averaging
from mrsimulator.simulator.sampling_scheme import ZCW_averaging
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

    sim.config.custom_sampling = ZCW_averaging(M=12)
    sim.run()
    spec_zcw = sim.methods[0].simulation.y[0].components[0].real

    sim.config.custom_sampling = STEP_averaging(N_alpha=51, N_beta=51)
    sim.run()
    spec_step = sim.methods[0].simulation.y[0].components[0].real

    np.testing.assert_almost_equal(spec_asg, spec_zcw, decimal=3)
    np.testing.assert_almost_equal(spec_asg, spec_step, decimal=3)
