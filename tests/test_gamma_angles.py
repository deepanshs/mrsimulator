"""Test gamma angles."""
import numpy as np
from mrsimulator import Simulator
from mrsimulator import Site
from mrsimulator import SpinSystem
from mrsimulator.method import Method
from mrsimulator.method import RotationEvent
from mrsimulator.method import SpectralDimension
from mrsimulator.method import SpectralEvent
from mrsimulator.method.lib import BlochDecaySpectrum

# import matplotlib.pyplot as plt


def setup_test(spin_system, volume="octant", sw=25000, n_gamma=500):
    mth_kwargs = {
        "channels": [spin_system.sites[0].isotope.symbol],
        "spectral_dimensions": [{"count": 1024, "spectral_width": sw}],
    }

    data = []
    for angle in [0, 54.735, 90]:
        method = BlochDecaySpectrum(
            rotor_angle=angle * np.pi / 180, rotor_frequency=0, **mth_kwargs
        )
        sim = Simulator(spin_systems=[spin_system], methods=[method])
        sim.config.integration_volume = volume
        sim.config.number_of_gamma_angles = 1 if angle == 0 else n_gamma
        sim.run(auto_switch=False)

        res = sim.methods[0].simulation.y[0].components[0]
        res /= res.max()
        data.append(res)

    # plt.plot(data[0])
    # plt.plot(data[1], "--")
    # plt.plot(data[2], "-.")
    # plt.show()

    np.testing.assert_almost_equal(data[0], data[1], decimal=2)
    np.testing.assert_almost_equal(data[0], data[2], decimal=1.6)


def test_csa_01():
    site = Site(isotope="13C", shielding_symmetric={"zeta": 50, "eta": 0.5})
    spin_system = SpinSystem(sites=[site])
    setup_test(spin_system, volume="octant", sw=2.5e4)


def test_quad_01():
    site = Site(
        isotope="27Al",
        shielding_symmetric={"zeta": 50, "eta": 0.5},
        quadrupolar={"Cq": 1e6, "eta": 0.1, "alpha": 1, "beta": 2, "gamma": 3},
    )
    spin_system = SpinSystem(sites=[site])
    setup_test(spin_system, volume="hemisphere", sw=1e6)


def test_2D():
    site_Ni = Site(
        isotope="2H",
        isotropic_chemical_shift=-97,  # in ppm
        shielding_symmetric=dict(
            zeta=-551,
            eta=0.12,
            alpha=62 * np.pi / 180,
            beta=114 * np.pi / 180,
            gamma=171 * np.pi / 180,
        ),
        quadrupolar=dict(Cq=77.2e3, eta=0.9),  # Cq in Hz
    )
    spin_system = SpinSystem(sites=[site_Ni])

    data = []
    for angle, n_gamma in zip([0, np.pi / 4], [1, 500]):
        shifting_d = Method(
            name="Shifting-d",
            channels=["2H"],
            magnetic_flux_density=9.395,  # in T
            rotor_frequency=0,  # in Hz
            rotor_angle=angle,  # in Hz
            spectral_dimensions=[
                SpectralDimension(
                    count=512,
                    spectral_width=2.5e5,  # in Hz
                    label="Quadrupolar frequency",
                    events=[
                        SpectralEvent(
                            transition_queries=[{"ch1": {"P": [-1]}}],
                            freq_contrib=["Quad1_2"],
                        ),
                        RotationEvent(),
                    ],
                ),
                SpectralDimension(
                    count=256,
                    spectral_width=2e5,  # in Hz
                    reference_offset=2e4,  # in Hz
                    label="Paramagnetic shift",
                    events=[
                        SpectralEvent(
                            transition_queries=[{"ch1": {"P": [-1]}}],
                            freq_contrib=["Shielding1_0", "Shielding1_2"],
                        )
                    ],
                ),
            ],
        )

        sim = Simulator(spin_systems=[spin_system], methods=[shifting_d])
        sim.config.integration_volume = "hemisphere"
        sim.config.number_of_gamma_angles = n_gamma
        sim.run(auto_switch=False)

        res = sim.methods[0].simulation.y[0].components[0]
        data.append(res / res.max())

    # _, ax = plt.subplots(1, 2)
    # ax[0].imshow(data[0].real)
    # ax[1].imshow(data[1].real)
    # plt.show()

    np.testing.assert_almost_equal(data[0], data[1], decimal=1.8)
