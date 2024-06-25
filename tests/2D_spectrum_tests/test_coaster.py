import numpy as np
from mrsimulator import Simulator
from mrsimulator import Site
from mrsimulator import SpinSystem
from mrsimulator.method import Method
from mrsimulator.method import MixingEvent
from mrsimulator.method import SpectralDimension
from mrsimulator.method import SpectralEvent
from mrsimulator.spin_system.tensors import SymmetricTensor


def setup_spin_system():
    site = Site(
        isotope="87Rb",
        isotropic_chemical_shift=-9,
        shielding_symmetric=SymmetricTensor(zeta=110, eta=0),
        quadrupolar=SymmetricTensor(Cq=3.5e6, eta=0.36, beta=70 * 3.14159 / 180),
    )
    return SpinSystem(sites=[site])


def process_spectrum(method):
    spin_system = setup_spin_system()

    sim = Simulator(spin_systems=[spin_system], methods=[method])
    sim.config.integration_volume = "hemisphere"
    sim.run()

    n_dim = len(method.spectral_dimensions)
    assert str(sim.methods[0].simulation.y[0].unit) == f"Hz^-{n_dim}"

    data = sim.methods[0].simulation.y[0].components[0]
    data /= data.max()
    return data


def coaster_simulation():
    coaster = Method(
        name="COASTER",
        channels=["87Rb"],
        magnetic_flux_density=9.4,  # in T
        rotor_angle=70.12 * np.pi / 180,  # in rads
        rotor_frequency=np.inf,
        spectral_dimensions=[
            SpectralDimension(
                count=512,
                spectral_width=4e4,  # in Hz
                reference_offset=-8e3,  # in Hz
                label="$\\omega_1$ (CSA)",
                events=[
                    SpectralEvent(transition_queries=[{"ch1": {"P": [3], "D": [0]}}]),
                    MixingEvent(ch1={"angle": np.pi * 109.5 / 180}),
                ],
            ),
            SpectralDimension(
                count=512,
                spectral_width=8e3,  # in Hz
                reference_offset=-4e3,  # in Hz
                label="$\\omega_2$ (Q)",
                events=[
                    SpectralEvent(transition_queries=[{"ch1": {"P": [-1], "D": [0]}}])
                ],
            ),
        ],
        affine_matrix=[[1, 0], [1 / 4, 3 / 4]],
    )
    return process_spectrum(coaster)


def csa_1D_projection():
    csa_only = Method(
        channels=["87Rb"],
        magnetic_flux_density=9.4,  # in T
        rotor_angle=70.12 * np.pi / 180,  # in rads
        rotor_frequency=np.inf,
        spectral_dimensions=[
            SpectralDimension(
                count=512,
                spectral_width=4e4,
                reference_offset=-8e3,
                events=[
                    SpectralEvent(
                        freq_contrib=["Quad2_0", "Shielding1_0", "Shielding1_2"],
                        transition_queries=[{"ch1": {"P": [3], "D": [0]}}],
                    ),
                ],
            )
        ],
    )
    return process_spectrum(csa_only)


def quad_1D_projection():
    quad_only = Method(
        channels=["87Rb"],
        magnetic_flux_density=9.4,  # in T
        rotor_angle=70.12 * np.pi / 180,  # in rads
        rotor_frequency=np.inf,
        spectral_dimensions=[
            SpectralDimension(
                count=512,
                spectral_width=8e3,
                reference_offset=-4e3,
                events=[
                    SpectralEvent(
                        fraction=1 / 4,
                        freq_contrib=["Quad2_0"],
                        transition_queries=[{"ch1": {"P": [3], "D": [0]}}],
                    ),
                    MixingEvent(ch1={"angle": np.pi * 109.5 / 180, "phase": 0}),
                    SpectralEvent(
                        fraction=3 / 4,
                        freq_contrib=["Quad2_0", "Quad2_2"],
                        transition_queries=[{"ch1": {"P": [-1], "D": [0]}}],
                    ),
                ],
            )
        ],
    )
    return process_spectrum(quad_only)


def test_01():
    coaster = coaster_simulation()
    csa_1D = csa_1D_projection()
    quad_1D = quad_1D_projection()

    csa_proj = coaster.sum(axis=1)
    quad_proj = coaster.sum(axis=0)

    csa_proj /= csa_proj.max()
    quad_proj /= quad_proj.max()

    np.testing.assert_almost_equal(csa_proj, csa_1D, decimal=2)
    np.testing.assert_almost_equal(quad_proj, quad_1D, decimal=1)
