"""Test amplitude for shift, reference offset, points, orientation averaging."""
import numpy as np
from mrsimulator import Simulator
from mrsimulator import Site
from mrsimulator import SpinSystem
from mrsimulator.method import Method
from mrsimulator.method import MixingEvent
from mrsimulator.method import SpectralDimension
from mrsimulator.method import SpectralEvent
from mrsimulator.method.lib import BlochDecaySpectrum
from mrsimulator.spin_system.tensors import SymmetricTensor


def pre_setup():
    site_1 = Site(isotope="13C", shielding_symmetric={"zeta": 50, "eta": 0.5})
    spin_system = SpinSystem(sites=[site_1])
    method = BlochDecaySpectrum(
        channels=["13C"], spectral_dimensions=[{"count": 1024, "spectral_width": 25000}]
    )

    sim = Simulator(spin_systems=[spin_system], methods=[method])
    return sim


def test_static_spinning_integral_amplitude():
    sim = pre_setup()
    sim.run()
    y_static = sim.methods[0].simulation.y[0].components[0].sum()

    sim.methods[0].spectral_dimensions[0].events[0].rotor_frequency = 1000  # in Hz
    sim.run()
    y_MAS = sim.methods[0].simulation.y[0].components[0].sum()

    sim.methods[0].spectral_dimensions[0].events[0].rotor_angle = 0.5  # in rad
    sim.run()
    y_VAS = sim.methods[0].simulation.y[0].components[0].sum()

    print(y_static, y_MAS, y_VAS)
    assert np.allclose(y_static, y_MAS), "Integral error from static to MAS."
    assert np.allclose(y_static, y_VAS), "Integral error from static to VAS."


def test_with_configuration_setting():
    sim = pre_setup()
    sim.config.integration_density = 64
    sim.run()
    y_static = sim.methods[0].simulation.y[0].components[0].sum()

    sim.config.integration_density = 128
    sim.run()
    y_static_1 = sim.methods[0].simulation.y[0].components[0].sum()

    sim.config.integration_volume = "sphere"
    sim.run()
    y_static_3 = sim.methods[0].simulation.y[0].components[0].sum()

    sim.config.integration_volume = "hemisphere"
    sim.run()
    y_static_2 = sim.methods[0].simulation.y[0].components[0].sum()

    sim.methods[0].spectral_dimensions[0].events[0].rotor_frequency = 1000  # in Hz
    sim.run()
    y_MAS = sim.methods[0].simulation.y[0].components[0].sum()

    sim.config.number_of_sidebands = 256
    sim.run()
    y_MAS_1 = sim.methods[0].simulation.y[0].components[0].sum()

    e = "Integral error from changing integration density."
    np.testing.assert_almost_equal(y_static, y_static_1, decimal=3, err_msg=e)

    e = "Integral error from changing integration volume octant to hemisphere."
    np.testing.assert_almost_equal(y_static, y_static_2, decimal=3, err_msg=e)

    e = "Integral error from changing integration volume from octant to sphere."
    np.testing.assert_almost_equal(y_static, y_static_3, decimal=3, err_msg=e)

    e = "Integral error from νr, integration density, integration volume."
    np.testing.assert_almost_equal(y_static, y_MAS, decimal=3, err_msg=e)

    e = "Integral error from number of sidebands."
    np.testing.assert_almost_equal(y_MAS, y_MAS_1, decimal=8, err_msg=e)


def test_number_of_points():
    sim = pre_setup()
    sim.methods[0].spectral_dimensions[0].count = 2048
    sim.run()
    y_static = sim.methods[0].simulation.y[0].components[0].sum()
    inc = sim.methods[0].simulation.x[0].increment.to("Hz").value

    sim.methods[0].spectral_dimensions[0].count = 4096
    sim.run()
    y_static_1 = sim.methods[0].simulation.y[0].components[0].sum()
    inc_1 = sim.methods[0].simulation.x[0].increment.to("Hz").value

    e = "Integral error from number_of_points"
    np.testing.assert_almost_equal(
        y_static * inc, y_static_1 * inc_1, decimal=8, err_msg=e
    )


def test_gamma_averaging():
    sim = pre_setup()
    sim.run()
    y_static = sim.methods[0].simulation.y[0].components[0].sum()

    sim.config.number_of_gamma_angles = 1
    sim.run(auto_switch=False)
    y_static_1 = sim.methods[0].simulation.y[0].components[0].sum()

    sim.config.number_of_gamma_angles = 500
    sim.run(auto_switch=False)
    y_static_2 = sim.methods[0].simulation.y[0].components[0].sum()

    e = "Integral error from number_of_gamma_angles"
    np.testing.assert_almost_equal(y_static, y_static_1, decimal=8, err_msg=e)

    e = "Integral error from number_of_gamma_angles"
    np.testing.assert_almost_equal(y_static, y_static_2, decimal=8, err_msg=e)


def test_complex_pathway_weight():
    """Create Method objects with complex pathway weights and test resulting complex
    spectral amplitudes.
    """
    # NOTE: This test can be cleaned up and/or moved to a different file
    site = Site(
        isotope="13C",
        isotropic_chemical_shift=-20,
        shielding_symmetric=SymmetricTensor(zeta=30, eta=0.3),
    )

    sys = SpinSystem(sites=[site])

    # Create different mixing phases
    n_phases = 32
    phases = np.asarray([n * np.pi / n_phases for n in range(n_phases)])

    # Create Methods emulating BlochDecaySpectrum, but with complex weights

    methods = [
        Method(
            channels=["13C"],
            rotor_frequency=0,
            rotor_angle=0,
            magnetic_flux_density=9.4,
            spectral_dimensions=[
                SpectralDimension(
                    count=4096,
                    spectral_width=1.5e4,  # 15 kHz
                    reference_offset=-2500,
                    events=[
                        # NOTE: Error thrown when a MixingEvent not sandwiched between
                        # SpectralEvent objects. Should this be addressed?
                        SpectralEvent(
                            transition_queries=[{"ch1": {"P": [+1]}}], fraction=0
                        ),  # Dummy spectral event
                        MixingEvent(ch1={"angle": np.pi, "phase": ph}),
                        SpectralEvent(
                            transition_queries=[{"ch1": {"P": [-1]}}], fraction=1
                        ),
                    ],
                )
            ],
        )
        for ph in phases
    ]

    pathway_weights = [mth.get_transition_pathways(sys)[0].weight for mth in methods]

    sim = Simulator(
        spin_systems=[sys],
        methods=methods,
    )
    sim.config.integration_density = 120

    sim.run()

    # Test complex scaling
    A = np.abs(sim.methods[0].simulation)
    for i, mth in enumerate(sim.methods):
        B = mth.simulation
        re_weight = np.real(pathway_weights[i])
        im_weight = np.imag(pathway_weights[i])

        assert B.real == A * re_weight
        assert B.imag == A * im_weight
