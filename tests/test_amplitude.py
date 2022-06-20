"""Test amplitude for shift, reference offset, points, orientation averaging."""
import numpy as np
from mrsimulator import Simulator
from mrsimulator import Site
from mrsimulator import SpinSystem
from mrsimulator.method.lib import BlochDecaySpectrum


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
    np.testing.assert_almost_equal(y_static, y_static_1, decimal=1, err_msg=e)

    e = "Integral error from changing integration volume."
    np.testing.assert_almost_equal(y_static, y_static_2, decimal=1, err_msg=e)

    e = "Integral error from νr, integration density, integration volume."
    np.testing.assert_almost_equal(y_static, y_MAS, decimal=1, err_msg=e)

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
