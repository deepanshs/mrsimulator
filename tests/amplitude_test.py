# -*- coding: utf-8 -*-
"""Test for shift and reference offset."""
import numpy as np
from mrsimulator import Isotopomer
from mrsimulator import Simulator
from mrsimulator import Site
from mrsimulator.methods import BlochDecayFT


def pre_setup():
    site_1 = Site(isotope="13C", shielding_symmetric={"zeta": 50, "eta": 0.5})
    isotopomer = Isotopomer(sites=[site_1])
    method = BlochDecayFT(
        isotope="13C", dimensions=[{"count": 1024, "spectral_width": 25000}]
    )

    sim = Simulator()
    sim.isotopomers.append(isotopomer)
    sim.methods += [method]
    return sim


def test_static_spinning_integral_amplitude():
    sim = pre_setup()
    sim.run()
    y_static = sim.methods[0].simulation.dependent_variables[0].components[0]

    sim.methods[0].spectral_dimensions[0].events[0].rotor_frequency = 1000  # in Hz
    sim.run()
    y_spinning_MAS = sim.methods[0].simulation.dependent_variables[0].components[0]

    sim.methods[0].spectral_dimensions[0].events[0].rotor_angle = 0.5  # in rad
    sim.run()
    y_spinning_VAS = sim.methods[0].simulation.dependent_variables[0].components[0]

    assert np.allclose(
        y_static.sum(), y_spinning_MAS.sum()
    ), "Integral error from static to MAS."
    assert np.allclose(
        y_static.sum(), y_spinning_VAS.sum()
    ), "Integral error from static to VAS."


def test_with_configuration_setting():
    sim = pre_setup()
    sim.config.integration_density = 64
    sim.run()
    y_static = sim.methods[0].simulation.dependent_variables[0].components[0]

    sim.config.integration_density = 128
    sim.run()
    y_static_1 = sim.methods[0].simulation.dependent_variables[0].components[0]

    sim.config.integration_volume = "hemisphere"
    sim.run()
    y_static_2 = sim.methods[0].simulation.dependent_variables[0].components[0]

    np.testing.assert_almost_equal(
        y_static.sum(),
        y_static_1.sum(),
        decimal=1,
        err_msg="Integral error from changing integration density.",
    )

    np.testing.assert_almost_equal(
        y_static.sum(),
        y_static_2.sum(),
        decimal=1,
        err_msg="Integral error from changing integration volume.",
    )

    sim.methods[0].spectral_dimensions[0].events[0].rotor_frequency = 1000  # in Hz
    sim.run()
    y_spinning_MAS = sim.methods[0].simulation.dependent_variables[0].components[0]

    np.testing.assert_almost_equal(
        y_static.sum(),
        y_spinning_MAS.sum(),
        decimal=1,
        err_msg="Integral error from νr, integration density, integration volume.",
    )

    sim.methods[0].spectral_dimensions[0].count = 256
    sim.run()
    y_spinning_MAS_1 = sim.methods[0].simulation.dependent_variables[0].components[0]

    np.testing.assert_almost_equal(
        y_spinning_MAS.sum(),
        y_spinning_MAS_1.sum(),
        decimal=8,
        err_msg="Integral error from number of sidebands.",
    )


def test_number_of_points():
    sim = pre_setup()
    sim.methods[0].spectral_dimensions[0].count = 2046
    sim.run()
    y_static = sim.methods[0].simulation.dependent_variables[0].components[0]

    sim.methods[0].spectral_dimensions[0].count = 4096
    sim.run()
    y_static_1 = sim.methods[0].simulation.dependent_variables[0].components[0]

    np.testing.assert_almost_equal(
        y_static.sum(),
        y_static_1.sum(),
        decimal=8,
        err_msg="Integral error from number_of_points",
    )
