# -*- coding: utf-8 -*-
"""Test for shift and reference offset."""
import numpy as np
from mrsimulator import Dimension
from mrsimulator import Isotopomer
from mrsimulator import Simulator
from mrsimulator import Site
from mrsimulator.methods import one_d_spectrum


def pre_setup():
    site_1 = Site(isotope="13C", shielding_symmetric={"zeta": 50, "eta": 0.5})
    isotopomer = Isotopomer(sites=[site_1])
    dimension = Dimension(isotope="13C", spectral_width=25000)

    sim = Simulator()
    sim.isotopomers.append(isotopomer)
    sim.dimensions.append(dimension)
    return sim


def test_static_spinning_integral_amplitude():
    sim = pre_setup()
    y_static = sim.run(method=one_d_spectrum)[1]

    sim.dimensions[0].rotor_frequency = 1000  # in Hz
    y_spinning_MAS = sim.run(method=one_d_spectrum)[1]

    sim.dimensions[0].rotor_angle = 0.5  # in rad
    y_spinning_VAS = sim.run(method=one_d_spectrum)[1]

    assert np.allclose(
        y_static.sum(), y_spinning_MAS.sum()
    ), "Integral error from static to MAS."
    assert np.allclose(
        y_static.sum(), y_spinning_VAS.sum()
    ), "Integral error from static to VAS."


def test_with_configuration_setting():
    sim = pre_setup()
    sim.config.integration_density = 64
    y_static = sim.run(method=one_d_spectrum)[1]

    sim.config.integration_density = 128
    y_static_1 = sim.run(method=one_d_spectrum)[1]

    sim.config.integration_volume = "hemisphere"
    y_static_2 = sim.run(method=one_d_spectrum)[1]

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

    sim.dimensions[0].rotor_frequency = 1000  # in Hz
    y_spinning_MAS = sim.run(method=one_d_spectrum)[1]

    np.testing.assert_almost_equal(
        y_static.sum(),
        y_spinning_MAS.sum(),
        decimal=1,
        err_msg="Integral error from νr, integration density, integration volume.",
    )

    sim.config.number_of_sidebands = 256
    y_spinning_MAS_1 = sim.run(method=one_d_spectrum)[1]

    np.testing.assert_almost_equal(
        y_spinning_MAS.sum(),
        y_spinning_MAS_1.sum(),
        decimal=8,
        err_msg="Integral error from number of sidebands.",
    )


def test_number_of_points():
    sim = pre_setup()
    sim.dimensions[0].number_of_points = 2046
    y_static = sim.run(method=one_d_spectrum)[1]

    sim.dimensions[0].number_of_points = 4096
    y_static_1 = sim.run(method=one_d_spectrum)[1]

    np.testing.assert_almost_equal(
        y_static.sum(),
        y_static_1.sum(),
        decimal=8,
        err_msg="Integral error from number_of_points",
    )
