import numpy as np
import pytest
from mrsimulator.utils.euler_angles import _add_two_euler_angles
from mrsimulator.utils.euler_angles import _euler_angles_to_angle_phase
from mrsimulator.utils.euler_angles import combine_euler_angles
from mrsimulator.utils.euler_angles import wrap_between_pi


__author__ = "Matthew Giammar"
__email__ = "giammar.7@osu.edu"


def test_wrap_between_pi():
    assert np.isclose(wrap_between_pi(0), 0.0)
    assert np.isclose(wrap_between_pi(2 * np.pi), 0.0)
    assert np.isclose(wrap_between_pi(-2 * np.pi), 0.0)
    assert np.isclose(wrap_between_pi(np.pi), np.pi)
    assert np.isclose(wrap_between_pi(-np.pi), np.pi)
    assert np.isclose(wrap_between_pi(-np.pi / 2), -np.pi / 2)
    assert np.isclose(wrap_between_pi(3 * np.pi / 2), -np.pi / 2)
    assert np.isclose(wrap_between_pi(5 * np.pi / 2), np.pi / 2)
    assert np.isclose(wrap_between_pi(-3 * np.pi / 2), np.pi / 2)
    assert np.isclose(wrap_between_pi(-5 * np.pi / 2), -np.pi / 2)


def test_add_euler_angles():
    # angle = 2pi/3, phase = 0
    a1 = np.pi / 2
    b1 = np.pi / 3
    g1 = -np.pi / 2
    a2 = np.pi / 2
    b2 = np.pi / 3
    g2 = -np.pi / 2

    assert np.allclose(
        _add_two_euler_angles(a1, b1, g1, a2, b2, g2),
        np.array([np.pi / 2, 2 * np.pi / 3, -np.pi / 2]),
    )

    # angle = 0, phase = undefined
    # Although phase undefined here, we decide to return phase of zero (alpha = pi/2)
    a1 = 0
    b1 = np.pi / 3
    g1 = 0
    a2 = 0
    b2 = -np.pi / 3
    g2 = 0

    assert np.allclose(
        _add_two_euler_angles(a1, b1, g1, a2, b2, g2),
        np.array([np.pi / 2, 0, -np.pi / 2]),
    )

    # angle = 0, phase = undefined
    # Although phase undefined here, we decide to return phase of zero (alpha = pi/2)
    a1 = 0
    b1 = np.pi / 3
    g1 = 0
    a2 = np.pi
    b2 = np.pi / 3
    g2 = -np.pi

    assert np.allclose(
        _add_two_euler_angles(a1, b1, g1, a2, b2, g2),
        np.array([np.pi / 2, 0, -np.pi / 2]),
    )

    # angle = pi, phase = pi / 5
    # When beta = 180 degrees, phase must be analytically derived
    a1 = np.pi / 5
    b1 = np.pi / 2
    g1 = -np.pi / 5
    a2 = np.pi / 5
    b2 = np.pi / 2
    g2 = -np.pi / 5

    assert np.allclose(
        _add_two_euler_angles(a1, b1, g1, a2, b2, g2),
        np.array([np.pi / 5, np.pi, -np.pi / 5]),
    )

    # Out of transverse plane rotation
    a1 = np.pi / 3
    b1 = np.pi / 2
    g1 = -np.pi / 3
    a2 = 0
    b2 = np.pi / 3
    g2 = 0
    assert np.allclose(
        _add_two_euler_angles(a1, b1, g1, a2, b2, g2),
        np.array([1.28976143, 2.01862872, -0.06440383]),
    )


def test_euler_angles_to_angle_phase():
    a = 0
    b = 1
    g = 1
    with pytest.raises(ValueError, match=".*Unable to convert.*"):
        _euler_angles_to_angle_phase(a, b, g)

    a = np.pi
    b = 1
    g = -np.pi
    assert np.allclose(_euler_angles_to_angle_phase(a, b, g), [1, -np.pi / 2])


def test_combine_euler_angles():
    angle_set_1 = [
        (np.pi / 2, np.pi / 4, -np.pi / 2),
        (-np.pi / 3, np.pi / 6, np.pi / 3),
        (np.pi / 2, np.pi / 2, np.pi / 2),
    ]

    # Check values tested against Mathematica
    alpha_should_be = 1.803605
    beta_should_be = 1.2346306
    gamma_should_be = 1.632287
    should_be = (alpha_should_be, beta_should_be, gamma_should_be)

    assert np.allclose(combine_euler_angles(angle_set_1), should_be)

    angle_set_2 = [
        (np.pi / 4, 0, np.pi / 4),
        (np.pi / 6, np.pi / 3, np.pi),
        (np.pi / 2, -np.pi / 4, -np.pi / 4),
        (0.25, 1, 0.5),
        (np.pi / 2, np.pi / 2, np.pi / 2),
    ]

    alpha_should_be = 2.601394
    beta_should_be = 1.0306734
    gamma_should_be = 2.029843
    should_be = (alpha_should_be, beta_should_be, gamma_should_be)

    assert np.allclose(combine_euler_angles(angle_set_2), should_be)
