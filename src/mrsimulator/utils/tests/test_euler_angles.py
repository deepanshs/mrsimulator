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


def test_add_two_euler_angles_sum_to_pi_rotation():
    """Two rotations of the same phase which sum to a pi rotation combine to the pi
    rotation of that phase, for a phase in any quadrant."""
    from mrsimulator.utils.euler_angles import _angle_phase_to_euler_angles
    from scipy.spatial.transform import Rotation

    for deg in [0, 30, 45, 90, 135, 180, 225, 270, 315, 350]:
        phase = np.deg2rad(deg)
        for angle in [np.pi / 2, np.pi / 3]:
            first = _angle_phase_to_euler_angles(angle, phase)
            second = _angle_phase_to_euler_angles(np.pi - angle, phase)
            out = _add_two_euler_angles(*first, *second)
            exact = Rotation.from_euler("zyz", first) * Rotation.from_euler(
                "zyz", second
            )
            res = Rotation.from_euler("zyz", out).as_matrix()
            assert np.allclose(res, exact.as_matrix(), atol=1e-12), (deg, angle)


def test_two_half_pi_mixing_events_equal_one_pi_mixing_event():
    """The pathway weight of two pi/2 mixing events of the same phase is the weight of
    one pi mixing event of that phase."""
    from mrsimulator import Method
    from mrsimulator import Site
    from mrsimulator import SpinSystem
    from mrsimulator.method import MixingEvent
    from mrsimulator.method import SpectralDimension
    from mrsimulator.method import SpectralEvent

    sys = SpinSystem(sites=[Site(isotope="1H")])

    def method(mix):
        events = [
            SpectralEvent(fraction=0.5, transition_queries=[{"ch1": {"P": [1]}}]),
            *mix,
            SpectralEvent(fraction=0.5, transition_queries=[{"ch1": {"P": [-1]}}]),
        ]
        return Method(
            channels=["1H"], spectral_dimensions=[SpectralDimension(events=events)]
        )

    for deg in [45, 135, 225, 315]:
        phase = np.deg2rad(deg)
        one = method([MixingEvent(ch1={"angle": np.pi, "phase": phase})])
        two = method([MixingEvent(ch1={"angle": np.pi / 2, "phase": phase})] * 2)
        w1 = one.get_transition_pathways(sys)[0].weight
        w2 = two.get_transition_pathways(sys)[0].weight
        assert np.allclose(w1, w2), deg
