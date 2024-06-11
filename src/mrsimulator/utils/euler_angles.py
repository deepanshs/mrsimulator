import numpy as np
from scipy.spatial.transform import Rotation


__author__ = "Matthew Giammar"
__email__ = "giammar.7@osu.edu"


TWO_PI = np.pi * 2


def wrap_between_pi(a: float) -> float:
    """Wraps the provided angle between (-pi and pi]"""
    a %= TWO_PI
    a -= np.sign(a) * TWO_PI if abs(a) > np.pi else 0
    return a


def combine_euler_angles(euler_angles: list) -> tuple:
    """Takes a list of Euler angles as (alpha, beta, gamma) tuples and converts all
    sequential rotations to a single set of Euler angles. Uses the ZYZ convention.

    Arguments:
        (list) euler_angles: List of Euler angles as (alpha, beta, gamma) tuples

    Returns:
        (tuple) Final set of Euler angles after applying all rotations.
    """
    alpha, beta, gamma = euler_angles[0]

    # Iterate over all Euler angles adding them to variables alpha, beta, gamma
    for elr_ang in euler_angles[1:]:
        alpha, beta, gamma = _add_two_euler_angles(alpha, beta, gamma, *elr_ang)

    return alpha, beta, gamma


def _angle_phase_to_euler_angles(angle: float, phase: float) -> tuple:
    """Convert the angle and phase of a mixing query to a set of Euler angles in
    the ZYZ convention. The returned angles will be constrained between (-pi, pi]

    Args:
        angle (float): Angle of mixing query between [0, 2pi)
        phase (float): Phase of mixing query between [0, 2pi)

    Returns:
        alpha, beta, gamma: Euler angles of the mixing query
    """
    # Wrap angle and phase between -pi and pi
    angle, phase = wrap_between_pi(angle), wrap_between_pi(phase)
    alpha = (np.pi / 2) - phase
    return wrap_between_pi(alpha), wrap_between_pi(angle), wrap_between_pi(-alpha)


def _euler_angles_to_angle_phase(alpha: float, beta: float, gamma: float) -> tuple:
    """Takes a set of Euler angles in the ZYZ convention and converts them to a
    mixing angle and phase. Provided alpha and gamma should be opposite of each other,
    otherwise, a ValueError is raised since the rotation vector does not lie in the XY
    plane.

    Args:
        alpha (float): First Euler angle in set
        beta (float): Second Euler angle in set
        gamma (float): Third Euler angle in set

    Returns:
        angle, phase: Angle and phase of the equivalent mixing query

    Raises:
        ValueError: Raised if alpha and gamma are not opposite of each other
    """
    if not np.isclose(alpha, -gamma):
        raise ValueError(
            "Unable to convert the provided Euler angles to an angle and phase"
        )

    phase = wrap_between_pi(gamma + np.pi / 2)

    return beta, phase


def _add_two_euler_angles(
    a1: float, b1: float, g1: float, a2: float, b2: float, g2: float
) -> tuple:
    """Adds two sets of euler angles -- (a1, b1, g1) and (a2, b2, g2) -- together.
    Also, check for edge cases where a gimbal lock would occur.

    If the result is the identity matrix, then beta = 0 and alpha, gamma are unbounded.
    As an arbitrary choice, alpha of pi/2 and gamma of -pi/2 are chosen.

    If the two phases are the same (i.e. a1 == a2 and g1 == g2) and b1 + b2 = pi, then
    scipy is unable to uniquely determine the last angle. In this case, the top left
    element of the rotation matrix equals -cos(2*phase) and alpha and gamma are computed
    from there.

    Arguments:
        a1 (float): alpha for first Euler angle
        b1 (float): beta for first Euler angle
        g1 (float): gamma for first Euler angle
        a2 (float): alpha for second Euler angle
        b2 (float): beta for second Euler angle
        g2 (float): gamma for second Euler angle

    Returns:
        alpha, beta, gamma: The resulting Euler angle
    """
    rot_1 = Rotation.from_euler("zyz", [a1, b1, g1])
    rot_2 = Rotation.from_euler("zyz", [a2, b2, g2])

    result = rot_1 * rot_2
    result_mat = result.as_matrix()

    # Check for identity matrix by comparing the diagonal to an array of ones
    if np.allclose(result_mat.diagonal(), [1.0, 1.0, 1.0]):
        return np.asarray([np.pi / 2, 0, -np.pi / 2])

    # Check if beta is 180 degrees and phase the same
    if np.isclose(result_mat[2][2], -1.0) and np.allclose([a1, g1], [a2, g2]):
        gamma = np.arccos(result_mat[0][0]) - np.pi
        gamma /= 2
        return np.asarray([-gamma, np.pi, gamma])

    return result.as_euler("zyz")
