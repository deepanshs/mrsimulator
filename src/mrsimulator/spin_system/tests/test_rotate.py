"""Tests for rotating different objects by a set of euler angles"""
import numpy as np
from mrsimulator.spin_system import SpinSystem
from mrsimulator.spin_system.coupling import Coupling
from mrsimulator.spin_system.site import Site
from mrsimulator.spin_system.tensors import SymmetricTensor
from mrsimulator.utils.euler_angles import combine_euler_angles


__author__ = "Matthew Giammar"
__email__ = "giammar.7@osu.edu"


angle_set_1 = [
    (np.pi / 2, np.pi / 4, -np.pi / 2),
    (np.pi / 3, np.pi / 6, np.pi / 3),
    (np.pi, -np.pi / 2, np.pi),
]


angle_set_2 = [
    (np.pi / 4, 0, np.pi / 4),
    (0.25, 1, 0.5),
    (np.pi, np.pi / 3, -np.pi),
]


def test_tensor_rotate():
    # Tensors initialized with all zero (alpha, beta, gamma)
    tens = SymmetricTensor(zeta=5, eta=0.2)
    should_be = combine_euler_angles(angle_set_1)

    tens.rotate(euler_angles=angle_set_1)

    assert np.allclose((tens.alpha, tens.beta, tens.gamma), should_be)

    tens = SymmetricTensor(zeta=5, eta=0.2)
    should_be = combine_euler_angles(angle_set_2)

    tens.rotate(euler_angles=angle_set_2)

    assert np.allclose((tens.alpha, tens.beta, tens.gamma), should_be)

    # Tensor with non-zero starting Euler angles
    tens = SymmetricTensor(
        zeta=5, eta=0.2, alpha=np.pi / 2, beta=np.pi / 2, gamma=np.pi / 2
    )
    should_be = combine_euler_angles(
        [(tens.alpha, tens.beta, tens.gamma)] + angle_set_1
    )

    tens.rotate(euler_angles=angle_set_1)

    assert np.allclose((tens.alpha, tens.beta, tens.gamma), should_be)


def test_site_rotate():
    site = Site(
        isotope="17O",
        shielding_symmetric={"zeta": 5, "eta": 0.2},
        quadrupolar={"Cq": 20, "eta": 0.5, "alpha": 0.5, "beta": 0.5, "gamma": 0.5},
    )

    site.rotate(angle_set_1)
    site_shield = (
        site.shielding_symmetric.alpha,
        site.shielding_symmetric.beta,
        site.shielding_symmetric.gamma,
    )
    site_quad = (
        site.quadrupolar.alpha,
        site.quadrupolar.beta,
        site.quadrupolar.gamma,
    )
    shield_should_be = combine_euler_angles(angle_set_1)
    quad_should_be = combine_euler_angles([(0.5, 0.5, 0.5)] + angle_set_1)

    assert np.allclose(site_shield, shield_should_be)
    assert np.allclose(site_quad, quad_should_be)


def test_coupling_rotate():
    coup = Coupling(
        site_index=[0, 1],
        j_symmetric={"zeta": 5, "eta": 0.2, "alpha": 0.5, "beta": 0.5, "gamma": 0.5},
    )

    coup.rotate(angle_set_2)
    j_symm = (coup.j_symmetric.alpha, coup.j_symmetric.beta, coup.j_symmetric.gamma)
    should_be = combine_euler_angles([(0.5, 0.5, 0.5)] + angle_set_2)

    assert np.allclose(should_be, j_symm)

    coup = Coupling(
        site_index=[0, 1], dipolar={"D": 100, "alpha": 1, "beta": 1, "gamma": 1}
    )

    coup.rotate(angle_set_2)
    dipolar = (coup.dipolar.alpha, coup.dipolar.beta, coup.dipolar.gamma)
    should_be = combine_euler_angles([(1, 1, 1)] + angle_set_2)

    assert np.allclose(should_be, dipolar)


def test_spin_system_rotate():
    site1 = Site(
        isotope="17O",
        shielding_symmetric={"zeta": 5, "eta": 0.2},
        quadrupolar={"Cq": 20, "eta": 0.5, "alpha": 0.5, "beta": 0.5, "gamma": 0.5},
    )
    site2 = Site(isotope="1H")

    coup = Coupling(
        site_index=[0, 1],
        j_symmetric={"zeta": 5, "eta": 0.2, "alpha": 1, "beta": 1, "gamma": 1},
    )

    sys = SpinSystem(sites=[site1, site2], couplings=[coup])

    sys.rotate(angle_set_1)

    shield_should_be = combine_euler_angles(angle_set_1)
    quad_should_be = combine_euler_angles([(0.5, 0.5, 0.5)] + angle_set_1)
    j_should_be = combine_euler_angles([(1, 1, 1)] + angle_set_1)
    site_shield = (
        sys.sites[0].shielding_symmetric.alpha,
        sys.sites[0].shielding_symmetric.beta,
        sys.sites[0].shielding_symmetric.gamma,
    )
    site_quad = (
        sys.sites[0].quadrupolar.alpha,
        sys.sites[0].quadrupolar.beta,
        sys.sites[0].quadrupolar.gamma,
    )
    coup_j = (
        sys.couplings[0].j_symmetric.alpha,
        sys.couplings[0].j_symmetric.beta,
        sys.couplings[0].j_symmetric.gamma,
    )

    assert np.allclose(shield_should_be, site_shield)
    assert np.allclose(quad_should_be, site_quad)
    assert np.allclose(j_should_be, coup_j)
