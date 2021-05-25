# -*- coding: utf-8 -*-
import numpy as np
from mrsimulator import Coupling
from mrsimulator import Site
from mrsimulator import SpinSystem
from mrsimulator.spin_system.split_spinsystems import new_systems_needed_matrix
from mrsimulator.spin_system.split_spinsystems import new_systems_needed_nosets


def setup_sites():
    A = Site(isotope="1H", isotropic_chemical_shift=0, name="a")
    B = Site(isotope="1H", isotropic_chemical_shift=2, name="b")
    C = Site(isotope="1H", isotropic_chemical_shift=4, name="c")
    D = Site(isotope="1H", isotropic_chemical_shift=6, name="d")
    E = Site(isotope="1H", isotropic_chemical_shift=8, name="e")
    F = Site(isotope="1H", isotropic_chemical_shift=10, name="f")
    return [A, B, C, D, E, F]


def setup_system_simplifiedSystem():

    A, B, C, D, E, F = setup_sites()
    sites = [A, B, C, D, E, F]
    AB_couple = Coupling(site_index=[0, 1], isotropic_j=10, name="AB")
    BC_couple = Coupling(site_index=[1, 2], isotropic_j=10, name="BC")
    AC_couple = Coupling(site_index=[0, 2], isotropic_j=10, name="AC")
    DF_couple = Coupling(site_index=[3, 5], isotropic_j=30, name="DF")
    couplings = [AB_couple, BC_couple, AC_couple, DF_couple]

    sys = SpinSystem(sites=sites, couplings=couplings, abundance=10)

    simplified_sys = [
        SpinSystem(sites=[E], abundance=10),
        SpinSystem(
            sites=[D, F],
            couplings=[Coupling(site_index=[0, 1], isotropic_j=30, name="DF")],
            abundance=10,
        ),
        SpinSystem(
            sites=[A, B, C], couplings=[AB_couple, AC_couple, BC_couple], abundance=10
        ),
    ]
    return sys, simplified_sys


def setup_uncoupled_system():
    A, B, C, D, E, F = setup_sites()
    sites = [A, B, C, D, E, F]
    sys = SpinSystem(sites=sites, abundance=25)

    simplified_sys = [SpinSystem(sites=[i], abundance=25) for i in sites]

    return sys, simplified_sys


def setup_somecoupled_system():
    A, B, C, D, E, F = setup_sites()
    sites = [A, B, C, D, E, F]

    AC_couple = Coupling(site_index=[0, 2], isotropic_j=20)
    DE_couple = Coupling(site_index=[3, 4], isotropic_j=40)

    sys = SpinSystem(sites=sites, couplings=[AC_couple, DE_couple], abundance=50)

    simplified_sys = [
        SpinSystem(sites=[B], abundance=50),
        SpinSystem(sites=[F], abundance=50),
        SpinSystem(
            sites=[A, C],
            couplings=[Coupling(site_index=[0, 1], isotropic_j=20)],
            abundance=50,
        ),
        SpinSystem(
            sites=[D, E],
            couplings=[Coupling(site_index=[0, 1], isotropic_j=40)],
            abundance=50,
        ),
    ]
    return sys, simplified_sys


def setup_partially_coupled_system():
    A, B, C, D, E, F = setup_sites()
    sites = [A, B, C, D, E, F]

    AC_couple = Coupling(site_index=[0, 2], isotropic_j=20)
    CE_couple = Coupling(site_index=[2, 4], isotropic_j=40)

    sys = SpinSystem(sites=sites, couplings=[AC_couple, CE_couple], abundance=50)

    simplified_sys = [
        SpinSystem(sites=[B], abundance=50),
        SpinSystem(sites=[D], abundance=50),
        SpinSystem(sites=[F], abundance=50),
        SpinSystem(
            sites=[A, C, E],
            couplings=[
                Coupling(site_index=[0, 1], isotropic_j=20),
                Coupling(site_index=[1, 2], isotropic_j=40),
            ],
            abundance=50,
        ),
    ]
    return sys, simplified_sys


def generic_test(sys, simplified_sys):
    sys_simplify = sys.simplify()
    assert len(sys_simplify) == len(simplified_sys)
    for system in sys_simplify:
        assert system in simplified_sys
        loc = simplified_sys.index(system)
        for site in system.sites:
            assert site in simplified_sys[loc].sites
        if system.couplings:
            for coupling in system.couplings:
                assert coupling in simplified_sys[loc].couplings
    # assert sys_simplify == simplified_sys


def test_simplify_1():
    sys, simplified_sys = setup_system_simplifiedSystem()
    generic_test(sys, simplified_sys)


def test_simplify_2():
    sys, simplified_sys = setup_uncoupled_system()
    generic_test(sys, simplified_sys)


def test_simplify_3():
    sys, simplified_sys = setup_somecoupled_system()
    generic_test(sys, simplified_sys)


def test_simplify_4():
    sys, simplified_sys = setup_partially_coupled_system()
    generic_test(sys, simplified_sys)


def test_new_systems_needed_matrix():
    test_matrix1 = np.array([[1, 1, 0], [1, 1, 0], [0, 0, 1]])
    test_matrix2 = np.array(
        [[1, 1, 0, 0], [1, 1, 5e-20, 0], [0, 0, 1, 1], [0, 0, 1, 1]]
    )
    sets1 = new_systems_needed_matrix(test_matrix1)
    sets2 = new_systems_needed_matrix(test_matrix2)
    assert sets1 == set({frozenset({0, 1}), frozenset({2})})
    assert sets2 == set({frozenset({0, 1}), frozenset({2, 3})})


def test_new_systems_needed_nosets():
    test_matrix1 = np.array([[1, 1, 0], [1, 1, 0], [0, 0, 1]])
    test_matrix2 = np.array(
        [[1, 1, 0, 0], [1, 1, 5e-20, 0], [0, 0, 1, 1], [0, 0, 1, 1]]
    )
    systems1 = new_systems_needed_nosets(test_matrix1)
    systems2 = new_systems_needed_nosets(test_matrix2)
    assert systems1 == [(0, 1), (2,)]
    assert systems2 == [(0, 1), (2, 3)]
