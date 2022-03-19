# -*- coding: utf-8 -*-
from mrsimulator import Coupling
from mrsimulator import Site
from mrsimulator import SpinSystem
from mrsimulator.spin_system.split_spinsystems import new_systems_needed_np

__author__ = "Alexis McCarthy"
__email__ = "mccarthy.677@osu.edu"


def setup_sites():
    A = Site(isotope="1H", isotropic_chemical_shift=0, name="a")
    B = Site(isotope="1H", isotropic_chemical_shift=2, name="b")
    C = Site(isotope="1H", isotropic_chemical_shift=4, name="c")
    D = Site(isotope="1H", isotropic_chemical_shift=6, name="d")
    E = Site(isotope="1H", isotropic_chemical_shift=8, name="e")
    F = Site(isotope="1H", isotropic_chemical_shift=10, name="f")
    return [A, B, C, D, E, F]


def setup_system_simplifiedSystem():
    """systems A-B-C, D-F, E"""
    sites = setup_sites()

    AB_couple = Coupling(site_index=[0, 1], isotropic_j=10, name="AB")
    BC_couple = Coupling(site_index=[1, 2], isotropic_j=10, name="BC")
    AC_couple = Coupling(site_index=[0, 2], isotropic_j=10, name="AC")
    DF_couple = Coupling(site_index=[3, 5], isotropic_j=30, name="DF")
    couplings = [AB_couple, BC_couple, AC_couple, DF_couple]
    sys = SpinSystem(sites=sites, couplings=couplings, abundance=10)

    A, B, C, D, E, F = sites
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
    """Systems A, B, Cc, D, E, F"""
    sites = setup_sites()
    sys = SpinSystem(sites=sites, abundance=25)

    simplified_sys = [SpinSystem(sites=[i], abundance=25) for i in sites]

    return sys, simplified_sys


def setup_somecoupled_system():
    """Systems A-C, B, D-E, F"""
    sites = setup_sites()

    AC_couple = Coupling(site_index=[0, 2], isotropic_j=20)
    DE_couple = Coupling(site_index=[3, 4], isotropic_j=40)
    sys = SpinSystem(sites=sites, couplings=[AC_couple, DE_couple], abundance=50)

    A, B, C, D, E, F = sites
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
    """Systems A-C-E, B, D, F"""
    sites = setup_sites()

    AC_couple = Coupling(site_index=[0, 2], isotropic_j=20)
    CE_couple = Coupling(site_index=[2, 4], isotropic_j=40)
    sys = SpinSystem(sites=sites, couplings=[AC_couple, CE_couple], abundance=50)

    A, B, C, D, E, F = sites
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


def setup_partially_coupled_system2():
    """Systems A-C-E-F, B-D"""
    sites = setup_sites()

    AC_couple = Coupling(site_index=[0, 2], isotropic_j=20)
    BD_couple = Coupling(site_index=[1, 3], isotropic_j=120)
    CE_couple = Coupling(site_index=[2, 4], isotropic_j=40)
    EF_couple = Coupling(site_index=[4, 5], isotropic_j=10)
    sys = SpinSystem(
        sites=sites,
        couplings=[AC_couple, BD_couple, CE_couple, EF_couple],
        abundance=50,
    )

    A, B, C, D, E, F = sites
    simplified_sys = [
        SpinSystem(
            sites=[A, C, E, F],
            couplings=[
                Coupling(site_index=[0, 1], isotropic_j=20),
                Coupling(site_index=[1, 2], isotropic_j=40),
                Coupling(site_index=[2, 3], isotropic_j=10),
            ],
            abundance=50,
        ),
        SpinSystem(
            sites=[B, D],
            couplings=[Coupling(site_index=[0, 1], isotropic_j=120)],
            abundance=50,
        ),
    ]
    return sys, simplified_sys


def setup_partially_coupled_system3():
    """Systems A-B-C-D-E-F"""
    sys, _ = setup_partially_coupled_system2()
    sys.couplings.insert(0, Coupling(site_index=[0, 1], isotropic_j=1020))
    return sys, [sys]


def generic_test(sys, simplified_sys, description):
    sys_simplify = sys.simplify()
    assert len(sys_simplify) == len(simplified_sys), description
    for system in sys_simplify:
        assert system in simplified_sys, description
        loc = simplified_sys.index(system)
        for site in system.sites:
            assert site in simplified_sys[loc].sites, description
        if system.couplings:
            for coupling in system.couplings:
                assert coupling in simplified_sys[loc].couplings, description


def test_builder():
    test_fn = [
        setup_system_simplifiedSystem,
        setup_uncoupled_system,
        setup_somecoupled_system,
        setup_partially_coupled_system,
        setup_partially_coupled_system2,
        setup_partially_coupled_system3,
    ]
    for fn in test_fn:
        sys, simplified_sys = fn()
        generic_test(sys, simplified_sys, str(fn))


def test_new_systems_needed_np():
    couplings1 = [[0, 1]]
    num_sites1 = 3
    couplings2 = [[0, 1], [2, 3]]
    num_sites2 = 4
    systems1 = new_systems_needed_np(couplings1, num_sites1)
    systems2 = new_systems_needed_np(couplings2, num_sites2)
    assert systems1 == [[0, 1], [2]]
    assert systems2 == [[0, 1], [2, 3]]
