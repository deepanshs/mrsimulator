# -*- coding: utf-8 -*-
from mrsimulator import Coupling
from mrsimulator import Site
from mrsimulator import SpinSystem


def setup_system_simplifiedSystem():
    A = Site(isotope="1H", isotropic_chemical_shift=0, name="a")
    B = Site(isotope="1H", isotropic_chemical_shift=2, name="b")
    C = Site(isotope="1H", isotropic_chemical_shift=4, name="c")
    D = Site(isotope="1H", isotropic_chemical_shift=6, name="d")
    E = Site(isotope="1H", isotropic_chemical_shift=8, name="e")
    F = Site(isotope="1H", isotropic_chemical_shift=10, name="f")
    sites = [A, B, C, D, E, F]

    AB_couple = Coupling(site_index=[0, 1], isotropic_j=10, name="AB")
    BC_couple = Coupling(site_index=[1, 2], isotropic_j=10, name="BC")
    AC_couple = Coupling(site_index=[0, 2], isotropic_j=10, name="AC")
    DF_couple = Coupling(site_index=[3, 5], isotropic_j=30, name="DF")
    couplings = [AB_couple, BC_couple, AC_couple, DF_couple]

    sys = SpinSystem(sites=sites, couplings=couplings, abundance=10)

    simplified_sys = [
        SpinSystem(
            sites=[A, B, C], couplings=[AB_couple, BC_couple, AC_couple], abundance=10
        ),
        SpinSystem(
            sites=[D, F],
            couplings=[Coupling(site_index=[0, 1], isotropic_j=10, name="DF")],
            abundance=10,
        ),
        SpinSystem(sites=[E], abundance=10),
    ]
    return sys, simplified_sys


def test_simplify_1():
    sys, simplified_sys = setup_system_simplifiedSystem()
    sys_simplify = sys.simplify()

    assert len(sys_simplify) == len(simplified_sys)
    assert sys_simplify == simplified_sys


# test_sys2 = SpinSystem(sites=sites, abundance=10)


# irr_sys1 = SpinSystem(
#     sites=[A, B, C], couplings=[AB_couple, BC_couple, AC_couple], abundance=10
# )
# irr_sys2 = SpinSystem(
#     sites=[D, F],
#     couplings=[Coupling(site_index=[0, 1], isotropic_j=10, name="DF")],
#     abundance=10,
# )
# irr_sys3 = SpinSystem(sites=[E], abundance=10)

# irr_sys = [irr_sys3, irr_sys2, irr_sys1]

# new_sys = split_spin_system(test_sys1)

# assert len(new_sys) == 3

# method_H = BlochDecaySpectrum(
#     channels=["1H"],
#     magnetic_flux_density=9.4,  # T
#     spectral_dimensions=[
#         {"count": 16000, "spectral_width": 5e3, "reference_offset": 2e3}  # in Hz
#     ],
# )


# sim = Simulator()
# sim.spin_systems = new_sys  # [test_sys1]  # irr_sys  #
# sim.methods = [method_H]
# sim.run()
# H_data = sim.methods[0].simulation

# normalized_H_data = H_data / (H_data.sum())
# plt.figure(figsize=(6, 4))  # set the figure size
# ax = plt.subplot(projection="csdm")
# ax.plot(
#     normalized_H_data.real,
#     color="black",
#     linewidth=0.5,
# )


# ax.invert_xaxis()
# plt.tight_layout()
# plt.show()
