# -*- coding: utf-8 -*-
import os
import timeit

import mrsimulator
import numpy as np
from mrsimulator import Simulator
from mrsimulator import Site
from mrsimulator import SpinSystem
from mrsimulator.methods import BlochDecayCentralTransitionSpectrum
from mrsimulator.methods import BlochDecaySpectrum

print(f"Benchmarking using mrsimulator version {mrsimulator.__version__}")
left_align = "Computation method"
right_align = "Average time"
size = os.get_terminal_size().columns - 15
print(f"{left_align:.<{size}}{right_align:>15}")


def generate_spin_half_spin_system(n=1000):
    iso = np.random.normal(loc=0.0, scale=10.0, size=n)
    zeta = np.random.normal(loc=50.0, scale=15.12, size=n)
    eta = np.random.normal(loc=0.25, scale=0.012, size=n)
    spin_systems = []
    for i, z, e in zip(iso, zeta, eta):
        site = Site(
            isotope="29Si",
            isotropic_chemical_shift=i,
            shielding_symmetric={"zeta": z, "eta": e},
        )
        spin_systems.append(SpinSystem(sites=[site]))
    return spin_systems


def generate_spin_half_int_quad_spin_system(n=1000):
    iso = np.random.normal(loc=0.0, scale=10.0, size=n)
    Cq = np.random.normal(loc=5.0e6, scale=1e5, size=n)
    eta = np.random.normal(loc=0.1, scale=0.012, size=n)
    spin_systems = []
    for i, c, e in zip(iso, Cq, eta):
        site = Site(
            isotope="17O", isotropic_chemical_shift=i, quadrupolar={"Cq": c, "eta": e}
        )
        spin_systems.append(SpinSystem(sites=[site]))
    return spin_systems


def generate_simulator(spin_systems, method):
    sim = Simulator()
    sim.spin_systems = spin_systems
    sim.methods = [method]
    return sim


def execute(sim):
    sim.run()


def CSA_static_method():
    return BlochDecaySpectrum(
        channels=["29Si"], spectral_dimensions=[dict(spectral_width=25000)]
    )


def CSA_MAS_method():
    return BlochDecaySpectrum(
        channels=["29Si"],
        rotor_frequency=5000,
        spectral_dimensions=[dict(spectral_width=25000)],
    )


def CSA_VAS_method():
    return BlochDecaySpectrum(
        channels=["29Si"],
        rotor_frequency=5000,
        rotor_angle=1.57079,
        spectral_dimensions=[dict(spectral_width=25000)],
    )


def quad_static_method():
    return BlochDecayCentralTransitionSpectrum(
        channels=["17O"], spectral_dimensions=[{"spectral_width": 50000}]
    )


def quad_MAS_method():
    return BlochDecayCentralTransitionSpectrum(
        channels=["17O"],
        rotor_frequency=14000,
        spectral_dimensions=[{"spectral_width": 50000}],
    )


def time_string(time):
    count = 0
    if time < 0.003:  # i.e 4 ms
        color = "\033[92m"
    if time >= 0.003 and time < 0.005:
        color = "\033[93m"
    if time >= 0.005:
        color = "\033[91m"

    while time < 1:
        time *= 1000  # ms
        count += 1
    if count == 0:
        return f"{time:.3f}  s", color
    if count == 1:
        return f"{time:.3f} ms", color
    if count == 2:
        return f"{time:.3f} µs", color


def blocks(n_blocks, fn_string, n):
    description = [
        "Static CSA spectrum",
        "Static quadrupolar (1st + 2nd order) spectrum",
        "CSA MAS sidebands spectrum",
        "CSA VAS spectrum",
        "MAS quadrupolar (1st + 2nd order) spectrum",
    ]
    for i in range(n_blocks):
        t = timeit.timeit(
            f"execute(sim[{i}])",
            setup=f"sim={fn_string}({n})",
            globals=globals(),
            number=1,
        )
        t /= n  # s
        t, color = time_string(t)
        end = "\033[0m"
        size = os.get_terminal_size().columns - 15
        left_align = f"Average time for {description[i]}"
        print(f"{end}{left_align:.<{size}}{color}{t:>15}{end}")


def level_n_CSA_static(n):
    # spin 1/2
    spin_sys = generate_spin_half_spin_system(n)

    method = CSA_static_method()  # static
    sim_csa_static = generate_simulator(spin_sys, method)

    method = CSA_MAS_method()  # mas
    sim_csa_mas = generate_simulator(spin_sys, method)

    method = CSA_VAS_method()  # vas
    sim_csa_vas = generate_simulator(spin_sys, method)

    # quad
    spin_sys = generate_spin_half_int_quad_spin_system(n)

    method = quad_static_method()  # static
    sim_quad_static = generate_simulator(spin_sys, method)

    method = quad_MAS_method()  # static
    sim_quad_mas = generate_simulator(spin_sys, method)

    return sim_csa_static, sim_quad_static, sim_csa_mas, sim_csa_vas, sim_quad_mas


class Benchmark:
    @staticmethod
    def l0():
        blocks(5, "level_n_CSA_static", 10)

    @staticmethod
    def l1():
        blocks(5, "level_n_CSA_static", 2000)

    @staticmethod
    def l2():
        blocks(5, "level_n_CSA_static", 10000)
