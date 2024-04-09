import os
import timeit

import mrsimulator.tests.tests as clib
import numpy as np
from mrsimulator import __version__
from mrsimulator import Simulator
from mrsimulator.method import Method
from mrsimulator.method.lib import BlochDecayCentralTransitionSpectrum
from mrsimulator.method.lib import BlochDecaySpectrum
from mrsimulator.utils.collection import single_site_system_generator

# import platform
# os_system = platform.system()


__author__ = ["Deepansh Srivastava", "Matthew D. Giammar"]
__email__ = ["srivastava.89@osu.edu", "giammar.7@buckeyemail.osu.edu"]


def generate_spin_half_spin_system(n=1000):
    iso = np.random.normal(loc=0.0, scale=10.0, size=n)
    zeta = np.random.normal(loc=50.0, scale=15.12, size=n)
    eta = np.random.normal(loc=0.25, scale=0.012, size=n)
    return single_site_system_generator(
        isotope="29Si",
        isotropic_chemical_shift=iso,
        shielding_symmetric={"zeta": zeta, "eta": eta},
    )


def generate_spin_half_int_quad_spin_system(n=1000):
    iso = np.random.normal(loc=0.0, scale=10.0, size=n)
    Cq = np.random.normal(loc=5.0e6, scale=1e5, size=n)
    eta = np.random.normal(loc=0.1, scale=0.012, size=n)
    return single_site_system_generator(
        isotope="17O",
        isotropic_chemical_shift=iso,
        quadrupolar={"Cq": Cq, "eta": eta},
    )


def generate_spin_half_int_csa_quad_spin_system(n=1000):
    iso = np.random.normal(loc=0.0, scale=10.0, size=n)
    zeta = np.random.normal(loc=150, scale=20, size=n)
    eta_z = np.random.normal(loc=0.3, scale=0.012, size=n)
    Cq = np.random.normal(loc=5.0e6, scale=5e5, size=n)
    eta_q = np.random.normal(loc=0.6, scale=0.02, size=n)
    beta = np.random.normal(loc=2.12, scale=0.1, size=n)
    return single_site_system_generator(
        isotope="17O",
        isotropic_chemical_shift=iso,
        shielding_symmetric={"zeta": zeta, "eta": eta_z},
        quadrupolar={"Cq": Cq, "eta": eta_q, "beta": beta},
    )


def generate_simulator(
    spin_systems, method, integration_volume="octant", number_of_sidebands=64
):
    sim = Simulator()
    sim.spin_systems = spin_systems
    sim.methods = [method]
    sim.config.integration_volume = integration_volume
    sim.config.number_of_sidebands = number_of_sidebands
    return sim


def execute(sim, n_jobs=1):
    sim.run(n_jobs=n_jobs)


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


def quad_static_2d_method():
    tq = [{"ch1": {"P": [-1], "D": [0]}}]
    to_rad = 3.14159 / 180
    return Method(
        channels=["17O"],
        magnetic_flux_density=4.2,  # in T
        spectral_dimensions=[
            {
                "count": 256,
                "spectral_width": 4e4,  # in Hz
                "reference_offset": -1e4,  # in Hz
                "events": [{"rotor_angle": 70.12 * to_rad, "transition_queries": tq}],
            },
            {
                "count": 512,
                "spectral_width": 5e4,  # in Hz
                "reference_offset": -5e3,  # in Hz
                "events": [{"rotor_angle": 54.74 * to_rad, "transition_queries": tq}],
            },
        ],
    )


def time_string(time):
    count = 0
    if time < 0.003:  # i.e 4 ms
        color = "\033[92m"
    if time >= 0.003 and time < 0.005:
        color = "\033[93m"
    if time >= 0.005:
        color = "\033[91m"

    # color = '' if os_system == 'Windows' else color

    while time < 1:
        time *= 1000  # ms
        count += 1
    if count == 0:
        return f"{time:.3f}  s", color
    if count == 1:
        return f"{time:.3f} ms", color
    if count == 2:
        return f"{time:.3f} µs", color


def terminal_start_setup():
    size = os.get_terminal_size().columns
    delimit = "-"
    print(f"{delimit:-<{size}}")
    left_align = "Computation method"
    right_align = "Average time"
    print(f"{left_align:<{size-15}}{right_align:>15}")
    print(f"{delimit:-<{size}}")


def terminal_end_setup(t, n, description):
    t /= n  # s
    t, color = time_string(t)
    end = "\033[0m"
    # end = '' if os_system == 'Windows' else "\033[0m"
    size = os.get_terminal_size().columns - 15
    print(f"{end}{description:.<{size}}{color}{t:>15}{end}")


def spectrum_blocks(n, level, n_jobs=1):
    description = [
        "Static CSA only spectrum",
        "Static quadrupolar (1st + 2nd order) only spectrum",
        "MAS CSA only sidebands spectrum",
        "VAS CSA only spectrum",
        "MAS quadrupolar (1st + 2nd order) only spectrum",
        "CSA-Quad 2D SAS spectrum (infinite speed spectrum)",
    ]
    print(f"\nLevel {level} results.")
    print("Average computation time for simulation one single-site spin system.")
    print(
        f"Reported value is the time per simulation averaged over {n} simulations "
        "generated for random tensor parameters."
    )
    terminal_start_setup()
    for i, des in enumerate(description):
        t = timeit.timeit(
            f"execute(sim[{i}], {n_jobs})",
            setup=f"sim=spectrum_simulation_benchmark({n})",
            globals=globals(),
            number=1,
        )
        terminal_end_setup(t, n, des)


def spectrum_simulation_benchmark(n):
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

    # quad csd 2d
    spin_sys = generate_spin_half_int_csa_quad_spin_system(n)

    method = quad_static_2d_method()
    sim_csd_quad_mas_2d = generate_simulator(
        spin_sys, method, integration_volume="hemisphere", number_of_sidebands=1
    )
    return (
        sim_csa_static,
        sim_quad_static,
        sim_csa_mas,
        sim_csa_vas,
        sim_quad_mas,
        sim_csd_quad_mas_2d,
    )


def interpolation_blocks(n, level):
    description = [
        "One-dimension interpolation (512 grid)",
        "Two-dimensional interpolation (512 x 512 grid)",
    ]
    print(f"\nLevel {level} results.")
    print("Average computation time for rendering triangles on an nD-grid.")
    print(f"Reported value is the time per render averaged over {10*n} triangles.")
    terminal_start_setup()
    for i, des in enumerate(description):
        t = timeit.timeit(
            f"interpolation_execute(tasks[{i}])",
            setup=f"tasks=interpolation_benchmark({10*n})",
            globals=globals(),
            number=1,
        )
        terminal_end_setup(t, 10 * n, des)


def interpolation_execute(args):
    fn, vertexes, amp = args
    if vertexes.ndim == 2:
        _ = [fn(list_, amp) for list_ in vertexes]
    if vertexes.ndim == 3:
        _ = [fn(list_[0], list_[1], amp) for list_ in vertexes]


def interpolation_benchmark(n):
    # 1D interpolation
    vertexes1 = (np.random.rand(n, 3) * 768) - 128
    amp1 = np.zeros(512)

    # 2D interpolation
    vertexes2 = np.random.rand(n, 2, 3) * 512
    amp2 = np.zeros((512, 512))

    return [
        [clib.triangle_interpolation1D, vertexes1, amp1],
        [clib.triangle_interpolation2D, vertexes2, amp2],
    ]


class Benchmark:
    @staticmethod
    def prep():
        print(f"Benchmarking using mrsimulator version {__version__}")

    @staticmethod
    def l0(n_jobs, interpolation, simulation):
        setup(10, 0, n_jobs, interpolation, simulation)

    @staticmethod
    def l1(n_jobs, interpolation, simulation):
        setup(2000, 1, n_jobs, interpolation, simulation)

    @staticmethod
    def l2(n_jobs, interpolation, simulation):
        setup(10000, 2, n_jobs, interpolation, simulation)


def setup(n, level, n_jobs, interpolation, simulation):
    if simulation:
        spectrum_blocks(n, level, n_jobs)
    if interpolation:
        interpolation_blocks(n, level)
