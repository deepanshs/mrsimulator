from mrsimulator.methods import one_d_spectrum
from mrsimulator import Simulator
import matplotlib.pyplot as plt
from mrsimulator.sandbox import _one_d_simulator

import numpy as np
from mrsimulator.python.simulator import simulator

from timeit import default_timer

import time

# time.sleep(5)
isotopomers = [
    # {
    #     "sites": [
    #         {
    #             "isotope_symbol": "1H",
    #             "isotropic_chemical_shift": "0 ppm",
    #             "shielding_symmetric": {"anisotropy": "15 kHz", "asymmetry": 0.25},
    #         },
    #         {
    #             "isotope_symbol": "1H",
    #             "isotropic_chemical_shift": "1 kHz",
    #             "shielding_symmetric": {"anisotropy": "10 kHz", "asymmetry": 0.5},
    #         },
    #         {
    #             "isotope_symbol": "1H",
    #             "isotropic_chemical_shift": "0 ppm",
    #             "shielding_symmetric": {"anisotropy": "12.3 kHz", "asymmetry": 0.21},
    #         },
    #         {
    #             "isotope_symbol": "1H",
    #             "isotropic_chemical_shift": "0 ppm",
    #             "shielding_symmetric": {"anisotropy": "15 kHz", "asymmetry": 0.25},
    #         },
    #         {
    #             "isotope_symbol": "1H",
    #             "isotropic_chemical_shift": "1 kHz",
    #             "shielding_symmetric": {"anisotropy": "10 kHz", "asymmetry": 0.5},
    #         },
    #         {
    #             "isotope_symbol": "1H",
    #             "isotropic_chemical_shift": "0 ppm",
    #             "shielding_symmetric": {"anisotropy": "12.3 kHz", "asymmetry": 0.21},
    #         },
    #         {
    #             "isotope_symbol": "1H",
    #             "isotropic_chemical_shift": "0 ppm",
    #             "shielding_symmetric": {"anisotropy": "15 kHz", "asymmetry": 0.25},
    #         },
    #         {
    #             "isotope_symbol": "1H",
    #             "isotropic_chemical_shift": "1 kHz",
    #             "shielding_symmetric": {"anisotropy": "10 kHz", "asymmetry": 0.5},
    #         },
    #         {
    #             "isotope_symbol": "1H",
    #             "isotropic_chemical_shift": "0 ppm",
    #             "shielding_symmetric": {"anisotropy": "12.3 kHz", "asymmetry": 0.21},
    #         },
    #         {
    #             "isotope_symbol": "1H",
    #             "isotropic_chemical_shift": "0 ppm",
    #             "shielding_symmetric": {"anisotropy": "15 kHz", "asymmetry": 0.25},
    #         },
    #         {
    #             "isotope_symbol": "1H",
    #             "isotropic_chemical_shift": "1 kHz",
    #             "shielding_symmetric": {"anisotropy": "10 kHz", "asymmetry": 0.5},
    #         },
    #         {
    #             "isotope_symbol": "1H",
    #             "isotropic_chemical_shift": "0 ppm",
    #             "shielding_symmetric": {"anisotropy": "12.3 kHz", "asymmetry": 0.21},
    #         },
    #     ],
    #     "abundance": "100 %",
    # },
    # {
    #     "sites": [
    #         {
    #             "isotope_symbol": "1H",
    #             "isotropic_chemical_shift": "10 ppm",
    #             "shielding_symmetric": {"anisotropy": "63.89 ppm", "asymmetry": 0.5},
    #         }
    #     ],
    #     "abundance": "100 %",
    # },
    {
        "sites": [
            {
                "isotope_symbol": "1H",
                "isotropic_chemical_shift": "0 ppm",
                "shielding_symmetric": {"anisotropy": "43.89 ppm", "asymmetry": 0.21},
            }
        ],
        "abundance": "100 %",
    }
]
spectrum = {
    "direct_dimension": {
        "magnetic_flux_density": "9.4 T",
        "rotor_frequency": "1e9 kHz",
        "rotor_angle": "0 deg",
        "number_of_points": 1024,
        "spectral_width": "125 kHz",
        "reference_offset": "0 kHz",
        "nucleus": "1H",
    }
}


# spectrum = {
#     "direct_dimension": {
#         "rotor_frequency": "0 kHz",
#         "number_of_points": 2048,
#         "spectral_width": "100 kHz",
#         "nucleus": "1H"
#     }
# }
# isotopomer = [
#     {
#         "sites": [
#             {
#                 "isotope_symbol": "1H",
#                 "isotropic_chemical_shift": "-120 Hz",
#                 "shielding_symmetric": {
#                     "anisotropy": "9.34 kHz",
#                     "asymmetry": 0.0
#                 }
#             }
#         ]
#     }
# ]

s = Simulator(isotopomers, spectrum)

# n = 1
# start = default_timer()
# [
#     simulator(s._spectrum_c, s._isotopomers_c, nt=90, number_of_sidebands=128)
#     for _ in range(n)
# ]
# print("python side time, py", (default_timer() - start) / float(n), " s")

# start = default_timer()
# [
#     s.run(
#         one_d_spectrum,
#         verbose=0,
#         geodesic_polyhedron_frequency=90,
#         number_of_sidebands=128,
#     )
#     for _ in range(n)
# ]
# print("python side time, c", (default_timer() - start) / float(n), " s")

# plt.plot(f, a)
# plt.show()

# print(s._spectrum_c)
# print(s._isotopomers_c)

# n = 1
# start = default_timer()
# f, a = simulator(s._spectrum_c, s._isotopomers_c, nt=90, number_of_sidebands=1)
# print("python side time, py", (default_timer() - start) / float(n), ' s')

# start = default_timer()
# freq, spec = s.run(
#     one_d_spectrum, geodesic_polyhedron_frequency=90, number_of_sidebands=1
# )
# print("python side time, c", (default_timer() - start) / float(n), " s")
# # assert np.allclose(spec, a)
# fig, ax = plt.subplots(2, 1)
# ax[0].plot(freq, spec / spec.max(), "r", label="c")
# ax[0].plot(f, a / a.max(), "b", label="py")
# ax[0].legend()
# ax[1].plot(f, (spec / spec.max()) - (a / a.max()))

# plt.show()

start = default_timer()
n = 1000
number_of_points = 256 * 2
f, a, time__ = _one_d_simulator(
    reference_offset=-128.0 * 2,
    increment=1,
    number_of_points=number_of_points,
    isotropic_chemical_shift=np.random.normal(0.0, 2.0, n),
    chemical_shift_anisotropy=np.random.normal(100.0, 60.0, n),
    chemical_shift_asymmetry=np.random.normal(0.5, 0.2, n),
    sample_rotation_frequency_in_Hz=0.0,
    rotor_angle=54.735,
)
# print('Execultion time in c ', time__)
print("Execution time in python", default_timer() - start)
plt.plot(f, a.reshape(n, number_of_points).sum(axis=0))
plt.show()

# start = default_timer()
# n = 1000
# number_of_points = 2048
# sw = 2500
# f, a, time__ = _one_d_simulator(
#     reference_offset=-sw / 2,
#     increment=sw / number_of_points,
#     number_of_points=number_of_points,
#     isotropic_chemical_shift=np.zeros(n),  # np.random.normal(0, 20, n),
#     chemical_shift_anisotropy=np.random.normal(500, 10, n),
#     chemical_shift_asymmetry=np.random.normal(0.75, 0.1, n),
#     # quadrupolar_coupling_constant=np.random.normal(5.0e6, 1e5, n),
#     # quadrupolar_asymmetry=np.random.normal(0.5, 0.1, n),
#     # larmor_frequency=67.4e6,
#     # spin_quantum_number=2.5,
#     sample_rotation_frequency=0.0,
#     rotor_angle=54.735,
#     number_of_sidebands=128,
# )
# print("Execution time in python", default_timer() - start)
# plt.plot(f, a.reshape(n, number_of_points).sum(axis=0))
# plt.show()
