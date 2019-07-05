from mrsimulator.methods import one_d_spectrum
from mrsimulator import Simulator
import matplotlib.pyplot as plt

# import numpy as np
from mrsimulator.python.simulator import simulator

from timeit import default_timer

# import time

# time.sleep(5)
isotopomers = [
    {
        "sites": [
            {
                "isotope_symbol": "1H",
                "isotropic_chemical_shift": "0 ppm",
                "shielding_symmetric": {"anisotropy": "15 kHz", "asymmetry": 0.25},
            },
            {
                "isotope_symbol": "1H",
                "isotropic_chemical_shift": "1 kHz",
                "shielding_symmetric": {"anisotropy": "10 kHz", "asymmetry": 0.5},
            },
            {
                "isotope_symbol": "1H",
                "isotropic_chemical_shift": "0 ppm",
                "shielding_symmetric": {"anisotropy": "12.3 kHz", "asymmetry": 0.21},
            },
            {
                "isotope_symbol": "1H",
                "isotropic_chemical_shift": "0 ppm",
                "shielding_symmetric": {"anisotropy": "15 kHz", "asymmetry": 0.25},
            },
            {
                "isotope_symbol": "1H",
                "isotropic_chemical_shift": "1 kHz",
                "shielding_symmetric": {"anisotropy": "10 kHz", "asymmetry": 0.5},
            },
            {
                "isotope_symbol": "1H",
                "isotropic_chemical_shift": "0 ppm",
                "shielding_symmetric": {"anisotropy": "12.3 kHz", "asymmetry": 0.21},
            },
            {
                "isotope_symbol": "1H",
                "isotropic_chemical_shift": "0 ppm",
                "shielding_symmetric": {"anisotropy": "15 kHz", "asymmetry": 0.25},
            },
            {
                "isotope_symbol": "1H",
                "isotropic_chemical_shift": "1 kHz",
                "shielding_symmetric": {"anisotropy": "10 kHz", "asymmetry": 0.5},
            },
            {
                "isotope_symbol": "1H",
                "isotropic_chemical_shift": "0 ppm",
                "shielding_symmetric": {"anisotropy": "12.3 kHz", "asymmetry": 0.21},
            },
            {
                "isotope_symbol": "1H",
                "isotropic_chemical_shift": "0 ppm",
                "shielding_symmetric": {"anisotropy": "15 kHz", "asymmetry": 0.25},
            },
            {
                "isotope_symbol": "1H",
                "isotropic_chemical_shift": "1 kHz",
                "shielding_symmetric": {"anisotropy": "10 kHz", "asymmetry": 0.5},
            },
            {
                "isotope_symbol": "1H",
                "isotropic_chemical_shift": "0 ppm",
                "shielding_symmetric": {"anisotropy": "12.3 kHz", "asymmetry": 0.21},
            },
        ],
        "abundance": "100 %",
    },
    {
        "sites": [
            {
                "isotope_symbol": "1H",
                "isotropic_chemical_shift": "10 ppm",
                "shielding_symmetric": {"anisotropy": "63.89 ppm", "asymmetry": 0.5},
            }
        ],
        "abundance": "100 %",
    },
    {
        "sites": [
            {
                "isotope_symbol": "1H",
                "isotropic_chemical_shift": "0 ppm",
                "shielding_symmetric": {"anisotropy": "43.89 ppm", "asymmetry": 0.21},
            }
        ],
        "abundance": "100 %",
    },
]
spectrum = {
    "direct_dimension": {
        "magnetic_flux_density": "9.4 T",
        "rotor_frequency": "0 kHz",
        "rotor_angle": "54.735 deg",
        "number_of_points": 2048,
        "spectral_width": "125 kHz",
        "reference_offset": "0 kHz",
        "nucleus": "1H",
    }
}
s = Simulator(isotopomers)
s.spectrum = spectrum

# n = 100
# start = default_timer()
# [
#     s.run(
#         one_d_spectrum,
#         verbose=0,
#         geodesic_polyhedron_frequency=90,
#         number_of_sidebands=128,
#     ) for _ in range(n)]
# print('python side time', (default_timer() - start) / float(n), ' s')
n = 1
start = default_timer()
f, a = s.run(
    one_d_spectrum,
    verbose=11,
    geodesic_polyhedron_frequency=48,
    number_of_sidebands=128,
)
print("python side time", (default_timer() - start) / float(n), " s")
# plt.plot(f, a)
# plt.show()

# print(s._spectrum_c)
# print(s._isotopomers_c)
n = 1

# start = default_timer()
# [simulator(s._spectrum_c, s._isotopomers_c) for _ in range(n)]
# print((default_timer() - start) / float(n), ' s')

start = default_timer()
freq, spec = simulator(s._spectrum_c, s._isotopomers_c, nt=48)
print("python side time", (default_timer() - start) / float(n), " s")
# assert np.allclose(spec, a)
fig, ax = plt.subplots(2, 1)
ax[0].plot(freq, spec / spec.max(), "r", label="py")
ax[0].plot(f, a / a.max(), "b", label="c")
ax[1].plot(f, (spec / spec.max()) - (a / a.max()))
plt.legend()
plt.show()
