# -*- coding: utf-8 -*-
import matplotlib.pyplot as plt
from mrsimulator import Simulator
from mrsimulator import SpinSystem
from mrsimulator.methods import BlochDecaySpectrum

spin_systems = [
    {
        "sites": [
            {
                "isotope": "1H",
                "isotropic_chemical_shift": "0 ppm",
                "shielding_symmetric": {"zeta": "13.89 ppm", "eta": 0.25},
            }
        ],
        "abundance": "100 %",
    }
]

method1 = {
    "channels": ["1H"],
    "magnetic_flux_density": "9.4 T",
    "rotor_frequency": "0 kHz",
    "rotor_angle": "54.735 deg",
    "spectral_dimensions": [
        {"count": 2048, "spectral_width": "25 kHz", "reference_offset": "0 Hz"}
    ],
}
method2 = {
    "channels": ["1H"],
    "magnetic_flux_density": "9.4 T",
    "rotor_frequency": "1 kHz",
    "rotor_angle": "54.735 deg",
    "spectral_dimensions": [
        {"count": 2048, "spectral_width": "25 kHz", "reference_offset": "0 Hz"}
    ],
}


sim = Simulator()
sim.spin_systems = [SpinSystem.parse_dict_with_units(item) for item in spin_systems]
sim.methods = [
    BlochDecaySpectrum.parse_dict_with_units(method1),
    BlochDecaySpectrum.parse_dict_with_units(method2),
]

sim.run()

freq1, amp1 = sim.methods[0].simulation.to_list()
freq2, amp2 = sim.methods[1].simulation.to_list()

fig, ax = plt.subplots(1, 2, figsize=(8, 3.5))
ax[0].plot(freq1, amp1, linewidth=1.0, color="k")
ax[0].set_xlabel(f"frequency ratio / {freq2.unit}")
ax[0].grid(color="gray", linestyle="--", linewidth=0.5, alpha=0.5)
ax[0].set_title("Static")
ax[1].plot(freq2, amp2, linewidth=1.0, color="k")
ax[1].set_xlabel(f"frequency ratio / {freq2.unit}")
ax[1].grid(color="gray", linestyle="--", linewidth=0.5, alpha=0.5)
ax[1].set_title("MAS")
plt.tight_layout(h_pad=0.1)
plt.show()
