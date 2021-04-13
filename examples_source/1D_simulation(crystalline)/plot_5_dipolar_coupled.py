#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
Coupled spin-1/2 (Static dipolar spectrum)
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

13C-1H static dipolar coupling simulation.
"""
# %%
import matplotlib.pyplot as plt
from mrsimulator import Simulator, SpinSystem
from mrsimulator.methods import BlochDecaySpectrum
import mrsimulator.signal_processing as sp
import mrsimulator.signal_processing.apodization as apo

# sphinx_gallery_thumbnail_number = 1

# %%
# **Spin Systems**
#
# Create a 13C-1H coupled spin system.
spin_system = SpinSystem(
    sites=[
        {"isotope": "13C", "isotropic_chemical_shift": 0.0},
        {"isotope": "1H", "isotropic_chemical_shift": 0.0},
    ],
    couplings=[{"site_index": [0, 1], "dipolar": {"D": -2e4}}],
)
# %%
# **Methods**
#
# Create a BlochDecaySpectrum method.
method = BlochDecaySpectrum(
    channels=["13C"],
    magnetic_flux_density=9.4,  # in T
    spectral_dimensions=[{"count": 2048, "spectral_width": 8.0e4}],
)

# %%
# **Simulator**
#
# Create the Simulator object and add the method and the spin system object.
sim = Simulator()
sim.spin_systems += [spin_system]  # add the spin system.
sim.methods += [method]  # add the method.
sim.run()

# %%
# **Post-Simulation Processing**
#
# Add post-simulation signal processing.
processor = sp.SignalProcessor(
    operations=[
        sp.IFFT(),
        apo.Exponential(FWHM="500 Hz"),
        sp.FFT(),
    ]
)
processed_data = processor.apply_operations(data=sim.methods[0].simulation)

# %%
# **Plot**
#
plt.figure(figsize=[4.5, 3.0])
ax = plt.subplot(projection="csdm")
ax.plot(processed_data.real, color="black", linewidth=1)
ax.invert_xaxis()
plt.tight_layout()
plt.show()
