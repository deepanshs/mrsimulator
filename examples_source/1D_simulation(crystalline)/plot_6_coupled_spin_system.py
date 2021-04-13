#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
Coupled spins 5/2-9/2 (Quad + J-coupling)
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

27Al-93Nb spin system spectrum.
"""
# %%
import matplotlib.pyplot as plt
from mrsimulator import Simulator, SpinSystem
from mrsimulator.methods import BlochDecayCTSpectrum
import mrsimulator.signal_processing as sp
import mrsimulator.signal_processing.apodization as apo

# sphinx_gallery_thumbnail_number = 1

# %%
# **Spin System**
#
# Create a 27Al-93Nb coupled spin system.
spin_system = SpinSystem(
    sites=[
        {
            "isotope": "27Al",
            "isotropic_chemical_shift": 0.0,  # in ppm
            "quadrupolar": {"Cq": 5.0e6, "eta": 0.0},  # Cq is in Hz
        },
        {
            "isotope": "93Nb",
            "isotropic_chemical_shift": 0.0,  # in ppm
        },
    ],
    couplings=[{"site_index": [0, 1], "isotropic_j": 200.0}],  # j-coupling in Hz
)

# %%
# **Method**
#
# Create a central transition selective Bloch decay spectrum method.
method = BlochDecayCTSpectrum(
    channels=["27Al"],
    magnetic_flux_density=9.4,  # in T
    rotor_frequency=5e3,  # in Hz
    spectral_dimensions=[
        {
            "count": 2048,
            "spectral_width": 4.0e4,  # in Hz
            "reference_offset": -2e3,  # in Hz
        }
    ],
)

# %%
# **Simulator**
#
# Create the Simulator object and add the method and the spin system object.
sim = Simulator()
sim.spin_systems += [spin_system]  # add the spin system
sim.methods += [method]  # add the method
sim.run()

# %%
# **Post-Simulation Processing**
#
# Add post-simulation signal processing.
processor = sp.SignalProcessor(
    operations=[
        sp.IFFT(),
        apo.Exponential(FWHM="30 Hz"),
        sp.FFT(),
    ]
)
processed_data = processor.apply_operations(data=sim.methods[0].simulation)

# %%
# **Plot**
#
# The plot of the simulation before signal processing.
plt.figure(figsize=[4.5, 3.0])
ax = plt.subplot(projection="csdm")
ax.plot(processed_data.real, color="black", linewidth=0.5)
ax.invert_xaxis()
plt.tight_layout()
plt.show()
