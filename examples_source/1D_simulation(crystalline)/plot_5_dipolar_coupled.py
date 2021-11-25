#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
Coupled spin-1/2 (Static dipolar spectrum)
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

¹³C-¹H static dipolar coupling simulation.
"""
# %%
import matplotlib.pyplot as plt

from mrsimulator import Simulator, SpinSystem, Site, Coupling
from mrsimulator.methods import BlochDecaySpectrum
from mrsimulator import signal_processing as sp
from mrsimulator.spin_system.tensors import SymmetricTensor

# sphinx_gallery_thumbnail_number = 1

# %%
# **Spin Systems**
#
# Create a 13C-1H coupled spin system.
spin_system = SpinSystem(
    sites=[
        Site(isotope="13C", isotropic_chemical_shift=0.0),
        Site(isotope="1H", isotropic_chemical_shift=0.0),
    ],
    couplings=[Coupling(site_index=[0, 1], dipolar=SymmetricTensor(D=-2e4))],
)
# %%
# **Methods**
#
# Create a BlochDecaySpectrum method.
method = BlochDecaySpectrum(
    channels=["13C"],
    magnetic_flux_density=9.4,  # in T
    spectral_dimensions=[dict(count=2048, spectral_width=8.0e4)],
)

# %%
# **Simulator**
#
# Create the Simulator object and add the method and the spin system object.
sim = Simulator()
sim.spin_systems = [spin_system]  # add the spin system.
sim.methods = [method]  # add the method.
sim.run()

# %%
# **Post-Simulation Processing**
#
# Add post-simulation signal processing.
processor = sp.SignalProcessor(
    operations=[
        sp.IFFT(),
        sp.apodization.Exponential(FWHM="500 Hz"),
        sp.FFT(),
    ]
)
processed_data = processor.apply_operations(data=sim.methods[0].simulation)

# %%
# **Plot**
#
plt.figure(figsize=(4.25, 3.0))
ax = plt.subplot(projection="csdm")
ax.plot(processed_data.real, color="black", linewidth=1)
ax.invert_xaxis()
plt.tight_layout()
plt.show()
