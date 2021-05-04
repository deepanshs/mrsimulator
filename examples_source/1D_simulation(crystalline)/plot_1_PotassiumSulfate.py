#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
Potassium Sulfate, 33S (I=3/2)
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

33S (I=3/2) quadrupolar spectrum simulation.
"""
# %%
# The following example is the :math:`^{33}\text{S}` NMR spectrum simulation of
# potassium sulfate (:math:`\text{K}_2\text{SO}_4`). The quadrupole tensor parameters
# for :math:`^{33}\text{S}` is obtained from Moudrakovski `et. al.` [#f3]_
import matplotlib.pyplot as plt

from mrsimulator import Simulator, SpinSystem, Site
from mrsimulator.methods import BlochDecayCTSpectrum
from mrsimulator import signal_processing as sp

# sphinx_gallery_thumbnail_number = 2

# %%
# **Step 1:** Create the spin system
site = Site(
    name="33S",
    isotope="33S",
    isotropic_chemical_shift=335.7,  # in ppm
    quadrupolar={"Cq": 0.959e6, "eta": 0.42},  # Cq is in Hz
)
spin_system = SpinSystem(sites=[site])

# %%
# **Step 2:** Create a central transition selective Bloch decay spectrum method.
method = BlochDecayCTSpectrum(
    channels=["33S"],
    magnetic_flux_density=21.14,  # in T
    rotor_frequency=14000,  # in Hz
    spectral_dimensions=[
        {
            "count": 2048,
            "spectral_width": 5000,  # in Hz
            "reference_offset": 22500,  # in Hz
            "label": r"$^{33}$S resonances",
        }
    ],
)

# %%
# **Step 3:** Create the Simulator object and add method and spin system objects.
sim = Simulator()
sim.spin_systems += [spin_system]  # add the spin system
sim.methods += [method]  # add the method

# %%
# **Step 4:** Simulate the spectrum.
sim.run()

# The plot of the simulation before signal processing.
plt.figure(figsize=(4.25, 3.0))
ax = plt.subplot(projection="csdm")
ax.plot(sim.methods[0].simulation.real, color="black", linewidth=1)
ax.invert_xaxis()
plt.tight_layout()
plt.show()

# %%
# **Step 5:** Add post-simulation signal processing.
processor = sp.SignalProcessor(
    operations=[sp.IFFT(), sp.apodization.Exponential(FWHM="10 Hz"), sp.FFT()]
)
processed_data = processor.apply_operations(data=sim.methods[0].simulation)

# The plot of the simulation after signal processing.
plt.figure(figsize=(4.25, 3.0))
ax = plt.subplot(projection="csdm")
ax.plot(processed_data.real, color="black", linewidth=1)
ax.invert_xaxis()
plt.tight_layout()
plt.show()

# %%
# .. [#f3] Moudrakovski, I., Lang, S., Patchkovskii, S., and Ripmeester, J. High field
#       :math:`^{33}\text{S}` solid state NMR and first-principles calculations in
#       potassium sulfates. J. Phys. Chem. A, 2010, **114**, *1*, 309–316.
#       `DOI: 10.1021/jp908206c <https://doi.org/10.1021/jp908206c>`_
