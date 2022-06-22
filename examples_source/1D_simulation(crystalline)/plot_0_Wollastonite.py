#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
Wollastonite, ²⁹Si (I=1/2)
^^^^^^^^^^^^^^^^^^^^^^^^^^

²⁹Si (I=1/2) spinning sideband simulation.
"""
# %%
# Wollastonite is a high-temperature calcium-silicate,
# :math:`\beta−\text{Ca}_3\text{Si}_3\text{O}_9`, with three distinct
# :math:`^{29}\text{Si}` sites. The :math:`^{29}\text{Si}` tensor parameters
# were obtained from Hansen `et al.` [#f1]_
import matplotlib.pyplot as plt

from mrsimulator import Simulator, SpinSystem, Site
from mrsimulator import signal_processor as sp
from mrsimulator.method.lib import BlochDecaySpectrum
from mrsimulator.spin_system.tensors import SymmetricTensor
from mrsimulator.method import SpectralDimension

# sphinx_gallery_thumbnail_number = 3

# %%
# **Step 1:** Create the sites.
Si29_1 = Site(
    isotope="29Si",
    isotropic_chemical_shift=-89.0,  # in ppm
    shielding_symmetric=SymmetricTensor(zeta=59.8, eta=0.62),  # zeta in ppm
)
Si29_2 = Site(
    isotope="29Si",
    isotropic_chemical_shift=-89.5,  # in ppm
    shielding_symmetric=SymmetricTensor(zeta=52.1, eta=0.68),  # zeta in ppm
)
Si29_3 = Site(
    isotope="29Si",
    isotropic_chemical_shift=-87.8,  # in ppm
    shielding_symmetric=SymmetricTensor(zeta=69.4, eta=0.60),  # zeta in ppm
)

# %%
# **Step 2:** Create the spin systems from these sites. Again, we create three
# single-site spin systems for better performance.
spin_systems = [
    SpinSystem(sites=[Si29_1]),
    SpinSystem(sites=[Si29_2]),
    SpinSystem(sites=[Si29_3]),
]

# %%
# **Step 3:** Create a Bloch decay spectrum method.
method = BlochDecaySpectrum(
    channels=["29Si"],
    magnetic_flux_density=14.1,  # in T
    rotor_frequency=1500,  # in Hz
    spectral_dimensions=[
        SpectralDimension(
            count=2048,
            spectral_width=25000,  # in Hz
            reference_offset=-10000,  # in Hz
            label=r"$^{29}$Si resonances",
        )
    ],
)

# A graphical representation of the method object.
plt.figure(figsize=(4, 2))
method.plot()
plt.show()

# %%
# **Step 4:** Create the Simulator object and add the method and spin system objects.
sim = Simulator()
sim.spin_systems = spin_systems  # add the spin systems
sim.methods = [method]  # add the method

# %%
# **Step 5:** Simulate the spectrum.
sim.run()

# The plot of the simulation before signal processing.
plt.figure(figsize=(4.25, 3.0))
ax = plt.subplot(projection="csdm")
ax.plot(sim.methods[0].simulation.real, color="black", linewidth=1)
ax.invert_xaxis()
plt.tight_layout()
plt.show()

# %%
# **Step 6:** Add post-simulation signal processing.
processor = sp.SignalProcessor(
    operations=[sp.IFFT(), sp.apodization.Exponential(FWHM="70 Hz"), sp.FFT()]
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
# .. [#f1] Hansen, M. R., Jakobsen, H. J., Skibsted, J., :math:`^{29}\text{Si}`
#       Chemical Shift Anisotropies in Calcium Silicates from High-Field
#       :math:`^{29}\text{Si}` MAS NMR Spectroscopy, Inorg. Chem. 2003,
#       **42**, *7*, 2368-2377.
#       `DOI: 10.1021/ic020647f <https://doi.org/10.1021/ic020647f>`_
