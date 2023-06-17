#!/usr/bin/env python
"""
Using Custom Isotopes
^^^^^^^^^^^^^^^^^^^^^

Simulating spectra of custom isotopes.
"""
# %%
# Mrsimulator includes the capability for simulating isotopes with user-defined
# attributes. In this example, we show how to define a custom isotopes and simulate
# spectra using custom isotopes.
from mrsimulator.spin_system.isotope import Isotope
from mrsimulator import Site, SpinSystem, Simulator
from mrsimulator.method import SpectralDimension
from mrsimulator.method.lib import BlochDecayCTSpectrum

import matplotlib.pyplot as plt

# sphinx_gallery_thumbnail_number = 1

# %%
# First, we register a new isotope symbol with custom attributes using the
# :py:meth:`~mrsimulator.spin_system.isotope.Isotope.register` method. The new symbol
# cannot match any real isotope symbols, and re-calling register on a custom symbol will
# update the stored attributes of that isotope.
Isotope.register(
    symbol="custom",
    spin_multiplicity=4,  # Spin of I=3/2
    gyromagnetic_ratio=12.3,
    quadrupole_moment=0.045,
)

# %%
# Create :py:class:`~mrsimulator.spin_system.Site` and
# :py:class:`~mrsimulator.method.Method` objects using the new isotope symbol. The
# syntax is the same as any other isotope.
site = Site(
    isotope="custom",
    isotropic_chemical_shift=-30,
    quadrupolar={"Cq": 1.3e6, "eta": 0.5},
)
sys = SpinSystem(sites=[site])

mth = BlochDecayCTSpectrum(
    channels=["custom"],
    spectral_dimensions=[
        SpectralDimension(spectral_width=12e3, reference_offset=-3000)
    ],
)

# %% Simulate and plot the spectrum
sim = Simulator(spin_systems=[sys], methods=[mth])
sim.run()

plt.figure(figsize=(4.25, 3))
plt.subplot(projection="csdm")
plt.plot(sim.methods[0].simulation.real, color="black", linewidth=1)
plt.tight_layout()
plt.show()
