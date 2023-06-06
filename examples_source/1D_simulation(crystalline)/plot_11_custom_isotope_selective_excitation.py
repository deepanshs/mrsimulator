#!/usr/bin/env python
"""
Selective Excitations using Custom Isotopes
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

Simulating spectra of custom isotopes
"""

# %%
# The isotope register method has the `copy_from` argument which can be used to create
# individually addressable channels which the same isotope attributes. This
# functionality can be used to emulate selective-excitation experiments by querying
# specific transitions on different sites
from mrsimulator.spin_system.isotope import Isotope
from mrsimulator import Site, SpinSystem, Coupling, Simulator, Method
from mrsimulator.method import SpectralDimension, SpectralEvent
from mrsimulator.method.lib import BlochDecaySpectrum, BlochDecayCTSpectrum
from mrsimulator.method.query import TransitionQuery
import mrsimulator.signal_processor as sp

import numpy as np
import matplotlib.pyplot as plt

# %%
# Below, we show how custom isotopes can be copied from existing isotopes, and how these
# custom isotopes can be used to simulate a selective-excitation spectrum of protons.
# First, two new isotopes -- ``"1H-a"`` and ``"1H-b"`` -- are registered from the known
# proton isotope by using the ``copy_from`` keyword argument.
Isotope.register("1H-a", copy_from="1H")
Isotope.register("1H-b", copy_from="1H")

# %%
# Next, we create two site objects using the previously registered ``"1H-a"`` and
# ``"1H-b"`` isotope symbols. The equivalent proton system (without custom isotopes) is
# also constructed as a comparison.
site_a = Site(isotope="1H", isotropic_chemical_shift=-0.5)
site_b = Site(isotope="1H", isotropic_chemical_shift=-2.0)
coupling_ab = Coupling(site_index=[0, 1], isotropic_j=48)
sys = SpinSystem(sites=[site_a, site_b], couplings=[coupling_ab])  # proton system

# Create 1D BlochDecaySpectrum for proton
spec_dim = SpectralDimension(count=10000, spectral_width=2000, reference_offset=-500)
mth = BlochDecaySpectrum(channels=["1H"], spectral_dimensions=[spec_dim])

# Create Simulator object for coupled proton system
sim = Simulator(spin_systems=[sys], methods=[mth])


# Create couple proton system from custom isotopes
site_a_custom = Site(isotope="1H-a", isotropic_chemical_shift=-0.5)
site_b_custom = Site(isotope="1H-b", isotropic_chemical_shift=-2.0)
sys_custom = SpinSystem(sites=[site_a_custom, site_b_custom], couplings=[coupling_ab])

# Create two BlochDecaySpectrum methods with custom isotope channels
mth_a = BlochDecaySpectrum(channels=["1H-a"], spectral_dimensions=[spec_dim])
mth_b = BlochDecaySpectrum(channels=["1H-b"], spectral_dimensions=[spec_dim])

# %%
# Create and run two simulator objects for the simulating the regular and selective
# excitation proton spectra.
sim = Simulator(spin_systems=[sys], methods=[mth])
sim_custom = Simulator(spin_systems=[sys_custom], methods=[mth_a, mth_b])

# Run both the simulations
sim.run()
sim_custom.run()

# %%
# Apply post-simulation signal processing and plot the spectra. Notice how the full
# proton spectrum (right) is the summation of the two individual proton spectra (left),
# and how the J-couplings are still present in the selective excitation spectra.
# Make processor for adding line broadening to spectra
processor = sp.SignalProcessor(
    operations=[sp.FFT(), sp.apodization.Exponential(FWHM="1 Hz"), sp.IFFT()]
)
spec_proton = processor.apply_operations(sim.methods[0].simulation).real
spec_a = processor.apply_operations(sim_custom.methods[0].simulation).real
spec_b = processor.apply_operations(sim_custom.methods[1].simulation).real

fig, ax = plt.subplots(1, 2, figsize=(6.5, 3), subplot_kw={"projection": "csdm"})

ax[0].plot(spec_proton)
ax[0].set_title("Full Proton Spectrum")
ax[1].plot(spec_a, label="1H-a")
ax[1].plot(spec_b, label="1H-b")
ax[1].legend()
ax[1].set_title("Individual proton spectra")

plt.tight_layout()
plt.show()
