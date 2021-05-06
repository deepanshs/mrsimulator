#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
1D PASS/MAT sideband order cross-section
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
"""
# %%
# This example illustrates the use of mrsimulator and LMFIT modules in fitting the
# sideband intensity profile across the isotropic chemical shift cross-section from a
# PASS/MAT dataset.
import csdmpy as cp
import matplotlib.pyplot as plt
from lmfit import Minimizer, report_fit

from mrsimulator import Simulator, SpinSystem, Site
from mrsimulator.methods import BlochDecaySpectrum
from mrsimulator import signal_processing as sp
from mrsimulator.utils import spectral_fitting as sf
from mrsimulator.utils import get_spectral_dimensions

# sphinx_gallery_thumbnail_number = 3

# %%
# Import the dataset
# ------------------
name = "https://sandbox.zenodo.org/record/745068/files/LHistidine_cross_section.csdf"
pass_cross_section = cp.load(name)

# standard deviation of noise from the dataset
sigma = 0.01285316

# For the spectral fitting, we only focus on the real part of the complex dataset.
pass_cross_section = pass_cross_section.real

# Convert the coordinates along each dimension from Hz to ppm.
_ = [item.to("ppm", "nmr_frequency_ratio") for item in pass_cross_section.dimensions]

# Normalize the spectrum.
max_amp = pass_cross_section.max()
pass_cross_section /= max_amp
sigma /= max_amp

# The plot of the dataset.
plt.figure(figsize=(4.25, 3.0))
ax = plt.subplot(projection="csdm")
ax.plot(pass_cross_section, "k", alpha=0.5)
ax.invert_xaxis()
plt.tight_layout()
plt.show()

# %%
# Create a fitting model
# ----------------------
# **Guess model**
#
# Create a guess list of spin systems. For fitting the sideband profile at an isotropic
# chemical shift cross-section from PASS/MAT datasets, set the isotropic_chemical_shift
# parameter of the site object as zero.
site = Site(
    isotope="13C",
    isotropic_chemical_shift=0,  #
    shielding_symmetric={"zeta": -70, "eta": 0.8},
)
spin_systems = [SpinSystem(sites=[site])]

# %%
# **Method**
#
# For the sideband-only cross-section, use the BlochDecaySpectrum method.

# Get the dimension information from the experiment.
spectral_dims = get_spectral_dimensions(pass_cross_section)

method = BlochDecaySpectrum(
    channels=["13C"],
    magnetic_flux_density=9.39,  # in T
    rotor_frequency=1500,  # in Hz
    spectral_dimensions=spectral_dims,
    experiment=pass_cross_section,  # also add the measurement to the method.
)

# Optimize the script by pre-setting the transition pathways for each spin system from
# the method.
for sys in spin_systems:
    sys.transition_pathways = method.get_transition_pathways(sys)

# %%
# **Guess Spectrum**

# Simulation
# ----------
sim = Simulator(spin_systems=spin_systems, methods=[method])
sim.run()

# Post Simulation Processing
# --------------------------
processor = sp.SignalProcessor(operations=[sp.Scale(factor=10)])
processed_data = processor.apply_operations(data=sim.methods[0].simulation).real

# Plot of the guess Spectrum
# --------------------------
plt.figure(figsize=(4.25, 3.0))
ax = plt.subplot(projection="csdm")
ax.plot(pass_cross_section, color="k", linewidth=1, label="Experiment")
ax.plot(processed_data, color="r", alpha=0.75, linewidth=1, label="guess spectrum")
plt.grid()
ax.invert_xaxis()
plt.legend()
plt.tight_layout()
plt.show()


# %%
# Least-squares minimization with LMFIT
# -------------------------------------
# First, create the fitting parameters.
# Use the :func:`~mrsimulator.utils.spectral_fitting.make_LMFIT_params` for a quick
# setup.
params = sf.make_LMFIT_params(sim, processor)

# Fix the value of the isotropic chemical shift to zero for pure anisotropic sideband
# amplitude simulation.
params["sys_0_site_0_isotropic_chemical_shift"].vary = False
print(params.pretty_print(columns=["value", "min", "max", "vary", "expr"]))

# %%
# Run the minimization using LMFIT
minner = Minimizer(sf.LMFIT_min_function, params, fcn_args=(sim, processor, sigma))
result = minner.minimize()
report_fit(result)

# %%
# The best fit solution
# ---------------------
best_fit = sf.bestfit(sim, processor)[0]
residuals = sf.residuals(sim, processor)[0]

# Plot the spectrum
plt.figure(figsize=(4.25, 3.0))
ax = plt.subplot(projection="csdm")
ax.plot(pass_cross_section, color="k", linewidth=1, label="Experiment")
ax.plot(best_fit, "r", alpha=0.75, linewidth=1, label="Best Fit")
ax.plot(residuals, alpha=0.75, linewidth=1, label="Residuals")
ax.invert_xaxis()
plt.grid()
plt.legend()
plt.tight_layout()
plt.show()
