#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
¹³C MAS NMR of Glycine (CSA) [1940 Hz] (TO REMOVE)
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
"""
# %%
# The following is a sideband least-squares fitting example of a
# :math:`^{13}\text{C}` MAS NMR spectrum of Glycine spinning at 1940 Hz.
# The following experimental dataset is a part of DMFIT [#f1]_ examples.
# We thank Dr. Dominique Massiot for sharing the dataset.
import csdmpy as cp
import matplotlib.pyplot as plt
from lmfit import Minimizer

from mrsimulator import Simulator, SpinSystem, Site
from mrsimulator.methods import BlochDecaySpectrum
from mrsimulator import signal_processing as sp
from mrsimulator.utils import spectral_fitting as sf
from mrsimulator.utils import get_spectral_dimensions

# sphinx_gallery_thumbnail_number = 3

# %%
# Import the dataset
# ------------------
host = "https://nmr.cemhti.cnrs-orleans.fr/Dmfit/Help/csdm/"
filename = "13C MAS 1940Hz - Glycine.csdf"
experiment = cp.load(host + filename)

# standard deviation of noise from the dataset
sigma = 3.890551

# For spectral fitting, we only focus on the real part of the complex dataset
experiment = experiment.real

# Convert the coordinates along each dimension from Hz to ppm.
_ = [item.to("ppm", "nmr_frequency_ratio") for item in experiment.dimensions]

# plot of the dataset.
plt.figure(figsize=(8, 4))
ax = plt.subplot(projection="csdm")
ax.plot(experiment, color="black", linewidth=0.5, label="Experiment")
ax.set_xlim(280, -10)
plt.grid()
plt.tight_layout()
plt.show()

# %%
# Create a fitting model
# ----------------------
# **Spin System**
C1 = Site(
    isotope="13C",
    isotropic_chemical_shift=176.0,  # in ppm
    shielding_symmetric={"zeta": 70, "eta": 0.6},  # zeta in Hz
)
C2 = Site(
    isotope="13C",
    isotropic_chemical_shift=43.0,  # in ppm
    shielding_symmetric={"zeta": 30, "eta": 0.5},  # zeta in Hz
)

spin_systems = [SpinSystem(sites=[C1], name="C1"), SpinSystem(sites=[C2], name="C2")]

# %%
# **Method**

# Get the spectral dimension parameters from the experiment.
spectral_dims = get_spectral_dimensions(experiment)

MAS = BlochDecaySpectrum(
    channels=["13C"],
    magnetic_flux_density=7.05,  # in T
    rotor_frequency=1940,  # in Hz
    spectral_dimensions=spectral_dims,
    experiment=experiment,  # experimental dataset
)

# Optimize the script by pre-setting the transition pathways for each spin system from
# the method.
for sys in spin_systems:
    sys.transition_pathways = MAS.get_transition_pathways(sys)

# %%
# **Guess Model Spectrum**

# Simulation
# ----------
sim = Simulator(spin_systems=spin_systems, methods=[MAS])
sim.config.decompose_spectrum = "spin_system"
sim.run()

# Post Simulation Processing
# --------------------------
processor = sp.SignalProcessor(
    operations=[
        sp.IFFT(),
        sp.apodization.Exponential(FWHM="20 Hz", dv_index=0),  # spin system 0
        sp.apodization.Exponential(FWHM="200 Hz", dv_index=1),  # spin system 1
        sp.FFT(),
        sp.Scale(factor=100),
    ]
)
processed_data = processor.apply_operations(data=sim.methods[0].simulation).real

# Plot of the guess Spectrum
# --------------------------
plt.figure(figsize=(8, 4))
ax = plt.subplot(projection="csdm")
ax.plot(experiment, color="black", linewidth=0.5, label="Experiment")
ax.plot(processed_data, linewidth=2, alpha=0.6)
ax.set_xlim(280, -10)
plt.grid()
plt.legend()
plt.tight_layout()
plt.show()


# %%
# Least-squares minimization with LMFIT
# -------------------------------------
# Use the :func:`~mrsimulator.utils.spectral_fitting.make_LMFIT_params` for a quick
# setup of the fitting parameters.
params = sf.make_LMFIT_params(sim, processor, include={"rotor_frequency"})
params["sys_1_site_0_shielding_symmetric_eta"].vary = False
print(params.pretty_print(columns=["value", "min", "max", "vary", "expr"]))

# %%
# **Solve the minimizer using LMFIT**
minner = Minimizer(sf.LMFIT_min_function, params, fcn_args=(sim, processor, sigma))
result = minner.minimize()
result

# %%
# The best fit solution
# ---------------------
best_fit = sf.bestfit(sim, processor)[0]
residuals = sf.residuals(sim, processor)[0]

# Plot the spectrum
plt.figure(figsize=(8, 4))
ax = plt.subplot(projection="csdm")
ax.plot(experiment, color="black", linewidth=0.5, label="Experiment")
ax.plot(residuals, color="gray", linewidth=0.5, label="Residual")
ax.plot(best_fit, linewidth=2, alpha=0.6)
ax.set_xlim(280, -10)
plt.grid()
plt.legend()
plt.tight_layout()
plt.show()

# %%
#
# .. [#f1] D.Massiot, F.Fayon, M.Capron, I.King, S.Le Calvé, B.Alonso, J.O.Durand,
#       B.Bujoli, Z.Gan, G.Hoatson, 'Modelling one and two-dimensional solid-state NMR
#       spectra.', Magn. Reson. Chem. **40** 70-76 (2002)
#       `DOI: 10.1002/mrc.984 <https://doi.org/10.1002/mrc.984>`_
