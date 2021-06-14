#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
²⁷Al MAS NMR of YAG (1st and 2nd order Quad)
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
"""
# %%
# The following is a quadrupolar lineshape fitting example for the 27Al MAS NMR of
# Yttrium aluminum garnet (YAG) crystal.
# The following experimental dataset is a part of DMFIT [#f1]_ examples.
# We thank Dr. Dominique Massiot for sharing the dataset.
import csdmpy as cp
import matplotlib.pyplot as plt
from lmfit import Minimizer

from mrsimulator import Simulator, Site, SpinSystem
from mrsimulator.methods import BlochDecaySpectrum
from mrsimulator import signal_processing as sp
from mrsimulator.utils import spectral_fitting as sf
from mrsimulator.utils import get_spectral_dimensions

# sphinx_gallery_thumbnail_number = 3

# %%
# Import the dataset
# ------------------
host = "https://nmr.cemhti.cnrs-orleans.fr/Dmfit/Help/csdm/"
filename = "27Al Quad MAS YAG 400MHz.csdf"
experiment = cp.load(host + filename)

# standard deviation of noise from the dataset
sigma = 0.4895381

# For spectral fitting, we only focus on the real part of the complex dataset
experiment = experiment.real

# Convert the coordinates along each dimension from Hz to ppm.
_ = [item.to("ppm", "nmr_frequency_ratio") for item in experiment.dimensions]

# plot of the dataset.
plt.figure(figsize=(8, 4))
ax = plt.subplot(projection="csdm")
ax.plot(experiment, color="black", linewidth=0.5, label="Experiment")
ax.set_xlim(1200, -1200)
plt.grid()
plt.tight_layout()
plt.show()

# %%
# Create a fitting model
# ----------------------
# **Guess model**
#
# Create a guess list of spin systems.
Al_1 = Site(
    isotope="27Al",
    isotropic_chemical_shift=76,  # in ppm
    quadrupolar={"Cq": 6e6, "eta": 0},  # Cq in Hz
)

Al_2 = Site(
    isotope="27Al",
    isotropic_chemical_shift=1,  # in ppm
    quadrupolar={"Cq": 5e5, "eta": 0.3},  # Cq in Hz
)
spin_systems = [
    SpinSystem(sites=[Al_1], name="AlO4"),
    SpinSystem(sites=[Al_2], name="AlO6"),
]


# %%
# **Method**

# Get the spectral dimension parameters from the experiment.
spectral_dims = get_spectral_dimensions(experiment)

MAS = BlochDecaySpectrum(
    channels=["27Al"],
    magnetic_flux_density=9.395,  # in T
    rotor_frequency=15250,  # in Hz
    spectral_dimensions=spectral_dims,
    experiment=experiment,  # add the measurement to the method.
)

# Optimize the script by pre-setting the transition pathways for each spin system from
# the method.
for sys in spin_systems:
    sys.transition_pathways = MAS.get_transition_pathways(sys)

# %%
# **Guess Spectrum**

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
        sp.apodization.Gaussian(FWHM="300 Hz"),
        sp.FFT(),
        sp.Scale(factor=50),
    ]
)
processed_data = processor.apply_operations(data=sim.methods[0].simulation).real

# Plot of the guess Spectrum
# --------------------------
plt.figure(figsize=(8, 4))
ax = plt.subplot(projection="csdm")
ax.plot(experiment, color="black", linewidth=0.5, label="Experiment")
ax.plot(processed_data, linewidth=2, alpha=0.6)
ax.set_xlim(1200, -1200)
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
ax.set_xlim(1200, -1200)
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
