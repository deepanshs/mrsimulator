#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
2H MAS NMR of Methionine
^^^^^^^^^^^^^^^^^^^^^^^^
"""
# %%
# The following is a least-squares fitting example of a :math:`^{2}\text{H}` MAS NMR
# spectrum of Methionine. The experimental dataset is a part of DMFIT [#f1]_ examples,
# and we acknowledge Dr. Dominique Massiot for sharing the dataset.
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
host = "https://nmr.cemhti.cnrs-orleans.fr/Dmfit/Help/csdm/"
filename = "2H%20methiodine%20MAS.csdf"
experiment = cp.load(host + filename)

# standard deviation of noise from the dataset
sigma = 0.3026282

# For spectral fitting, we only focus on the real part of the complex dataset
experiment = experiment.real

# Convert the coordinates along each dimension from Hz to ppm.
_ = [item.to("ppm", "nmr_frequency_ratio") for item in experiment.dimensions]

# plot of the dataset.
plt.figure(figsize=(4.25, 3.0))
ax = plt.subplot(projection="csdm")
ax.plot(experiment, "k", alpha=0.5)
ax.set_xlim(600, -700)
plt.grid()
plt.tight_layout()
plt.show()

# %%
# Create a fitting model
# ----------------------
# **Spin System**
H_2 = Site(
    isotope="2H",
    isotropic_chemical_shift=-57.0,  # in ppm,
    quadrupolar={"Cq": 3e4, "eta": 0},  # Cq in Hz
)

spin_systems = [SpinSystem(sites=[H_2])]

# %%
# **Method**

# Get the spectral dimension paramters from the experiment.
spectral_dims = get_spectral_dimensions(experiment)

method = BlochDecaySpectrum(
    channels=["2H"],
    magnetic_flux_density=9.4,  # in T
    rotor_frequency=4517,  # in Hz
    spectral_dimensions=spectral_dims,
    experiment=experiment,  # experimental dataset
)

# Optimize the script by pre-setting the transition pathways for each spin system from
# the method.
for sys in spin_systems:
    sys.transition_pathways = method.get_transition_pathways(sys)

# %%
# **Guess Model Spectrum**

# Simulation
# ----------
sim = Simulator(spin_systems=spin_systems, methods=[method])
sim.run()

# Post Simulation Processing
# --------------------------
processor = sp.SignalProcessor(
    operations=[
        sp.IFFT(),
        sp.apodization.Exponential(FWHM="130 Hz"),
        sp.FFT(),
        sp.Scale(factor=100),
    ]
)
processed_data = processor.apply_operations(data=sim.methods[0].simulation).real

# Plot of the guess Spectrum
# --------------------------
plt.figure(figsize=(4.25, 3.0))
ax = plt.subplot(projection="csdm")
ax.plot(experiment, "k", linewidth=1, label="Experiment")
ax.plot(processed_data, "r", alpha=0.75, linewidth=1, label="guess spectrum")
ax.set_xlim(600, -700)
plt.grid()
plt.legend()
plt.tight_layout()
plt.show()


# %%
# Least-squares minimization with LMFIT
# -------------------------------------
# Use the :func:`~mrsimulator.utils.spectral_fitting.make_LMFIT_params` for a quick
# setup of the fitting parameters.
params = sf.make_LMFIT_params(sim, processor)
print(params.pretty_print(columns=["value", "min", "max", "vary", "expr"]))

# %%
# **Solve the minimizer using LMFIT**
minner = Minimizer(sf.LMFIT_min_function, params, fcn_args=(sim, processor, sigma))
result = minner.minimize(method="powell")
report_fit(result)

# %%
# The best fit solution
# ---------------------
best_fit = sf.bestfit(sim, processor)[0]
residuals = sf.residuals(sim, processor)[0]

# Plot the spectrum
plt.figure(figsize=(4.25, 3.0))
ax = plt.subplot(projection="csdm")
ax.plot(experiment, "k", linewidth=1, label="Experiment")
ax.plot(best_fit, "r", alpha=0.75, linewidth=1, label="Best Fit")
ax.plot(residuals, alpha=0.75, linewidth=1, label="Residual")
ax.set_xlim(600, -700)
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
