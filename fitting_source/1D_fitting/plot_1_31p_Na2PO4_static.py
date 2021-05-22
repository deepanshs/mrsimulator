#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
31P static NMR of crystalline Na2PO4 (CSA)
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
"""
# %%
# The following is a CSA static least-squares fitting example of a
# :math:`^{31}\text{P}` MAS NMR spectrum of :math:`\text{Na}_{2}\text{PO}_{4}`.
# The following experimental dataset is a part of DMFIT [#f1]_ examples, and we
# acknowledge Dr. Dominique Massiot for sharing the dataset.
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
filename = "31P Phophonate Static.csdf"
experiment = cp.load(host + filename)

# standard deviation of noise from the dataset
sigma = 3.258224

# For spectral fitting, we only focus on the real part of the complex dataset
experiment = experiment.real

# Convert the coordinates along each dimension from Hz to ppm.
_ = [item.to("ppm", "nmr_frequency_ratio") for item in experiment.dimensions]

# plot of the dataset.
plt.figure(figsize=(4.25, 3.0))
ax = plt.subplot(projection="csdm")
ax.plot(experiment, color="black", linewidth=0.5, label="Experiment")
ax.set_xlim(200, -200)
plt.grid()
plt.tight_layout()
plt.show()

# %%
# Create a fitting model
# ----------------------
# **Spin System**
P_31 = Site(
    isotope="31P",
    isotropic_chemical_shift=5.0,  # in ppm,
    shielding_symmetric={"zeta": -80, "eta": 0.5},  # zeta in Hz
)

spin_systems = [SpinSystem(sites=[P_31])]

# %%
# **Method**

# Get the spectral dimension parameters from the experiment.
spectral_dims = get_spectral_dimensions(experiment)

static1D = BlochDecaySpectrum(
    channels=["31P"],
    magnetic_flux_density=9.395,  # in T
    rotor_frequency=0,  # in Hz
    spectral_dimensions=spectral_dims,
    experiment=experiment,  # experimental dataset
)

# Optimize the script by pre-setting the transition pathways for each spin system from
# the method.
for sys in spin_systems:
    sys.transition_pathways = static1D.get_transition_pathways(sys)

# %%
# **Guess Model Spectrum**

# Simulation
# ----------
sim = Simulator(spin_systems=spin_systems, methods=[static1D])
sim.run()

# Post Simulation Processing
# --------------------------
processor = sp.SignalProcessor(
    operations=[
        sp.IFFT(),
        sp.apodization.Gaussian(FWHM="3000 Hz"),
        sp.FFT(),
        sp.Scale(factor=4000),
    ]
)
processed_data = processor.apply_operations(data=sim.methods[0].simulation).real

# Plot of the guess Spectrum
# --------------------------
plt.figure(figsize=(4.25, 3.0))
ax = plt.subplot(projection="csdm")
ax.plot(experiment, color="black", linewidth=0.5, label="Experiment")
ax.plot(processed_data, linewidth=2, alpha=0.6, label="Guess Spectrum")
ax.set_xlim(200, -200)
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
params.pop("sys_0_abundance")
print(params.pretty_print(columns=["value", "min", "max", "vary", "expr"]))

# %%
# **Solve the minimizer using LMFIT**
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
ax.plot(experiment, color="black", linewidth=0.5, label="Experiment")
ax.plot(residuals, color="gray", linewidth=0.5, label="Residual")
ax.plot(best_fit, linewidth=2, alpha=0.6, label="Best Fit")
ax.set_xlim(200, -200)
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
