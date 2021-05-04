#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
87Rb 2D 3QMAS NMR of RbNO3
^^^^^^^^^^^^^^^^^^^^^^^^^^
"""
# %%
# The following is a 3QMAS fitting example for :math:`\text{RbNO}_3`. The dataset was
# acquired and shared by Brendan Wilson.
import numpy as np
import csdmpy as cp
import matplotlib.pyplot as plt
from lmfit import Minimizer, report_fit

from mrsimulator import Simulator
from mrsimulator.methods import ThreeQ_VAS
from mrsimulator import signal_processing as sp
from mrsimulator.utils import get_spectral_dimensions
from mrsimulator.utils.collection import single_site_system_generator
from mrsimulator.utils.spectral_fitting import LMFIT_min_function, make_LMFIT_params

# sphinx_gallery_thumbnail_number = 3

# %%
# Import the dataset
# ------------------
filename = "https://sandbox.zenodo.org/record/784924/files/RbNO3_MQMAS.csdf"
experiment = cp.load(filename)

# standard deviation of noise from the dataset
sigma = 2.119073

# For spectral fitting, we only focus on the real part of the complex dataset
experiment = experiment.real

# Convert the coordinates along each dimension from Hz to ppm.
_ = [item.to("ppm", "nmr_frequency_ratio") for item in experiment.dimensions]

# Normalize the spectrum
max_amp = experiment.max()

# plot of the dataset.
levels = (np.arange(14) + 1) * max_amp / 15  # contours are drawn at these levels.
options = dict(levels=levels, alpha=0.75, linewidths=0.5)  # plot options

plt.figure(figsize=(4.25, 3.0))
ax = plt.subplot(projection="csdm")
ax.contour(experiment, colors="k", **options)
ax.set_xlim(-20, -50)
ax.set_ylim(-45, -65)
plt.tight_layout()
plt.show()

# %%
# Create a fitting model
# ----------------------
# **Guess model**
#
# Create a guess list of spin systems.

shifts = [-26.4, -28.5, -32]  # in ppm
Cq = [1.65e6, 1.9e6, 1.7e6]  # in  Hz
eta = [0.3, 0.8, 0.4]
abundance = [33.33, 33.33, 33.33]  # in %

spin_systems = single_site_system_generator(
    isotopes="87Rb",
    isotropic_chemical_shifts=shifts,
    quadrupolar={"Cq": Cq, "eta": eta},
    abundance=abundance,
)


# %%
# **Method**

# Create the 3QMAS method.
# Get the spectral dimension paramters from the experiment.
spectral_dims = get_spectral_dimensions(experiment)

method = ThreeQ_VAS(
    channels=["87Rb"],
    magnetic_flux_density=9.4,  # in T
    spectral_dimensions=spectral_dims,
    experiment=experiment,  # add the measurement to the method.
)

# Optimize the script by pre-setting the transition pathways for each spin system from
# the das method.
for sys in spin_systems:
    sys.transition_pathways = method.get_transition_pathways(sys)

# %%
# **Guess Spectrum**

# Simulation
# ----------
sim = Simulator()
sim.spin_systems = spin_systems  # add the spin systems
sim.methods = [method]  # add the method.
sim.run()

# Post Simulation Processing
# --------------------------
processor = sp.SignalProcessor(
    operations=[
        # Gaussian convolution along both dimensions.
        sp.IFFT(dim_index=(0, 1)),
        sp.apodization.Gaussian(FWHM="0.08 kHz", dim_index=0),
        sp.apodization.Gaussian(FWHM="0.1 kHz", dim_index=1),
        sp.FFT(dim_index=(0, 1)),
        sp.Scale(factor=50),
    ]
)
processed_data = processor.apply_operations(data=sim.methods[0].simulation).real

# Plot of the guess Spectrum
# --------------------------
plt.figure(figsize=(4.25, 3.0))
ax = plt.subplot(projection="csdm")
ax.contour(experiment, colors="k", **options)
ax.contour(processed_data, colors="r", linestyles="--", **options)
ax.set_xlim(-20, -50)
ax.set_ylim(-45, -65)
plt.tight_layout()
plt.show()


# %%
# Least-squares minimization with LMFIT
# -------------------------------------
# Use the :func:`~mrsimulator.utils.spectral_fitting.make_LMFIT_params` for a quick
# setup of the fitting parameters.
params = make_LMFIT_params(sim, processor)
print(params.pretty_print(columns=["value", "min", "max", "vary", "expr"]))

# %%
# **Solve the minimizer using LMFIT**
minner = Minimizer(LMFIT_min_function, params, fcn_args=(sim, processor, sigma))
result = minner.minimize()
report_fit(result)

# %%
# The best fit solution
# ---------------------
sim.run()
processed_data = processor.apply_operations(data=sim.methods[0].simulation).real

# Plot the spectrum
plt.figure(figsize=(4.25, 3.0))
ax = plt.subplot(projection="csdm")
ax.contour(experiment, colors="k", **options)
ax.contour(processed_data, colors="r", linestyles="--", **options)
ax.set_xlim(-20, -50)
ax.set_ylim(-45, -65)
plt.grid()
plt.tight_layout()
plt.show()
