#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
13C 2D MAT NMR of L-Histidine
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
"""
# %%
# The following is an illustration for fitting 2D MAT/PASS datasets. The example dataset
# is a :math:`^{13}\text{C}` 2D MAT spectrum of L-Histidine from Walder `et al.` [#f1]_
import numpy as np
import csdmpy as cp
import matplotlib.pyplot as plt
from lmfit import Minimizer, report_fit

from mrsimulator import Simulator
from mrsimulator.methods import SSB2D
from mrsimulator import signal_processing as sp
from mrsimulator.utils import spectral_fitting as sf
from mrsimulator.utils import get_spectral_dimensions
from mrsimulator.utils.collection import single_site_system_generator

# sphinx_gallery_thumbnail_number = 3

# %%
# Import the dataset
# ------------------
filename = "https://sandbox.zenodo.org/record/814455/files/1H13C_CPPASS_LHistidine.csdf"
mat_data = cp.load(filename)

# standard deviation of noise from the dataset
sigma = 0.4192854

# For the spectral fitting, we only focus on the real part of the complex dataset.
mat_data = mat_data.real

# Convert the coordinates along each dimension from Hz to ppm.
_ = [item.to("ppm", "nmr_frequency_ratio") for item in mat_data.dimensions]

# %%
# When using the SSB2D method, ensure the horizontal dimension of the dataset is the
# isotropic dimension. Here, we apply an appropriate transpose operation to the dataset.
mat_data = mat_data.T  # transpose

# plot of the dataset.
max_amp = mat_data.max()
levels = (np.arange(24) + 1) * max_amp / 25  # contours are drawn at these levels.
options = dict(levels=levels, alpha=0.75, linewidths=0.5)  # plot options

plt.figure(figsize=(4.25, 3.0))
ax = plt.subplot(projection="csdm")
ax.contour(mat_data, colors="k", **options)
ax.set_xlim(200, 10)
ax.invert_yaxis()
plt.grid()
plt.tight_layout()
plt.show()


# %%
# Create a fitting model
# ----------------------
# **Guess model**
#
# Create a guess list of spin systems.

shifts = [120, 128, 135, 175, 55, 25]  # in ppm
zeta = [-70, -65, -60, -60, -10, -10]  # in ppm
eta = [0.8, 0.4, 0.9, 0.3, 0.0, 0.0]

spin_systems = single_site_system_generator(
    isotopes="13C",
    isotropic_chemical_shifts=shifts,
    shielding_symmetric={"zeta": zeta, "eta": eta},
    abundance=100 / 6,
)

# %%
# **Method**

# Create the DAS method.
# Get the spectral dimension paramters from the experiment.
spectral_dims = get_spectral_dimensions(mat_data)

ssb = SSB2D(
    channels=["13C"],
    magnetic_flux_density=9.395,  # in T
    rotor_frequency=1500,  # in Hz
    spectral_dimensions=spectral_dims,
    experiment=mat_data,  # add the measurement to the method.
)

# Optimize the script by pre-setting the transition pathways for each spin system from
# the method.
for sys in spin_systems:
    sys.transition_pathways = ssb.get_transition_pathways(sys)

# %%
# **Guess Spectrum**

# Simulation
# ----------
sim = Simulator()
sim.spin_systems = spin_systems  # add the spin systems
sim.methods = [ssb]  # add the method
sim.run()

# Post Simulation Processing
# --------------------------
processor = sp.SignalProcessor(
    operations=[
        # Lorentzian convolution along the isotropic dimensions.
        sp.FFT(axis=0),
        sp.apodization.Exponential(FWHM="50 Hz"),
        sp.IFFT(axis=0),
        sp.Scale(factor=60),
    ]
)
processed_data = processor.apply_operations(data=sim.methods[0].simulation).real

# Plot of the guess Spectrum
# --------------------------
plt.figure(figsize=(4.25, 3.0))
ax = plt.subplot(projection="csdm")
ax.contour(mat_data, colors="k", **options)
ax.contour(processed_data, colors="r", linestyles="--", **options)
ax.set_xlim(200, 10)
plt.grid()
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
result = minner.minimize()
report_fit(result)


# %%
# The best fit solution
# ---------------------
best_fit = sf.bestfit(sim, processor)[0]

# Plot of the best fit solution
plt.figure(figsize=(4.25, 3.0))
ax = plt.subplot(projection="csdm")
ax.contour(mat_data, colors="k", **options)
ax.contour(best_fit, colors="r", linestyles="--", **options)
ax.set_xlim(200, 10)
plt.grid()
plt.tight_layout()
plt.show()

# %%
# .. [#f1] B. J. Walder, K. K. Dey, D. C. Kaseman, J. H. Baltisberger, and P. J.
#       Grandinetti, Sideband separation experiments in NMR with phase incremented
#       echo train acquisition, J. Phys. Chem. 2013, **138**, 174203-1-12.
#       `DOI: 10.1063/1.4803142 <https://doi.org/10.1063/1.4803142>`_
