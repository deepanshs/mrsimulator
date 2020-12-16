#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
13C 2D MAT NMR of L-Histidine
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
"""
# %%
# The following is an illustration for fitting 2D MAT/PASS datasets. The example dataset
# is a :math:`^{13}\text{C}` 2D MAT spectrum of L-Histidine from Walder `et. al.` [#f1]_
import numpy as np
import csdmpy as cp
import matplotlib as mpl
import matplotlib.pyplot as plt
import mrsimulator.signal_processing as sp
import mrsimulator.signal_processing.apodization as apo
from mrsimulator import Simulator
from mrsimulator.methods import SSB2D
from mrsimulator.utils import get_spectral_dimensions
from mrsimulator.utils.collection import single_site_system_generator
from mrsimulator.utils.spectral_fitting import LMFIT_min_function, make_LMFIT_params
from lmfit import Minimizer, report_fit

# global plot configuration
mpl.rcParams["figure.figsize"] = [4.5, 3.0]
mpl.rcParams["lines.linewidth"] = 0.5
mpl.rcParams["grid.linestyle"] = "--"
# sphinx_gallery_thumbnail_number = 3

# %%
# Import the dataset
# ------------------
filename = "https://sandbox.zenodo.org/record/687656/files/1H13C_CPPASS_LHistidine.csdf"
mat_data = cp.load(filename)

# For the spectral fitting, we only focus on the real part of the complex dataset.
mat_data = mat_data.real

# Convert the coordinates along each dimension from Hz to ppm.
_ = [item.to("ppm", "nmr_frequency_ratio") for item in mat_data.dimensions]

# Normalize the spectrum
mat_data /= mat_data.max()

# %%
# When using the SSB2D method, ensure the horizontal dimension of the dataset is the
# isotropic dimension. Here, we apply an appropriate transpose operation to the dataset.
mat_data = mat_data.T  # transpose

# plot of the dataset.
levels = (np.arange(10) + 0.3) / 15  # contours are drawn at these levels.
ax = plt.subplot(projection="csdm")
cb = ax.contour(mat_data, colors="k", levels=levels, alpha=0.75)
plt.colorbar(cb)
ax.set_xlim(200, 10)
ax.invert_yaxis()
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
    magnetic_flux_density=9.4,  # in T
    rotor_frequency=1500,  # in Hz
    spectral_dimensions=spectral_dims,
    experiment=mat_data,  # add the measurement to the method.
)

# Optimize the script by pre-setting the transition pathways for each spin system from
# the das method.
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
        apo.Exponential(FWHM="50 Hz"),
        sp.IFFT(axis=0),
        sp.Scale(factor=0.6),
    ]
)
processed_data = processor.apply_operations(data=sim.methods[0].simulation).real

# Plot of the guess Spectrum
# --------------------------
ax = plt.subplot(projection="csdm")
cb = ax.contour(mat_data, colors="k", levels=levels, alpha=0.75)
ax.contour(processed_data, colors="r", linestyles="--", levels=levels, alpha=0.75)
plt.colorbar(cb)
ax.set_xlim(200, 10)
plt.grid()
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
minner = Minimizer(LMFIT_min_function, params, fcn_args=(sim, processor))
result = minner.minimize()
report_fit(result)


# %%
# The best fit solution
# ---------------------
sim.run()
processed_data = processor.apply_operations(data=sim.methods[0].simulation).real

# Plot of the best fit solution
ax = plt.subplot(projection="csdm")
cb = ax.contour(mat_data, colors="k", levels=levels, alpha=0.75)
ax.contour(processed_data, colors="r", linestyles="--", levels=levels, alpha=0.75)
plt.colorbar(cb)
ax.set_xlim(200, 10)
plt.grid()
plt.tight_layout()
plt.show()

# %%
# .. [#f1] B. J. Walder, K. K. Dey, D. C. Kaseman, J. H. Baltisberger, and P. J.
#       Grandinetti, Sideband separation experiments in NMR with phase incremented
#       echo train acquisition, J. Phys. Chem. 2013, **138**, 174203-1-12.
#       `DOI: 10.1063/1.4803142 <https://doi.org/10.1063/1.4803142>`_
