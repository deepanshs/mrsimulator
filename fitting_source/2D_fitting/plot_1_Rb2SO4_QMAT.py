#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
⁸⁷Rb 2D QMAT NMR of Rb₂SO₄
^^^^^^^^^^^^^^^^^^^^^^^^^^
"""
# %%
# The following is an illustration for fitting 2D QMAT/QPASS datasets. The example
# dataset is a :math:`^{87}\text{Rb}` 2D QMAT spectrum of :math:`\text{Rb}_2\text{SO}_4`
# from Walder `et al.` [#f1]_
import numpy as np
import csdmpy as cp
import matplotlib.pyplot as plt
from lmfit import Minimizer

from mrsimulator import Simulator, SpinSystem, Site
from mrsimulator.method.lib import SSB2D
from mrsimulator import signal_processor as sp
from mrsimulator.utils import spectral_fitting as sf
from mrsimulator.utils import get_spectral_dimensions
from mrsimulator.spin_system.tensors import SymmetricTensor

# sphinx_gallery_thumbnail_number = 3

# %%
# Import the dataset
# ------------------
filename = "https://ssnmr.org/sites/default/files/mrsimulator/Rb2SO4_QMAT.csdf"
qmat_dataset = cp.load(filename)

# standard deviation of noise from the dataset
sigma = 6.530634

# For the spectral fitting, we only focus on the real part of the complex dataset.
qmat_dataset = qmat_dataset.real

# Convert the coordinates along each dimension from Hz to ppm.
_ = [item.to("ppm", "nmr_frequency_ratio") for item in qmat_dataset.dimensions]

# plot of the dataset.
max_amp = qmat_dataset.max()
levels = (np.arange(31) + 0.15) * max_amp / 32  # contours are drawn at these levels.
options = dict(levels=levels, alpha=1, linewidths=0.5)  # plot options

plt.figure(figsize=(8, 3.5))
ax = plt.subplot(projection="csdm")
ax.contour(qmat_dataset.T, colors="k", **options)
ax.set_xlim(200, -200)
ax.set_ylim(75, -120)
plt.grid()
plt.tight_layout()
plt.show()

# %%
# Create a fitting model
# ----------------------
# **Guess model**
#
# Create a guess list of spin systems.
Rb_1 = Site(
    isotope="87Rb",
    isotropic_chemical_shift=16,  # in ppm
    quadrupolar=SymmetricTensor(Cq=5.3e6, eta=0.1),  # Cq in Hz
)
Rb_2 = Site(
    isotope="87Rb",
    isotropic_chemical_shift=40,  # in ppm
    quadrupolar=SymmetricTensor(Cq=2.2e6, eta=0.95),  # Cq in Hz
)

spin_systems = [SpinSystem(sites=[s]) for s in [Rb_1, Rb_2]]

# %%
# **Method**
#
# Create the SSB2D method.

# Get the spectral dimension parameters from the experiment.
spectral_dims = get_spectral_dimensions(qmat_dataset)

PASS = SSB2D(
    channels=["87Rb"],
    magnetic_flux_density=9.395,  # in T
    rotor_frequency=2604,  # in Hz
    spectral_dimensions=spectral_dims,
    experiment=qmat_dataset,  # add the measurement to the method.
)

# Optimize the script by pre-setting the transition pathways for each spin system from
# the method.
for sys in spin_systems:
    sys.transition_pathways = PASS.get_transition_pathways(sys)

# %%
# **Guess Spectrum**

# Simulation
# ----------
sim = Simulator(spin_systems=spin_systems, methods=[PASS])
sim.run()

# Post Simulation Processing
# --------------------------
processor = sp.SignalProcessor(
    operations=[
        # Lorentzian convolution along the isotropic dimensions.
        sp.FFT(dim_index=0),
        sp.apodization.Gaussian(FWHM="100 Hz"),
        sp.IFFT(dim_index=0),
        sp.Scale(factor=1e4),
    ]
)
processed_dataset = processor.apply_operations(dataset=sim.methods[0].simulation).real

# Plot of the guess Spectrum
# --------------------------
plt.figure(figsize=(8, 3.5))
ax = plt.subplot(projection="csdm")
ax.contour(qmat_dataset.T, colors="k", **options)
ax.contour(processed_dataset.T, colors="r", linestyles="--", **options)
ax.set_xlim(200, -200)
ax.set_ylim(75, -120)
plt.grid()
plt.tight_layout()
plt.show()


# %%
# Least-squares minimization with LMFIT
# -------------------------------------
# Use the :func:`~mrsimulator.utils.spectral_fitting.make_LMFIT_params` for a quick
# setup of the fitting parameters.
params = sf.make_LMFIT_params(sim, processor)
params["SP_0_operation_1_Gaussian_FWHM"].min = 0
print(params.pretty_print(columns=["value", "min", "max", "vary", "expr"]))

# %%
# **Solve the minimizer using LMFIT**
minner = Minimizer(sf.LMFIT_min_function, params, fcn_args=(sim, processor, sigma))
result = minner.minimize()
result


# %%
# The best fit solution
# ---------------------
best_fit = sf.bestfit(sim, processor)[0].real

# Plot of the best fit solution
plt.figure(figsize=(8, 3.5))
ax = plt.subplot(projection="csdm")
ax.contour(qmat_dataset.T, colors="k", **options)
ax.contour(best_fit.T, colors="r", linestyles="--", **options)
ax.set_xlim(200, -200)
ax.set_ylim(75, -120)
plt.grid()
plt.tight_layout()
plt.show()

# %%
# .. [#f1] B. J. Walder, K. K. Dey, D. C. Kaseman, J. H. Baltisberger, and P. J.
#       Grandinetti, Sideband separation experiments in NMR with phase incremented
#       echo train acquisition, J. Phys. Chem. 2013, **138**, 174203-1-12.
#       `DOI: 10.1063/1.4803142 <https://doi.org/10.1063/1.4803142>`_
