#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
CoCl₂.2D₂O, ²H (I=1) Shifting-d echo
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

²H (I=1) 2D NMR CSA-Quad 1st order correlation spectrum.
"""
# %%
# The following is an example of fitting static shifting-*d* echo NMR correlation
# spectrum of :math:`\text{NiCl}_2\cdot 2\text{D}_2\text{O}` crystalline solid. The
# spectrum used here is from Walder `et al.` [#f1]_.
import numpy as np
import csdmpy as cp
import matplotlib.pyplot as plt
from lmfit import Minimizer, report_fit

from mrsimulator import Simulator, Site, SpinSystem
from mrsimulator.methods import Method2D
from mrsimulator import signal_processing as sp
from mrsimulator.utils import spectral_fitting as sf
from mrsimulator.utils import get_spectral_dimensions

# sphinx_gallery_thumbnail_number = 3

# %%
# Import the dataset
# ------------------
filename = "https://sandbox.zenodo.org/record/830903/files/CoCl2.2D2O.csdf"
experiment = cp.load(filename)

# standard deviation of noise from the dataset
sigma = 11.11578

# For spectral fitting, we only focus on the real part of the complex dataset
experiment = experiment.real

# Convert the coordinates along each dimension from Hz to ppm.
_ = [item.to("ppm", "nmr_frequency_ratio") for item in experiment.dimensions]

# plot of the dataset.
max_amp = experiment.max()
levels = (np.arange(24) + 1) * max_amp / 25  # contours are drawn at these levels.
options = dict(levels=levels, linewidths=0.5)  # plot options

plt.figure(figsize=(4.25, 3.0))
ax = plt.subplot(projection="csdm")
ax.contour(experiment, colors="k", **options)
ax.set_xlim(2500, -1500)
ax.set_ylim(2000, -2200)
plt.grid()
plt.tight_layout()
plt.show()

# %%
# Create a fitting model
# ----------------------
# **Guess model**
#
# Create a guess list of spin systems.
site = Site(
    isotope="2H",
    isotropic_chemical_shift=200,  # in ppm
    shielding_symmetric={
        "zeta": -1300,  # in ppm
        "eta": 0.2,
        "alpha": np.pi,  # in rads
        "beta": np.pi / 2,  # in rads
        "gamma": np.pi / 2,  # in rads
    },
    quadrupolar={"Cq": 110e3, "eta": 0.83},  # Cq in Hz
)

spin_systems = [SpinSystem(sites=[site])]

# %%
# **Method**
#
# Use the generic 2D method, `Method2D`, to generate a shifting-d echo method.

# Get the spectral dimension parameters from the experiment.
spectral_dims = get_spectral_dimensions(experiment)

shifting_d = Method2D(
    channels=["2H"],
    magnetic_flux_density=9.395,  # in T
    spectral_dimensions=[
        {
            **spectral_dims[0],
            "label": "Quadrupolar frequency",
            "events": [
                {
                    "rotor_frequency": 0,
                    "transition_query": {"P": [-1]},
                    "freq_contrib": ["Quad1_2"],
                }
            ],
        },
        {
            **spectral_dims[1],
            "label": "Paramagnetic shift",
            "events": [
                {
                    "rotor_frequency": 0,
                    "transition_query": {"P": [-1]},
                    "freq_contrib": ["Shielding1_0", "Shielding1_2"],
                }
            ],
        },
    ],
    experiment=experiment,  # also add the measurement to the method.
)

# Optimize the script by pre-setting the transition pathways for each spin system from
# the method.
for sys in spin_systems:
    sys.transition_pathways = shifting_d.get_transition_pathways(sys)

# %%
# **Guess Spectrum**

# Simulation
# ----------
sim = Simulator(spin_systems=spin_systems, methods=[shifting_d])
sim.config.integration_volume = "hemisphere"
sim.run()

# Post Simulation Processing
# --------------------------
processor = sp.SignalProcessor(
    operations=[
        # Gaussian convolution along both dimensions.
        sp.IFFT(dim_index=0),
        sp.apodization.Gaussian(FWHM="10 kHz", dim_index=0),  # along dimension 0
        sp.FFT(dim_index=0),
        sp.Scale(factor=5e8),
    ]
)
processed_data = processor.apply_operations(data=sim.methods[0].simulation).real

# Plot of the guess Spectrum
# --------------------------
plt.figure(figsize=(4.25, 3.0))
ax = plt.subplot(projection="csdm")
ax.contour(experiment, colors="k", **options)
ax.contour(processed_data, colors="r", linestyles="--", **options)
ax.set_xlim(2500, -1500)
ax.set_ylim(2000, -2200)
plt.grid()
plt.tight_layout()
plt.show()

# %%
# Least-squares minimization with LMFIT
# -------------------------------------
# Use the :func:`~mrsimulator.utils.spectral_fitting.make_LMFIT_params` for a quick
# setup of the fitting parameters.
params = sf.make_LMFIT_params(sim, processor)
params["sys_0_site_0_shielding_symmetric_alpha"].vary = False
params["sys_0_site_0_shielding_symmetric_beta"].vary = False
params["sys_0_site_0_shielding_symmetric_gamma"].vary = False
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

# Plot the spectrum
plt.figure(figsize=(4.25, 3.0))
ax = plt.subplot(projection="csdm")
ax.contour(experiment, colors="k", **options)
ax.contour(best_fit, colors="r", linestyles="--", **options)
ax.set_xlim(2500, -1500)
ax.set_ylim(2000, -2200)
plt.grid()
plt.tight_layout()
plt.show()


# %%
# Image plots with residuals
# --------------------------
residuals = sf.residuals(sim, processor)[0]

fig, ax = plt.subplots(
    1, 3, sharey=True, figsize=(10, 3.0), subplot_kw={"projection": "csdm"}
)
vmax, vmin = experiment.max(), experiment.min()
for i, dat in enumerate([experiment, best_fit, residuals]):
    ax[i].imshow(dat, aspect="auto", cmap="gist_ncar_r", vmax=vmax, vmin=vmin)
    ax[i].set_xlim(2500, -1500)
ax[0].set_ylim(2000, -2200)
plt.tight_layout()
plt.show()
# %%
# .. [#f1] Walder B.J, Patterson A.M., Baltisberger J.H, and Grandinetti P.J
#       Hydrogen motional disorder in crystalline iron group chloride dihydrates
#       spectroscopy, J. Chem. Phys. (2018)  **149**, 084503.
#       `DOI: 10.1063/1.5037151 <https://doi.org/10.1063/1.5037151>`_
