#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
NiCl₂.2D₂O, ²H (I=1) Shifting-d echo
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
from lmfit import Minimizer

from mrsimulator import Simulator, Site, SpinSystem
from mrsimulator import signal_processing as sp
from mrsimulator.utils import spectral_fitting as sf
from mrsimulator.utils import get_spectral_dimensions
from mrsimulator.spin_system.tensors import SymmetricTensor
from mrsimulator.method import Method, SpectralDimension, SpectralEvent, MixingEvent

# sphinx_gallery_thumbnail_number = 3

# %%
# Import the dataset
# ------------------
filename = "https://ssnmr.org/sites/default/files/mrsimulator/NiCl2.2D2O.csdf"
experiment = cp.load(filename)

# standard deviation of noise from the dataset
sigma = 7.500

# For spectral fitting, we only focus on the real part of the complex dataset
experiment = experiment.real

# Convert the coordinates along each dimension from Hz to ppm.
_ = [item.to("ppm", "nmr_frequency_ratio") for item in experiment.dimensions]

# plot of the dataset.
max_amp = experiment.max()
levels = (np.arange(29) + 1) * max_amp / 30  # contours are drawn at these levels.
options = dict(levels=levels, linewidths=0.5)  # plot options

plt.figure(figsize=(4.25, 3.0))
ax = plt.subplot(projection="csdm")
ax.contour(experiment, colors="k", **options)
ax.set_xlim(1000, -1000)
ax.set_ylim(1500, -1500)
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
    isotropic_chemical_shift=-90,  # in ppm
    shielding_symmetric=SymmetricTensor(
        zeta=-610,  # in ppm
        eta=0.15,
        alpha=0.7,  # in rads
        beta=2.0,  # in rads
        gamma=3.0,  # in rads
    ),
    quadrupolar=SymmetricTensor(Cq=75.2e3, eta=0.9),  # Cq in Hz
)

spin_systems = [SpinSystem(sites=[site])]

# %%
# **Method**
#
# Use the generic method, `Method`, to generate a shifting-d echo method. The
# reported shifting-d 2D sequence is a correlation of the shielding frequencies to the
# first-order quadrupolar frequencies. Here, we create a correlation method using the
# :attr:`~mrsimulator.method.event.freq_contrib` attribute, which acts as a switch
# for including the frequency contributions from interaction during the event.
#
# In the following method, we assign the ``["Quad1_2"]`` and
# ``["Shielding1_0", "Shielding1_2"]`` as the value to the ``freq_contrib`` key. The
# *Quad1_2* is an enumeration for selecting the first-order second-rank quadrupolar
# frequency contributions. *Shielding1_0* and *Shielding1_2* are enumerations for
# the first-order shielding with zeroth and second-rank tensor contributions,
# respectively. See :ref:`freq_contrib_api` for details.

# Get the spectral dimension parameters from the experiment.
spectral_dims = get_spectral_dimensions(experiment)

shifting_d = Method(
    channels=["2H"],
    magnetic_flux_density=9.395,  # in T
    rotor_frequency=0,
    spectral_dimensions=[
        SpectralDimension(
            **spectral_dims[0],
            label="Quadrupolar frequency",
            events=[
                SpectralEvent(
                    transition_query=[{"ch1": {"P": [-1]}}],
                    freq_contrib=["Quad1_2"],
                ),
                MixingEvent(query="NoMixing"),
            ],
        ),
        SpectralDimension(
            **spectral_dims[1],
            label="Paramagnetic shift",
            events=[
                SpectralEvent(
                    transition_query=[{"ch1": {"P": [-1]}}],
                    freq_contrib=["Shielding1_0", "Shielding1_2"],
                )
            ],
        ),
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
        sp.IFFT(dim_index=(0, 1)),
        sp.apodization.Gaussian(FWHM="5 kHz", dim_index=0),  # along dimension 0
        sp.apodization.Gaussian(FWHM="5 kHz", dim_index=1),  # along dimension 1
        sp.FFT(dim_index=(0, 1)),
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
ax.set_xlim(1000, -1000)
ax.set_ylim(1500, -1500)
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
result

# %%
# The best fit solution
# ---------------------
best_fit = sf.bestfit(sim, processor)[0].real

# Plot the spectrum
plt.figure(figsize=(4.25, 3.0))
ax = plt.subplot(projection="csdm")
ax.contour(experiment, colors="k", **options)
ax.contour(best_fit, colors="r", linestyles="--", **options)
ax.set_xlim(1000, -1000)
ax.set_ylim(1500, -1500)
plt.grid()
plt.tight_layout()
plt.show()

# %%
# Image plots with residuals
# --------------------------
residuals = sf.residuals(sim, processor)[0].real

fig, ax = plt.subplots(
    1, 3, sharey=True, figsize=(10, 3.0), subplot_kw={"projection": "csdm"}
)
vmax, vmin = experiment.max(), experiment.min()
for i, dat in enumerate([experiment, best_fit, residuals]):
    ax[i].imshow(dat, aspect="auto", cmap="gist_ncar_r", vmax=vmax, vmin=vmin)
    ax[i].set_xlim(1000, -1000)
ax[0].set_ylim(1500, -1500)
plt.tight_layout()
plt.show()
# %%
# .. [#f1] Walder B.J, Patterson A.M., Baltisberger J.H, and Grandinetti P.J
#       Hydrogen motional disorder in crystalline iron group chloride dihydrates
#       spectroscopy, J. Chem. Phys. (2018)  **149**, 084503.
#       `DOI: 10.1063/1.5037151 <https://doi.org/10.1063/1.5037151>`_
