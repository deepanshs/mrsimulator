#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
29Si 1D MAS spinning sideband (Xonotlite)
=========================================
"""
# %%
# The following is an example for submitting the NMR tensor parameters to mpcontribs.
# We use the :math:`^{29}\text{Si}` 1D MAS NMR spectrum of Xonotlite crystal by Hansen
# et al. [#f1]_ for demonstration.
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
filename = "https://sandbox.zenodo.org/record/744498/files/xonotlite.csdf"
experiment = cp.load(filename).real

# standard deviation of noise from the dataset
sigma = 2.819601

# Convert the coordinates along each dimension from Hz to ppm.
_ = [item.to("ppm", "nmr_frequency_ratio") for item in experiment.dimensions]

# Plot of the synthetic dataset.
plt.figure(figsize=(4.25, 3.0))
ax = plt.subplot(projection="csdm")
ax.plot(experiment, "k", alpha=0.5)
ax.invert_xaxis()
plt.tight_layout()
plt.show()

# %%
# Create a fitting model
# ----------------------
# **Guess model**
#
# Create a guess list of spin systems. There are three crystallographic
# :math:`^{29}\text{Si}` sites in Xonotlite.
s1 = Site(
    isotope="29Si",
    isotropic_chemical_shift=-97.17,  # in ppm,
    shielding_symmetric={"zeta": 35.0, "eta": 0.0},  # zeta in ppm
)
s2 = Site(
    isotope="29Si",
    isotropic_chemical_shift=-86.3,  # in ppm,
    shielding_symmetric={"zeta": 50.0, "eta": 0.5},  # zeta in ppm
)
s3 = Site(
    isotope="29Si",
    isotropic_chemical_shift=-87.2,  # in ppm,
    shielding_symmetric={"zeta": 44.0, "eta": 0.5},  # zeta in ppm
)
spin_systems = [
    SpinSystem(name="Q3", sites=[s1], abundance=25),
    SpinSystem(name="Q2 (1)", sites=[s2], abundance=75 / 2),
    SpinSystem(name="Q2 (2)", sites=[s3], abundance=75 / 2),
]

# %%
# **Method**

# Get the spectral dimension paramters from the experiment.
spectral_dims = get_spectral_dimensions(experiment)

method = BlochDecaySpectrum(
    channels=["29Si"],
    magnetic_flux_density=14.1,  # in T
    rotor_frequency=1800.0,  # in Hz
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
sim = Simulator(spin_systems=spin_systems, methods=[method])
sim.run()

# Post Simulation Processing
# --------------------------
processor = sp.SignalProcessor(
    operations=[
        sp.IFFT(),  # inverse FFT to convert frequency based spectrum to time domain.
        sp.apodization.Exponential(FWHM="50 Hz"),  # apodization of time domain signal.
        sp.FFT(),  # forward FFT to convert time domain signal to frequency spectrum.
        sp.Scale(factor=500),  # scale the frequency spectrum.
    ]
)
processed_data = processor.apply_operations(data=sim.methods[0].simulation).real

# Plot of the guess Spectrum
# --------------------------
plt.figure(figsize=(4.25, 3.0))
ax = plt.subplot(projection="csdm")
ax.plot(experiment, "k", linewidth=1, label="Experiment")
ax.plot(processed_data, "r", alpha=0.75, linewidth=1, label="guess spectrum")
ax.invert_xaxis()
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

params.pop("sys_0_abundance")
params.pop("sys_1_abundance")
params.pop("sys_2_abundance")
params["sys_0_site_0_shielding_symmetric_eta"].vary = False
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
plt.figure(figsize=(4.25, 3.0))
ax = plt.subplot(projection="csdm")
ax.plot(experiment, "k", linewidth=1, label="Experiment")
ax.plot(best_fit, "r", alpha=0.75, linewidth=1, label="Best Fit")
ax.plot(residuals, alpha=0.75, linewidth=1, label="Residuals")
ax.invert_xaxis()
plt.xlabel("$^{29}$Si frequency / ppm")
plt.grid()
plt.legend()
plt.tight_layout()
plt.show()


# %%
# Submitting data to MPContribs
# -----------------------------
#
# To contribute to MPContribs, we need to export the mrsimulator objects to a list of
# mp-compatible data dictionaries. At present, MPContribs only support data contribution
# on a per NMR site basis and, therefore, we only generate mp contributions for
# uncoupled spin systems. Use the ``mrsimulator.contribs`` module to create data
# dictionaries as follows.

# %%
# from mrsimulator.contribs import mpcontribs_export
# from pprint import pprint

# mp_project = "lsdi_nmr_exp_test"  # this should be your mpcontribs project name
# cards = mpcontribs_export(
#     sim,
#     [processor],
#     project=mp_project,
#     identifier="Ca6Si6O17(OH)2",
#     exp_dict={
#         "90degreePulseLength": "6 µs",
#         "relaxationDelay": "8 s",
#         "numberOfScans": 7224,
#         "referenceCompound": "TMS",
#     },
# )
# print("Number of contributions", len(cards))
# pprint(cards[0]["data"])

# %%
# Here, ``cards`` hold a list of mp-data dictionaries. In this example, it corresponds
# to three---the number of uncoupled spin systems.
# To submit contributions, use the mpcontribs client as shown below.

# from mpcontribs.client import Client
#
# client = Client(<YOUR-API-KEY>)  # uses MPCONTRIBS_API_KEY envvar.
# client.submit_contributions(cards)

# %%
# .. [#f1] Hansen, M. R., Jakobsen, H. J., Skibsted, J., :math:`^{29}\text{Si}`
#       Chemical Shift Anisotropies in Calcium Silicates from High-Field
#       :math:`^{29}\text{Si}` MAS NMR Spectroscopy, Inorg. Chem. 2003,
#       **42**, *7*, 2368-2377.
#       `DOI: 10.1021/ic020647f <https://doi.org/10.1021/ic020647f>`_
