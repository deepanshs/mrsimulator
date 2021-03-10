#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
29Si 1D MAS spinning sideband (Xonotlite)
=========================================
"""
# %%
# The following is a spinning sideband fitting example for :math:`^{29}\text{Si}` 1D
# MAS NMR spectrum of Xonotlite crystal, acquired by by Hansen et al. [#f1]_
import csdmpy as cp
import matplotlib as mpl
import matplotlib.pyplot as plt
import mrsimulator.signal_processing as sp
import mrsimulator.signal_processing.apodization as apo
from mrsimulator import Simulator, SpinSystem, Site
from mrsimulator.methods import BlochDecaySpectrum
from lmfit import Minimizer, fit_report
from mrsimulator.utils import get_spectral_dimensions
from mrsimulator.utils.spectral_fitting import LMFIT_min_function, make_LMFIT_params

font = {"size": 9}
mpl.rc("font", **font)
mpl.rcParams["figure.figsize"] = [4.5, 3.0]
mpl.rcParams["grid.linestyle"] = "--"
# sphinx_gallery_thumbnail_number = 3

# %%
# Import the dataset
# ------------------
# Use the `csdmpy <https://csdmpy.readthedocs.io/en/stable/index.html>`_
# module to load the synthetic dataset as a CSDM object.
filename = "xonotlite.csdf"
exp_data = cp.load(filename).real

# standard deviation of noise from the dataset
sigma = 2.819601

# convert the dimension coordinates from Hz to ppm
exp_data.dimensions[0].to("ppm", "nmr_frequency_ratio")

# Normalize the spectrum
max_amp = exp_data.max()
exp_data /= max_amp
sigma /= max_amp

# Plot of the synthetic dataset.
ax = plt.subplot(projection="csdm")
ax.plot(exp_data, "k", alpha=0.5)
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
    shielding_symmetric={"zeta": 33, "eta": 0.01},  # zeta in ppm
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
spectral_dims = get_spectral_dimensions(exp_data)

method = BlochDecaySpectrum(
    channels=["29Si"],
    magnetic_flux_density=14.1,  # in T
    rotor_frequency=1800.0,  # in Hz
    spectral_dimensions=spectral_dims,
    experiment=exp_data,  # add the measurement to the method.
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
sim.spin_systems = spin_systems
sim.methods = [method]
sim.run()

# Post Simulation Processing
# --------------------------
processor = sp.SignalProcessor(
    operations=[
        sp.IFFT(),  # inverse FFT to convert frequency based spectrum to time domain.
        apo.Exponential(FWHM="50 Hz"),  # apodization of time domain signal.
        sp.FFT(),  # forward FFT to convert time domain signal to frequency spectrum.
        sp.Scale(factor=0.6),  # scale the frequency spectrum.
    ]
)
processed_data = processor.apply_operations(data=sim.methods[0].simulation).real

# Plot of the guess Spectrum
# --------------------------
ax = plt.subplot(projection="csdm")
ax.plot(exp_data, "k", linewidth=1, label="Experiment")
ax.plot(processed_data, "r", alpha=0.5, linewidth=2.5, label="guess spectrum")
ax.invert_xaxis()
plt.legend()
plt.grid()
plt.tight_layout()
plt.show()

# %%
# Least-squares minimization with LMFIT
# -------------------------------------
# Use the :func:`~mrsimulator.utils.spectral_fitting.make_LMFIT_params` for a quick
# setup of the fitting parameters.
params = make_LMFIT_params(sim, processor)

params.pop("sys_0_abundance")
params.pop("sys_1_abundance")
params.pop("sys_2_abundance")
params["sys_0_site_0_shielding_symmetric_eta"].vary = False
print(params.pretty_print(columns=["value", "min", "max", "vary", "expr"]))

# %%
# **Solve the minimizer using LMFIT**
minner = Minimizer(LMFIT_min_function, params, fcn_args=(sim, processor, sigma))
result = minner.minimize()
print(fit_report(result))

# %%
# The best fit solution
# ---------------------
sim.run()
processed_data = processor.apply_operations(data=sim.methods[0].simulation).real

# Plot the spectrum
ax = plt.subplot(projection="csdm")
plt.plot(exp_data, "k", linewidth=1, label="Experiment")
plt.plot(processed_data, "r", alpha=0.5, linewidth=2.5, label="Best Fit")
plt.xlabel("$^{17}$O frequency / ppm")
plt.legend()
plt.grid()
plt.tight_layout()
plt.show()


# %%
# Mpcontribs export
# -----------------
#
# Export the site data as Mpcontribs card.
from mrsimulator.contribs import mpcontribs_export
from pprint import pprint

cards = mpcontribs_export(
    sim,
    project="lsdi_nmr_exp_test",
    identifier="Ca6Si6O17(OH)2",
    exp_dict={
        "90degreePulseLength": "6 µs",
        "relaxationDelay": "8 s",
        "numberOfScans": 7224,
        "referenceCompound": "TMS",
    },
)
pprint(cards[0])

# %%
# .. [#f1] Hansen, M. R., Jakobsen, H. J., Skibsted, J., :math:`^{29}\text{Si}`
#       Chemical Shift Anisotropies in Calcium Silicates from High-Field
#       :math:`^{29}\text{Si}` MAS NMR Spectroscopy, Inorg. Chem. 2003,
#       **42**, *7*, 2368-2377.
#       `DOI: 10.1021/ic020647f <https://doi.org/10.1021/ic020647f>`_
