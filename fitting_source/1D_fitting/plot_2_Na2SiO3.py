#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
17O MAS NMR of crystalline Na2SiO3 (2nd order quad)
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
"""
# %%
# In this example, we illustrate the use of the mrsimulator objects to
#
# - create a quadrupolar fitting model using Simulator and SignalProcessor objects,
# - use the fitting model to perform a least-squares analysis, and
# - extract the fitting parameters from the model.
#
# We use the `LMFIT <https://lmfit.github.io/lmfit-py/>`_ library to fit the spectrum.
# The following example shows the least-squares fitting procedure applied to the
# :math:`^{17}\text{O}` MAS NMR spectrum of :math:`\text{Na}_{2}\text{SiO}_{3}` [#f5]_.
# The dataset was shared by Dr. Philip Grandinetti.
#
# Start by importing the relevant modules.
import csdmpy as cp
import matplotlib.pyplot as plt
from lmfit import Minimizer, report_fit

from mrsimulator import Simulator, SpinSystem, Site
from mrsimulator.methods import BlochDecayCTSpectrum
from mrsimulator import signal_processing as sp
from mrsimulator.utils import spectral_fitting as sf
from mrsimulator.utils import get_spectral_dimensions

# sphinx_gallery_thumbnail_number = 3

# %%
# Import the dataset
# ------------------
#
# Import the experimental data. We use dataset file serialized with the CSDM
# file-format, using the
# `csdmpy <https://csdmpy.readthedocs.io/en/stable/index.html>`_ module.
filename = "https://sandbox.zenodo.org/record/814455/files/Na2SiO3_O17.csdf"
experiment = cp.load(filename)

# standard deviation of noise from the dataset
sigma = 1.931335

# For spectral fitting, we only focus on the real part of the complex dataset
experiment = experiment.real

# Convert the dimension coordinates from Hz to ppm.
experiment.x[0].to("ppm", "nmr_frequency_ratio")

# plot of the dataset.
plt.figure(figsize=(4.25, 3.0))
ax = plt.subplot(projection="csdm")
ax.plot(experiment, "k", alpha=0.5)
ax.set_xlim(100, -50)
plt.grid()
plt.tight_layout()
plt.show()

# %%
# Create a fitting model
# ----------------------
#
# A fitting model is a composite of ``Simulator`` and ``SignalProcessor`` objects.
#
# **Step 1:** Create initial guess sites and spin systems
O17_1 = Site(
    isotope="17O",
    isotropic_chemical_shift=60.0,  # in ppm,
    quadrupolar={"Cq": 4.2e6, "eta": 0.5},  # Cq in Hz
)

O17_2 = Site(
    isotope="17O",
    isotropic_chemical_shift=40.0,  # in ppm,
    quadrupolar={"Cq": 2.4e6, "eta": 0},  # Cq in Hz
)

spin_systems = [SpinSystem(sites=[s], abundance=50) for s in [O17_1, O17_2]]

# %%
# **Step 2:** Create the method object. Create an appropriate method object that closely
# resembles the technique used in acquiring the experimental data. The attribute values
# of this method must meet the experimental conditions, including the acquisition
# channels, the magnetic flux density, rotor angle, rotor frequency, and the
# spectral/spectroscopic dimension.
#
# In the following example, we set up a central transition selective Bloch decay
# spectrum method where the spectral/spectroscopic dimension information, i.e., count,
# spectral_width, and the reference_offset, is extracted from the CSDM dimension
# metadata using the :func:`~mrsimulator.utils.get_spectral_dimensions` utility
# function. The remaining attribute values are set to the experimental conditions.

# get the count, spectral_width, and reference_offset information from the experiment.
spectral_dims = get_spectral_dimensions(experiment)

method = BlochDecayCTSpectrum(
    channels=["17O"],
    magnetic_flux_density=9.395,  # in T
    rotor_frequency=14000,  # in Hz
    spectral_dimensions=spectral_dims,
    experiment=experiment,  # experimental dataset
)

# A method object queries every spin system for a list of transition pathways that are
# relevant for the given method. Since the method and the number of spin systems remain
# the same during the least-squares fit, a one-time query is sufficient. To avoid
# querying for the transition pathways at every iteration in a least-squares fitting,
# evaluate the transition pathways once and store it as follows
for sys in spin_systems:
    sys.transition_pathways = method.get_transition_pathways(sys)

# %%
# **Step 3:** Create the Simulator object and add the method and spin system objects.
sim = Simulator(spin_systems=spin_systems, methods=[method])
sim.run()

# %%
# **Step 4:** Create a SignalProcessor class object and apply the post-simulation
# signal processing operations.
processor = sp.SignalProcessor(
    operations=[
        sp.IFFT(),
        sp.apodization.Gaussian(FWHM="100 Hz"),
        sp.FFT(),
        sp.Scale(factor=200.0),
    ]
)
processed_data = processor.apply_operations(data=sim.methods[0].simulation).real

# %%
# **Step 5:** The plot of the data and the guess spectrum.
plt.figure(figsize=(4.25, 3.0))
ax = plt.subplot(projection="csdm")
ax.plot(experiment, "k", linewidth=1, label="Experiment")
ax.plot(processed_data, "r", alpha=0.75, linewidth=1, label="guess spectrum")
ax.set_xlim(100, -50)
plt.legend()
plt.grid()
plt.tight_layout()
plt.show()


# %%
# Least-squares minimization with LMFIT
# -------------------------------------
#
# Once you have a fitting model, you need to create the list of parameters to use in the
# least-squares fitting. For this, you may use the
# `Parameters <https://lmfit.github.io/lmfit-py/parameters.html>`_ class from *LMFIT*,
# as described in the previous example.
# Here, we make use of a utility function,
# :func:`~mrsimulator.utils.spectral_fitting.make_LMFIT_params`, that considerably
# simplifies the LMFIT parameters generation process.
#
# **Step 6:** Create a list of parameters.
params = sf.make_LMFIT_params(sim, processor)

# %%
# The `make_LMFIT_params` parses the instances of the ``Simulator`` and the
# ``PostSimulator`` objects for parameters and returns an LMFIT `Parameters` object.
#
# **Customize the Parameters:**
# You may customize the parameters list, ``params``, as desired. Here, we remove the
# abundance of the two spin systems and constrain it to the initial value of 50% each.
params.pop("sys_0_abundance")
params.pop("sys_1_abundance")
print(params.pretty_print(columns=["value", "min", "max", "vary", "expr"]))

# %%
# **Step 7:** Perform least-squares minimization. For the user's convenience, we also
# provide a utility function,
# :func:`~mrsimulator.utils.spectral_fitting.LMFIT_min_function`, for evaluating the
# difference vector between the simulation and experiment, based on
# the parameters update. You may use this function directly as the argument of the
# LMFIT Minimizer class, as follows,
minner = Minimizer(sf.LMFIT_min_function, params, fcn_args=(sim, processor, sigma))
result = minner.minimize()
report_fit(result)

# %%
# **Step 8:** The plot of the fit and the measurement data.

# Best fit spectrum
best_fit = sf.bestfit(sim, processor)[0]
residuals = sf.residuals(sim, processor)[0]

plt.figure(figsize=(4.25, 3.0))
ax = plt.subplot(projection="csdm")
ax.plot(experiment, "k", linewidth=1, label="Experiment")
ax.plot(best_fit, "r", alpha=0.75, linewidth=1, label="Best Fit")
ax.plot(residuals, alpha=0.75, linewidth=1, label="Residuals")
ax.set_xlabel("$^{17}$O frequency / ppm")
ax.set_xlim(100, -50)
plt.legend()
plt.grid()
plt.tight_layout()
plt.show()

# %%
#
# .. [#f5] T. M. Clark, P. Florian, J. F. Stebbins, and P. J. Grandinetti,
#       An :math:`^{17}\text{O}` NMR Investigation of Crystalline Sodium Metasilicate:
#       Implications for the Determination of Local Structure in Alkali Silicates,
#       J. Phys. Chem. B. 2001, **105**, 12257-12265.
#       `DOI: 10.1021/jp011289p  <https://doi.org/10.1021/jp011289p>`_
