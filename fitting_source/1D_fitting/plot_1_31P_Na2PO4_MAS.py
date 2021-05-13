#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
31P MAS NMR of crystalline Na2PO4 (CSA)
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
"""
# %%
# In this example, we illustrate the use of the mrsimulator objects to
#
# - create a CSA fitting model using Simulator and SignalProcessor objects,
# - use the fitting model to perform a least-squares analysis, and
# - extract the fitting parameters from the model.
#
# We use the `LMFIT <https://lmfit.github.io/lmfit-py/>`_ library to fit the spectrum.
# The following example shows the least-squares fitting procedure applied to the
# :math:`^{31}\text{P}` MAS NMR spectrum of :math:`\text{Na}_{2}\text{PO}_{4}`.
# The following experimental dataset is a part of DMFIT [#f1]_ examples, and we
# acknowledge Dr. Dominique Massiot for sharing the dataset.
#
# Start by importing the relevant modules.
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
#
# Import the experimental data. We use dataset file serialized with the CSDM
# file-format, using the
# `csdmpy <https://csdmpy.readthedocs.io/en/stable/index.html>`_ module.
host = "https://nmr.cemhti.cnrs-orleans.fr/Dmfit/Help/csdm/"
filename = "31P Phosphate 6kHz.csdf"
experiment = cp.load(host + filename)

# standard deviation of noise from the dataset
sigma = 1.523217

# For spectral fitting, we only focus on the real part of the complex dataset
experiment = experiment.real

# Convert the dimension coordinates from Hz to ppm.
experiment.x[0].to("ppm", "nmr_frequency_ratio")

# plot of the dataset.
plt.figure(figsize=(4.25, 3.0))
ax = plt.subplot(projection="csdm")
ax.plot(experiment, color="black", linewidth=0.5, label="Experiment")
ax.set_xlim(150, -150)
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
P_31 = Site(
    isotope="31P",
    isotropic_chemical_shift=5.0,  # in ppm,
    shielding_symmetric={"zeta": -80, "eta": 0.5},  # zeta in Hz
)

spin_systems = [SpinSystem(sites=[P_31])]

# %%
# **Step 2:** Create the method object. Create an appropriate method object that closely
# resembles the technique used in acquiring the experimental data. The attribute values
# of this method must meet the experimental conditions, including the acquisition
# channels, the magnetic flux density, rotor angle, rotor frequency, and the
# spectral/spectroscopic dimension.
#
# In the following example, we set up a Bloch decay spectrum method where the
# spectral/spectroscopic dimension information, i.e., count, spectral_width, and the
# reference_offset, is extracted from the CSDM dimension metadata using the
# :func:`~mrsimulator.utils.get_spectral_dimensions` utility function. The remaining
# attribute values are set to the experimental conditions.
#

# get the count, spectral_width, and reference_offset information from the experiment.
spectral_dims = get_spectral_dimensions(experiment)

method = BlochDecaySpectrum(
    channels=["31P"],
    magnetic_flux_density=9.395,  # in T
    rotor_frequency=6000,  # in Hz
    spectral_dimensions=spectral_dims,
    experiment=experiment,  # experimental dataset
)

# A method object queries every spin system for a list of transition pathways that are
# relevant to the given method. Since the method and the number of spin systems remains
# unchanged during the least-squares analysis, a one-time query is sufficient. To avoid
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
        sp.apodization.Exponential(FWHM="0.3 kHz"),
        sp.FFT(),
        sp.Scale(factor=300),
        sp.baseline.ConstantOffset(offset=-2),
    ]
)
processed_data = processor.apply_operations(data=sim.methods[0].simulation).real

# %%
# **Step 5:** The plot of the data and the guess spectrum.
plt.figure(figsize=(4.25, 3.0))
ax = plt.subplot(projection="csdm")
ax.plot(experiment, color="black", linewidth=0.5, label="Experiment")
ax.plot(processed_data, linewidth=2, alpha=0.6, label="Guess Spectrum")
ax.set_xlim(150, -150)
plt.legend()
plt.grid()
plt.tight_layout()
plt.show()


# %%
# Least-squares minimization with LMFIT
# -------------------------------------
#
# Once you have a fitting model, you need to create the list of parameters to use in the
# least-squares minimization. For this, you may use the
# `Parameters <https://lmfit.github.io/lmfit-py/parameters.html>`_ class from *LMFIT*,
# as described in the previous example.
# Here, we make use of a utility function,
# :func:`~mrsimulator.utils.spectral_fitting.make_LMFIT_params`, to simplifies the
# LMFIT parameters generation process. By default, the function only creates parameters
# from the SpinSystem and SignalProcessor objects. Often, in spectrum with sidebands,
# spinning speed may not be accurately known; and is, therefore, included as a fitting
# parameter. To include a keyword from the method object, use the *include* argument
# of the function, as follows,
#
# **Step 6:** Create a list of parameters.
params = sf.make_LMFIT_params(sim, processor, include={"rotor_frequency"})

# %%
# The `make_LMFIT_params` parses the instances of the ``Simulator`` and the
# ``PostSimulator`` objects for parameters and returns a LMFIT `Parameters` object.
#
# **Customize the Parameters:**
# You may customize the parameters list, ``params``, as desired. Here, we remove the
# abundance parameter.
params.pop("sys_0_abundance")
print(params.pretty_print(columns=["value", "min", "max", "vary", "expr"]))

# %%
# **Step 7:** Perform the least-squares minimization. For the user's convenience, we
# also provide a utility function,
# :func:`~mrsimulator.utils.spectral_fitting.LMFIT_min_function`, for evaluating the
# difference vector between the simulation and experiment, based on
# the parameters update. You may use this function directly as the argument of the
# LMFIT Minimizer class, as follows,
minner = Minimizer(sf.LMFIT_min_function, params, fcn_args=(sim, processor, sigma))
result = minner.minimize()
report_fit(result)

# %%
# **Step 8:** The plot of the fit, measurement, and residuals.

# Best fit spectrum
best_fit = sf.bestfit(sim, processor)[0]
residuals = sf.residuals(sim, processor)[0]

plt.figure(figsize=(4.25, 3.0))
ax = plt.subplot(projection="csdm")
ax.plot(experiment, color="black", linewidth=0.5, label="Experiment")
ax.plot(residuals, color="gray", linewidth=0.5, label="Residual")
ax.plot(best_fit, linewidth=2, alpha=0.6, label="Best Fit")
ax.set_xlabel(r"$^{31}$P frequency / ppm")
ax.set_xlim(150, -150)
plt.legend()
plt.grid()
plt.tight_layout()
plt.show()

# %%
#
# .. [#f1] D.Massiot, F.Fayon, M.Capron, I.King, S.Le Calvé, B.Alonso, J.O.Durand,
#       B.Bujoli, Z.Gan, G.Hoatson, 'Modelling one and two-dimensional solid-state NMR
#       spectra.', Magn. Reson. Chem. **40** 70-76 (2002)
#       `DOI: 10.1002/mrc.984 <https://doi.org/10.1002/mrc.984>`_
