#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
Fitting Crystalline Sodium Metasilicate
=======================================
"""
# %%
# In this example, we illustrate the use of the mrsimulator objects to
#
# - create a spin system fitting model,
# - use the fitting model to perform a least-squares fit on the experimental, and
# - extract the tensor parameters of the spin-system model.
#
# We will be using the `LMFIT <https://lmfit.github.io/lmfit-py/>`_ methods to
# establish fitting parameters and fit the spectrum. The following example illustrates
# the least-squares fitting on a :math:`^{17}\text{O}` measurement of
# :math:`\text{Na}_{2}\text{SiO}_{3}` [#f5]_.
#
# We will begin by importing relevant modules and presetting figure style and layout.
import csdmpy as cp
import matplotlib as mpl
import matplotlib.pyplot as plt
import mrsimulator.signal_processing as sp
import mrsimulator.signal_processing.apodization as apo
from lmfit import Minimizer
from lmfit import report_fit
from mrsimulator import Simulator
from mrsimulator import Site
from mrsimulator import SpinSystem
from mrsimulator.methods import BlochDecayCentralTransitionSpectrum
from mrsimulator.utils.spectral_fitting import LMFIT_min_function
from mrsimulator.utils.spectral_fitting import make_LMFIT_parameters

font = {"size": 9}
mpl.rc("font", **font)
mpl.rcParams["figure.figsize"] = [4.25, 3.0]
# sphinx_gallery_thumbnail_number = 3

# %%
# Import the dataset
# ------------------
#
# Import the experimental data. In this example, we will import the dataset file
# serialized with the CSDM file-format, using the
# `csdmpy <https://csdmpy.readthedocs.io/en/latest/index.html>`_ module.
filename = "https://osu.box.com/shared/static/kfgt0jxgy93srsye9pofdnoha6qy58qf.csdf"
oxygen_experiment = cp.load(filename)

# For spectral fitting, we only focus on the real part of the complex dataset
oxygen_experiment = oxygen_experiment.real

# Convert the dimension coordinates from Hz to ppm.
oxygen_experiment.dimensions[0].to("ppm", "nmr_frequency_ratio")

# Normalize the spectrum
oxygen_experiment /= oxygen_experiment.max()

# plot of the dataset.
ax = plt.subplot(projection="csdm")
ax.plot(oxygen_experiment, color="black", linewidth=1)
ax.set_xlim(-50, 100)
ax.invert_xaxis()
plt.tight_layout()
plt.show()

# %%
# Create a fitting model
# ----------------------
#
# Next, we will create a ``simulator`` object that we use to fit the spectrum. We will
# start by creating the guess ``SpinSystem`` objects.
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
    quadrupolar={"Cq": 2.4e6, "eta": 0.5},  # Cq in Hz
)

system_object = [SpinSystem(sites=[s]) for s in [O17_1, O17_2]]  # from the above code

# %%
# **Step 2:** Create the method object. Note, when performing the least-squares fit, you
# must create an appropriate method object which matches the method used in acquiring
# the experimental data. The attribute values of this method must match the
# exact conditions under which the experiment was acquired. This including the
# acquisition channels, the magnetic flux density, rotor angle, rotor frequency, and
# the spectral/spectroscopic dimension. In the following example, we set up a central
# transition selective Bloch decay spectrum method, where we obtain the
# spectral/spectroscopic information from the metadata of the CSDM dimension. The
# remaining attribute values are set to the experimental conditions.

# get the count, increment, and coordinates_offset info from experiment CSDM dimension.
count = oxygen_experiment.dimensions[0].count
increment = oxygen_experiment.dimensions[0].increment.to("Hz").value
offset = oxygen_experiment.dimensions[0].coordinates_offset.to("Hz").value

method = BlochDecayCentralTransitionSpectrum(
    channels=["17O"],
    magnetic_flux_density=9.4,  # in T
    rotor_frequency=14000,  # in Hz
    spectral_dimensions=[
        {
            "count": count,
            "spectral_width": count * increment,  # in Hz
            "reference_offset": offset,  # in Hz
        }
    ],
)

# %%
# Assign the experimental dataset to the ``experiment`` attribute of the above method.
method.experiment = oxygen_experiment

# %%
# **Step 3:** Create the Simulator object and add the method and spin system objects.
sim = Simulator()
sim.spin_systems = system_object
sim.methods = [method]

# %%
# **Step 4:** Simulate the spectrum.
for iso in sim.spin_systems:
    # To avoid querying for transition pathways at every iteration, evaluate it once
    # and save the transition pathways information following
    iso.transition_pathways = method.get_transition_pathways(iso).tolist()

sim.run()

# %%
# **Step 5:** Create the SignalProcessor class object and apply the post-simulation
# signal processing operations.
processor = sp.SignalProcessor(
    operations=[
        sp.IFFT(),
        apo.Exponential(FWHM="40 Hz"),
        sp.FFT(),
        sp.Scale(factor=0.6),
    ]
)
processed_data = processor.apply_operations(data=sim.methods[0].simulation)

# %%
# **Step 6:** The plot of initial guess simulation (black) along with the experiment
# (red) is shown below.
ax = plt.subplot(projection="csdm")
ax.plot(processed_data.real, color="black", linewidth=1, label="guess spectrum")
ax.plot(oxygen_experiment, c="r", linewidth=1.5, alpha=0.5, label="experiment")
ax.set_xlim(-50, 100)
ax.invert_xaxis()
plt.legend()
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
# :func:`~mrsimulator.utils.spectral_fitting.make_LMFIT_parameters`, that considerably
# simplifies the LMFIT parameters generation process.
#
# **Step 7:** Create a list of parameters.
params = make_LMFIT_parameters(sim, processor)

# %%
# The `make_LMFIT_parameters` parses the instances of the ``Simulator`` and the
# ``PostSimulator`` objects for parameters and returns an LMFIT `Parameters` object.
#
# **Customize the Parameters:**
# You may customize the parameters list, ``params``, as desired. Here, we add a
# constraint on the fit by fixing the site abundances for the spin systems at index
# 1 and 2 to 50%.
params["sys_0_abundance"].vary = False  # fix the abundance
params.pretty_print()

# %%
# **Step 8:** Perform least-squares minimization. For the user's convenience, we also
# provide a utility function,
# :func:`~mrsimulator.utils.spectral_fitting.LMFIT_min_function`, for evaluating the
# difference vector between the simulation and experiment, based on
# the parameters update. You may use this function directly as the argument of the
# LMFIT Minimizer class, as follows,
minner = Minimizer(LMFIT_min_function, params, fcn_args=(sim, processor))
result = minner.minimize()
report_fit(result)

# %%
# **Step 9:** The plot of the fit, measurement and the residuals is shown below.
plt.figsize = (4, 3)
x, y_data = oxygen_experiment.to_list()
residual = result.residual
plt.plot(x, y_data, label="Spectrum")
plt.plot(x, y_data - residual, "r", alpha=0.5, label="Fit")
plt.plot(x, residual, alpha=0.5, label="Residual")

plt.xlabel("$^{17}$O frequency / ppm")
plt.xlim(-50, 100)
plt.gca().invert_xaxis()
plt.grid(which="major", axis="both", linestyle="--")
plt.legend()
plt.tight_layout()
plt.show()

# %%
#
# .. [#f5] T. M. Clark, P. Florian, J. F. Stebbins, and P. J. Grandinetti,
#       An :math:`^{17}\text{O}` NMR Investigation of Crystalline Sodium Metasilicate:
#       Implications for the Determination of Local Structure in Alkali Silicates,
#       J. Phys. Chem. B. 2001, **105**, 12257-12265.
#       `DOI: 10.1021/jp011289p  <https://doi.org/10.1021/jp011289p>`_
