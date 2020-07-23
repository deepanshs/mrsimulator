#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
Fitting Crystalline Sodium Metasilicate
=======================================
.. sectionauthor:: Maxwell C. Venetos <maxvenetos@gmail.com>
"""
# %%
# Often, after obtaining an NMR measurement we must fit tensors to our data so we can
# obtain the tensor parameters. In this example, we will illustrate the use of the
# *mrsimulator* method to simulate the experimental spectrum and fit the simulation to
# the data allowing us to extract the tensor parameters for our spin systems. We will
# be using the `LMFIT <https://lmfit.github.io/lmfit-py/>`_ methods to establish
# fitting parameters and fit the spectrum. The following examples will show fitting with
# two synthetic :math:`^{29}\text{Si}` spectra--cuspidine and wollastonite--as well as
# the measurements from an :math:`^{17}\text{O}` experiment on
# :math:`\text{Na}_{2}\text{SiO}_{3}`. The *mrsimulator* library and data make use of
# CSDM compliant files. In this example we will fit a simulation to an experimentally
# obtained :math:`^{17}\text{O}` spectrum.
# We use the :math:`^{17}\text{O}` tensor information from Grandinetti `et. al.` [#f5]_
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
from mrsimulator.spectral_fitting import LMFIT_min_function
from mrsimulator.spectral_fitting import make_LMFIT_parameters

font = {"size": 9}
mpl.rc("font", **font)
mpl.rcParams["figure.figsize"] = [4.25, 3.0]
# sphinx_gallery_thumbnail_number = 3

# %%
# Import the dataset
# ------------------
#
# **Step 0**
# Import the experimental data. In this example, we will import the data file serialized
# with the CSDM file-format. We use the
# `csdmpy <https://csdmpy.readthedocs.io/en/latest/index.html>`_ module to load the
# dataset.
filename = "https://osu.box.com/shared/static/kfgt0jxgy93srsye9pofdnoha6qy58qf.csdf"
oxygen_experiment = cp.load(filename)

# For spectral fitting, we only focus on the real part of the complex dataset
oxygen_experiment = oxygen_experiment.real

# Convert the dimension coordinates from Hz to ppm.
oxygen_experiment.dimensions[0].to("ppm", "nmr_frequency_ratio")

# plot of the dataset.
ax = plt.subplot(projection="csdm")
ax.plot(oxygen_experiment, color="black", linewidth=1)
ax.invert_xaxis()
plt.tight_layout()
plt.show()

# %%
# Create a "fitting model"
# ------------------------
#
# Next, we will want to create a ``simulator`` object that we will use to fit to our
# spectrum. We will start by creating a guess ``SpinSystem`` objects.
#
# **Step 1** Create the guess sites.
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

# %%
# **Step 2** Create the spin systems for the guess sites.
system_object = [SpinSystem(sites=[s]) for s in [O17_1, O17_2]]  # from the above code

# %%
# **Step 3** Create the method object. Note, when fitting, you must create an
# appropriate method object which matches the method used in acquiring the experimental
# data. The attribute values of this method must be set to match the exact conditions
# under which the experiment was acquired. This including the acquisition channels, the
# magnetic flux density, rotor angle, rotor frequency, and the spectral/spectroscopic
# dimension. In the following example, we set up a central transition selective Bloch
# decay spectrum method, where the spectral/spectroscopic information is obtained from
# the dimension metadata of the dimension. The remaining attribute values are set to
# the experimental conditions.

# get the count, increment, and coordinates_offset info from the experiment dimension.
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
# Now we add the experimental data to the ``experiment`` attribute of the above method.
method.experiment = oxygen_experiment

# %%
# **Step 4** Create the Simulator object and add the method and spin-system objects.
sim = Simulator()
sim.spin_systems = system_object
sim.methods = [method]

# %%
# **Step 5** Simulate the spectrum.
for iso in sim.spin_systems:
    # To avoid querying at every iteration we will save the relevant transition pathways
    iso.transition_pathways = method.get_transition_pathways(iso).tolist()
sim.run()

# %%
# **Step 6** Create the SignalProcessor class object and add to it the list of
# operations that will be applied to the simulated spectrum. The generic list of
# operations in NMR includes the line-broadening function. In the following, we add a
# scaling and a Lorentzian line-broadening function. Here, the Lorentzian
# line-broadening is defined as an exponential apodization operation sandwiched
# between two Fourier transformations.
factor = oxygen_experiment.max()
op_list = [
    sp.IFFT(),
    apo.Exponential(Lambda=100),
    sp.FFT(),
    sp.Scale(factor=factor),
]
post_sim = sp.SignalProcessor(operations=op_list)

# %%
# **Step 7** Process and plot the spectrum.
processed_data = post_sim.apply_operations(data=sim.methods[0].simulation)

ax = plt.subplot(projection="csdm")
ax.plot(processed_data.real, color="black", linewidth=1)
ax.invert_xaxis()
plt.tight_layout()
plt.show()

# %%
# Least-squares minimization with LMFIT
# -------------------------------------
#
# Once we have our simulation we must create our list of parameters to use in our
# fitting. We will be using the
# `Parameters <https://lmfit.github.io/lmfit-py/parameters.html>`_ class from *LMFIT*.
#
# In each experiment the number of spin systems and sites present as well as their
# attributes may vary. To simplify the parameter list creation we will use the
# :func:`~mrsimulator.spectral_fitting.make_LMFIT_parameters`
#
# **Step 8** Create a list of parameters to vary during fitting.
params = make_LMFIT_parameters(sim, post_sim)

# %%
# **Customize the Parameters:**
# You may customize the parameters list, ``params``, as desired. Here, we add a
# constraint on the fit by fixing the site abundances for the spin systems at index
# 1 and 2 to 50%.
params["ISO_0_abundance"].vary = False  # fix the abundance
params.pretty_print()

# %%
# **Step 9** Perform the minimization.
minner = Minimizer(LMFIT_min_function, params, fcn_args=(sim, post_sim))
result = minner.minimize()
report_fit(result)

# %%
# **Step 10** Plot the fitted spectrum.
plt.figsize = (4, 3)
residual = oxygen_experiment.copy()
residual[:] = result.residual
plt.plot(*oxygen_experiment.to_list(), label="Spectrum")
plt.plot(*(oxygen_experiment - residual).to_list(), "r", alpha=0.5, label="Fit")
plt.plot(*residual.to_list(), alpha=0.5, label="Residual")

plt.xlabel("$^{17}$O frequency / ppm")
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
