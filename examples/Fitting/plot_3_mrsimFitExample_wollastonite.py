#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
Fitting to multiple spectra.
^^^^^^^^^^^^^^^^^^^^^^^^^^^^
"""
# sphinx_gallery_thumbnail_number = 3
#%%
# Often, after obtaining an NMR measurement we must fit tensors to our data so we can
# obtain the tensor parameters. In this example, we will illustrate the use of the *mrsimulator*
# method to simulate the experimental spectrum and fit the simulation to the data allowing us to
# extract the tensor parameters for our isotopomers. We will be using the `LMFIT <https://lmfit.github.io/lmfit-py/>`_
# methods to establish fitting parameters and fit the spectrum. The following examples will show fitting with
# two synthetic :math:`^{29}\text{Si}` spectra--cuspidine and wollastonite--as well as the
# measurements from an :math:`^{17}\text{O}` experiment on :math:`\text{Na}_{2}\text{SiO}_{3}`.
# The *mrsimulator* library and data make use of CSDM compliant files.
# In this example we will be creating a synthetic spectrum of cuspidine from reported tensor
# parameters and then fit a simulation to the spectrum to demonstrate a simple fitting procedure.
# The :math:`^{29}\text{Si}` tensor parameters were obtained from Hansen et. al. [#f1]_
#
# We will begin by importing *matplotlib* and establishing figure size.
import matplotlib.pylab as pylab
import matplotlib.pyplot as plt

params = {"figure.figsize": (4.5, 3), "font.size": 9}
pylab.rcParams.update(params)

#%%
# Next we will import `csdmpy <https://csdmpy.readthedocs.io/en/latest/index.html>`_ and loading the data files.
import csdmpy as cp

filename_1500Hz = (
    "https://osu.box.com/shared/static/222qlg06fuxanw6vleiybqsm7vpxxhd0.csdf"
)
filename_6000Hz = (
    "https://osu.box.com/shared/static/remqy02f9oikd5dx6ilus3puzlpsfnow.csdf"
)

wollastonite_1 = cp.load(filename_1500Hz).real
wollastonite_2 = cp.load(filename_6000Hz).real


wollastonite_1.dimensions[0].to("ppm", "nmr_frequency_ratio")
wollastonite_2.dimensions[0].to("ppm", "nmr_frequency_ratio")

x1, y1 = wollastonite_1.to_list()
x2, y2 = wollastonite_2.to_list()

fig, ax = plt.subplots(2, 1, figsize=(6, 7))
ax[0].plot(x1, y1)
ax[0].set_xlabel("$^{29}$Si frequency / ppm")
ax[0].set_xlim(x1.value.max(), x1.value.min())
ax[0].grid(color="gray", linestyle="--", linewidth=0.5, alpha=0.5)

ax[1].plot(x2, y2)
ax[1].set_xlabel("$^{29}$Si frequency / ppm")
ax[1].set_xlim(x2.value.max(), x2.value.min())
ax[1].grid(color="gray", linestyle="--", linewidth=0.5, alpha=0.5)

plt.tight_layout()
plt.show()
#%%
# In order to to fit a simulation to the data we will need to establish a ``Simulator`` object. We will
# use approximate initial parameters to generate our simulation:

#%%


from mrsimulator import Simulator, Isotopomer, Site
from mrsimulator.methods import BlochDecaySpectrum
from mrsimulator.post_simulation import PostSimulator

sim = Simulator()

method_1 = BlochDecaySpectrum(
    channels=["29Si"],
    magnetic_flux_density=14.1,  # in T
    rotor_frequency=1500,  # in Hz
    spectral_dimensions=[
        {"count": 2046, "spectral_width": 25000, "reference_offset": -10000}
    ],
)

method_2 = BlochDecaySpectrum(
    channels=["29Si"],
    magnetic_flux_density=14.1,  # in T
    rotor_frequency=6000,  # in Hz
    spectral_dimensions=[
        {"count": 2046, "spectral_width": 25000, "reference_offset": -10000}
    ],
)

PS = PostSimulator(
    scale=1, apodization=[{"args": [50], "function": "Lorentzian", "dimension": 0}]
)

sim.methods = [method_1, method_2]
sim.methods[0].post_simulation = PS
sim.methods[1].post_simulation = PS

sim.methods[0].experiment = wollastonite_1
sim.methods[1].experiment = wollastonite_2

#%%
# Next, we will need a list of parameters that will be used in the fit. the *LMFIT* library allows us to create
# a list of parameters rather easily using the `Parameters <https://lmfit.github.io/lmfit-py/parameters.html>`_ class.
# We have created a function to parse the
# ``simulator`` object for available parameters and construct an *LMFIT* ``Parameter`` object which is shown in the next two
# examples on fitting. Here, however, we will construct the parameter list explicitly to demonstrate how the parameters
# are created.

#%%


S29_1 = Site(
    isotope="29Si",
    isotropic_chemical_shift=-89.0,
    shielding_symmetric={"zeta": 60, "eta": 0.6},
)
S29_2 = Site(
    isotope="29Si",
    isotropic_chemical_shift=-90,
    shielding_symmetric={"zeta": 50, "eta": 0.7},
)
S29_3 = Site(
    isotope="29Si",
    isotropic_chemical_shift=-88,
    shielding_symmetric={"zeta": 70, "eta": 0.60},
)

isotopomers = [Isotopomer(sites=[site]) for site in [S29_1, S29_2, S29_3]]

sim.isotopomers = isotopomers

sim.run()

sim.methods[0].simulation.dimensions[0].to("ppm", "nmr_frequency_ratio")
sim.methods[1].simulation.dimensions[0].to("ppm", "nmr_frequency_ratio")

sim_x_1, sim_y_1 = sim.methods[0].simulation.to_list()
sim_x_2, sim_y_2 = sim.methods[1].simulation.to_list()

fig, ax = plt.subplots(2, 1, figsize=(6, 7))
ax[0].plot(sim_x_1, sim_y_1)
ax[0].set_xlabel("$^{29}$Si frequency / ppm")
ax[0].set_xlim(sim_x_1.value.max(), sim_x_1.value.min())
ax[0].grid(color="gray", linestyle="--", linewidth=0.5, alpha=0.5)

ax[1].plot(sim_x_2, sim_y_2)
ax[1].set_xlabel("$^{29}$Si frequency / ppm")
ax[1].set_xlim(sim_x_2.value.max(), sim_x_2.value.min())
ax[1].grid(color="gray", linestyle="--", linewidth=0.5, alpha=0.5)

plt.tight_layout()
plt.show()

#%%
# Next, we will need a list of parameters that will be used in the fit. the *LMFIT* library allows us to create
# a list of parameters rather easily using the `Parameters <https://lmfit.github.io/lmfit-py/parameters.html>`_ class.
# We have created a function to parse the
# ``simulator`` object for available parameters and construct an *LMFIT* ``Parameter`` object which is shown in the next two
# examples on fitting. Here, however, we will construct the parameter list explicitly to demonstrate how the parameters
# are created.

#%%


from lmfit import Minimizer, report_fit
from mrsimulator import spectral_fitting

params = spectral_fitting.make_fitting_parameters(sim)
params.pretty_print()


#%%
# We will next set up an error function that will update the simulation throughout the minimization.
# We will construct a simple function here to demonstrate the *LMFIT* library, however, the next examples
# will showcase a fitting function provided in the *mrsimulator* library which automates the process.
minner = Minimizer(spectral_fitting.min_function, params, fcn_args=(sim, "Lorentzian"))
result = minner.minimize()
report_fit(result)


#%%
# After the fit, we can plot the new simulated spectrum.


residual_1 = sim.methods[0].experiment.copy()
residual_2 = sim.methods[1].experiment.copy()

residual_1[:] = result.residual[0:2046]
residual_2[:] = result.residual[2046:]

fig, ax = plt.subplots(2, 1, figsize=(6, 7))
ax[0].plot(*sim.methods[0].experiment.to_list(), label="Spectrum")
ax[0].plot(
    *(sim.methods[0].experiment - residual_1).to_list(), "r", alpha=0.5, label="Fit"
)
ax[0].plot(*residual_1.to_list(), alpha=0.5, label="Residual")

ax[0].set_xlabel("$^{29}$Si frequency / ppm")
ax[0].set_xlim(sim_x_1.value.max(), sim_x_1.value.min())
ax[0].grid(color="gray", linestyle="--", linewidth=0.5, alpha=0.5)

ax[1].plot(*sim.methods[1].experiment.to_list(), label="Spectrum")
ax[1].plot(
    *(sim.methods[1].experiment - residual_2).to_list(), "r", alpha=0.5, label="Fit"
)
ax[1].plot(*residual_2.to_list(), alpha=0.5, label="Residual")

ax[1].set_xlabel("$^{29}$Si frequency / ppm")
ax[1].set_xlim(sim_x_2.value.max(), sim_x_2.value.min())
ax[1].grid(color="gray", linestyle="--", linewidth=0.5, alpha=0.5)

plt.tight_layout()
plt.show()

#%%
# .. [#f1] Hansen, M. R., Jakobsen, H. J., Skibsted, J., :math:`^{29}\text{Si}`
#       Chemical Shift Anisotropies in Calcium Silicates from High-Field
#       :math:`^{29}\text{Si}` MAS NMR Spectroscopy, Inorg. Chem. 2003,
#       **42**, *7*, 2368-2377.
#       `DOI: 10.1021/ic020647f <https://doi.org/10.1021/ic020647f>`_


# %%
