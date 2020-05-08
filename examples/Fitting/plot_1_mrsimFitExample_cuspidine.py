#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
Fitting Cusipidine.
^^^^^^^^^^^^^^^^^^^
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
# Next we will import `csdmpy <https://csdmpy.readthedocs.io/en/latest/index.html>`_ and loading the data file.
import csdmpy as cp
import numpy as np

synthetic_experiment = cp.load("synthetic_cuspidine_test.csdf")
synthetic_experiment.dimensions[0].to("ppm", "nmr_frequency_ratio")

x1, y1 = synthetic_experiment.to_list()

plt.plot(x1, y1)
plt.xlabel("$^{29}$Si frequency / ppm")
plt.xlim(x1.value.max(), x1.value.min())
plt.grid(color="gray", linestyle="--", linewidth=0.5, alpha=0.5)
plt.tight_layout()
plt.show()

#%%
# In order to to fit a simulation to the data we will need to establish a ``Simulator`` object. We will
# use approximate initial parameters to generate our simulation:

#%%


from mrsimulator import Simulator, Isotopomer, Site
from mrsimulator.methods import BlochDecayFT

S29 = Site(
    isotope="29Si",
    isotropic_chemical_shift=-80.0,
    shielding_symmetric={"zeta": -60, "eta": 0.6},
)

method = BlochDecayFT(
    channel="29Si",
    magnetic_flux_density=7.1,
    rotor_frequency=780,
    dimensions=[{"count": 2046, "spectral_width": 25000, "reference_offset": -5000}],
)

sim = Simulator()
sim.isotopomers += [Isotopomer(name="Si29", sites=[S29], abundance=100)]
sim.methods += [method]
sim.run()
sim.methods[0].simulation.dimensions[0].to("ppm", "nmr_frequency_ratio")

x, y = sim.methods[0].simulation.to_list()

plt.plot(x, y)
plt.xlabel("$^{29}$Si frequency / ppm")
plt.xlim(x.value.max(), x.value.min())
plt.grid(color="gray", linestyle="--", linewidth=0.5, alpha=0.5)
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


from lmfit import Minimizer, Parameters, report_fit

params = Parameters()

params.add(name="iso", value=-80)
params.add(name="eta", value=0.6, min=0, max=1)
params.add(name="zeta", value=-60)
params.add(name="sigma", value=200)
params.add(name="factor", value=1)


#%%
# We will next set up an error function that will update the simulation throughout the minimization.
# We will construct a simple function here to demonstrate the *LMFIT* library, however, the next examples
# will showcase a fitting function provided in the *mrsimulator* library which automates the process.
from mrsimulator.apodization import Apodization


def test_function(params, data, sim):
    values = params.valuesdict()

    intensity = data.dependent_variables[0].components[0].real

    # Here, we update simulation parameters iso, eta, and zeta for the site object
    site = sim.isotopomers[0].sites[0]
    site.isotropic_chemical_shift = values["iso"]
    site.shielding_symmetric.eta = values["eta"]
    site.shielding_symmetric.zeta = values["zeta"]

    # here we run the simulation
    sim.run()
    # here we apodize the signal to simulate line broadening
    y = sim.apodize(Apodization.Lorentzian, sigma=values["sigma"])

    y_factored = y * values["factor"]

    return intensity - y_factored


#%%
# With the synthetic data, simulation, and the parameters we are ready to perform the fit. To fit, we use
# the *LMFIT* `Minimizer <https://lmfit.github.io/lmfit-py/fitting.html>`_ class.

#%%

minner = Minimizer(test_function, params, fcn_args=(synthetic_experiment, sim))
result = minner.minimize()
report_fit(result)

#%%
# After the fit, we can plot the new simulated spectrum.


plt.figsize = (4, 3)
residual = synthetic_experiment.copy()
residual[:] = result.residual
plt.plot(*synthetic_experiment.to_list(), label="Spectrum")
plt.plot(*(synthetic_experiment - residual).to_list(), "r", alpha=0.5, label="Fit")
plt.plot(*residual.to_list(), alpha=0.5, label="Residual")

plt.xlim(x1.value.max(), x1.value.min())
plt.xlabel("Frequency / Hz")
plt.grid(which="major", axis="both", linestyle="--")
plt.legend()

plt.tight_layout()
plt.show()

#%%
# .. [#f1] Hansen, M. R., Jakobsen, H. J., Skibsted, J., :math:`^{29}\text{Si}`
#       Chemical Shift Anisotropies in Calcium Silicates from High-Field
#       :math:`^{29}\text{Si}` MAS NMR Spectroscopy, Inorg. Chem. 2003,
#       **42**, *7*, 2368-2377.
#       `DOI: 10.1021/ic020647f <https://doi.org/10.1021/ic020647f>`_


# %%
