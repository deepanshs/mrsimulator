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
# extract the tensor parameters for our isotopomers. We will be using the *LMFIT* methods to
# establish fitting parameters and fit the spectrum. The following examples will show fitting with
# two synthetic :math:`^{29}\text{Si}` spectra--cuspidine and wollastonite--as well as the
# measurements from an :math:`^{17}\text{O}` experiment on :math:`\text{Na}_{2}\text{SiO}_{3}`.
# The *mrsimulator* library and data make use of CSDM compliant files.
# In this example we will be creating a synthetic spectrum of cuspidine from reported tensor
# parameters and then fit a simulation to the spectrum to demonstrate a simple fitting procedure.
#
# Begin by importing *csdmpy* and loading the data file. We will also be establishing *matplotlib*
# figure size here.
import csdmpy as cp
import matplotlib.pylab as pylab

params = {"figure.figsize": (4.5, 3), "font.size": 9}
pylab.rcParams.update(params)

synthetic_experiment = cp.load("synthetic_cuspidine_test.csdf")

cp.plot(synthetic_experiment)

#%%
# In order to to fit a simulation to the data we will need to establish a ``simulation`` object. We will
# use approximate initial parameters to generate our simulation:

#%%


from mrsimulator import Simulator, Isotopomer, Site, Dimension
from mrsimulator.methods import one_d_spectrum

S29 = Site(
    isotope="29Si",
    isotropic_chemical_shift=-80.0,
    shielding_symmetric={"zeta": -60, "eta": 0.6},
)

dimension = Dimension(
    isotope="29Si",
    magnetic_flux_density=7.1,  # in T
    number_of_points=2046,
    spectral_width=25000,  # in Hz
    reference_offset=-10000,  # in Hz
    rotor_frequency=780,  # in Hz
)

sim = Simulator()
sim.isotopomers += [Isotopomer(name="Si29", sites=[S29], abundance=100)]
sim.dimensions += [dimension]
sim.run(method=one_d_spectrum)
sim_data = sim.as_csdm_object()

cp.plot(sim_data)

#%%
# Next, we will need a list of parameters that will be used in the fit. the *LMFIT* library allows us to create
# a list of parameters rather easily using the ``Parameters()`` class. We have created a function to parse the
# simulation object for available parameters and construct an *LMFIT* ``Parameter`` object which is shown in the next two
# examples on fitting. Here, however, we will construct the parameter list explicitly to demonstrate how the parameters
# are created.
#
# One thing to note is that the names of our parameters must correspond to their addresses within the simulation object
# in order to update the simulation during the fit. The *LMFIT* library does not allow for the use of special characters
# such as "\[", "\]", or "." so our current workaround is converting the special characters to their corresponding HTML
# character code numbers and converting back to the special character when updating the simulation.

#%%


from lmfit import Minimizer, Parameters, report_fit

params = Parameters()

params.add(name="iso", value=-80)
params.add(name="eta", value=0.6, min=0, max=1)
params.add(name="zeta", value=-60)
params.add(name="sigma", value=0)
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
    sim.run(one_d_spectrum)
    # here we apodize the signal to simulate line broadening
    y = sim.apodize(Apodization.Lorentzian, sigma=values["sigma"])

    y_factored = y * values["factor"]

    return intensity - y_factored


#%%
# With the synthetic data, simulation, and the parameters we are ready to perform the fit. To fit, we use
# the *LMFIT* ``Minimizer`` class.

#%%

minner = Minimizer(test_function, params, fcn_args=(synthetic_experiment, sim))
result = minner.minimize()
report_fit(result)

#%%
# After the fit, we can plot the new simulated spectrum using the *matplotlib* library.


import matplotlib.pyplot as plt


plt.figsize = (4, 3)
residual = synthetic_experiment.copy()
residual[:] = result.residual
plt.plot(*synthetic_experiment.to_list(), label="Spectrum")
plt.plot(*(synthetic_experiment - residual).to_list(), "r", alpha=0.5, label="Fit")
plt.plot(*residual.to_list(), alpha=0.5, label="Residual")

plt.xlabel("Frequency / Hz")
plt.grid(which="major", axis="both", linestyle="--")
plt.legend()

plt.tight_layout()
plt.show()
