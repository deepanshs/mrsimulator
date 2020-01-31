#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
Fitting a simulation to a synthetic spectrum - Wollastonite.
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
"""
#%%
# Often, after obtaining an NMR measurement we must fit tensors to our data so we can obtain the tensor parameters. In this example, we will illustrate the use of the *mrsimulator* method to simulate the experimental spectrum and fit the simulation to the data allowing us to extract the tensor parameters for our isotopomers. We will be using the *LMFIT* methods to establish fitting parameters and fit the spectrum. The following examples will show the measurements from an :math:`^{17}\text{O}` experiment on :math:`\text{Na}_{2}\text{SiO}_{3}` and a :math:`^{29}\text{Si}` experiment on wollastonite. The *mrsimulator* library and data make use of CSDM compliant files.
#
# In this example we will be creating a synthetic spectrum of wollastonite and attempting to fit a simulation to the synthetic spectrum. We will create the synthetic spectrum using tensor parameters from wollastonite and ensure that our fitting method is able to fit the spectrum.
#
# As before with :math:`\text{Na}_{2}\text{SiO}_{3}`, we will be creating simulation objects
#%%
import matplotlib.pylab as pylab

params = {"figure.figsize": (4.5, 3)}
pylab.rcParams.update(params)
import matplotlib
import matplotlib.pyplot as plt
import csdmpy as cp
from mrsimulator import Simulator, Isotopomer, Site, Dimension
from mrsimulator import SymmetricTensor as st
from mrsimulator.methods import one_d_spectrum

S29_1 = Site(
    isotope="29Si",
    isotropic_chemical_shift=-89.0,
    shielding_symmetric={"zeta": 59.8, "eta": 0.62},
)
S29_2 = Site(
    isotope="29Si",
    isotropic_chemical_shift=-89.5,
    shielding_symmetric={"zeta": 52.1, "eta": 0.68},
)
S29_3 = Site(
    isotope="29Si",
    isotropic_chemical_shift=-87.8,
    shielding_symmetric={"zeta": 69.4, "eta": 0.60},
)

dimension = Dimension(
    isotope="29Si",
    magnetic_flux_density=14.1,  # in T
    number_of_points=2046,
    spectral_width=25000,  # in Hz
    reference_offset=-10000,  # in Hz
    rotor_frequency=0,  # in Hz
)


isotopomers = [Isotopomer(sites=[site]) for site in [S29_1, S29_2, S29_3]]

synth_wollastonite = Simulator()
synth_wollastonite.isotopomers += isotopomers
synth_wollastonite.dimensions += [dimension]

synth_wollastonite.run(method=one_d_spectrum)
synth_data = synth_wollastonite.as_csdm_object()
cp.plot(synth_data)
plt.show()

#%%
# In order to create a synthetic spectrum we will add noise and line broadening to the
# spectrum. The line broadening is performed using the Apodization class which is a
# member of a simulator object.

#%%


from mrsimulator.apodization import Apodization
import numpy as np

synth_data = synth_wollastonite.as_csdm_object()

x = synth_data.dimensions[0]
y = synth_data.dependent_variables[0].components[0]
synth_data.dependent_variables[0].components[0] = synth_wollastonite.apodize(
    Apodization.Lorentzian, sigma=200
)

noise = np.random.normal(-0.025, 0.025, y.shape)
synthetic = y + noise
synth_data.dependent_variables[0].components[0] = synthetic

cp.plot(synth_data)

#%%
# We have just created a synthetic experimental data set and now we must create the simulation to fit to this synthetic data. We will use another initial guess to generate our simulation:

#%%


S29_1 = Site(
    isotope="29Si",
    isotropic_chemical_shift=-91.0,
    shielding_symmetric={"zeta": 60, "eta": 0.6},
)
S29_2 = Site(
    isotope="29Si",
    isotropic_chemical_shift=-90,
    shielding_symmetric={"zeta": 50, "eta": 0.70},
)
S29_3 = Site(
    isotope="29Si",
    isotropic_chemical_shift=-88,
    shielding_symmetric={"zeta": 70, "eta": 0.55},
)

isotopomers = [Isotopomer(sites=[site]) for site in [S29_1, S29_2, S29_3]]

sim = Simulator()
sim.isotopomers += isotopomers
sim.dimensions += [dimension]
sim.run(method=one_d_spectrum)
sim_data = sim.as_csdm_object()

cp.plot(sim_data)

#%%
# Again, we will need a list of parameters to change during the fitting:

#%%


from lmfit import Minimizer, Parameters, report_fit
from mrsimulator import spectral_fitting

params = Parameters()
params = spectral_fitting.make_fitting_parameters(sim)
params.add(name="sigma", value=0)
params.add(
    name="factor", value=sim_data.dependent_variables[0].components[0].max(), min=0
)

#%%
# With the synthetic data, simulation, and the parameters we are ready to perform the fit

#%%


minner = Minimizer(
    spectral_fitting.fcn2min,
    params,
    fcn_args=(synth_data.dependent_variables[0].components[0].real, sim, "Lorentzian"),
)
result = minner.minimize()
report_fit(result)


#%%


sim.run(method=one_d_spectrum)
sim_data = sim.as_csdm_object()
values = params.valuesdict()


x_simulated = sim_data.dimensions[0]
y_simulated = sim_data.dependent_variables[0].components[0]
sim_data.dependent_variables[0].components[0] = spectral_fitting.line_broadening(
    x_simulated, y_simulated, result.params["sigma"], 0
).real

plt.plot(
    x.coordinates, synth_data.dependent_variables[0].components[0], label="Spectrum"
)
plt.plot(
    x_simulated.coordinates,
    sim_data.dependent_variables[0].components[0] * result.params["factor"],
    "r",
    alpha=0.5,
    label="Fit",
)
plt.plot(
    x.coordinates,
    synth_data.dependent_variables[0].components[0].real
    - sim_data.dependent_variables[0].components[0] * result.params["factor"],
    alpha=0.5,
    label="Residual",
)
plt.grid(which="major", axis="both", linestyle="--")
plt.xlabel("Frequency / Hz")
plt.legend()

plt.show()
