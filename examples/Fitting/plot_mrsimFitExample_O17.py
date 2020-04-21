#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
Fitting Sodium Silicate.
^^^^^^^^^^^^^^^^^^^^^^^^
"""
# sphinx_gallery_thumbnail_number = 3
#%%
# Often, after obtaining an NMR measurement we must fit tensors to our data so we can
# obtain the tensor parameters. In this example, we will illustrate the use of the
# *mrsimulator* method to simulate the experimental spectrum and fit the simulation to
# the data allowing us to extract the tensor parameters for our isotopomers. We will
# be using the *LMFIT* methods to establish fitting parameters and fit the spectrum.
# The following examples will show the measurements from an :math:`^{17}\text{O}`
# experiment on :math:`\text{Na}_{2}\text{SiO}_{3}` and a :math:`^{29}\text{Si}`
# experiment on wollastonite. The *mrsimulator* library and data make use of CSDM
# compliant files.
#
# Import the *csdmpy* model and load the dataset. We will also be importing *matplotlib* to
# adjust the plotting output.
import csdmpy as cp
import matplotlib.pylab as pylab
import numpy as np

params = {"figure.figsize": (4.5, 3), "font.size": 9}
pylab.rcParams.update(params)


filename = "https://osu.box.com/shared/static/kfgt0jxgy93srsye9pofdnoha6qy58qf.csdf"
oxygen_experiment = cp.load(filename).real

cp.plot(oxygen_experiment)

#%%
# Next, we will want to create a simulation object that we will use to fit to our
# spectrum. We will need to import the necessary libraries for the *mrsimulator*
# methods. We will then create isotopomer objects. While the tensor parameters we
# initialize are not exact, they serve as a good initial guess.

#%%


from mrsimulator import Simulator, Isotopomer, Site, Dimension
from mrsimulator import SymmetricTensor as st

sim = Simulator()
O17_1 = Site(
    isotope="17O",
    isotropic_chemical_shift=60.0,
    quadrupolar={"Cq": 4200000, "eta": 0.5},
)
O17_2 = Site(
    isotope="17O", isotropic_chemical_shift=40, quadrupolar={"Cq": 2400000, "eta": 0.18}
)

isotopomers = [Isotopomer(sites=[site]) for site in [O17_1, O17_2]]

dimension = Dimension(
    isotope="17O",
    magnetic_flux_density=9.4,  # in T
    number_of_points=oxygen_experiment.dimensions[0].count,  # 2046,
    spectral_width=oxygen_experiment.dimensions[0].count
    * oxygen_experiment.dimensions[0].increment.to("Hz").value,  # 20000, # in Hz
    reference_offset=oxygen_experiment.dimensions[0]
    .coordinates_offset.to("Hz")
    .value,  # 1220.7031, # in Hz
    rotor_frequency=14000,  # in Hz
)

sim.isotopomers += isotopomers
sim.dimensions += [dimension]
sim.run()
sim_data = sim.as_csdm_object()

cp.plot(sim_data)

#%%
# Once we have our simulation we must create our list of parameters to use in our
# fitting. We will be using the Parameters() class from LMFIT. For more information on
# the LMFIT library please visit [LMFIT documentation](https://lmfit.github.io/lmfit-py/).

#%%

from lmfit import Minimizer, Parameters, report_fit

params = Parameters()

#%%
# In each experiment the number of isotopomers and sites present as well as their
# attributes may vary. To simplify the parameter list creation we will use the
# *make_fitting_parameters()* method from *mrsimulator.spectral_fitting*

#%%


from mrsimulator import spectral_fitting

params = spectral_fitting.make_fitting_parameters(sim)
params.add(
    name="sigma", value=oxygen_experiment.dimensions[0].increment.to("Hz").value, min=0
)
params.add(
    name="factor",
    value=oxygen_experiment.dependent_variables[0].components[0].max().real,
    min=0,
)

#%%
# With an experimental spectrum, a simulaton, and a list of parameters we are now
# ready to perform a fit. This fit will be performed using the *LMFIT* library as
# well as our *fcn2min* function from *mrsimulator.spectral_fitting*. The arguments
# for *fcn2min* are the intensities from the experimental data and the simulation
# CSDM object. Reporting the results of the fit will give us our tensor parameters.
#
# The naming convention of the variables may seem a bit odd at first. LMFIT
# parameter names do not accept special characters such as "\[," "\], "or "."
# in the names of the parmeters, so in order to maintain the addresses of each
# parameter in the simulator object we converted the square brakets and periods
# to their HTML numbers. For example, isotopomers\[0\].sites\[0\].isotropic_chemical_shift
# is converted to isotopomers9109346sites9109346isotropic_chemical_shift, and during
# any updates to the simualtion object the name is converted back to the address before
# execution of the assignments.

#%%


minner = Minimizer(
    spectral_fitting.fcn2min,
    params,
    fcn_args=(
        oxygen_experiment.dependent_variables[0].components[0].real,
        sim,
        "Lorentzian",
    ),
)
result = minner.minimize()
report_fit(result)

#%%
# Next, we can compare the fit to the experimental data using the *matplotlib* library:

#%%


import matplotlib
import matplotlib.pyplot as plt
from mrsimulator.apodization import Apodization

sim.run()
sim_data = sim.as_csdm_object()
values = params.valuesdict()

x_simulated = sim_data.dimensions[0]
y_simulated = sim_data.dependent_variables[0].components[0]
sim_data.dependent_variables[0].components[0] = sim.apodize(
    Apodization.Lorentzian, sigma=result.params["sigma"]
)


plt.plot(
    oxygen_experiment.dimensions[0].coordinates,
    oxygen_experiment.dependent_variables[0].components[0].real,
    label="Experimental",
)
plt.plot(
    x_simulated.coordinates,
    sim_data.dependent_variables[0].components[0] * result.params["factor"],
    "r",
    alpha=0.5,
    label="Fit",
)
plt.plot(
    x_simulated.coordinates,
    oxygen_experiment.dependent_variables[0].components[0]
    - sim_data.dependent_variables[0].components[0] * result.params["factor"],
    alpha=0.5,
    label="Residual",
)
plt.xlabel("Frequency / Hz")
plt.legend()
plt.show()
